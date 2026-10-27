"""Reproducible whole-observation distribution references for epsilon budgets.

Sample present, unflagged finite complex values independently by column,
configuration, correlation and auto/cross class. SplitMix64 is a stateless
sampling hash, not decoder randomness. Small strata use the minimum-priority
sample rather than an absolute-valued error floor. No source writes.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from .error_budget import NORMALIZED_DEFAULTS
from .inspection import _open, CORRELATIONS

CONFIG_COLUMNS = ('DATA_DESC_ID', 'FIELD_ID', 'ARRAY_ID',
                  'PROCESSOR_ID', 'STATE_ID', 'FEED1', 'FEED2')
SCHEMA = 'varquanms.distribution-references/1'
RECORD = np.dtype([('index', '<u8'), ('value', '<c8')])


def key(context):
    return json.dumps(context, sort_keys=True, separators=(',', ':'))


def table_context(t):
    if not t.nrows():
        raise ValueError('Distribution references require nonempty MAIN')
    return {c:int(t.getcell(c, 0)) for c in CONFIG_COLUMNS if c in t.colnames()}


def configuration_runs(t):
    names = [c for c in CONFIG_COLUMNS if c in t.colnames()]
    changed = np.zeros(t.nrows()-1, dtype=bool)
    for name in names:
        a = np.asarray(t.getcol(name))
        if a.ndim != 1 or a.dtype.kind not in 'iu':
            raise ValueError('Configuration identities must be integer scalars')
        changed |= a[1:] != a[:-1]
    starts = np.r_[0, np.flatnonzero(changed)+1]
    stops = np.r_[starts[1:], t.nrows()]
    if len(starts) > 4096:
        raise ValueError('More than 4096 configuration runs are unsupported')
    identities = {name:np.asarray(t.getcol(name))[starts] for name in names}
    return [(int(a),int(b),{c:int(identities[c][i]) for c in names})
            for i,(a,b) in enumerate(zip(starts,stops))]


def priorities(indices, stream):
    x = np.asarray(indices, dtype=np.uint64) ^ np.uint64(stream)
    with np.errstate(over='ignore'):
        x = x + np.uint64(0x9E3779B97F4A7C15)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def _bytes(indices, values):
    records = np.empty(len(indices), dtype=RECORD)
    records['index'] = indices; records['value'] = values
    return records.tobytes()


class Sampler:
    def __init__(self, fraction, minimum, stream):
        self.fraction, self.minimum, self.stream = fraction, minimum, stream
        self.count = self.selected = 0
        self.mean = np.zeros(2); self.m2 = np.zeros(2)
        self.digest = hashlib.sha256()
        self.ids = np.empty(0,dtype='<u8'); self.ranks = self.ids.copy()
        self.values = np.empty(0,dtype='<c8')

    def add(self, indices, values):
        if not len(values): return
        self.count += len(values)
        ranks = priorities(indices, self.stream)
        selected = np.ones(len(values),bool) if self.fraction == 1 else ranks < int(self.fraction*2**64)
        a = values[selected]
        if len(a):
            self.digest.update(_bytes(indices[selected],a))
            x = np.column_stack((a.real,a.imag)).astype(np.float64)
            mean = x.mean(axis=0); m2 = np.sum((x-mean)**2,axis=0)
            delta = mean-self.mean; total = self.selected+len(a)
            self.m2 += m2+delta**2*(self.selected*len(a)/total)
            self.mean += delta*(len(a)/total); self.selected = total
        # Keep a bounded minimum-size, uniform priority sample for small strata.
        if len(values) > self.minimum:
            keep = np.argpartition(ranks,self.minimum-1)[:self.minimum]
            ranks, indices, values = ranks[keep], indices[keep], values[keep]
        ranks = np.r_[self.ranks,ranks]; indices = np.r_[self.ids,indices]; values = np.r_[self.values,values]
        keep = np.lexsort((indices,ranks))[:self.minimum]
        self.ranks,self.ids,self.values = ranks[keep],indices[keep],values[keep]

    def finish(self):
        minimum = self.selected < min(self.minimum,self.count)
        count = len(self.values) if minimum else self.selected
        if minimum:
            order = np.argsort(self.ids); a = self.values[order]
            x = np.column_stack((a.real,a.imag)).astype(np.float64)
            mean = x.mean(axis=0); var = np.var(x,axis=0)
            digest = hashlib.sha256(_bytes(self.ids[order],a)).hexdigest()
        else:
            mean = self.mean; var = self.m2/count if count else np.zeros(2)
            digest = self.digest.hexdigest()
        return dict(eligible_samples=self.count, fraction_selected_samples=self.selected,
            sampled_samples=count, actual_fraction=count/self.count if self.count else 0.,
            minimum_sample_used=minimum, sample_sha256=digest,
            mean_real=float(mean[0]),mean_imag=float(mean[1]),
            std_real=float(np.sqrt(var[0])) if count else None,
            std_imag=float(np.sqrt(var[1])) if count else None,
            component_std=float(np.sqrt(var.mean())) if count else None)


def sample_references(source, config):
    if config.normalized_bounds is None:
        raise ValueError('Reference sampling requires normalized_bounds')
    started = time.perf_counter()
    source = Path(source)
    policy = dict(NORMALIZED_DEFAULTS, **config.normalized_bounds)
    groups, samplers = {}, {}
    with _open(source/'DATA_DESCRIPTION') as dd, _open(source/'POLARIZATION') as pol:
        products = {i:[CORRELATIONS[int(c)] for c in pol.getcell('CORR_TYPE',int(dd.getcell('POLARIZATION_ID',i)))]
                    for i in range(dd.nrows())}
    with _open(source) as t:
        if not t.nrows(): raise ValueError('Empty MS')
        runs = configuration_runs(t)
        for start,stop,context in runs:
            group_key = key(context); names = products[context['DATA_DESC_ID']]
            groups.setdefault(group_key,dict(context=context,columns={}))
            for lo in range(start,stop,config.io_chunk_rows):
                count = min(stop-lo,config.io_chunk_rows)
                flags = np.asarray(t.getcol('FLAG',lo,count),bool)
                if 'FLAG_ROW' in t.colnames():
                    flags |= np.asarray(t.getcol('FLAG_ROW',lo,count),bool)[:,None,None]
                auto = t.getcol('ANTENNA1',lo,count)==t.getcol('ANTENNA2',lo,count)
                for column in config.columns:
                    a = np.asarray(t.getcol(column,lo,count))
                    if a.dtype != np.dtype('complex64') or a.ndim != 3 or a.shape != flags.shape or a.shape[2] != len(names):
                        raise ValueError('Sampling requires regular complex64 row/channel/correlation cells within each configuration')
                    if stop*a.shape[1] >= 2**64: raise ValueError('Sampling index overflow')
                    indices = np.arange(lo*a.shape[1],(lo+count)*a.shape[1],dtype='<u8').reshape(a.shape[:2])
                    for p,name in enumerate(names):
                        valid = ~flags[:,:,p] & np.isfinite(a[:,:,p])
                        for label,rows in [('auto',auto),('cross',~auto)]:
                            identity = (group_key,column,name,label)
                            if identity not in samplers:
                                stream = int.from_bytes(hashlib.sha256(json.dumps([policy['sampling_seed'],*identity],separators=(',',':')).encode()).digest()[:8],'little')
                                samplers[identity] = Sampler(policy['sampling_fraction'],policy['minimum_samples'],stream)
                            mask = valid & rows[:,None]
                            samplers[identity].add(indices[mask],a[:,:,p][mask])
        source_rows = t.nrows()
    for (group_key,column,name,label),sampler in samplers.items():
        groups[group_key]['columns'].setdefault(column,{}).setdefault(name,{})[label] = sampler.finish()
    return dict(schema=SCHEMA,source_rows=source_rows,groups=groups,
        sampling=dict(algorithm='splitmix64-row-channel-product-v1',fraction=policy['sampling_fraction'],
                      seed=policy['sampling_seed'],minimum_samples=policy['minimum_samples']),
        definition='sqrt((population_var(real)+population_var(imag))/2); float64 merged moments; unflagged finite original complex samples; not thermal noise',
        elapsed_seconds=time.perf_counter()-started)


def column_references(t, column, references):
    if references is None or references.get('schema') != SCHEMA:
        raise ValueError('Missing distribution reference inventory')
    context = key(table_context(t))
    try:
        result = references['groups'][context]['columns'][column]
    except KeyError as exc:
        raise ValueError('Missing configuration/column reference') from exc
    for classes in result.values():
        for item in classes.values():
            scale = item['component_std']
            if scale is not None and (not math.isfinite(scale) or scale < 0):
                raise ValueError('Invalid sampled component std')
    return result


def validate_reference_manifest(manifest):
    """Check normalized-policy provenance before accepting an archive."""
    configuration = manifest.get('config',{})
    if configuration.get('normalized_bounds') is None:
        return
    from .config import Config
    Config(**configuration).validate()
    policy = dict(NORMALIZED_DEFAULTS,**configuration['normalized_bounds'])
    references = manifest.get('distribution_references')
    if not isinstance(references,dict) or references.get('schema') != SCHEMA:
        raise ValueError('Missing or unsupported distribution references')
    expected = dict(algorithm='splitmix64-row-channel-product-v1',fraction=policy['sampling_fraction'],
                    seed=policy['sampling_seed'],minimum_samples=policy['minimum_samples'])
    if references.get('sampling') != expected:
        raise ValueError('Distribution sampling policy mismatch')
    if type(references.get('source_rows')) is not int or references['source_rows'] < manifest['source']['main_row_count']:
        raise ValueError('Invalid distribution reference source row count')
    groups = references.get('groups')
    if not isinstance(groups,dict) or not groups:
        raise ValueError('Empty distribution reference groups')
    for identity,group in groups.items():
        if key(group['context']) != identity:
            raise ValueError('Distribution reference identity mismatch')
        if set(group['columns']) != set(configuration['columns']):
            raise ValueError('Distribution reference columns mismatch')
        for products in group['columns'].values():
            for classes in products.values():
                if set(classes) != {'auto','cross'}:
                    raise ValueError('Missing distribution reference baseline class')
                for item in classes.values():
                    total,selected,count = (item[k] for k in ('eligible_samples','fraction_selected_samples','sampled_samples'))
                    if any(type(n) is not int or n<0 for n in (total,selected,count)) or max(selected,count)>total or (total>0 and count==0):
                        raise ValueError('Invalid distribution sample count')
                    if item['actual_fraction'] != (count/total if total else 0.):
                        raise ValueError('Invalid distribution sampled fraction')
                    if len(item['sample_sha256']) != 64 or any(c not in '0123456789abcdef' for c in item['sample_sha256']):
                        raise ValueError('Invalid distribution sample checksum')
                    for name in ('std_real','std_imag','component_std'):
                        value=item[name]
                        if ((count==0 and value is not None) or
                            (count>0 and (type(value) not in (float,int) or not math.isfinite(value) or value<0))):
                            raise ValueError('Invalid distribution reference std')
                    if count and not math.isclose(item['component_std'],math.hypot(item['std_real'],item['std_imag'])/math.sqrt(2),rel_tol=1e-14,abs_tol=0.):
                        raise ValueError('Distribution component std convention differs')
    for table in manifest['tables']:
        for name,column in table['columns'].items():
            if column.get('encoding') != 'quantized': continue
            identity=column.get('distribution_reference_key')
            if identity not in groups or column.get('distribution_references') != groups[identity]['columns'].get(name):
                raise ValueError('Visibility distribution references are missing or inconsistent')
            if set(column['distribution_references']) != set(column['product_names']):
                raise ValueError('Distribution correlation identities differ')
