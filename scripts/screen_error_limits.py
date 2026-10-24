"""Matched-block absolute/relative/hybrid sweep and rate-matched full-test plan.

Thresholds for full imaging are selected using block bytes only, before any
candidate image is examined. All candidates share the same codec and groups.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import json
from pathlib import Path
import time

import numpy as np

from varquanms.codec import encode_block
from varquanms.config import Config
from varquanms.error_budget import component_signal_power
from varquanms.inspection import _open, CORRELATIONS
from varquanms.ordering import all_orders
from varquanms.partitioned import context_runs


def load_blocks(ms, base, groups=64):
    _, runs = context_runs(ms)
    jobs = []
    for integration in sorted({0, len(runs)//2, len(runs)-1}):
        start, stop = runs[integration]
        with _open(ms) as t:
            data = t.getcol('CORRECTED_DATA', start, stop-start)
            flags = t.getcol('FLAG', start, stop-start) | t.getcol('FLAG_ROW', start, stop-start)[:,None,None]
            pairs = np.column_stack([t.getcol(c, start, stop-start) for c in ('ANTENNA1','ANTENNA2')])
            uvw = t.getcol('UVW', start, stop-start)
            ddids = np.unique(t.getcol('DATA_DESC_ID', start, stop-start))
            from varquanms.references import CONFIG_COLUMNS, key
            configuration_key = key({c:int(t.getcell(c,start)) for c in CONFIG_COLUMNS if c in t.colnames()})
        if len(ddids) != 1 or data.shape[1] != 48:
            raise ValueError('This study requires one 48-channel configuration per time context')
        with _open(ms/'DATA_DESCRIPTION') as t: polid = t.getcell('POLARIZATION_ID', int(ddids[0]))
        with _open(ms/'POLARIZATION') as t: names = [CORRELATIONS[int(v)] for v in t.getcell('CORR_TYPE', int(polid))]
        data[flags] = 0j
        orders, _, _ = all_orders(data, flags, pairs, uvw, base)
        for product, order in enumerate(orders):
            bounds = order['group_boundaries']
            for index in np.unique(np.linspace(0,len(bounds)-2,min(groups,len(bounds)-1)).astype(int)):
                rows = order['permutation'][bounds[index]:bounds[index+1]]
                for lo in range(0,data.shape[1],base.channel_segment):
                    hi=min(data.shape[1],lo+base.channel_segment)
                    a, f = data[rows,lo:hi,product], flags[rows,lo:hi,product]
                    jobs.append(dict(integration=integration, product=names[product], configuration_key=configuration_key,
                        group=int(index),channel_lo=lo,channel_hi=hi,
                        auto=bool(np.all(pairs[rows,0] == pairs[rows,1])),
                        median_baseline_m=float(np.median(np.linalg.norm(uvw[rows],axis=1))),
                        signal_component_power=component_signal_power(a,f), values=a, flags=f))
    return jobs


def configuration(base, family, parameter, pivot_power):
    # One parameter controls an RMS ceiling and matching max/bias/row guards.
    # Relative parameter is a fraction; hybrid parameter is its absolute MSE cap.
    options = dict(max_component_mse=None, max_complex_error=None,
                   max_complex_bias=None, max_baseline_component_mse=None,
                   max_p99_error=None, relative_bounds=None)
    if family in ('absolute','hybrid'):
        options.update(max_component_mse=parameter, max_complex_error=float(10*np.sqrt(parameter)),
                       max_complex_bias=float(np.sqrt(parameter)), max_baseline_component_mse=2*parameter)
    if family in ('relative','hybrid'):
        r = parameter if family == 'relative' else np.sqrt(parameter/pivot_power)
        options['relative_bounds'] = dict(reference='block-rms-v1',max_component_rmse=float(r),
            max_complex_error=float(10*r),max_complex_bias=float(r),max_baseline_component_rmse=float(np.sqrt(2)*r))
    cfg = replace(base, **options)
    cfg.validate()
    return cfg


def evaluate(jobs, config):
    stamp = time.perf_counter()
    def run(job):
        b = encode_block(job['values'],job['flags'],config,config.bounds(job['product']))
        return dict(bytes=len(b.blob),samples=b.decoded.size,bits=b.bits,
                    gamma=b.parameters.get('gamma'),exceptions=b.exceptions,metrics=b.metrics)
    with ThreadPoolExecutor(max_workers=8) as pool: blocks = list(pool.map(run,jobs))
    n = sum(b['metrics']['valid_samples'] for b in blocks)
    error = sum(b['metrics']['component_mse']*b['metrics']['valid_samples'] for b in blocks)
    signal = sum(j['signal_component_power']*b['metrics']['valid_samples'] for j,b in zip(jobs,blocks))
    relative = [b['metrics']['component_mse']/j['signal_component_power']
                for j,b in zip(jobs,blocks) if j['signal_component_power'] > 0]
    return dict(config=config.to_dict(),payload_bytes=sum(b['bytes'] for b in blocks),
                component_mse=error/max(1,n),global_relative_rmse=float(np.sqrt(error/signal)) if signal else 0.,
                mean_block_relative_mse=float(np.mean(relative)) if relative else 0.,
                block_bits=dict(Counter(str(b['bits']) if b['bits'] else 'float32' for b in blocks)),
                elapsed_seconds=time.perf_counter()-stamp,blocks=blocks)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('ms',type=Path); p.add_argument('output',type=Path)
    args = p.parse_args(); args.output.mkdir(exist_ok=False)
    base = replace(Config.load('configs/calibrated-gamma-ratio.json'),workers=24)
    jobs = load_blocks(args.ms,base)
    # Common median of nonzero cross-correlation block powers, all products.
    pivot = float(np.median([j['signal_component_power'] for j in jobs if not j['auto'] and j['signal_component_power'] > 0]))
    cache = {}
    def measure(family, value):
        key = (family,float(value))
        if key not in cache:
            r = evaluate(jobs,configuration(base,family,float(value),pivot))
            r.update(family=family,parameter=float(value))
            cache[key] = r
            print(json.dumps({k:v for k,v in r.items() if k not in ('blocks','config')}),flush=True)
        return cache[key]
    for family, values in [('absolute',[1,16,64,256,1024,4096,16384,65536]),
                           ('relative',[.0005,.001,.002,.004,.008,.016,.032,.064,.128]),
                           ('hybrid',[16,64,256,1024,4096,16384])]:
        for value in values: measure(family,value)
    selected=[]
    for budget in (64,256,4096):
        target = measure('absolute',budget)['payload_bytes']
        for family in ('absolute','relative','hybrid'):
            if family == 'absolute': chosen = measure(family,budget)
            else:
                low,high = ((1e-5,1.) if family == 'relative' else (1e-2,1e7))
                # Bisection on byte counts only, retaining closest actual point.
                options=[]
                for _ in range(12):
                    middle = np.sqrt(low*high)
                    candidate = measure(family,middle); options.append(candidate)
                    if candidate['payload_bytes'] > target: low=middle
                    else: high=middle
                chosen=min(options,key=lambda r:abs(r['payload_bytes']-target))
            label=f'{family}-t{budget}'
            cfgpath=args.output/(label+'.json'); cfgpath.write_text(json.dumps(chosen['config'],indent=2)+'\n')
            selected.append(dict(name=label,family=family,parameter=chosen['parameter'],
                target_absolute_mse=budget,target_screen_bytes=target,screen_bytes=chosen['payload_bytes'],
                screen_rate_difference_fraction=chosen['payload_bytes']/target-1,config=chosen['config']))
    report=dict(schema='varquanms.error-limit-screen/1',source=str(args.ms),pivot_component_power=pivot,
        blocks=[{k:v for k,v in j.items() if k not in ('values','flags')} for j in jobs],
        experiments=list(cache.values()),selected=selected,
        selection='Three absolute targets; 12-step log bisection matches block payload bytes for relative and hybrid. No images used to choose thresholds.',
        limitations=['Block-relative RMS includes coherent signal; not samplewise fractional error or thermal-noise accuracy.',
                     'Screen bytes include headers, exceptions and secondary compression, not complete archive overhead.',
                     'Same 768 blocks from integrations 0,7,14 for every trial; all full tests use all 15 integrations.'])
    (args.output/'screen.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (args.output/'plan.json').write_text(json.dumps(selected,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__': main()
