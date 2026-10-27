"""Measure experimental MX payloads; HDF5 holds indexed binary block records.

This writes a visibility-component experiment, not a standalone MS archive.
Optional MS export requires the original source for exact auxiliary information.
"""

import argparse
from collections import Counter
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import platform
import time

from casacore.tables import table
import h5py
import numpy as np

from vamp.archive import compare_ms
from vamp.codec import encode_block as uniform_encode, decode_block as uniform_decode, error_metrics, within_bounds
from vamp.config import Config
from vamp.experimental_mx import encode_batch, decode_block, decode_component
from vamp.inspection import inspect_ms
from vamp.ordering import all_orders
from vamp.visibility import _flags, quality


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(4*1024**2):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ms', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--variant', choices=['unit-mega', 'max-mega', 'unit-micro', 'max-micro', 'uniform'], required=True)
    parser.add_argument('--sample-groups', type=int, default=0, help='evenly spaced groups per product; zero evaluates the entire MS')
    parser.add_argument('--batch', type=int, default=128)
    parser.add_argument('--export-ms', action='store_true')
    parser.add_argument('--ordering', choices=['amplitude', 'class-original', 'length', 'length-orientation'], default='amplitude')
    parser.add_argument('--mse', type=float, help='override component MSE; per-baseline ceiling becomes twice this')
    args = parser.parse_args()
    if args.output.exists() or args.batch < 1 or args.sample_groups < 0 or args.export_ms and args.sample_groups:
        raise ValueError('Use a new output directory; export requires the full observation')
    source = args.ms.resolve()
    if args.output.resolve().is_relative_to(source):
        raise ValueError('Output must be outside source')
    args.output.mkdir(parents=True)
    started = time.perf_counter()
    config = Config.load(Path(__file__).resolve().parents[1]/'configs/performance-adaptive.json')
    config = replace(config, baseline_group=32, channel_segment=24, ordering=args.ordering, workers=1)
    if args.mse is not None:
        config=replace(config,max_component_mse=args.mse,max_baseline_component_mse=2*args.mse)
    config.validate()
    bounds = config.bounds()
    report = inspect_ms(source, ['DATA'])
    if report['layout']['initial_regular_profile_issues'] or report['layout']['organization'] != 'channels_inside_cells':
        raise ValueError('MX experiment requires the inspected regular single-configuration layout')
    names = report['layout']['configurations'][0]['correlation_names']
    with table(str(source), readonly=True, ack=False) as t:
        values = t.getcol('DATA')
        if values.dtype != np.complex64:
            raise ValueError('MX experiment requires complex64')
        flags = _flags(t)
        pairs = np.column_stack([t.getcol('ANTENNA1'), t.getcol('ANTENNA2')])
        uvw = t.getcol('UVW')
    read_seconds = time.perf_counter()-started
    stamp = time.perf_counter()
    orders, _, _ = all_orders(values, flags, pairs, uvw, config)
    jobs = []
    for p, order in enumerate(orders):
        groups = list(zip(order['group_boundaries'][:-1], order['group_boundaries'][1:]))
        selected = range(len(groups)) if not args.sample_groups else np.unique(np.linspace(0, len(groups)-1, min(args.sample_groups,len(groups)), dtype=int))
        for g in selected:
            start, stop = groups[g]
            for lo in range(0, values.shape[1], 24):
                jobs.append((p, int(start), int(stop), lo, min(lo+24, values.shape[1])))
    order_seconds = time.perf_counter()-stamp
    reconstructed = values.copy()
    path = args.output/'visibility.mx.h5.partial'
    counts, diagnostics, account = Counter(), Counter(), Counter()
    offsets, records = [0], []
    checked_samples, error_sum, largest_error, largest_row_mse, largest_bias = 0, 0., 0., 0., 0.
    timings = Counter()
    with h5py.File(path, 'w') as h5:
        h5.attrs.update(format='varquanms.experimental-mx-component/1', complete=False,
                        contains_full_ms_metadata=False, variant=args.variant, source=source.name,
                        full_observation=not args.sample_groups, shape=values.shape)
        h5.create_dataset('configuration', data=json.dumps(config.to_dict()).encode())
        payload = h5.create_dataset('payload', (0,), maxshape=(None,), dtype='u1', chunks=(1024**2,), fletcher32=True)
        for p, order in enumerate(orders):
            h5.create_dataset(f'permutations/{p}', data=order['permutation'], compression='gzip', compression_opts=1, fletcher32=True)
        for start in range(0, len(jobs), args.batch):
            batch_jobs = jobs[start:start+args.batch]
            stamp = time.perf_counter()
            channels = batch_jobs[0][4]-batch_jobs[0][3]
            if any(hi-lo != channels for _,_,_,lo,hi in batch_jobs):
                raise ValueError('Use --batch 1 for heterogeneous edge channel counts')
            batch = np.zeros((len(batch_jobs), 32, channels), np.complex64)
            masks = np.ones(batch.shape, bool)
            lengths = []
            for i,(p,a,b,lo,hi) in enumerate(batch_jobs):
                rows = orders[p]['permutation'][a:b]
                batch[i,:len(rows)] = values[rows,lo:hi,p]
                masks[i,:len(rows)] = flags[rows,lo:hi,p]
                lengths.append(len(rows))
            timings['gathering'] += time.perf_counter()-stamp
            stamp = time.perf_counter()
            if args.variant == 'uniform':
                encoded = [uniform_encode(v[:n], m[:n], config, bounds) for v,m,n in zip(batch,masks,lengths)]
            else:
                outer, granularity = args.variant.split('-')
                encoded, diag = encode_batch(batch,masks,lengths,bounds,outer_mode=outer,granularity=granularity)
                diagnostics.update(diag)
            timings['candidate_selection_packing_compression'] += time.perf_counter()-stamp
            stamp = time.perf_counter()
            blobs = []
            for i, (record, job) in enumerate(zip(encoded,batch_jobs)):
                p,a,b,lo,hi = job
                rows = orders[p]['permutation'][a:b]
                decoded = (uniform_decode(record.blob,(b-a,hi-lo)) if args.variant == 'uniform' else decode_block(record.blob))
                metrics = error_metrics(batch[i,:b-a],decoded,~masks[i,:b-a],percentile=False)
                if not within_bounds(metrics,bounds):
                    raise ValueError('Serialized block failed source error bounds')
                protected = masks[i,:b-a] | ~np.isfinite(batch[i,:b-a]) | ((batch[i,:b-a].real == 0)&np.signbit(batch[i,:b-a].real)) | ((batch[i,:b-a].imag == 0)&np.signbit(batch[i,:b-a].imag))
                if decoded[protected].tobytes() != batch[i,:b-a][protected].tobytes():
                    raise ValueError('Protected source bits changed')
                reconstructed[rows,lo:hi,p] = decoded
                checked_samples += metrics['valid_samples']
                error_sum += metrics['component_mse']*metrics['valid_samples']
                largest_error = max(largest_error,metrics['max_complex_error'])
                largest_row_mse = max(largest_row_mse,metrics['max_baseline_component_mse'])
                largest_bias = max(largest_bias,metrics['complex_bias_abs'])
                if args.variant == 'uniform':
                    counts[str(record.bits)] += (b-a)*(hi-lo)*2
                    account.update(record.accounting)
                else:
                    counts.update(record.formats); account.update(record.accounting)
                blobs.append(record.blob); offsets.append(offsets[-1]+len(record.blob)); records.append(job)
            timings['decode_and_source_checks'] += time.perf_counter()-stamp
            stamp = time.perf_counter()
            payload.resize((offsets[-1],))
            payload[offsets[-len(blobs)-1]:] = np.frombuffer(b''.join(blobs),np.uint8)
            timings['payload_writing'] += time.perf_counter()-stamp
            if start//args.batch % 20 == 0:
                print(json.dumps({'completed_blocks':min(start+args.batch,len(jobs)), 'total_blocks':len(jobs)}),flush=True)
        h5.create_dataset('offsets',data=np.asarray(offsets,dtype='<u8'),fletcher32=True)
        h5.create_dataset('records',data=np.asarray(records,dtype='<u4'),fletcher32=True)
        if not args.sample_groups:
            stamp=time.perf_counter()
            summary, arrays = quality(values,reconstructed,flags,names,workers=1)
            for name, metrics in summary['by_product'].items():
                if not within_bounds(metrics,bounds):
                    raise ValueError('Whole-product bounds failed: '+name)
            for name,array in arrays.items():
                h5.create_dataset('quality/'+name,data=array,compression='gzip',compression_opts=1,fletcher32=True)
            timings['whole_product_quality'] += time.perf_counter()-stamp
        else:
            summary = None
        h5.attrs['complete'] = True
        h5.flush()
    path.rename(args.output/'visibility.mx.h5')
    component_seconds = time.perf_counter()-started
    if not args.sample_groups:
        stamp=time.perf_counter()
        independent=decode_component(args.output/'visibility.mx.h5')
        if independent.tobytes() != reconstructed.tobytes():
            raise ValueError('Independent file decode differs from encode-time decode')
        del independent
        timings['separate_component_read_decode'] = time.perf_counter()-stamp
    record = {'schema':'vamp.mx-evaluation/1','status':'complete','variant':args.variant,
              'source':source.name,'source_shape':list(values.shape),'component_only':True,
              'sample_groups_per_product':args.sample_groups,'selection':'evenly spaced amplitude ranks including extremes' if args.sample_groups else 'all source samples',
              'ordering':args.ordering,'block_count':len(jobs),'complex_samples':sum((b-a)*(hi-lo) for p,a,b,lo,hi in jobs),
              'bounds':bounds,'configuration':config.to_dict(),'component_file_bytes':(args.output/'visibility.mx.h5').stat().st_size,
              'configuration_role':'Grouping/error bounds and uniform control; MX element formats are defined by variant, not configuration.bits/codebook',
              'serialized_block_bytes':offsets[-1],'block_offsets_bytes':8*len(offsets),'record_map_bytes':20*len(jobs),
              'accounting_before_joint_compression':dict(account),'real_component_format_counts':dict(counts),
              'fixed_format_passing_blocks':dict(diagnostics),
              'measured_errors':{'component_mse':error_sum/max(1,checked_samples),'max_complex_error':largest_error,
                                 'max_within_block_baseline_mse':largest_row_mse,'max_block_bias':largest_bias},
              'quality':summary,'timings_seconds':{'inspection_and_read':read_seconds,'ordering':order_seconds,**timings,'complete_component_experiment':component_seconds},
              'timing_scope':'one warm/unspecified-cache shared-host process, one worker; includes inspection, selection, packing, writing, independent per-block decode/checks; excludes separate file decode, export and imaging; diagnostic timing with possible documentation/synthetic-test work, not comparable to production encode-only timing',
              'environment':{'python':platform.python_version(),'numpy':np.__version__,'cpu_count':os.cpu_count()},
              'component_sha256':digest(args.output/'visibility.mx.h5')}
    if args.export_ms:
        stamp=time.perf_counter()
        restored=args.output/'restored.ms'
        with table(str(source),readonly=True,ack=False) as t:
            copied=t.copy(str(restored),deep=True,valuecopy=True); copied.close()
        with table(str(restored),readonly=False,ack=False) as t:
            t.putcol('DATA',reconstructed); t.flush()
        record['export_comparison']=compare_ms(source,restored,approximate_columns=('DATA',))
        with table(str(restored),readonly=True,ack=False) as t:
            if t.getcol('DATA').tobytes()!=reconstructed.tobytes():
                raise ValueError('Exported DATA differs from decoded block values')
        record['timings_seconds']['ms_export_and_exact_comparison']=time.perf_counter()-stamp
    (args.output/'report.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:record[k] for k in ['variant','component_file_bytes','serialized_block_bytes','measured_errors','timings_seconds']},indent=2))


if __name__=='__main__':
    main()
