"""Deterministic block screen; estimates payload rates, not whole-archive sizes."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import time

from casacore.tables import table
import numpy as np

from varquanms.codec import encode_block
from varquanms.config import Config
from varquanms.ordering import all_orders
from varquanms.partitioned import context_runs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ms', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    _, runs = context_runs(args.ms)
    start, stop = runs[0]
    with table(str(args.ms), ack=False, readonly=True, lockoptions='usernoread') as t:
        values = t.getcol('CORRECTED_DATA', startrow=start, nrow=stop-start)
        flags = t.getcol('FLAG', startrow=start, nrow=stop-start)
        flags |= t.getcol('FLAG_ROW', startrow=start, nrow=stop-start)[:, None, None]
        pairs = np.column_stack([t.getcol(c, startrow=start, nrow=stop-start) for c in ('ANTENNA1', 'ANTENNA2')])
        uvw = t.getcol('UVW', startrow=start, nrow=stop-start)
    with table(str(args.ms / 'POLARIZATION'), ack=False, readonly=True, lockoptions='usernoread') as t:
        types = t.getcell('CORR_TYPE', 0)
    labels = {9:'XX', 10:'XY', 11:'YX', 12:'YY'}
    names = [labels[int(v)] for v in types]
    values[flags] = 0
    base = Config.load('configs/calibrated-uniform.json').to_dict()
    methods = [
        ('baseline', dict(ordering='original', baseline_group=1, channel_segment=24)),
        ('original32', dict(ordering='class-original', baseline_group=32, channel_segment=24)),
        ('amplitude32', dict(ordering='amplitude', baseline_group=32, channel_segment=24)),
        ('buckets32', dict(ordering='buckets', baseline_group=32, channel_segment=24)),
        ('geometry32', dict(ordering='length-orientation', baseline_group=32, channel_segment=24)),
        ('amplitude128', dict(ordering='amplitude', baseline_group=128, channel_segment=48)),
        ('power32', dict(ordering='amplitude', baseline_group=32, channel_segment=24,
                         codebook='power-midrise-v1', power_center='median', power_gammas=[1., .5])),
    ]
    results = []
    for method, options in methods:
        cfg = Config.load(None, **{**base, **options})
        orders, _, _ = all_orders(values, flags, pairs, uvw, cfg)
        jobs = []
        for p, order in enumerate(orders):
            boundaries = order['group_boundaries']
            selected = np.unique(np.linspace(0, len(boundaries)-2, min(128, len(boundaries)-1)).astype(int))
            for g in selected:
                rows = order['permutation'][boundaries[g]:boundaries[g+1]]
                for lo in range(0, values.shape[1], cfg.channel_segment):
                    hi = min(lo+cfg.channel_segment, values.shape[1])
                    jobs.append((p, int(g), lo, values[rows,lo:hi,p], flags[rows,lo:hi,p]))
        for mse in (.01, 1., 64., 256.):
            name = f'{method}-mse{mse:g}'
            config = Config.load(None, **{**base, **options, 'max_component_mse':mse,
                                         'max_baseline_component_mse':2*mse,
                                         'max_complex_error':float(10*np.sqrt(mse)), 'max_complex_bias':float(np.sqrt(mse))})
            (args.output / (name+'.json')).write_text(json.dumps(config.to_dict(), indent=2)+'\n')
            stamp = time.perf_counter()
            def run(job):
                p,g,lo,a,f = job
                block = encode_block(a,f,config,config.bounds(names[p]),stream=(0,p,g,lo//config.channel_segment))
                return len(block.blob), block.decoded.size, block.bits, block.metrics, block.accounting
            with ThreadPoolExecutor(max_workers=8) as pool:
                encoded = list(pool.map(run,jobs))
            samples = sum(e[1] for e in encoded)
            nvalid = sum(e[3]['valid_samples'] for e in encoded)
            histogram = Counter(str(e[2]) if e[2] else 'float32' for e in encoded)
            cost = Counter()
            for e in encoded: cost.update(e[4])
            result = dict(name=name, config=config.to_dict(), sampled_blocks=len(encoded), sampled_samples=samples,
                          payload_bytes=sum(e[0] for e in encoded), bytes_per_sample=sum(e[0] for e in encoded)/samples,
                          component_mse=sum(e[3]['component_mse']*e[3]['valid_samples'] for e in encoded)/max(1,nvalid),
                          max_complex_error=max(e[3]['max_complex_error'] for e in encoded), block_bits=dict(histogram),
                          accounting=dict(cost), elapsed_seconds=time.perf_counter()-stamp)
            results.append(result)
            print(json.dumps({k:v for k,v in result.items() if k not in ('config','accounting')}), flush=True)
    report = dict(schema='varquanms.calibrated-screen/1', source=str(args.ms), source_rows=[start,stop],
                  sampling='Up to 128 evenly spaced quantization groups per product; all frequency segments. Deterministic, not random or a whole-MS estimate.',
                  limitations='Different groupings sample different values. Rates exclude metadata, maps, quality arrays and HDF5 overhead. Use full encodes and images for selection.',
                  experiments=results)
    (args.output/'screen.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
