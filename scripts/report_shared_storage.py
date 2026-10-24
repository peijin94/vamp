"""Record exact shared-storage measurements separately from quantization."""
import json
import os
from pathlib import Path
import platform

import casacore
import h5py
import numpy as np

from varquanms.container import Reader
from varquanms.imaging import file_sha256
from varquanms.partitioned import partition_reader
from varquanms.shared import load_shared

ROOT = Path('results/calibrated-sharing')


def collect():
    baseline = json.loads(Path('reports/calibrated-template-tests.json').read_text())
    experiments = {e['name']:e for e in baseline['experiments']}
    results = []
    for label, filename, source in [
        ('lossless', 'lossless-compaction.json', Path('results/calibrated-20261023/lossless-calib.bundle')),
        ('amplitude256', 'compaction.json', Path('results/calibrated-template/amplitude256.bundle')),
        ('sixbit2048', 'sixbit2048-compaction.json', Path('results/calibrated-template/sixbit2048.bundle')),
    ]:
        record = json.loads((ROOT/filename).read_text())
        manifest = load_shared(record['archive'])
        original = json.loads((source/'manifest.json').read_text())
        arrays = 0
        for a,b in zip(original['partitions'],manifest['partitions']):
            with Reader(source/a['archive'],verify=False) as old, partition_reader(record['archive'],b,verify=False) as new:
                if old.h5['manifest'][...].tobytes() != new.h5['manifest'][...].tobytes():
                    raise ValueError('Scientific manifest changed during compaction')
                if old.entries != new.entries:
                    raise ValueError('Logical array inventory changed')
                arrays += len(new.entries)
        old = experiments[label]
        record.update(name=label, all_original_partition_manifests_bitwise_equal=True,
                      verified_array_references=arrays, file_sha256=file_sha256(record['archive']),
                      image_comparison=old['comparison']['products']['image-image.fits'],
                      image_comparison_note='Previous imaging of unchanged decoder inputs; no new CLEAN run',
                      saved_bytes=old['archive']['archive_bytes']-record['archive_bytes'])
        results.append(record)
    a = experiments['amplitude256']['archive']
    phases = dict(partition_copy=a['sum_partition_copy_seconds'], inspection=a['sum_inspection_seconds'],
                  encoding=a['sum_encode_seconds_excluding_inspection'])
    phases['verification_integrity_and_other'] = a['elapsed_seconds_including_copy_inspection_verification']-sum(phases.values())
    export = json.loads((ROOT/'export-verification.json').read_text())
    identical = json.loads((ROOT/'identical-export.json').read_text())
    assert export['status']=='passed' and identical['status']=='passed'
    assert identical['all_supported_values_and_metadata_bitwise_equal']
    assert identical['original_source_import_checksums_match']
    return dict(schema='varquanms.shared-storage-experiment/1', experiments=results,
                recommended_profile='amplitude256', original_encoding_phase_seconds=phases,
                amplitude256_export_verification=export, amplitude256_matches_previously_imaged_export=identical,
                tests=dict(local_suite_passed=95), hardware=baseline['hardware'],
                software=dict(python=platform.python_version(), numpy=np.__version__, h5py=h5py.__version__,
                              hdf5=h5py.version.hdf5_version, python_casacore=casacore.__version__),
                limitations=['Compaction is an additional operation on existing archives, not an MS-to-archive speedup.',
                             'Shared host, caches not evicted; some validation and compaction jobs overlapped.',
                             'All three shared archives passed integrity and manifest equality checks; full MS export equality was additionally checked for amplitude256.',
                             'No new scientific fidelity claim: unchanged Stokes I image inputs, no additional quantization; Q/U/V and held-out data still untested.'])


def plot(result):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/varquanms-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    entries = result['experiments']
    x = np.arange(len(entries))
    fig, ax = plt.subplots(figsize=(8,4.8), constrained_layout=True)
    for offset, values, label, color in [
        (-.18, [e['storage_compaction']['source_bundle_bytes']/1e6 for e in entries], 'Directory bundle', '#888888'),
        (.18, [e['archive_bytes']/1e6 for e in entries], 'Shared HDF5 file', '#1f77b4'),
    ]:
        bars = ax.bar(x+offset, values, width=.36, label=label, color=color)
        ax.bar_label(bars, labels=[f'{v:.1f}' for v in values], padding=3)
    ax.set(xticks=x, xticklabels=['Lossless visibility', 'Amplitude 128×48', 'Six-bit experiment'],
           ylabel='Whole archive (decimal MB)', ylim=(0,1200),
           title='Exact storage sharing: unchanged visibilities and imaging inputs\nSame calibrated-only preservation policy for every bar')
    ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True); ax.legend()
    fig.savefig('devnotes/figs/calibrated-shared-size.png',dpi=180); plt.close(fig)


if __name__ == '__main__':
    result = collect()
    Path('reports/calibrated-shared-storage.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    plot(result)
    from report_calibrated_tests import write_note
    write_note(json.loads(Path('reports/calibrated-template-tests.json').read_text()))
