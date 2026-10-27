"""Encode a regular calibrated MS once with casacore's Dysco storage manager.

Run inside a casacore environment that has DyscoStMan. The input must already
use the calibrated-only policy (flagged visibility values zeroed). The output
keeps the original row metadata, flags, weights and subtables, and stores the
calibrated visibility column as Dysco-compressed DATA. It never overwrites an
input or an existing output. This is an experiment, not an archive format.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import numpy as np
from casacore.tables import makearrcoldesc, table


def regular_layout(source):
    """Fail before encoding if Dysco's documented regular-time-block rule fails."""
    times = np.asarray(source.getcol('TIME'))
    if not len(times) or not np.all(np.isfinite(times)):
        raise ValueError('TIME must contain finite records')
    edges = np.r_[0, np.flatnonzero(times[1:] != times[:-1]) + 1, len(times)]
    sizes = np.diff(edges)
    if len(sizes) < 2 or not np.all(sizes == sizes[0]) or np.any(np.diff(times[edges[:-1]]) <= 0):
        raise ValueError('Dysco requires at least two ordered, equally sized time blocks')
    identities = ('ANTENNA1', 'ANTENNA2', 'DATA_DESC_ID', 'FIELD_ID')
    count = int(sizes[0])
    first = {name: source.getcol(name, 0, count) for name in identities}
    for start in edges[1:-1]:
        for name, expected in first.items():
            if not np.array_equal(expected, source.getcol(name, int(start), count)):
                raise ValueError(f'Dysco time-block ordering differs in {name}')
    return {'time_blocks': len(sizes), 'rows_per_time_block': count,
            'identity_columns': list(identities)}


def inventory(directory, sequence):
    files, data_files = {}, {}
    prefix = f'table.f{sequence}'
    for path in sorted(directory.rglob('*')):
        if path.is_symlink():
            raise ValueError('Unexpected symlink in completed MS')
        if path.is_file():
            stat = path.stat()
            name = str(path.relative_to(directory))
            files[name] = {'bytes': stat.st_size, 'allocated_bytes': stat.st_blocks * 512}
            if path.parent == directory and (path.name == prefix or path.name.startswith(prefix + '_')):
                data_files[name] = files[name]
    if not data_files:
        raise ValueError('Dysco storage files were not found')
    return {'file_count': len(files), 'total_bytes': sum(x['bytes'] for x in files.values()),
            'allocated_bytes': sum(x['allocated_bytes'] for x in files.values()),
            'dysco_data_manager_files': data_files,
            'dysco_data_manager_bytes': sum(x['bytes'] for x in data_files.values()),
            'note': 'Dysco DATA-manager files include code payload, scales and manager headers; generic table metadata is in the complete-MS total.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--bits', type=int, choices=(2, 3, 4, 6, 8, 10, 12, 16), required=True)
    p.add_argument('--normalization', choices=('AF', 'RF'), required=True)
    p.add_argument('--truncation', type=float, choices=(1.5, 2.5, 3.5), default=2.5)
    p.add_argument('--flagged-input', choices=('zero', 'nan'), default='zero',
                   help='Use NaN for Dysco normalisation or zero as an explicit control')
    p.add_argument('--chunk-rows', type=int, default=2048)
    args = p.parse_args()
    source_path, output = args.input.resolve(), args.output.resolve()
    partial = output.with_name(output.name + '.partial')
    if not source_path.is_dir() or not output.parent.is_dir() or args.chunk_rows < 1:
        raise ValueError('Invalid input, output parent or row chunk')
    if output.is_relative_to(source_path) or any(x.exists() for x in (output, partial, args.report)):
        raise ValueError('Use a new output/report outside the input MS')
    started = time.perf_counter()
    report = {'schema': 'vamp.dysco-encode/1', 'status': 'running',
              'input': str(source_path), 'output': str(output), 'partial_output': str(partial),
              'bits_per_real_component': args.bits, 'normalization': args.normalization,
              'distribution': 'TruncatedGaussian', 'distribution_truncation': args.truncation,
              'input_visibility_column': 'CORRECTED_DATA', 'output_visibility_column': 'DATA',
              'flagged_input_values': 'complex NaN' if args.flagged_input == 'nan' else 'positive complex zero',
              'flag_columns': 'preserved',
              'weights': 'WEIGHT stored exactly; no WEIGHT_SPECTRUM in this input',
              'source_raw_DATA_and_MODEL_DATA': 'omitted',
              'spec_weightBitCount': 12,
              'note': 'The 12-bit weight setting is not applied to a column: only DATA is bound to DyscoStMan.'}
    args.report.write_text(json.dumps(report, indent=2) + '\n')
    try:
        with table(str(source_path), readonly=True, ack=False, lockoptions='usernoread') as source:
            if 'CORRECTED_DATA' not in source.colnames() or 'WEIGHT_SPECTRUM' in source.colnames():
                raise ValueError('This benchmark expects CORRECTED_DATA and row WEIGHT only')
            report['layout'] = regular_layout(source)
            cell = np.asarray(source.getcell('CORRECTED_DATA', 0))
            if cell.dtype != np.complex64 or cell.ndim != 2:
                raise ValueError('Visibility cells must be two-dimensional complex64')
            report['cell_shape'] = list(cell.shape)
            retained = [c for c in source.colnames() if c not in ('DATA', 'CORRECTED_DATA', 'MODEL_DATA')]
            selected = source.query(columns=','.join(retained))
            stage = time.perf_counter()
            copied = selected.copy(str(partial), deep=True, valuecopy=True)
            copied.close()
            selected.close()
            report['metadata_copy_seconds'] = time.perf_counter() - stage
            out = table(str(partial), readonly=False, ack=False)
            # The input's private markers describe a different physical column
            # arrangement and a bitwise DATA copy; they would be false here.
            removed = []
            for key in ('VARQUANMS_CALIBRATED_ONLY', 'VARQUANMS_DERIVED_DATA', 'VARQUANMS_WEIGHT_POLICY'):
                if key in out.getkeywords():
                    out.removekeyword(key)
                    removed.append(key)
            out.putkeyword('DYSCO_COMPARISON', {'source_column': 'CORRECTED_DATA',
                           'stored_column': 'DATA', 'lossy': True,
                           'flagged_input_values': report['flagged_input_values'],
                           'original_raw_DATA_and_MODEL_DATA': 'omitted'})
            report['replaced_private_keywords'] = removed
            descriptor = makearrcoldesc('DATA', 0j, shape=list(cell.shape), options=1,
                                        valuetype='complex',
                                        comment=source.getcoldesc('CORRECTED_DATA').get('comment', ''),
                                        keywords=source.getcolkeywords('CORRECTED_DATA'))
            spec = {'dataBitCount': args.bits, 'weightBitCount': 12,
                    'distribution': 'TruncatedGaussian', 'normalization': args.normalization,
                    'distributionTruncation': args.truncation}
            out.addcols(descriptor, {'TYPE': 'DyscoStMan', 'NAME': 'DyscoData', 'SPEC': spec})
            manager = out.getdminfo('DATA')
            if manager['TYPE'] != 'DyscoStMan':
                raise ValueError('DATA is not bound to DyscoStMan')
            report['manager'] = manager
            stage = time.perf_counter()
            for start in range(0, source.nrows(), args.chunk_rows):
                count = min(args.chunk_rows, source.nrows()-start)
                values = np.asarray(source.getcol('CORRECTED_DATA', start, count))
                flags = source.getcol('FLAG', start, count) | source.getcol('FLAG_ROW', start, count)[:, None, None]
                if values.dtype != np.complex64 or values.shape != (count, *cell.shape):
                    raise ValueError('Heterogeneous visibility shape or dtype')
                if not np.isfinite(values[~flags]).all():
                    raise ValueError('Unflagged nonfinite values are not supported by this benchmark')
                if np.any(values[flags].view('u4') != 0):
                    raise ValueError('Input flagged values are not positive complex zero')
                if args.flagged_input == 'nan':
                    # Dysco's published convention ignores flagged samples by
                    # encoding them as NaN, rather than including zeros in its
                    # scale statistics. Flags themselves remain unchanged.
                    values = values.copy()
                    values[flags] = complex(float('nan'), float('nan'))
                out.putcol('DATA', values, start, count)
            out.flush()
            report['putcol_and_flush_seconds'] = time.perf_counter() - stage
            stage = time.perf_counter()
            out.close()  # Dysco can finish substantial asynchronous work here.
            report['close_seconds'] = time.perf_counter() - stage
            report['rows'] = source.nrows()
        with table(str(partial), readonly=True, ack=False, lockoptions='usernoread') as check:
            if check.nrows() != report['rows'] or check.getdminfo('DATA')['TYPE'] != 'DyscoStMan':
                raise ValueError('Completed output did not reopen with Dysco DATA')
            sequence = int(check.getdminfo('DATA')['SEQNR'])
        report['storage'] = inventory(partial, sequence)
        report['work_seconds'] = time.perf_counter() - started
        usage = resource.getrusage(resource.RUSAGE_SELF)
        report.update(peak_rss_bytes=usage.ru_maxrss * 1024, cpu_user_seconds=usage.ru_utime,
                      cpu_system_seconds=usage.ru_stime)
        os.rename(partial, output)
        report['status'] = 'complete'
    except BaseException as error:
        report.update(status='failed', error=f'{type(error).__name__}: {error}',
                      work_seconds=time.perf_counter()-started)
        args.report.write_text(json.dumps(report, indent=2, default=lambda x: x.tolist()) + '\n')
        raise
    args.report.write_text(json.dumps(report, indent=2, default=lambda x: x.tolist(), allow_nan=False) + '\n')
    print(json.dumps({'status': report['status'], 'bits': args.bits, 'normalization': args.normalization,
                      'bytes': report['storage']['total_bytes'], 'work_seconds': report['work_seconds']}), flush=True)


if __name__ == '__main__':
    main()
