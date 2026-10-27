"""Check one-pass Dysco output against the matched calibrated input."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from vamp.archive import discover_tables
from vamp.inspection import _open
from vamp.metadata import casacore_keywords, equal, logical_descriptor, relocate
from run_wsclean import snapshot


def flagged_value_counts(values, policy):
    """Count marker mismatches and, for zero markers, noncanonical zero bits."""
    if policy == 'complex NaN':
        return int(np.count_nonzero(~(np.isnan(values.real) & np.isnan(values.imag)))), 0
    if policy == 'positive complex zero':
        return int(np.count_nonzero(values != 0)), int(np.count_nonzero(values.view('u4')))
    raise ValueError('Unknown flagged visibility policy')


def verify(source, decoded):
    source, decoded = Path(source).resolve(), Path(decoded).resolve()
    started = time.perf_counter()
    a_graph = {s['path'] for s in discover_tables(source).values()}
    b_graph = {s['path'] for s in discover_tables(decoded).values()}
    if a_graph != b_graph:
        raise ValueError('Subtable graphs differ')
    exact = []
    for relative in sorted(a_graph):
        with _open(source / relative) as a, _open(decoded / relative) as b:
            if a.nrows() != b.nrows() or not equal(a.info(), b.info()):
                raise ValueError(f'Table row count or information differs: {relative}')
            main = relative == '.'
            names = [n for n in a.colnames() if not main or n not in ('DATA','CORRECTED_DATA','MODEL_DATA')]
            expected = set(names) | ({'DATA'} if main else set())
            if expected != set(b.colnames()):
                raise ValueError(f'Unexpected columns in {relative}')
            ak, bk = dict(a.getkeywords()), dict(b.getkeywords())
            if main:
                for name in ('VARQUANMS_CALIBRATED_ONLY','VARQUANMS_DERIVED_DATA','VARQUANMS_WEIGHT_POLICY'):
                    ak.pop(name, None)
                marker = bk.pop('DYSCO_COMPARISON', None)
                if not isinstance(marker, dict) or marker.get('lossy') is not True or marker.get('stored_column') != 'DATA':
                    raise ValueError('Dysco provenance marker missing or invalid')
                flagged_policy = marker.get('flagged_input_values')
                if flagged_policy not in ('positive complex zero', 'complex NaN'):
                    raise ValueError('Unknown flagged visibility policy')
                if b.getdminfo('DATA')['TYPE'] != 'DyscoStMan':
                    raise ValueError('DATA is not Dysco-compressed')
                desc_a, desc_b = a.getcoldesc('CORRECTED_DATA'), b.getcoldesc('DATA')
                for key in ('valueType','ndim','keywords'):
                    if not equal(desc_a[key], desc_b[key]):
                        raise ValueError(f'DATA descriptor differs in {key}')
            if not equal(relocate(casacore_keywords(ak),source),relocate(casacore_keywords(bk),decoded)):
                raise ValueError(f'Table keywords differ: {relative}')
            for column in names:
                da = logical_descriptor({column: a.getcoldesc(column)})
                db = logical_descriptor({column: b.getcoldesc(column)})
                if not equal(relocate(da,source),relocate(db,decoded)):
                    raise ValueError(f'Column descriptor differs: {relative}/{column}')
                for start in range(0, a.nrows(), 2048):
                    count = min(2048,a.nrows()-start)
                    try:
                        x = np.asarray(a.getcol(column,start,count))
                        y = np.asarray(b.getcol(column,start,count))
                    except RuntimeError:
                        for row in range(start,start+count):
                            defined_a, defined_b = a.iscelldefined(column,row), b.iscelldefined(column,row)
                            if defined_a != defined_b or (defined_a and not equal(a.getcell(column,row),b.getcell(column,row))):
                                raise ValueError(f'Exact value differs: {relative}/{column}, row {row}')
                    else:
                        if not equal(x,y):
                            raise ValueError(f'Exact values differ: {relative}/{column}, rows {start}:{start+count}')
                exact.append(f'{relative}/{column}')
    flags_total = flagged_changed = flagged_noncanonical_zero_components = unflagged_nonfinite = 0
    with _open(source) as a, _open(decoded) as b:
        for start in range(0,a.nrows(),2048):
            count = min(2048,a.nrows()-start)
            flags = a.getcol('FLAG',start,count) | a.getcol('FLAG_ROW',start,count)[:,None,None]
            values = b.getcol('DATA',start,count)
            flags_total += int(flags.sum())
            discarded = values[flags]
            bad, noncanonical = flagged_value_counts(discarded, flagged_policy)
            flagged_changed += bad
            flagged_noncanonical_zero_components += noncanonical
            unflagged_nonfinite += int(np.count_nonzero(~np.isfinite(values[~flags])))
    if flagged_changed or unflagged_nonfinite:
        raise ValueError('Dysco output did not preserve the declared discarded-value marker or produced unflagged nonfinite values')
    return {'schema':'vamp.dysco-verification/1','status':'passed',
            'source':str(source),'decoded':str(decoded),'table_count':len(a_graph),
            'decoded_sha256':snapshot(decoded),
            'exact_columns':exact,'exact_column_count':len(exact),
            'flagged_samples':flags_total,'flagged_input_policy':flagged_policy,
            'flagged_values_policy_mismatch':flagged_changed,
            'flagged_values_not_zero':flagged_changed if flagged_policy == 'positive complex zero' else None,
            'flagged_values_not_nan':flagged_changed if flagged_policy == 'complex NaN' else None,
            'flagged_noncanonical_zero_components':flagged_noncanonical_zero_components if flagged_policy == 'positive complex zero' else None,
            'unflagged_nonfinite_values':unflagged_nonfinite,
            'metadata_scope':'Table graph, row counts, table keywords, table info, exact column descriptors and values; storage-manager assignments excluded.',
            'intentional_differences':['Only calibrated visibility DATA is stored; original DATA/MODEL_DATA and duplicate CORRECTED_DATA omitted.',
                'DATA uses a direct fixed-shape Dysco descriptor; the actual selected input cells are complex64 48×4.',
                'Three vamp-specific policy keywords are replaced by an explicit DYSCO_COMPARISON provenance keyword.'],
            'elapsed_seconds':time.perf_counter()-started}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path)
    p.add_argument('decoded',type=Path)
    p.add_argument('--report',type=Path,required=True)
    args=p.parse_args()
    if args.report.exists():
        raise ValueError('Report must be new')
    result=verify(args.source,args.decoded)
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','table_count','exact_column_count','flagged_samples','elapsed_seconds')}),flush=True)


if __name__=='__main__':
    main()
