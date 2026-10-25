"""Experimental ordered partition bundle CLI; see devnotes for support limits."""
import argparse
import json
from pathlib import Path

from varquanms.cli import _write_new
from varquanms.config import Config
from varquanms.partitioned import decode_partitioned, encode_partitioned


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("encode", "decode", "compact"))
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--keep-calib-only", action="store_true")
    parser.add_argument("--partition-workers", type=int, default=1)
    data_export = parser.add_mutually_exclusive_group()
    data_export.add_argument("--ms-v2-data", action="store_true", default=None)
    data_export.add_argument("--no-ms-v2-data", dest="ms_v2_data", action="store_false", default=None)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists() or not args.report.parent.is_dir():
        raise ValueError("Report must be a new file in an existing directory")
    if args.operation != "encode" and (args.config or args.keep_calib_only or args.partition_workers != 1):
        raise ValueError("Decode and compaction use only the configuration stored in the bundle")
    if args.operation != "decode" and args.ms_v2_data is not None:
        raise ValueError('--ms-v2-data is an export option')
    overrides = {"keep_calib_only": True, "columns": ["CORRECTED_DATA"]} if args.keep_calib_only else {}
    if args.operation == 'compact':
        from varquanms.shared import compact_bundle
        result = compact_bundle(args.source, args.output)
    elif args.operation == 'encode':
        result = encode_partitioned(args.source, args.output, Config.load(args.config, **overrides),
                                    partition_workers=args.partition_workers)
    else:
        result = decode_partitioned(args.source, args.output, ms_v2_data=args.ms_v2_data)
    _write_new(args.report, json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
