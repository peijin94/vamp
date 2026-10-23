"""Strict command-line entry points for inspection, archival and reconstruction."""

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

from . import __version__


def _check_output(path, source):
    path = Path(path).resolve()
    if path == source or path.is_relative_to(source):
        raise ValueError(f"Inspection output cannot be inside the source MS: {path}")
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    if not path.parent.is_dir():
        raise ValueError(f"Output parent does not exist: {path.parent}")
    return path


def _write_new(path, contents):
    # link() publishes only a fully written file and refuses overwrite atomically.
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     prefix=path.name + ".", suffix=".partial", delete=False) as f:
        temporary = Path(f.name)
        try:
            f.write(contents)
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="visarchive", description="Measurement Set archive research tools")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    inspect = sub.add_parser("inspect", help="scan a frozen MS read-only; write a JSON inventory")
    inspect.add_argument("input", type=Path)
    inspect.add_argument("--report", required=True, type=Path)
    inspect.add_argument("--summary", type=Path, help="also save the human-readable summary")
    inspect.add_argument("--columns", nargs="+", default=["DATA"], help="visibility columns proposed for compression (default DATA)")
    inspect.add_argument("--chunk-rows", type=int, default=1024)
    encode_parser = sub.add_parser("encode", help="create an independently decodable HDF5 archive")
    encode_parser.add_argument("input", type=Path)
    encode_parser.add_argument("output", type=Path)
    encode_parser.add_argument("--config", type=Path)
    encode_parser.add_argument("--workers", type=int, help="bounded CPU worker count; default 1 or config value")
    encode_parser.add_argument("--burn-flag", "--burn_flag", action="store_true", default=None,
                               help="destructive: discard MAIN flags; export defined flag cells as all-false (no reflagging)")
    encode_parser.add_argument("--keep-calib-only", "--keep_calib_only", action="store_true", default=None,
                               help="keep only CORRECTED_DATA; omit DATA/FLOAT_DATA/MODEL_DATA and zero flagged values; flags exact, weights default to 4 bits")
    encode_parser.add_argument("--weight-bits", type=int, choices=(4, 32), help="calibrated-only weights: packed 4-bit codebooks (default) or exact float32")
    encode_parser.add_argument("--weight-max-relative-error", type=float, help="calibrated-only weights: maximum relative error (default 0.01); failing chunks use float32")
    encode_parser.add_argument("--report", type=Path)
    decode_parser = sub.add_parser("decode", help="restore a new MS without overwriting an existing path")
    decode_parser.add_argument("input", type=Path)
    decode_parser.add_argument("output", type=Path)
    decode_parser.add_argument("--report", type=Path)
    data_export = decode_parser.add_mutually_exclusive_group()
    data_export.add_argument("--ms-v2-data", action="store_true", default=None, help="default for calibrated-only archives: add DATA as a copy of decoded CORRECTED_DATA; does not recover raw DATA")
    data_export.add_argument("--no-ms-v2-data", dest="ms_v2_data", action="store_false", default=None, help="export a CORRECTED_DATA-only projection (not compatible with the CASA imager)")
    verify_parser = sub.add_parser("verify", help="verify integrity and compare original/exported values")
    verify_parser.add_argument("input", type=Path)
    verify_parser.add_argument("--source", type=Path)
    verify_parser.add_argument("--restored", type=Path)
    verify_parser.add_argument("--report", type=Path)
    benchmark_parser = sub.add_parser("benchmark", help="run a strict experiment matrix in fresh processes")
    benchmark_parser.add_argument("input", type=Path)
    benchmark_parser.add_argument("--matrix", type=Path, required=True)
    benchmark_parser.add_argument("--output", type=Path, required=True)
    benchmark_parser.add_argument("--experiments", nargs="+", help="run only these named matrix entries")
    args = parser.parse_args(argv)
    try:
        if args.command == "benchmark":
            from .benchmark import run
            result = run(args.input, args.matrix, args.output, args.experiments)
            print(json.dumps({"status": result["status"], "experiments": len(result["experiments"]), "results": str(args.output / "results.json")}, indent=2))
            return 0
        if args.command != "inspect":
            from .archive import decode, encode, verify
            from .config import Config
            protected_source = args.source if args.command == "verify" and args.source else args.input
            report_path = _check_output(args.report, protected_source.resolve()) if args.report else None
            if args.command == "encode":
                overrides = {} if args.burn_flag is None else {"burn_flag": True}
                if args.workers is not None:
                    overrides["workers"] = args.workers
                if args.keep_calib_only:
                    overrides.update(keep_calib_only=True, columns=["CORRECTED_DATA"])
                if args.weight_bits is not None:
                    overrides["weight_bits"] = args.weight_bits
                if args.weight_max_relative_error is not None:
                    overrides["weight_max_relative_error"] = args.weight_max_relative_error
                configuration = Config.load(args.config, **overrides)
                if not configuration.keep_calib_only and (args.weight_bits is not None or args.weight_max_relative_error is not None):
                    raise ValueError("Weight precision options require --keep-calib-only")
                if configuration.burn_flag:
                    print("WARNING: burn_flag discards MAIN FLAG/FLAG_CATEGORY/FLAG_ROW values. "
                          "Export resets them to false; no flagging algorithm is run.", file=sys.stderr)
                result = encode(args.input, args.output, configuration)
            elif args.command == "decode":
                result = decode(args.input, args.output, ms_v2_data=args.ms_v2_data)
            else:
                result = verify(args.input, args.source, args.restored)
            serialized = json.dumps(result, indent=2, allow_nan=False) + "\n"
            if report_path:
                _write_new(report_path, serialized)
            print(serialized, end="")
            return 0
        if args.chunk_rows < 1:
            raise ValueError("--chunk-rows must be positive")
        source = args.input.resolve()
        report_path = _check_output(args.report, source)
        summary_path = _check_output(args.summary, source) if args.summary else None
        if summary_path == report_path:
            raise ValueError("Report and summary paths must be different")
        from .inspection import human_summary, inspect_ms
        report = inspect_ms(source, args.columns, args.chunk_rows)
        summary = human_summary(report)
        _write_new(report_path, json.dumps(report, indent=2, allow_nan=False) + "\n")
        if summary_path:
            _write_new(summary_path, summary)
        print(summary, end="")
        return 0
    except (ValueError, OSError, RuntimeError, TypeError, ImportError) as error:
        print(f"visarchive: {error}", file=sys.stderr)
        return 1
