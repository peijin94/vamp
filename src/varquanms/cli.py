"""Command-line entry point for the first (inspection) milestone."""

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
    parser = argparse.ArgumentParser(prog="visarchive", description="Measurement Set archive research tools (inspection milestone)")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    inspect = sub.add_parser("inspect", help="scan a frozen MS read-only; write a JSON inventory")
    inspect.add_argument("input", type=Path)
    inspect.add_argument("--report", required=True, type=Path)
    inspect.add_argument("--summary", type=Path, help="also save the human-readable summary")
    inspect.add_argument("--columns", nargs="+", default=["DATA"], help="visibility columns proposed for compression (default DATA)")
    inspect.add_argument("--chunk-rows", type=int, default=1024)
    args = parser.parse_args(argv)
    try:
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
