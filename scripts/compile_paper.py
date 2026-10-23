"""Compile the paper with Pandoc and Typst; equations and vector figures stay typeset."""

import argparse
import datetime
import json
from pathlib import Path
import subprocess
import tempfile

import pypandoc
import typst

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "docs/paper"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PAPER / "varquantMS.pdf")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="varquanms-paper-") as temporary:
        intermediate = Path(temporary) / "paper.typ"
        subprocess.run([
            pypandoc.get_pandoc_path(), str(PAPER / "varquantMS.md"),
            "--from=markdown-implicit_figures", "--to=typst", "--standalone",
            "--syntax-highlighting=none", "--fail-if-warnings", "--wrap=none",
            "--template=" + str(PAPER / "pdf-template.typ"),
            "--lua-filter=" + str(PAPER / "pdf-filter.lua"),
            "--output=" + str(intermediate),
        ], check=True, cwd=ROOT)
        result = typst.compile(intermediate.read_bytes(), root=str(ROOT),
                               timestamp=datetime.datetime(2026, 10, 23, tzinfo=datetime.timezone.utc))
        if not result.startswith(b"%PDF-"):
            raise ValueError("Compiler did not produce a PDF")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        # Publish only after the complete document compiles successfully.
        with tempfile.NamedTemporaryFile(dir=args.output.parent, suffix=".partial", delete=False) as stream:
            staged = Path(stream.name)
            stream.write(result)
        staged.replace(args.output)
    print(json.dumps({"pdf": str(args.output), "bytes": len(result),
                      "pandoc": str(pypandoc.get_pandoc_version()), "typst": typst.__version__}, indent=2))


if __name__ == "__main__":
    main()
