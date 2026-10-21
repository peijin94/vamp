"""Verify a local sample tar and every extracted regular file against its source.

Ownership and filesystem allocation are deliberately excluded from byte identity.
This is copy verification, not archive round-trip or CASA compatibility validation.
"""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile


def digest(stream):
    result = hashlib.sha256()
    for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
        result.update(block)
    return result.hexdigest()


def verify(source_tar, local_tar, extracted_parent):
    with source_tar.open("rb") as f:
        original_hash = digest(f)
    with local_tar.open("rb") as f:
        copy_hash = digest(f)
    if original_hash != copy_hash:
        raise ValueError("Source and local tar SHA-256 differ")
    inventory = []
    seen = set()
    with tarfile.open(local_tar, "r:") as archive:
        for member in archive:
            relative = PurePosixPath(member.name)
            if relative.is_absolute() or ".." in relative.parts or not relative.parts:
                raise ValueError(f"Unsafe tar member: {member.name}")
            if member.name in seen:
                raise ValueError(f"Duplicate tar member: {member.name}")
            seen.add(member.name)
            path = extracted_parent.joinpath(*relative.parts)
            if path.is_symlink() or not path.resolve().is_relative_to(extracted_parent.resolve()):
                raise ValueError(f"Unsafe extracted path: {path}")
            if member.isdir():
                if not path.is_dir():
                    raise ValueError(f"Missing extracted directory: {path}")
                continue
            if not member.isfile():
                raise ValueError(f"Unsupported tar member type: {member.name}")
            if path.stat().st_size != member.size:
                raise ValueError(f"File length mismatch: {path}")
            with archive.extractfile(member) as stream:
                expected = digest(stream)
            with path.open("rb") as stream:
                actual = digest(stream)
            if expected != actual:
                raise ValueError(f"Extracted file SHA-256 mismatch: {path}")
            inventory.append({"path": member.name, "bytes": member.size, "sha256": actual})
    roots = {PurePosixPath(item["path"]).parts[0] for item in inventory}
    actual_files = {str(p.relative_to(extracted_parent))
                    for root in roots for p in (extracted_parent / root).rglob("*") if p.is_file()}
    if actual_files != {item["path"] for item in inventory}:
        raise ValueError("Extracted directory has missing or unexpected files")
    return {"schema": "varquanms.sample-copy/1", "status": "verified",
            "source_tar": str(source_tar.resolve()), "local_tar": str(local_tar),
            "extracted_parent": str(extracted_parent), "tar_sha256": copy_hash,
            "tar_bytes": local_tar.stat().st_size, "regular_file_count": len(inventory),
            "file_bytes": sum(item["bytes"] for item in inventory), "files": inventory,
            "note": "Source/local tar SHA-256 and every extracted regular-file SHA-256 match. Source has only been read. No archive codec or scientific-fidelity check is implied."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_tar", type=Path)
    parser.add_argument("local_tar", type=Path)
    parser.add_argument("extracted_parent", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    try:
        from varquanms.cli import _check_output, _write_new
        if args.report.resolve().is_relative_to(args.extracted_parent.resolve()):
            raise ValueError("Report must be outside the extracted data directory")
        report_path = _check_output(args.report, args.source_tar.resolve())
        report = verify(args.source_tar, args.local_tar, args.extracted_parent)
        _write_new(report_path, json.dumps(report, indent=2) + "\n")
        print(f"Verified tar SHA-256 {report['tar_sha256']} and {report['regular_file_count']} extracted files ({report['file_bytes']:,} bytes)")
    except (ValueError, OSError, tarfile.TarError) as error:
        print(f"copy verification: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
