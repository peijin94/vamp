#!/usr/bin/env bash
# Prepare a reviewable initial Git commit without writing the workspace .git.
set -euo pipefail

project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
bundle_path="$project_dir/vamp-initial.bundle"
if [[ -e "$bundle_path" ]]; then
    printf 'Refusing to overwrite %s\n' "$bundle_path" >&2
    exit 1
fi

staging_dir=$(mktemp -d /tmp/vamp-git.XXXXXX)
for item in .gitignore README.md pyproject.toml src tests scripts docs devnotes configs .github; do
    cp -a -- "$project_dir/$item" "$staging_dir/"
done
mkdir -p "$staging_dir/reports"
for report in "$project_dir"/reports/*; do
    if [[ -f "$report" ]]; then
        cp -p -- "$report" "$staging_dir/reports/"
    fi
done

git -C "$staging_dir" init --initial-branch=main
git -C "$staging_dir" add -- .gitignore README.md pyproject.toml src tests scripts docs devnotes reports configs .github
if git -C "$staging_dir" ls-files | LC_ALL=C grep -Eq '(^data/|\.ms(/|$)|\.ms\.tar$|^\.(aws|codex|agents)/)'; then
    printf 'Unexpected observation or environment files in staged inventory\n' >&2
    exit 1
fi
git -C "$staging_dir" diff --cached --check
git -C "$staging_dir" commit --message 'Initialize vamp inspection and archive experiment specification'
git -C "$staging_dir" bundle create "$bundle_path.partial" main HEAD
git -C "$staging_dir" bundle verify "$bundle_path.partial"
# Hard-link publication refuses an existing target and only exposes a complete bundle.
ln -- "$bundle_path.partial" "$bundle_path"
unlink -- "$bundle_path.partial"
printf 'Bundle: %s\nPrepared repository: %s\n' "$bundle_path" "$staging_dir"
sha256sum -- "$bundle_path"
