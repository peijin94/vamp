# Rename to vamp

The repository, distribution, Python package, CLI and paper are named `vamp`.
The paper title is *Visibility-Adaptive Mixed Precision for Radio Interferometric
Measurement Sets*. The former `visarchive` command and `varquanms` Python
namespace are compatibility aliases; new code should use `vamp`.

Version-1 archives retain their established `varquanms.*` on-disk schema and
format IDs and their `VARQUANMS_*` MS provenance keywords. Those identifiers
are part of the stored format, not the project name, and changing them during a
rename would make existing archives or their reconstructed metadata harder to
interpret. New benchmark configurations and non-archive report schemas use
`vamp.*`; the benchmark also accepts the old configuration identifier.
Historical experiment reports are retained byte-for-byte, including original
software labels and checksums. Renderers accept both old and new report labels.

## Validation

An existing 282,825,653-byte calibrated archive produced before the rename was
decoded by vamp without a source MS. The reconstructed 931,920-row MS passed
146 checks against the original source, including exact UVW, flags, weights,
sigma and metadata, plus a bitwise comparison of the reconstructed calibrated
column against the archive decoder. The source-file snapshot matched the
pre-rename snapshot.

The earlier and newly reconstructed exports were then compared directly:
all **147 column checks were bitwise exact**, including `DATA` and
`CORRECTED_DATA`. This shows that the rename did not change decoded values
for this archive. The [machine-readable validation](../reports/vamp-rename-validation.json)
records the archive SHA-256 and individual checks. The calibrated-only archive
already discarded flagged visibility values under its stated policy.

The Dysco image, zoom and bit-distribution panels were regenerated from the
same saved FITS images and measurements with the new labels. No compression,
CLEAN imaging, calibration or scientific threshold was rerun or changed for
this rename. Historical measurements and their recorded checksums remain intact.

```bash
python -m pip install -e '.[validation]'
vamp --help
python -m unittest discover -s tests -v
PYTHONPATH=src:scripts python scripts/build_dysco_paper_comparison.py render
PYTHONPATH=src:scripts python scripts/update_dysco_paper.py
python scripts/compile_paper.py
```
