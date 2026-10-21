# varquanMS

Adaptive visibility archive experiments for single-integration CASA Measurement Sets.

Status: **working research prototype**. Read-only inspection, lossless HDF5
archives, new-MS reconstruction, fixed/adaptive uniform 2/4/6/8-bit quantization,
exact exceptions, float32 fallback, baseline ordering and reproducible benchmarks
are implemented. The current suite has 36 passing tests. This is not a validated
default archival profile for arbitrary future science.

The experiment specification is in [docs/design.md](docs/design.md). The
codec/storage contract is in [docs/format.md](docs/format.md). Actual measurements
and the limits of the checks performed are in
[docs/experiment-results.md](docs/experiment-results.md).

The real MS passed exact lossless reconstruction. The initial lossless archive
is 331.2 MB versus 409.6 MB of source files (1.24×). One engineering adaptive
profile produced 171.2 MB (2.39×); it uses explicit absolute error budgets in
stored DATA units, not a validated thermal-noise model. Sorting and per-baseline
controls are recorded separately, including their complete byte costs.

## Initial sample

One MS was copied from
`/nas8/lwa/event-MS/20260904BigType35pol/59MHz/19/20260904_190008_59MHz.ms.tar`.
The source remains untouched. The local tar is under `data/raw/`, and the
extracted MS is `data/working/20260904_190008_59MHz.ms/`. These large files are
excluded from Git. Copy integrity is recorded in
[reports/sample-copy-verification.json](reports/sample-copy-verification.json).

The [inspection summary](reports/initial-inspection.txt) and
[machine-readable report](reports/initial-inspection.json) describe the actual input:

These initial reports predate codec implementation; their capability/status
notes are historical. A fresh `inspect` output reflects the current prototype.

| Property | Observed value |
| --- | --- |
| MAIN rows | 62,128 |
| Organization | One spectral window; channels inside each DATA cell |
| DATA cell | 192 channels × 4 products, complex64 |
| Product order | XX, XY, YX, YY |
| Antennas / autocorrelation rows | 352 / 352 |
| DATA logical bytes | 381,714,432 |
| MS physical file lengths | 409,594,838 bytes |
| TIME / TIME_CENTROID | One distinct value each; they differ |
| UVW | float64 metres, ITRF; 61,777 distinct bitwise vectors |
| Duplicate row keys / flagged DATA / nonfinite DATA | 0 / 0 / 0 |

All columns in this sample share a StandardStMan storage manager. Physical
file lengths cannot be attributed exactly to each column. UVW is stored once
per row and has no 192-fold channel-record repetition in this input. Its
deduplication benefit must be measured after including a reconstruction map.

## Usage

Use Python 3.10 or newer with NumPy, python-casacore and h5py. Install the optional
`plot` extra for the standalone scientific plotting scripts.

```bash
python -m pip install -e .
visarchive inspect INPUT.ms --report inspection.json --summary inspection.txt
visarchive encode INPUT.ms exact.h5 --config configs/lossless.json --report encode.json
visarchive decode exact.h5 RESTORED.ms
visarchive verify exact.h5 --source INPUT.ms --restored RESTORED.ms --report verification.json
python -m unittest discover -s tests -v
```

The already-installed local environment can run without package installation:

```bash
PYTHONPATH=src /data1/users/pz47/miniconda3/bin/python -m varquanms inspect \
  data/working/20260904_190008_59MHz.ms \
  --report reports/new-inspection.json --summary reports/new-inspection.txt
PYTHONPATH=src /data1/users/pz47/miniconda3/bin/python -m unittest discover -s tests -v
```

Inspection defaults to selecting `DATA` for compression experiments.
Use `--columns DATA CORRECTED_DATA` to select other present visibility columns.
Inspection itself does not compress or change them. The tool scans every defined cell,
reports all descriptors/keywords/storage managers and subtable values, and
records numeric special values. Numeric metadata arrays include dtype, shape,
and exact base64 bytes. `values` fields are readable previews.

Reports are created atomically, refuse overwrite, and cannot be written inside
the input MS. Unknown options fail. External subtable references and malformed
required references fail clearly. Heterogeneous cells, mixed contexts, duplicates,
and single-channel records are inspected. Lossless mode preserves these layouts
through original rows and explicit cells. The initial quantized profile rejects
them when they need additional partitions/maps. The layout classification is not archival acceptance.
Use a frozen snapshot: read-only `usernoread` access avoids lock-file writes but
does not protect against concurrent writers in another process.

Copy verification can be repeated independently:

```bash
PYTHONPATH=src python scripts/verify_sample_copy.py \
  /nas8/lwa/event-MS/20260904BigType35pol/59MHz/19/20260904_190008_59MHz.ms.tar \
  data/raw/20260904_190008_59MHz.ms.tar data/working \
  --report reports/new-copy-verification.json
```

## Quantization and ordering

`configs/engineering-adaptive.json` is an explicit engineering experiment, not
an endorsed scientific tolerance. The default configuration is lossless. Lossy
mode requires MSE and maximum-error criteria; per-baseline and per-product bounds,
bias/p99 checks, and a user-supplied validated noise model are configurable.
Weights and sigmas are never treated as noise estimates automatically.

```bash
visarchive encode INPUT.ms experiment.h5 --config configs/engineering-adaptive.json
visarchive decode experiment.h5 APPROXIMATE.ms
visarchive verify experiment.h5 --source INPUT.ms --restored APPROXIMATE.ms --report errors.json
visarchive benchmark INPUT.ms --matrix configs/benchmark-smoke.json --output results/my-run \
  --experiments lossless-no-dedup per-baseline-8bit grouped-amplitude grouped-buckets
```

Uniform blocks share a stored float32 real/imag step, with no offset. Codes are
truly packed, including six-bit values. Flags, nonfinite values, signed zeros and
outliers retain exact complex64 bits. Every candidate is decoded and checked;
the smallest passing encoded block competes with exact float32 storage.
Reconstructed float32 values remain approximations for quantized samples.

Orderings include original, baseline length, stable amplitude and logarithmic
buckets, with independent or shared product order. Spectra remain intact. Actual
permutations/bucket boundaries and original row presence are preserved. Quantized
mode currently requires one integration and one regular configuration. It holds
selected columns in RAM (configurable 2 GiB per-column limit), rather than offering
an out-of-core encoder.

Selected archive columns currently must be DATA, CORRECTED_DATA or MODEL_DATA.
FLOAT_DATA, when present alongside these columns, is preserved as an exact
auxiliary column. FLOAT_DATA-only inputs are not supported by the archive profile.

Reports count payloads, indexes, permutations, flagged/exception values, scales,
quality arrays, metadata and container overhead. Before/after secondary-compression
counts are distinguished. Errors include product/baseline/channel summaries and
error covariance across products and adjacent channels. Unknown configuration
options, codebooks and rounding modes fail. Distribution-adapted codebooks,
subtractive dithering, offsets, scale quantization and Dysco comparison remain
future experiments.

## Verification and CASA scope

Lossless tests cover both channel organizations, missing/duplicate records,
heterogeneous shapes, mixed fields, antenna order, optional metadata, special bits
and multiple I/O chunks. Codec tests cover every width, partial blocks, corruption,
stochastic streams, sorting, exact escape and independent block reads. Literal
fixtures and an independent standard-library block decoder are included.

The real lossless and adaptive exports were reopened with python-casacore and
CASA 5.4.1. A 64×64 natural-weight MFS Stokes I dirty-image workflow was run on
both. CASA warned about outdated IERS/leap-second data and missing OVRO-LWA
observatory metadata. The image comparison is a smoke test, not precision
astrometry, frequency conversion, calibration or weak-polarization acceptance.
See [image comparison](reports/image-comparison.json).

CASA imaging writes HISTORY even when scratch model columns are disabled.
`scripts/check_casa.py` therefore creates a separate working copy. Pristine
verified exports are `results/lossless/verified.ms` and
`results/adaptive/verified.ms`; experiment outputs are excluded from Git.

The source tar and all 65 extracted files were checksum-verified again after
codec tests, confirming the source copy stayed untouched.

## Publishing

The intended repository is public and named `varquanMS`. It could not be
created from the initial execution environment: GitHub DNS was unavailable,
GitHub CLI authentication failed, the GitHub connector was disabled by an
administrator, and the workspace `.git` directory was read-only. A Git bundle
is prepared separately so the completed prototype can be imported on a writable,
authenticated host. See [docs/publishing.md](docs/publishing.md).
