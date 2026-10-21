# varquanMS

Adaptive visibility archive experiments for single-integration CASA Measurement Sets.

Status: **inspection milestone**. A read-only inspection CLI, synthetic MS fixtures,
and the first real-data inspection are available. Encoding, MS reconstruction,
quantization, compression benchmarks, and scientific-fidelity validation remain
planned. No compression ratio or scientific-compatibility claim is established.

The experiment specification is in [docs/design.md](docs/design.md). The first
required next step is a lossless archive and exact reconstruction test, before
introducing low-bit quantization.

## Initial sample

One MS was copied from
`/nas8/lwa/event-MS/20260904BigType35pol/59MHz/19/20260904_190008_59MHz.ms.tar`.
The source remains untouched. The local tar is under `data/raw/`, and the
extracted MS is `data/working/20260904_190008_59MHz.ms/`. These large files are
excluded from Git. Copy integrity is recorded in
[reports/sample-copy-verification.json](reports/sample-copy-verification.json).

The [inspection summary](reports/initial-inspection.txt) and
[machine-readable report](reports/initial-inspection.json) describe the actual input:

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

Use Python 3.10 or newer with NumPy and python-casacore. HDF5 is the current
prototype container candidate, available through the optional `archive` extra;
its selection remains subject to the lossless implementation.

```bash
python -m pip install -e .
visarchive inspect INPUT.ms --report inspection.json --summary inspection.txt
python -m unittest discover -s tests -v
```

The already-installed local environment can run without package installation:

```bash
PYTHONPATH=src /data1/users/pz47/miniconda3/bin/python -m varquanms inspect \
  data/working/20260904_190008_59MHz.ms \
  --report reports/new-inspection.json --summary reports/new-inspection.txt
PYTHONPATH=src /data1/users/pz47/miniconda3/bin/python -m unittest discover -s tests -v
```

Inspection defaults to selecting `DATA` for the future compression experiment.
Use `--columns DATA CORRECTED_DATA` to select other present visibility columns.
Selection does not compress or change them. The tool scans every defined cell,
reports all descriptors/keywords/storage managers and subtable values, and
records numeric special values. Numeric metadata arrays include dtype, shape,
and exact base64 bytes. `values` fields are readable previews.

Reports are created atomically, refuse overwrite, and cannot be written inside
the input MS. Unknown options fail. External subtable references and malformed
required references fail clearly. Heterogeneous cells, mixed contexts, duplicates,
and single-channel records are inspected and marked when they require explicit
future reconstruction maps. The layout classification is not archival acceptance.
Use a frozen snapshot: read-only `usernoread` access avoids lock-file writes but
does not protect against concurrent writers in another process.

Copy verification can be repeated independently:

```bash
PYTHONPATH=src python scripts/verify_sample_copy.py \
  /nas8/lwa/event-MS/20260904BigType35pol/59MHz/19/20260904_190008_59MHz.ms.tar \
  data/raw/20260904_190008_59MHz.ms.tar data/working \
  --report reports/new-copy-verification.json
```

## Publishing

The intended repository is public and named `varquanMS`. It could not be
created from the initial execution environment: GitHub DNS was unavailable,
GitHub CLI authentication failed, the GitHub connector was disabled by an
administrator, and the workspace `.git` directory was read-only. A Git bundle
is prepared separately so the initial commit can be imported on a writable,
authenticated host. See [docs/publishing.md](docs/publishing.md).
