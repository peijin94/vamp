# Archive and codec format, version 1

The prototype uses scientific schema `varquanms.archive/1` in an ordinary HDF5
file. Root attributes identify `varquanms.hdf5`, schema version 1, a complete
flag and the manifest SHA-256. The `manifest` dataset contains strict UTF-8 JSON.
No pickle, evaluation, network access, model or downloaded codebook is used.
Decode with `visarchive decode ARCHIVE OUTPUT.ms`; NumPy, h5py and python-casacore
are the only runtime dependencies. These packages must already be installed.

## Scientific contents

The manifest contains the source layout, original table descriptors and storage
assignments, table info, configuration, software/provenance, timings, table/column
reconstruction specifications and a dataset inventory. All stored numeric
payloads are explicitly little endian. The inventory records path, dtype,
shape, SHA-256, logical bytes, allocated HDF5 dataset bytes and byte category.
Whole-file bytes minus dataset allocations count all container/manifest overhead.

HDF5 numeric array chunks target 4 MiB; a single cell is limited to 64 MiB.
Chunk Fletcher32 checksums supplement inventory SHA-256. Visibility I/O streams
use independently compressed quantization blocks in 1 MiB HDF5 chunks, with an
offset index and CRC32 for each quantization block. Reading one block never
requires decoding an earlier block. Manifest SHA-256 detects accidental changes;
it is not an authentication signature. External/soft HDF5 links are rejected.

Each MS table keeps its original row count and column order. Exact columns use
dense arrays, a constant cell, a bitwise dictionary with uint32 row map, or
explicit typed cell records. Undefined cells and defined empty cells differ.
Dictionary selection currently uses logical size before conventional compression;
the benchmark reports actual physical results. UVW is never regenerated or
merged by tolerance. Source UVW stays in metres with its original measure.

The regular quantized axis convention is `[original row/baseline, channel,
correlation]`. A product-specific permutation stores `perm[j]`, the original row
at archived position j. Restore by `restored[perm,:,product] = decoded`.
Identity permutations are implicit, and shared product order uses one stored
permutation. Group boundary arrays define exact partial edges and auto/cross,
bucket and no-valid-data boundaries. Geometry/antenna grouping also respects its
configured radial/angular/endpoint boundaries. Each product records the grouping
definition and a group-size histogram; the stored permutation and boundaries,
not a re-run of the sort algorithm, determine reconstruction.
Frequency segments have exact final edges.
Original correlation and frequency order come from unchanged MS metadata.
The row-presence rule is identity to the original MS; no absent rows are padded.

Metadata uses tagged JSON records, arrays, strings, lists and scalar values.
Numeric arrays include dtype, shape and base64 C-order bits. Numeric scalar kinds
and special bits are retained. Table references are relative tokens rebased to
the new output table graph. Export materializes columns using StandardStMan;
storage layout is outside exactness. Table/column logical descriptors, keywords,
units, references, subtable relations and supported auxiliary values remain exact,
except for the explicitly enabled destructive flag profile below.

### Optional destructive MAIN flags (`burned-flags-v1`)

Configuration `burn_flag: true` replaces MAIN `FLAG`, `FLAG_CATEGORY` and
`FLAG_ROW` values with reconstruction-only column specifications. Each contains
`encoding: "burned-flags-v1"`, `row_count`, `dtype: "|b1"`, and either a common
`cell_shape` or `shape_chunks`. A shape chunk records `start`, `stop` and an
inventoried JSON dataset `path`; its per-row values are dimension lists or null
for undefined cells. No original Boolean payload is stored. Decoding fills
defined cells with false and leaves undefined cells undefined. Empty defined
cells retain their shape. Subtable flags are outside this policy and stay exact.
Readers reject a mismatch between the opt-in and column encodings; older readers
that lack this encoding fail rather than silently recover an exact flag mask.

The archive config and provenance and all reports disclose the flag policy.
Export adds MAIN keyword `VARQUANMS_FLAG_POLICY` with the encoding identifier,
affected columns, reconstruction description and `reflagging_performed: false`.
Verification permits exactly this metadata addition and checks generated flags
against the declared reconstruction, not source flag-value equality. Other
metadata and weights remain exact. A lossless visibility mode plus burned flags
is not a lossless MS round trip. Quantization still protects source-flagged
visibility values and computes metrics using source flags; only the masks are
discarded. See [development notes](../devnotes/burn-flag.md) for use and risks.

## Uniform block header

Each block starts with the 16-byte little-endian structure `<BBHIIf`:

| Field | Type | Meaning |
| --- | --- | --- |
| bits | uint8 | 0 for exact complex64; otherwise 2, 4, 6 or 8 per component |
| compression | uint8 | 0 uncompressed body; 1 zlib-wrapped DEFLATE body |
| reserved | uint16 | Must be zero |
| samples | uint32 | Number of complex samples, 1 through 262,144 |
| exceptions | uint32 | Number of exact exception samples |
| step | float32 | Nonnegative stored uniform step, or zero for an exact block |

The configuration's `gzip` selection uses the HDF5 gzip filter for ordinary
arrays and Python zlib for independently compressed block bodies. A block uses
secondary compression only when it shortens its body. Body decoding is bounded
by the exact size derived from the header; truncation and trailing bytes fail.

For b bits, let `M = 2**(b-1)-1`. Codes 0 through 2M reconstruct signed integer
levels -M through M. The final code `2**b-1` is reserved and rejected. Thus the
codebook is exactly symmetric and represents zero, at the cost of one unused
bit pattern. The immutable codec identifier is `uniform-symmetric-v1`.

Decode each real and imaginary component separately in this exact order:

```text
signed = int16(code) - M
component = float32(signed) * stored_float32_step   # one float32 multiplication
```

Assemble complex64 by interleaving the resulting float32 component bits; do not
use complex arithmetic that could change signed zeros. There is no offset.
Nonfinite reconstruction outside the exception stream is rejected. For a zero
step, encoded ordinary zeros use code M. Signed zeros that would change bits are
exact exceptions. Float32 escape blocks contain the original little-endian
complex64 bits and require no exception stream.

The quantized raw body contains:

1. `ceil(2 * samples * bits / 8)` packed-code bytes.
2. `exceptions` uint32 positions, strictly increasing within the flattened
   baseline-major, channel-minor complex sample domain.
3. `exceptions` original complex64 values, preserving both component bit patterns.

Code order is real, imaginary, real, imaginary. Codes and bytes are least
significant bit first. Unused high bits in the last byte must be zero. Six-bit
codes really occupy six bits. Exceptions replace whole complex samples after
ordinary reconstruction. Flags/FLAG_ROW, NaN/Inf, signed zeros and out-of-range
samples bypass scale estimation as applicable and retain their original bits.
Ordinary representable zeros need no exception.

## Encoder choices and error domain

Scale candidates use the maximum absolute real/imag component of each eligible
sample and configured quantiles of that statistic. Step is the cutoff divided
by M, rounded to float32 and advanced to the next float32 when necessary to cover
the cutoff. One step serves both components. Outliers become exact exceptions;
there is no silent saturation. Candidates are evaluated using the actual decoder
and stored step. Total serialized header/body bytes plus the fixed index cost
determine selection; exact storage competes on the same basis. HDF5 allocation,
permutations, counts and quality-array overhead are also reported in file totals.

Nearest mode divides in float64 and uses ties-to-even rounding. The optional
`stochastic-pcg64-v1` mode chooses between adjacent uniform integer levels using
their interpolation probability. PCG64 uses NumPy SeedSequence words
`[seed_low32, seed_high32, selected_column_index, product, baseline_group_index,
frequency_segment_index, bits, scale_quantile_index, component]`. Seeds are
configured, not searched for better errors. Codes contain the outcome; decoding
needs no randomness. This is not subtractive dithering. Floating arithmetic,
exceptions and finite support qualify any theoretical unbiasedness statement.

MSE/bias/percentile checks use unflagged finite samples. Flagged, nonfinite and
signed-zero sample values are checked for exact bit restoration independently.
Checks include per-block component MSE, complex maximum/bias, optional p99, and
per-baseline MSE. Whole-product metrics and arrays by baseline/frequency are
recomputed before publication. Verification recomputes the same statistics from
the source and checks every block. Original weights are never adjusted.

Quality arrays include baseline MSE/max error, channel MSE, centered complex error
covariance across products and adjacent channels, paired counts, amplitude medians
and valid counts. All their bytes count toward the visibility archive totals.
Code/exception categories before joint secondary compression are exact logical
counts; physical bytes of a jointly compressed body cannot be attributed uniquely
to its ingredients. Reports make this distinction explicit.

## Support boundaries and durable reconstruction

With default flag preservation, lossless mode preserves both channel organizations, heterogeneous cells,
duplicates, missing row records and mixed fields/configurations through original
table rows and explicit cells. Quantized mode currently requires one integration,
one data description, one distinguishing row context, uniform complex64 visibility
cells, full spectra inside cells and no duplicate row keys. Unsupported cases
fail clearly. External subtables, unreferenced nested tables, non-table files,
record-valued columns and nonempty private table keywords are unsupported.

Source snapshots must be frozen. The prototype holds selected quantized columns
in RAM, with a configurable size limit (2 GiB per column by default); it is not
an out-of-core visibility encoder. Error budgets are required for lossy modes.
The provided engineering budgets do not validate thermal noise or science.
Nonuniform codebooks, complex offsets, subtractive dithering, scale quantization,
predictors, BDA, implicit Stokes conversion and a virtual MS storage manager are
not implemented; unsupported codebook/rounding names are rejected.

Incomplete archives carry `.partial` filenames and `complete=false`. An archive
is published only after metadata and inventory finalization. Export creates a
new staging directory, restores values and all table references, reopens tables,
then renames the directory. Existing outputs and writes inside the source are
rejected. Physical casacore files need not match the source storage layout.

Literal conformance bytes are in `tests/fixtures/codec-v1.json`. The independent
standard-library block decoder in `scripts/reference_block_decoder.py` is tested
against the NumPy decoder, including exact special-value exceptions. The full
reference MS reconstruction path lives in `src/varquanms/archive.py`.
