# Archive and codec format, version 1

The prototype uses scientific schema `varquanms.archive/1` in an ordinary HDF5
file. Root attributes identify `varquanms.hdf5`, schema version 1, a complete
flag and the manifest SHA-256. The `manifest` dataset contains strict UTF-8 JSON.
No pickle, evaluation, network access, model or downloaded codebook is used.
Decode with `visarchive decode ARCHIVE OUTPUT.ms`; NumPy, h5py and python-casacore
are the only runtime dependencies. These packages must already be installed.

## Encoder error limits

Absolute limits retain their existing units and whole-product checks. Optional
`relative_bounds` uses `reference: "block-rms-v1"`. Its reference is
`S = sqrt(mean(real(V)^2 + imag(V)^2)/2)` over unflagged finite samples in the
original quantization block, evaluated in float64. The signal center is included;
this is neither residual spread nor a noise estimate. Available fractional limits
are `max_component_rmse`, `max_baseline_component_rmse`, `max_complex_error`,
`max_complex_bias` and `max_p99_error`. RMSE fractions become `(fraction*S)^2`
MSE limits; the other fractions become `fraction*S`. The baseline reference is
the same block S, not a separate per-baseline denominator. If omitted, the
baseline RMSE fraction defaults to the component RMSE fraction.

Every configured absolute and relative limit must pass: the effective threshold
is their minimum. Relative-only mode requires component RMSE and maximum-error
limits. Zero-reference blocks get zero relative limits, without a hidden floor;
flagged and nonfinite values never determine S. Overflow fails explicitly. This
is a block-normalized criterion, not a bound on each sample's fractional error.
Source verification recomputes each block's reference and measured errors.
Relative criteria are block-local; absolute criteria also retain their existing
whole-product checks. Stored configuration fully records the policy. Payload,
decoder arithmetic, exceptions and exact fallback are unchanged.

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
| bits | uint8 | 0 for exact complex64; otherwise 4, 8, 12 or 16 per component (legacy 2/6 supported) |
| compression | uint8 | 0 uncompressed body; 1 zlib-wrapped DEFLATE body |
| mapping | uint16 | 0: `uniform-symmetric-v1`; 1: `uniform-midrise-v1`; 2: `power-midrise-v1`; exact blocks require 0 |
| samples | uint32 | Number of complex samples, 1 through 262,144 |
| exceptions | uint32 | Number of exact exception samples |
| step | float32 | Nonnegative stored uniform step, or zero for an exact block |

The configuration's `gzip` selection uses the HDF5 gzip filter for ordinary
arrays and zlib for independently compressed block bodies. A block uses
secondary compression only when it shortens its body. Body decoding is bounded
by the exact size derived from the header; truncation and trailing bytes fail.

For b bits, let `M = 2**(b-1)-1`. Codes 0 through 2M reconstruct signed integer
levels -M through M. The final code `2**b-1` is reserved and rejected. Thus the
codebook is exactly symmetric and represents zero, at the cost of one unused
bit pattern. The immutable codec identifier is `uniform-symmetric-v1`.
Mapping 0 is the original header's zero reserved field. Existing archives and
their numerical interpretation are unchanged. Old decoders reject mapping 1.
The default candidate ladder is 4/8/12/16 bits per component. Legacy 2/6-bit
archives and explicitly selected legacy configurations remain supported. Adding
12/16-bit widths does not change existing codebooks, block headers or packed
bytes at older widths; older decoders reject unsupported widths explicitly.

Decode each real and imaginary component separately in this exact order:

```text
signed = int32(code) - M
component = float32(signed) * stored_float32_step   # one float32 multiplication
```

Assemble complex64 by interleaving the resulting float32 component bits; do not
use complex arithmetic that could change signed zeros. There is no offset.
Nonfinite reconstruction outside the exception stream is rejected. For a zero
step, encoded ordinary zeros use code M. Signed zeros that would change bits are
exact exceptions. Float32 escape blocks contain the original little-endian
complex64 bits and require no exception stream.

### Optional full-range unsigned mapping

`uniform-midrise-v1` uses mapping 1 and every unsigned code, including all ones:

```text
N = 2**b
C = (N - 1) / 2                       # exact half-integer, not an integer
step = ceil_float32(cutoff / C)       # upward float32 rounding, not integer ceil
code = rint_even(float64(x) / float64(step) + C)
coefficient = float32(code) - float32(C)
component = coefficient * stored_float32_step  # separate float32 multiply
```

Codes lie in `[0,N-1]`; reconstruction levels span `[-C*step,+C*step]`.
Values outside this interval use exact exceptions, with no silent clipping.
The stored step serves both components and the existing actual-byte selection,
error checks, and exact fallback still apply. The manifest codebook must agree
with every quantized block's mapping; exact blocks use mapping 0 under either
codebook. Unknown mappings fail explicitly.

Nonzero-step levels straddle zero. With nearest-even, ordinary positive zero
encodes as `N/2` and reconstructs as `+step/2` for supported widths. This is an
explicit opt-in loss of ordinary zero, which can bias zero-rich data. Negative
zero samples, flagged samples and nonfinite values remain exact exceptions.
Zero-step blocks encode ordinary zeros as `N/2` and reconstruct positive zero.
Stochastic rounding operates on the shifted code coordinate above with the same
independent component streams as the symmetric codec; no decoder RNG is needed.
Exception placeholders also use `N/2` before their original bits are restored.

All packing, byte order, maximum block size and exception rules below are shared.
There is no extra sign bit, exponent field, or separately stored offset. At a
common cutoff the step ratio to the default codec is `(N-2)/(N-1)`, a small
improvement at high precision. This codec is optional, not a default replacement.

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
Twelve-bit codes really occupy twelve bits (three bytes per complex sample);
sixteen-bit codes occupy two little-endian bytes per component. Sixteen-bit
quantization is scaled integer quantization, not IEEE float16 or lossless storage.

## Encoder choices and error domain

### Median-centered power-law blocks

`power-midrise-v1` uses mapping 2. The base header's final float32 stores the
nonnegative residual **range t**, not a quantization step. Three additional
little-endian float32 values immediately follow the base header, uncompressed:
`center_real, center_imag, gamma`. Total header size is **28 bytes**. The packed
codes/exception body has the same layout and compression identifiers as above,
starting at byte 28. Exact blocks still use the original 16-byte mapping-0 header.
All unsigned codes `0,...,2**b-1` are valid. Unknown/nonfinite centers, negative or
nonfinite t, unsupported gamma, malformed lengths and nonfinite decoded values
are rejected.

For each quantization block (one correlation product), the default center is the
componentwise median over eligible samples: separate real/imaginary medians of
unflagged, finite, non-negative-zero-protected samples. Medians are calculated in
float64 and stored as float32 **before** residual/scale estimation. `power_center:
"zero"` is a control. The range candidates are max/.995 (or configured quantiles)
of `max(abs(real residual),abs(imag residual))`, rounded upward to float32. A
sample beyond ±t in either component is stored as a complete exact exception.
No implicit clipping occurs. An empty eligible domain uses center=(0,0), t=0.

Write `M=2**b-1`. In exact arithmetic the signed compander is
`y=sign(x-center)*abs((x-center)/t)**gamma`, with unsigned coordinate
`(y+1)*M/2`. Gamma below 1 increases level density near the center; it is not
statistical skewness and does not favor one sign. Initially gamma is a unique
configured subset of **[1, 0.5, 0.25]**, stored exactly as float32. Other powers
fail explicitly. The default candidate set is [1, 0.5]. Its inverse requires
only fixed products and no external codebooks or implementation-dependent pow.

Decode each component in exactly this order, IEEE round-to-nearest/even:

```text
u = float64(2*integer_code - M) / float64(M)
a = abs(u)
if gamma != 1: a = float64(a*a)
if gamma == 0.25: a = float64(a*a)
r = float64(copysign(a,u) * float64(stored_t))
x_hat = float32(float64(r + float64(stored_component_center)))
```

Multiplication/addition are separate, not fused. Only the final result is rounded
to float32. t=0 encodes the constant component centers with code `2**(b-1)`.
For t=0 the decoder copies each stored center's float32 bits directly, bypassing
the arithmetic above (including preserving a negative-zero center).
There is no exact center level at nonzero t in ideal arithmetic; its two nearest
levels lie at `center ± t/M**(1/gamma)`. Final float32 rounding can collapse
nearby levels. Ordinary positive-zero values may change, as with the opt-in
midrise codec; source negative zeros, flags and nonfinite bits remain protected.

The encoder forms the companded coordinate in float64 using zero, one or two
square roots. It evaluates both bracketing codes with the actual decoder and
selects the nearest **original-value-space** result; distance ties choose an
even code. Thus `nearest-even` denotes nearest reconstructed value here, not
rounding the companded coordinate. This is deliberate: transformed-space rounding
does not minimize visibility-domain squared error. Stochastic rounding is not
implemented for this codec and is rejected before writing an archive.

Search order is ascending bits, configured scale quantiles, then configured
gammas. Actual compressed header+body size decides among candidates satisfying
all existing measured bounds; the first wins exact byte ties. Exact float32
always competes. Product summaries record gamma/block counts, all 28-byte headers,
exception bytes and remaining archive categories. The public decoder uses only
archive contents and these operations. Defaults remain the original codebook;
select `configs/power-adaptive.json` for this experiment.

Optional encoder-only `power_gamma_policy` values `max-std-v1` and
`max-std-inlier-v1` narrow the gamma candidates; the default `search` retains
the configured search above. For residual pairs after subtracting stored
component centers (medians by default),
define `sigma=sqrt((var(real)+var(imag))/2)` using float64 population variances,
and `R=max(abs(component residual))/sigma`. Exclude the same protected samples
used for scale estimation. The inlier variant additionally excludes samples
whose maximum component residual exceeds the candidate's stored float32 radius.
This variant computes a separate R for each radius, before the bit-width search.

With `power_ratio_reference=R0` (default 4), choose gamma 1 when `R<=R0`,
1/2 when `R0<R<=2*R0`, and 1/4 otherwise. Empty or zero-variance populations
choose 1. If the configured `power_gammas` excludes that choice, select the
nearest allowed value in log2 space, preferring the larger gamma at a tie.
R is a signal-spread statistic, not a thermal-noise estimate or a scientific
error bound. All reconstruction checks, exception handling and exact fallback
still apply. The chosen gamma is already stored in the block header; decoder
bytes, version and arithmetic are unchanged. These options are recorded in the
archive config and rejected for other codebooks when a ratio policy is selected.

### Uniform encoder choices

Scale candidates use the maximum absolute real/imag component of each eligible
sample and configured quantiles of that statistic. Step is the cutoff divided
by M (or C for midrise), rounded to float32 and advanced to the next float32 when necessary to cover
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

### Optional CPU acceleration and quality arithmetic

An optional C extension accelerates the encoder; it does not change codec IDs,
bit packing, exception semantics or decoder dependencies. It uses float64
division, ties-to-even integer rounding and the prescribed float32 reconstruction,
with fast-math and floating contraction disabled. The complete candidate search
uses the same zlib bodies and encoded-byte comparison. Conservative native
prechecks avoid expensive compression of failing candidates. Every winner is
decoded by the portable decoder and checked strictly with NumPy. A candidate
near a numeric bound that fails that check triggers the full reference search.
Stochastic rounding and configured p99 bounds retain the reference selection path.

`workers` (1–64, default 1) controls bounded CPU pools. Block submission order,
PRNG stream IDs and output order are deterministic. Only the calling thread uses
casacore and HDF5; worker queues hold at most twice the configured worker count.
Worker elapsed-time sums overlap and must not be added to pipeline wall time.
Native `*_cpu` counters measure thread CPU time, separately from elapsed time.

Quality summaries marked `calculation: "native-moments-v2"` use fixed 1024-row
partitions, ordered float64 moment reductions, paired masks and exact protected-bit
checks. Worker count does not change these reductions. This changes floating
reduction order from the older NumPy implementation, so stored quality numbers
may differ in their last bits while the decoded visibilities remain identical.
Archives without this marker retain the original `numpy-v1` calculation.
Verification with the native implementation compares its quality outputs exactly.
Portable verification without the extension recomputes NumPy moments with relative
tolerance `1e-10` and absolute tolerance `1e-12` times the error variance (or its
square root for amplitude/bias quantities). Counts and protected values remain
exact, all source metadata/coordinates remain exact, and reconstruction-error
thresholds are still enforced without relaxed bounds. Inventory checksums also
cover every quality array. Unknown calculation IDs are rejected.

The writer caches two 4 MiB HDF5 chunks so partial row writes do not repeatedly
compress a chunk larger than HDF5's cache. The archive still includes chunk,
dataset and manifest integrity checks and fsync before final publication.

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
The optional median/power codebook supplies a complex offset and nonuniform
levels. Subtractive dithering, scale quantization, predictors, BDA, implicit
Stokes conversion and a virtual MS storage manager are not implemented;
unsupported codebook/rounding names are rejected.

Incomplete archives carry `.partial` filenames and `complete=false`. An archive
is published only after metadata and inventory finalization. Export creates a
new staging directory, restores values and all table references, reopens tables,
then renames the directory. Existing outputs and writes inside the source are
rejected. Physical casacore files need not match the source storage layout.

Literal conformance bytes are in `tests/fixtures/codec-v1.json`. The independent
standard-library block decoder in `scripts/reference_block_decoder.py` is tested
against the NumPy decoder, including exact special-value exceptions. The full
reference MS reconstruction path lives in `src/varquanms/archive.py`.

## Calibrated-only extension and ordered partition bundles

The explicit `keep_calib_only` configuration/CLI option omits MAIN DATA,
FLOAT_DATA and MODEL_DATA, selects CORRECTED_DATA, and replaces values selected
by FLAG OR FLAG_ROW with positive complex zero before encoding. Row/cell shape,
presence, flag bits and original UVW are retained. The original flagged values
cannot be recovered. This option is incompatible with `burn_flag`; the default
archive profile is unchanged. Exported descriptors contain the
`VARQUANMS_CALIBRATED_ONLY` and `VARQUANMS_WEIGHT_POLICY` provenance records.

By default, decoding this profile adds a DATA column as an exact copy of decoded
CORRECTED_DATA, with a `VARQUANMS_DERIVED_DATA` record stating that original raw
DATA was not restored. This supplies the nominal data column required by MS v2
and the installed CASA imager. It adds no DATA payload to the archive.
`--no-ms-v2-data` requests the limited CORRECTED_DATA-only projection explicitly;
WSClean can image it with `-data-column CORRECTED_DATA`. Verification checks both
the declared source projection and bitwise equality of the derived DATA copy.

For this profile, WEIGHT and WEIGHT_SPECTRUM use `weight-codebook4-v1` by default.
Each column specification records dtype `<f4`, cell shape, total row count and
contiguous independently decoded chunk records. A `codebook4` record contains
`start`, `stop`, the exact scalar `count`, and inventory paths to a float32
codebook (1–16 entries) and a packed byte array. Scalars are flattened in C order;
the first code occupies the low nibble. Unused high-nibble padding must be zero.
Decoding is the exact float32 lookup `book[code]`, with no decoder fitting or
floating arithmetic. Codes outside the codebook, inconsistent counts and
missing/overlapping chunks are errors. The alternative `float32` chunk kind
references an exact dense float32 array with the expected row/cell shape.

For at most 16 distinct original float32 bit patterns, the encoder stores those
patterns, ordered by their unsigned little-endian uint32 values. Otherwise it
tries 16 stored float32 representatives uniformly spanning the finite nonnegative
range. It chooses the nearest stored representative, with ties to the even code.
The default maximum relative error per positive weight is .01; zero and special
weights must stay exact, and positive weights must remain positive. Failure uses
float32. `weight_bits=32` selects the existing exact column encodings;
`weight_max_relative_error=0` forbids any numerical weight error. SIGMA and
SIGMA_SPECTRUM remain exact. Quality records include achieved error extrema,
lossy/exact/fallback counts and packed byte counts. HDF5 inventory/chunk checksums
cover codebooks and codes. This is a numerical storage policy, not a validated
noise or sensitivity model.

The experimental `varquanms.ordered-partitions/1` directory bundle stores an
ordered list of ordinary archives with contiguous `[start,stop)` source row
ranges, archive lengths, SHA-256 hashes and checksummed verification reports.
`manifest.sha256` covers `manifest.json`; the complete directory is published
only when every partition succeeds. Each partition repeats the full exact
subtable graph. Decode checks coverage, row counts, archive configurations,
logical descriptors, subtable content hashes and every archive inventory before
assembling the rows. It needs no original MS. The wrapper supports at most 4,096
contiguous context runs, and each run must satisfy the ordinary archive preflight.
Repeated subtable bytes and all outer manifest/report overhead count toward
bundle size. No time or channel averaging occurs.

Detailed policy, commands and calibration caveats are in
[the calibrated imaging development note](../devnotes/calibrated-imaging.md).

### Shared partition storage (`varquanms.shared-partitions/1`)

`python scripts/timeslot_archive.py compact INPUT.bundle OUTPUT.h5 --report report.json`
losslessly repacks an already verified directory bundle into one HDF5 file.
Decode either representation with the same script's `decode` operation.
This is a storage transformation, with no second quantization pass.

Root attributes are `format=varquanms.shared-partitions`, `schema_version=1`,
`complete`, and `manifest_sha256`. The root `manifest` is a gzip-compressed
uint8 UTF-8 JSON array, at most 32 MiB. It retains the ordered bundle configuration,
row ranges and original source hashes, adds storage-compaction provenance, and
identifies `partitions/pNNNN` and `verification/vNNNN` for each partition.
The verification datasets preserve the original JSON bytes and SHA-256.

Each partition group contains the unchanged scientific manifest, archive
attributes and relative payload paths of an ordinary `varquanms.archive/1`.
Arrays with equal roles, shapes, dtypes, verified SHA-256 and bitwise contents
share one HDF5 dataset through hard links. Numerical equality alone is insufficient:
signed zeros and NaN payloads remain distinct. No external/soft links, virtual
datasets, external dataset storage or aliased/cyclic groups are supported.
Later partitions and their quantization blocks remain directly addressable.

The full reader validates links, row coverage, configurations, subtable/schema
agreement, both levels of manifests, verification reports and every logical
array before creating an MS. Missing or corrupt members fail before publication.
Source bundles remain untouched; construction uses a `.partial` file and a new
final destination. NumPy, h5py and python-casacore suffice for independent export.

Original per-partition byte counts remain provenance, not the shared file's
physical accounting. `shared_summary` counts each unique HDF5 object once;
all remaining file bytes count as container, links, manifests and verification
overhead. Compaction time is additional to the original MS encoding time.
