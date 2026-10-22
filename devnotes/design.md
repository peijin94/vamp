# Adaptive visibility archive: experiment specification

Status: implementation proposal. Inspection has been run on the initial sample.
The first lossless/adaptive measurements are in `experiment-results.md`. This
specification still describes the broader planned experiment, and its future
steps must not be mistaken for implemented or validated capabilities.

## Objective and hard requirements

Build an independently decodable archive that expands into a valid CASA MS v2.
The immediate experiment uses one integration, roughly 62,000 baselines,
192 channels, and four complex64 products. The actual sample is documented in
the inspection reports; its layout takes precedence over these estimates.

- Preserve UVW exactly as original float64 metre values and preserve measure
  references. Never regenerate it from antenna positions or time. Share vectors
  only when their bytes match; retain geometry-record and row mappings.
- Preserve time and frequency resolution, ordered antenna pairs, feed/field/
  data-description/scan/array/observation/processor/state context, original row,
  channel, product order, duplicates, and row presence. Missing rows and flagged
  present rows are different. Do not conjugate or swap antennas to canonicalize.
- Read correlation conventions from metadata. Support linear products and
  Stokes without implicit conversion. Keep incompatible configurations in
  separate partitions; reject unimplemented layouts before finalizing archives.
- Preserve flags, weights, sigma, metadata, auxiliary columns and flagged sample
  values losslessly in the first profile. Preserve nonfinite and signed-zero
  bits. Store original weights unchanged; quantization errors are separate.
- Quantization is explicitly lossy. Reconstructed float32 values do not regain
  discarded precision. Include an exact complex64 escape block whenever a
  candidate fails checks or costs more than lossless storage.
- Export to a new MS, never overwrite the source. Mark incomplete output clearly.
  Exactness concerns supported values/metadata, not casacore physical file layout.

## Primary sources and prior art

Read before codec implementation:

1. [MS v2 specification](https://casacore.github.io/casacore-notes/229.html).
2. [Offringa (2016), Dysco](https://arxiv.org/abs/1609.02019),
   [full text](https://arxiv.org/html/1609.02019v1).
3. [Casacore table and storage-manager overview](https://casacore.github.io/casacore-notes/255.html).
4. [QLoRA](https://arxiv.org/abs/2305.14314).
5. [XRADIO/MS v4 overview](https://xradio.readthedocs.io/en/latest/measurement_set/overview.html).

Dysco already combines normalization, nonlinear quantization and dithering,
including factorized antenna/frequency scales. Investigate adaptive precision,
grouping, exceptions, explicit error accounting and archival reconstruction;
do not claim Dysco's ingredients as novel. Its published observations provide
benchmarks, not a fidelity promise for this sample. NF4 is an experimental
codebook candidate; LLM task accuracy does not demonstrate scientific fidelity.
MS v4 motivates labeled partitioned arrays but is not an exact MS v2 replacement.

## Inspection and support boundaries

Record every column descriptor, dtype, cell shape, units, measure reference,
table/column keyword, storage-manager assignment and optional column. Scan all
rows for distinct times, identity context, autos, flags/weights, nonfinite and
signed-zero values, all-flagged regions, duplicates and heterogeneous shapes.
Inventory spectral frequencies, widths, references, CORR_TYPE and CORR_PRODUCT.
Report measured file lengths/allocations and logical array bytes separately.

The first archive target is one integration, one regular data-description
configuration, full spectra inside DATA cells, with original row ordering.
An inspector may describe other layouts, but encoding them requires explicit
support. In particular 192 single-channel records/windows can repeat UVW while
192 channels within one DATA cell do not. Casacore may already compress metadata.

Inspection currently reports duplicate rows by full distinguishing row keys
and verifies equality of candidate duplicates across all MAIN columns. It
counts missing triangular pairs relative to the ANTENNA table only when pair
orientation is consistently triangular. It does not infer an observation's
intended baseline inventory or channels absent from spectral-window metadata.

## Canonical schema and reconstruction

```text
visibility[baseline,channel,correlation] : complex64
flags[baseline,channel,correlation]      : bool
baseline_pairs[baseline,2]              : integer
uvw[geometry_record,3]                  : original float64
time[time_record]                       : original dtype
frequency[channel]                     : original dtype
```

Retain an extensible, initially length-one time axis. Partition mixed contexts.
Use compact declarative row/cell reconstruction rules for regular layouts,
otherwise explicit maps or a clear rejection. Preserve distinct time centroids,
phase directions, numerical UVW differences and all original subtable relations.
Frequency-scaled wavelength coordinates are an imaging operation, never a
reconstruction step. Keep row presence independent from FLAG/FLAG_ROW.

Scientific schema and container are separate. Evaluate installed HDF5 for the
single-file prototype; consider Zarr for chunk/object access. Logical groups:
manifest, ms_metadata, subtables, coordinates, reconstruction, visibility_chunks,
exceptions, flags_weights, auxiliary_columns, quality and integrity.
Use portable little-endian dtypes, bounded independently addressable chunks,
immutable codec identifiers, checksums and a complete inventory. No executable
pickle, remote codebook dependency or custom container. Include a reference
decoder and small conformance fixtures.

## Precision, ordering and blocks

Each complex64 sample is 8 bytes. A b-bit codec uses b bits for each component:
2/4/6/8 bits correspond to 4/8/12/16 complex payload bits, ideally 16/8/5.33/4×
smaller before metadata. For exactly 62,000 × 192 × 4 samples, DATA is
380,928,000 bytes. A uint32 permutation of 62,000 baselines costs 248,000 bytes;
four cost 992,000 bytes. A permutation per complex sample would cost
190,464,000 bytes and is outside the initial design.

Keep each baseline spectrum intact. Initially separate products and auto/cross
classes. Compute robust `median(hypot(real,imag))` in float64 over valid channels;
exclude flags/nonfinite values, retain valid counts and exact original samples.
Use deterministic special groups for no-valid-data baselines. Compare original,
baseline-length, stable amplitude and factor-of-two logarithmic bucket order.
Define a separate zero bucket, reference scale and deterministic ties. Persist
actual permutations and bucket boundaries. Compare shared vs per-product order;
frequency-specific order is later and must pay for its additional mappings.

```python
# perm[j] is the original baseline at archived position j.
archived = original[perm, :, product]
restored[perm, :, product] = decoded
```

Reordering only helps by changing scale/codebook sharing or subsequent encoding.
The mandatory control is per-baseline segmented scaling without sorting.
Compare at matched measured error and include full mapping costs.

I/O chunks and quantization blocks are different. First sweep baseline groups
1/8/32/128 and frequency segments 16/32/64/192, one product per block, including
exact partial edges. Never cross configuration or auto/cross boundaries. Begin
with independent float32 scales, no offset and one shared real/imag scale.
Simulate the stored arithmetic/precision when evaluating candidates and reject
nonfinite/overflow reconstructions. Handle zero blocks explicitly.

Later compare a stored complex median or documented mean offset:
`V_hat = mu + scale * (codebook[q_real] + i*codebook[q_imag])`.
Offsets must justify their bytes and shared errors. Do not quantize scales in
version one. Predictor/sky-residual codecs are later and must be self-contained.

## Codec and adaptive allocation

Implement fully specified uniform quantization first, then independently test
symmetric distribution-adapted codebooks, including NF4 if appropriate. No
initial amplitude/phase path. Specify codebook entries, per-component bits,
nearest rounding/ties, scale meaning, decoder arithmetic/order, exception and
clipping rules, byte/bit order, real/imag order, padding and decoded count.
Actually pack 6-bit codes into six bits. Reject truncated/corrupt chunks.

Try 2, 4, 6 and 8 bits plus exact float32 components. Choose the smallest actual
encoded candidate satisfying every configured bound, including scales, indexes,
exceptions, permutations, flagged-value streams, padding and secondary encoding.
Compare isolated RFI exceptions with promotion, splitting and lossless blocks.
Never silently saturate an unflagged value. Preserve special values by original
complex64 bits. A representable ordinary zero needs no exception; preserve a
signed zero when reconstruction would change its bits.

Test stochastic rounding after the deterministic baseline. For adjacent levels
lo/hi, choose hi with `(x-lo)/(hi-lo)` in-range probability; clipping invalidates
the unbiasedness assumption. Record PRNG, deterministic seed and independent
block/component/product streams. Do not optimize seeds for measured error.
Stored stochastic codes need no decoder randomness; subtractive dithering does.
Record recompression as additional lossy provenance.

## Error budgets and science

```text
e = V_hat - V
component_mse = mean(abs(e)**2) / 2
complex_bias = mean(e)
```

Report real/imag error, maximum absolute complex error, high percentiles, and
checks by baseline/product/frequency as well as whole blocks. Evaluate actual
stored-scale reconstruction. Higher bits alone do not guarantee scientific
acceptability. Reject or escape a block whenever any configured check fails.

Only use `component_mse <= epsilon * sigma_th**2` with a separately validated
per-real-component thermal-noise variance and recorded calibration/weight/sigma
conventions. Under independent unbiased additive errors, RMS rises by
`sqrt(1+epsilon)`; a 1% variance budget is an example, not an approved default.
Median amplitude, spectral spread, solar emission, fringes or RFI are not thermal
noise estimates. One integration cannot establish temporal noise. Without a
validated model, report absolute bounds and no thermal-noise guarantee.

Keep weights exact; optional derived effective weights require explicit
assumptions. Conservative archive profiles must record every threshold, including
absolute/leakage constraints for weak Q/U/V and cross-hands independently of
strong I or parallel products. Flagged-data deletion is a separate destructive
profile and is excluded initially.

## Implementation order and acceptance

1. Inspect actual source and declare layout/support boundaries. (Initial sample
   completed; synthetic inspection fixtures supplied.)
2. Build lossless normalized archive and deterministic new-MS reconstruction. (Implemented.)
3. Prove exact supported metadata, UVW, flags, weights, values and special bits. (Fixtures and initial source/export passed.)
4. Add uniform fixed-bit codec, packing, exceptions, float32 escape, no sorting. (Implemented.)
5. Establish the per-baseline segmented scaling control. (Implemented; initial control recorded.)
6. Add amplitude sorting/buckets and grouped scales with invertible permutations. (Implemented.)
7. Add adaptive byte-aware selection and full error/byte accounting. (Implemented.)
8. Compare stochastic/dither and distribution-adapted codebooks. (Stochastic rounding implemented/tested; broader comparisons pending.)
9. Benchmark Dysco where available; record any nonequivalent or unsupported mode.
10. Validate science products and document measured limits before default advice.

Do not prioritize scale quantization, learned predictors, sample-level sorting,
GPU kernels, a custom container or a virtual MS storage manager.

Implemented CLI: `encode`, `verify`, `decode`, `benchmark` and `inspect`.
Configuration must reject unknown options and include
columns, ordering/statistics, group/segment sizes, bits, codebook, rounding/seed,
absolute/noise bounds, exception policy and exact fallback. Do not add placeholder
commands that generate apparently complete archives.

Acceptance requires source-vs-export bitwise comparison of every declared exact
column and metadata, restored row/channel/product order and subtable relations,
visibility errors matching encode-time reports, independently readable chunks,
and a self-contained decoder. Reopen using python-casacore and CASA when available;
run representative imaging/calibration before claiming compatibility. Unavailable
checks must remain distinct from passed checks.

## Benchmark matrix and final report

Mandatory comparisons: exact archive; original-order per-baseline segmented
4/6/8-bit codec; original-order grouped fixed/adaptive codec; per-product amplitude
sort and buckets at identical groups/segments; Dysco with recorded settings.
Compare total bytes versus measured error, plus visibility-only payload ratios.
Count metadata, flags, weights, permutations, exceptions and all other bytes.
Measure secondary compression rather than assuming noise-like codes compress.

Time input reading, statistics, sorting, gathering, quantization, packing,
secondary compression, writing, decoding and MS export separately. Record
hardware, threads, software, peak RAM and cold/warm cache method. Read rows
efficiently before gathering baseline IDs in memory; no sorted-baseline disk
read loop or presorted intermediate MS is needed.

Scientific checks: bias/error distributions by product/baseline/frequency,
cross-product and adjacent-channel covariance, image differences/dynamic range,
calibration solution changes, spectra/narrow features, solar peak flux and weak
polarization where relevant. Use separate multiple-frame data to assess transient
light curves; one frame cannot validate temporal behavior.

Fixtures must cover both channel organizations, missing/duplicate records, mixed
fields, antenna orientation, partial blocks, zeros/all-flagged blocks, RFI,
nonfinite/signed-zero bits, very different product amplitudes and exact escape.
Test every packing width, corruption/truncation and independent block reads when
those components exist.

The final report must answer whether grouping beats per-baseline scaling,
sorting beats coarse buckets, which block sizes work, the 2/4/6/8/lossless
selection frequencies, what dominates bytes and which tested science products
are preserved. No conclusion favoring sorting or four bits is predetermined.
