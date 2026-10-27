# Calibrated-only profile and 15-integration imaging experiment

23 October 2026. The input was already on `inti`; no SSH transfer was needed.
The supplied tar was copied into `data/raw/inti-20261023/` and extracted into
`data/working/lwa-15slot-flat.ms`. Both tar copies and all 108 extracted files
match the supplied SHA-256 manifest. The imported source is never modified.

## Inspected input

The MS contains 931,920 rows: 15 contiguous integrations of 62,128 ordered
antenna pairs. There are **48 channels**, with products XX, XY, YX, YY in that
order. Channel centres span 63.965576–68.463623 MHz with 95.703125 kHz spacing.
Each of DATA, CORRECTED_DATA and MODEL_DATA contains 1,431,429,120 logical bytes.
The complete source occupies 4,428,203,949 bytes. There are 5,280 autocorrelation
rows, no duplicate row keys, and 926,641 distinct original float64 UVW vectors.
UVW remains in metres; its recorded reference is ITRF. TIME/TIME_CENTROID retain
their original float64 bits and UTC reference.

FLAG contains 41,752,612 true elements out of 178,928,640 samples; 175,350 rows
also have FLAG_ROW set. The combined flag mask has the same sample count.
The selected CORRECTED_DATA column contains no nonfinite samples.
HISTORY records CASA `applycal`; imaging explicitly selects CORRECTED_DATA.
This establishes the supplied calibration history, not the accuracy of the
calibration, flux scale, polarization response or a thermal-noise model.

See [inspection](../reports/calibrated-inspection.json),
[summary](../reports/calibrated-inspection.txt) and
[import checksums](../reports/calibrated-import.json). Large quality-statistic
subtables are fully scanned and hashed, with metadata and byte inventories in
the report; their individual cells are not expanded into millions of JSON
records. Archive preservation of those subtables is unaffected.

## Explicit calibrated-only policy

The user authorized this separate destructive profile and explicitly confirmed
removal of MODEL_DATA. Default archive behavior remains unchanged.

`--keep-calib-only` selects CORRECTED_DATA and omits DATA, FLOAT_DATA and MODEL_DATA
when present. Samples selected by `FLAG OR FLAG_ROW` lose their original complex
values and reconstruct as positive complex zero. The array shape and row/channel/
correlation presence are retained. Flags themselves, including FLAG_CATEGORY,
remain exact. Zero placeholders may occupy compressed blocks/indexes; no original
flagged visibility payload is retained. Unflagged nonfinite values and signed
zeros retain the ordinary exact-exception policy.

The option cannot be combined with `--burn-flag`: this profile needs the original
flags to identify discarded values and keep them excluded from imaging.
The exported MS carries `VARQUANMS_CALIBRATED_ONLY` and `VARQUANMS_WEIGHT_POLICY`
keywords. Verification applies this declared projection to the source and checks
all remaining columns, rather than claiming an exact original-MS round trip.

The CORRECTED_DATA-only intermediate works with WSClean when that column is named
explicitly, but it is not the nominal MS v2 data-column organization: the
[MS v2 specification](https://casacore.github.io/casacore-notes/229.html) calls for
DATA, FLOAT_DATA or LAG_DATA. Installed CASA 5.4.1's imager rejected the intermediate
because DATA/FLOAT_DATA was absent. Following the user's explicit choice, decoding
this profile now **adds DATA by default as a bitwise copy of decoded CORRECTED_DATA**.
It never stores or reconstructs the original raw DATA. A `VARQUANMS_DERIVED_DATA`
keyword records the distinction. Compatibility export increases the exported MS
size, not the archive size. Verification checks that derived DATA exactly matches
the decoded calibrated column. `--no-ms-v2-data` explicitly requests the limited
CORRECTED_DATA-only projection; `--ms-v2-data` states the default explicitly.

## Four-bit weights

Within this profile, WEIGHT and optional WEIGHT_SPECTRUM default to `weight_bits=4`.
SIGMA and SIGMA_SPECTRUM remain lossless. Each bounded row chunk first tries an
exact dictionary of up to 16 **float32 bit patterns**. Codes really occupy four
bits each, low nibble first; an odd final count has zero padding in the high
nibble. The float32 codebook is stored with the chunk, so decoding is only a
lookup and does not require a fitting algorithm.

When more than 16 patterns exist, the candidate has 16 uniformly spaced float32
representatives spanning the finite nonnegative values. Choose the nearest
stored representative in float64 arithmetic, resolving equal distances to the
even code. A candidate must preserve zero weights and positive-weight signs and
satisfy, for every positive finite weight,

```text
abs(w_hat - w) / w <= weight_max_relative_error
```

The default bound is 0.01, an engineering setting, not scientific acceptance.
Chunks that fail use exact float32. Nonfinite, negative and signed-zero patterns
are either represented in the exact dictionary or take the exact fallback.
`--weight-max-relative-error 0` requires exact values; `--weight-bits 32` selects
the existing exact array path. Error extrema, lossy/exact/fallback chunk counts,
codebook bytes, packed bytes and container overhead remain auditable.

This input has only three distinct weights, approximately 3.9591246e-8,
3.9799620e-8 and 4.0007993e-8. It therefore needs no weight approximation at four
bits. The full WEIGHT column has 3,727,680 float32 elements (14,910,720 logical
bytes); its code payload is 1,863,840 bytes before secondary compression and
codebook/container overhead. These weights are not assumed to be a validated
thermal-noise estimate.

## Partition bundles

`scripts/timeslot_archive.py` is an experimental wrapper around the existing
single-context archives. It splits contiguous runs on exact bits of all supported
row-context columns, including TIME, centroid, field, DDID, scan, observation,
feeds and interval/exposure. It does not use the antenna pair as the complete
row identity. Each partition independently passes the ordinary inspection and
archive checks. Source ordering is restored by explicit `[start, stop)` row
rules; no time/channel averaging or UVW regeneration occurs.

A bundle contains independently decodable HDF5 archives, verified row ranges,
member SHA-256 values and a checksummed manifest. Every archive repeats the
original subtables, including unused FIELD/OBSERVATION rows and quality tables.
This repetition is included in total bytes. The decoder reconstructs the full MS
using only bundle contents; the source is used later only for verification.
Failed or interrupted outputs retain a `.partial` name. Interleaved contexts
that would require more than 4,096 runs are rejected; this is not arbitrary-MS
partition support. Unknown private casacore keywords still fail; standard
hypercolumn declarations are restored exactly for retained columns.

## Reproduction and comparison policy

```bash
# Single integration:
vamp encode INPUT.ms OUTPUT.h5 --keep-calib-only \
  --config configs/calibrated-power.json --report encode.json
vamp decode OUTPUT.h5 RESTORED.ms
vamp verify OUTPUT.h5 --source INPUT.ms --restored RESTORED.ms --report verify.json

# All 15 integrations, without averaging:
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/timeslot_archive.py encode INPUT.ms NEW.bundle \
  --config configs/calibrated-power.json --report bundle-encode.json
PYTHONPATH=src python scripts/validate_partition_export.py NEW.bundle INPUT.ms \
  RESTORED.ms --report export-verification.json

PYTHONPATH=src python scripts/run_wsclean.py RESTORED.ms --data-column CORRECTED_DATA \
  --output NEW-IMAGES --size 2048 --scale 4amin --threads 1 --niter 10000 \
  --beam-from REFERENCE-IMAGES/image-image.fits
```

The matrix contains calibrated-only lossless visibility storage, uniform-symmetric
quantization, and median/power quantization. All three use identical column,
flagged-value and four-bit weight policies. Quantized configurations use per-product
amplitude ordering, 128 baselines × 48 channels × one correlation, 4/8/12/16-bit
components, max/.995 range candidates, and 24 workers. Absolute bounds in stored
CORRECTED_DATA units are component MSE 1, per-baseline component MSE 2, maximum
complex error 10 and maximum complex bias 1. No thermal-noise budget is asserted.

The original-input reference uses WSClean 3.7 from the pinned installed image,
2048² pixels, 4-arcmin cells, Briggs robust 0, Stokes I, multiscale CLEAN, mgain .8,
automask 5 and threshold 1. It reached the 10,000-iteration cap. Its restoring
beam is 643.149452 × 463.096951 arcsec, PA 43.436339 degrees. Candidate images use
that same beam. All imaging runs operate on private MS copies, with source
checksums checked before and after.

Whole-source size reduction combines the explicitly discarded columns/values
with compression. The lossless calibrated-only control is needed to isolate
the incremental effect of visibility quantization. Equal configured error ceilings
are not equal achieved error. Bundle timings include copying, inspecting and
verifying each partition; encoding-only times are reported separately. The three
codec runs share this host and execute concurrently, so their times are workload
observations, not isolated latency benchmarks or evidence of a 15-second full-MS
conversion. Joint imaging of 15 integrations does not validate light curves.

## Completed measurements

All 15 partitions were encoded, independently reconstructed and checked against
the full input. The original 108 source files still match the supplied manifest.
The lossless projection passes 146 exact column checks; quantized projections
pass 145, and reconstructed CORRECTED_DATA matches the archived decoder output
bitwise. The default-compatible power export additionally verifies DATA against
decoded CORRECTED_DATA. All 90 regression tests pass.

| Profile | Whole bundle bytes | Component MSE | Total encoding time, excluding inspection (s) | Complete bundle, including copy/inspection/verification (s) |
| --- | ---: | ---: | ---: | ---: |
| Calibrated-only lossless visibility | 1,079,628,707 | 0 | 90.29 | 318.19 |
| Uniform symmetric | 634,554,109 | 0.672812 | 132.68 | 369.43 |
| Median/power | 634,658,559 | 0.764378 | 146.76 | 396.51 |

Four-bit weights are exact in all 915 chunks, with no fallback. Their stored
codebooks and codes occupy 16,932 + 49,796 = **66,728 bytes**, before container
overhead. Original SIGMA values and flag arrays remain exact. Four-bit weight
approximation on a wider distribution is therefore not scientifically tested
by this dataset.

The complete calibrated-only lossless image control is pixel-identical to the
original input in all five WSClean products. For the two lossy visibility runs:

| Product / metric | Uniform | Median/power |
| --- | ---: | ---: |
| Dirty relative RMS difference | 0.007089% | 0.007537% |
| Restored relative RMS difference | 5.1252% | 6.8909% |
| Model relative RMS difference | 5.6456% | 6.6382% |
| Residual relative RMS difference | 8.7828% | 12.8460% |
| PSF difference | 0 | 0 |
| Restored difference RMS (Jy/beam) | 0.098868 | 0.132928 |
| Restored maximum absolute difference (Jy/beam) | 13.6999 | 26.6953 |

The large relative differences after capped, automatically masked multiscale
CLEAN must not be hidden by the small dirty-image error. The original-input
reference RMS is 1.929042 Jy/beam in the restored image, not a thermal-noise
measurement. These tests establish neither equivalent deconvolved science
products nor an advantage for the power candidate.

The compatible power export reopens in CASA 5.4.1 and passes a 64×64,
600-arcsec-cell natural-weight Stokes I dirty-image smoke workflow on
CORRECTED_DATA. The installed binary requires
`CASAPATH='/common/casa/casa-release-5.4.1-31.el6 linux'` to find its existing
observatory data. The earlier CORRECTED_DATA-only projection failed the CASA
imager's DATA-column requirement; the derived calibrated DATA copy resolves
that issue. This smoke workflow is not a CASA 6, calibration or flux-scale test.

Evidence: [complete experiment](../reports/calibrated-evaluation.json),
[lossless images](../reports/wsclean-calibrated-lossless.json),
[uniform images](../reports/wsclean-calibrated-uniform.json),
[power images](../reports/wsclean-calibrated-power.json),
[compatible export](../reports/calibrated-power-compatible-export.json), and
[CASA workflow](../reports/casa-calibrated-power.json). The
[paper](../docs/paper/vamp.md) contains equations, the profile diagram,
rate–distortion plot, actual image panels and explicitly clipped diagnostic
difference maps. Full-resolution FITS and archives remain in the ignored
`results/calibrated-20261023/` directory.
