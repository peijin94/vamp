# Median-centered signed power-law quantization

23 October 2026. Implemented as the optional `power-midrise-v1` codec. The default
remains `uniform-symmetric-v1`; nothing is silently converted to the new mapping.

## Mathematical validation

For a complex quantization block, use two component medians:
`C_mid = median(real V) + i median(imag V)`. A median of amplitudes would not
supply the signed real/imaginary offsets needed here. Exclude source-flagged,
nonfinite and negative-zero-protected samples; calculate in float64 and store
the two medians as float32 before evaluating any candidate residuals.

Let `r=x-c` for one component, and t be an upward-rounded float32 range candidate
derived from the maximum or configured quantile of the largest component
residual in each eligible complex sample. The proposed idea becomes

```text
M = 2**b - 1
y = sign(r) * abs(r/t)**gamma
q_continuous = M/2 * (1+y)
u_q = (2*q-M)/M
x_hat(q) = c + t*sign(u_q)*abs(u_q)**(1/gamma)
q = 0,...,M
```

This uses exactly b unsigned bits per component. **Gamma below one increases
level density near the median.** Gamma one is uniform; gamma above one would
decrease central density. This parameter does not measure statistical skewness:
the residual codebook remains symmetric. Median subtraction can remove a
coherent offset, but is not a noise estimate or a prediction of scientific
fidelity.

The initial exact gamma options are 1, 1/2 and 1/4; [1, 1/2] is the candidate
default. Their inverse powers use fixed binary64 products, followed by scaling,
addition of the stored component median, and one float32 rounding. No generic
pow or downloaded table/model is required. Unsupported powers, including .75,
fail explicitly in this version. Extension to other powers requires a new
specified reproducible decoder, not silently approximating a requested gamma.

The encoder evaluates the two bracketing codes and chooses the nearest decoded
value in the **original visibility domain**, with even-code ties. Rounding the
continuous transformed coordinate alone does not minimize original-domain
squared error. This distinction matters for nonlinear levels. A gamma candidate
is accepted only after reconstruction with the stored parameters and all existing
error checks. Actual encoded bytes decide among passing candidates; exact
complex64 always competes. There is no distribution-based heuristic assigning
bits or gamma without evaluating the data.

This is companding, not a new quantization principle. Relevant primary references:
[AudioLabs quantization](https://www.audiolabs-erlangen.de/resources/MIR/FMP/C2/C2S2_DigitalSignalQuantization.html)
and [Offringa's nonlinear visibility quantization](https://arxiv.org/html/1609.02019v1).
These sources do not establish fidelity for this observation or these settings.

## Archival contract

The new mapping uses the existing header with mapping ID 2 and stores t in its
last float32 field. A 12-byte extension stores `center_real, center_imag, gamma`,
making a 28-byte header. M is derived from the stored bit width, so no separate
integer M is necessary. At 1,944 quantized blocks this adds 23,328 header bytes.
The container stays HDF5; blocks remain independently readable and checksummed.

The range and gamma are shared between real/imaginary components, while their
medians are separate. Each 128×192×1 block in the measured configuration gets its
own parameters. All configured bits, scale quantiles and gamma choices compete;
none are silently clipped. Exact exceptions still protect flags, negative-zero
samples, nonfinite payloads and outliers. Constant blocks use t=0 and reconstruct
their centers exactly. With nonzero t, there is no center level in ideal
arithmetic: the closest residual magnitudes are `t/M**(1/gamma)`; float32 rounding
can collapse nearby levels. Ordinary positive-zero loss remains explicit.

Mapping 0/1 bytes retain their original meaning. Old decoders reject mapping 2.
Both the portable NumPy implementation and native encoder agree with the
independent standard-library block decoder. Nearest reconstructed-value rounding
is supported; stochastic rounding for this mapping is explicitly rejected.

All 79 tests pass: literal fixtures, every packed width, monotonic levels,
original-domain nearest-level selection, stored-median precision, constants,
special values, RFI, exact fallback, invalid parameters, truncation, native versus
portable selection (including near-boundary constraints), old codecs and complete
MS export verification. Full-source scientific conclusions remain separate.

## Experiment design

The same inspected 62,128-row, 192-channel, four-product MS is used throughout.
Preserve flags, weights, coordinates and all supported auxiliary columns exactly.
Use amplitude ordering per product, 128×192×1 blocks, 4/8/12/16 bits, maximum/.995
range candidates, gzip level 1, and 24 workers. Original engineering ceilings are
component MSE 1e-7, per-baseline MSE 2e-7, maximum complex error .003 and complex
bias .0003, with no validated noise model.

Separate zero-center/uniform, median-center/uniform, zero-center/adaptive-gamma
and median-center/adaptive-gamma controls. The zero/uniform new-format control
also isolates its header and decoder arithmetic from the older midrise codec.
Test the original symmetric codec at the original bounds and at 1.5e-7/3e-7
and 1.4e-7/2.8e-7 MSE ceilings to compare achieved rate/distortion. Other bounds
and preservation policies stay fixed. Equal configured bounds alone are not an
equal-error comparison.

The original combined implementation measured 15.98/15.67/15.56 seconds in fresh
processes. Native profiling attributed about 71 thread-CPU seconds to candidate
quantization and 88 to secondary compression. Bounding a candidate's compressed
output to the size that could still beat the current winner reduced compression
work to about 62 thread-CPU seconds without changing selected bytes. A full output
buffer makes that candidate a loser; it does not publish a truncated stream.
This optimization preserves complete candidate selection and is tested against
the unbounded portable reference.

The final three process times are **19.852, 14.404 and 14.436 s**, median **14.436 s**.
The first run exceeded the 15-second target; the result is not a guarantee that
every run finishes within 15 seconds. Its archive finalization took about 2.55 s
versus approximately .5 s in the initial run; it also had a larger process-level
remainder. The host/cache/filesystem are shared, with no isolation or cache
eviction, so causes beyond the measured regions are not established. No imaging
overlaps any timed encode. Reports retain all repetitions, including the outlier.

## Measured byte/error comparison

| Variant | Whole archive MB | Aggregate component MSE |
| --- | ---: | ---: |
| Original symmetric, original bounds | 106.731363 | 2.811694574e-8 |
| New header, zero center, gamma 1 | 106.732296 | 2.805425551e-8 |
| Median center, gamma 1 | 106.733067 | 2.782965100e-8 |
| Zero center, adaptive gamma | 104.656037 | 4.343666348e-8 |
| Median center, adaptive gamma | 104.656745 | 4.336359607e-8 |
| Symmetric, MSE ceilings 1.5e-7 / 3e-7 | 101.484000 | 4.379775723e-8 |
| Symmetric, MSE ceilings 1.4e-7 / 2.8e-7 | 102.533698 | 4.150311716e-8 |

All 12 full-MS encodes pass independent source verification. The median-plus-power
profile selects 40/1,552/348/4 blocks at 4/8/12/16 bits and no exact blocks on this
source; exact escape remains tested and available. It selects gamma 1 in 1,449
blocks and gamma 1/2 in 495, out of 1,944 blocks. Gamma 1/4 is supported and tested,
but excluded from this full-MS candidate grid. Repetition and compression-search
acceleration preserve every stored scientific dataset hash.

At identical configured ceilings, the combined profile is about 1.94% smaller
than the original symmetric archive, with about 54.23% higher aggregate MSE.
Median-only centering lowers MSE slightly but does not reduce whole-file bytes.
The extra gamma flexibility mostly improves serialized code compressibility;
it does not produce a large shift toward fewer stored bits on this source.
Offsets can help individual coherent/autocorrelation blocks without dominating
whole-observation storage.

The 1.4e-7 uniform control is smaller **and** has lower aggregate MSE than the
combined power profile. This does not prove uniform is superior for every
constraint: that control has different block/per-baseline ceilings, and product,
maximum-error and imaging comparisons need not follow the aggregate-MSE ranking.
It does rule out claiming an unconditional aggregate rate–distortion improvement
from the smaller file at the original ceilings. The codec remains optional.
For the cross-correlation rows, the combined profile has XY/YX component MSEs
of approximately 5.1143e-8 and 5.1202e-8. The original symmetric whole-product
values are approximately 2.6603e-8 and 2.6034e-8. Weak-product accuracy therefore
must be assessed per product; clustering code levels near a block median is not
by itself a weak-polarization preservation guarantee. Autocorrelation-only XX/YY
imaginary biases remain about 1.0394e-4 / 5.7905e-5 in the new reconstruction.

## Reconstruction and imaging

The combined archive exports a new MS and passes independent reconstruction
checks. All 132 supported non-DATA column entries and metadata are exact;
quantized DATA equals the archive reconstruction and its measured encode-time
error statistics. Source and exported row/channel/correlation order are retained.
All 65 original source file SHA-256 hashes are unchanged. The exported MS reopens
with python-casacore and is consumed by the pinned WSClean container. A separate
new CASA application check was not run for this codec.

WSClean uses 2048² pixels, 4 arcmin/pixel, Briggs 0, multiscale CLEAN, Stokes I,
one worker and the original reference restoring beam. Relative RMS differences
over the valid celestial footprint are **0.231471% dirty / 1.147691% restored**,
with 10.956040% model and 2.570448% residual differences; the PSF is unchanged.
The run takes 248.01 s separately from encoding and reaches the 10,000-component
cap, so convergence is not claimed. It leaves its input MS unchanged.

The earlier strict symmetric profile gives a 0.7137% restored-image difference.
The power codec has no demonstrated imaging advantage here. The relaxed uniform
controls were not imaged; lower aggregate visibility MSE is not a guarantee of
better nonlinear CLEAN agreement. This is a cross-correlation imaging check,
not acceptance of autocorrelation offsets, weak polarization, calibration
solutions, calibrated solar flux, narrow spectral features or temporal behavior.

## Reproduction

```bash
visarchive encode INPUT.ms NEW.h5 --config configs/power-adaptive.json
visarchive decode NEW.h5 RESTORED.ms
visarchive verify NEW.h5 --source INPUT.ms --restored RESTORED.ms
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 visarchive benchmark INPUT.ms \
  --matrix configs/benchmark-power.json --output results/new-power-controls
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 visarchive benchmark INPUT.ms \
  --matrix configs/benchmark-power-final.json --output results/new-power-timings
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 visarchive benchmark INPUT.ms \
  --matrix configs/benchmark-power-fidelity.json --output results/new-power-fidelity
python scripts/report_power.py --figures-only
```

Detailed results: [machine-readable study](../reports/power-evaluation.json),
[image comparison](../reports/wsclean-power.json), [format](../docs/format.md),
and [paper](../docs/paper/varquantMS.md). Large archives, MSs and FITS products
remain local under `results/power-20261023/` and are excluded from Git.
