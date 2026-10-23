# Full-range unsigned quantization

23 October 2026. Continue the scaled unsigned-integer route; retain MX only as
the documented comparison experiment. The default codebook remains
`uniform-symmetric-v1`. The optional new codec is `uniform-midrise-v1`.

## Formula correction

In the suggested expression, `M=2**b-1` and `q=0,...,2*M-1` give `2**(b+1)-2`
possible codes. They require **b+1 total component bits**, not b: 30 codes need
5 bits. The old codec already packs unsigned integers and uses `2**b-1` distinct
reconstruction values. It leaves one unused code, not an unused sign bit.

For exactly b bits, use `N=2**b`, `C=(N-1)/2`, `delta=ceil_float32(t/C)`, and
`q=rint_even(x/delta+C)` in `0,...,N-1`. Decode by float32 subtraction of C and
float32 multiplication by the stored step. `ceil_float32` means rounding upward
to a representable float32, **not rounding to an integer**. For example, a
cutoff of 0.01 at 4 bits needs a step near 0.001333, not 1.

At fixed cutoff, the step ratio is `(N-2)/(N-1)`. The high-resolution uniform
error approximation predicts MSE reductions of 12.89%, 0.783%, 0.0488% and
0.00305% at 4/8/12/16 bits. Actual structured visibility errors can improve or
worsen; scale rounding and exception choices are evaluated explicitly.

The full range has no zero level. Nearest-even maps ordinary positive zero to
`+delta/2`. This is authorized ordinary-zero loss, not an exception for every
zero. All-zero blocks still use delta=0 and preserve positive zero. Negative-zero
samples, source-flagged values and nonfinite bit patterns remain exact. In a
zero-rich block, positive zero errors can accumulate into a coherent bias.
Stochastic rounding is supported with independent component streams, but changes
zero bias into zero-component variance; it is not used in these measurements.

## Implementation and compatibility

The existing 16-byte block header's uint16 mapping field is 0 for the original
codec, 1 for midrise. No additional header or container is needed. Exact blocks
always use mapping 0, under either manifest codebook. The manifest and quantized
block mapping must agree. Unknown mappings fail. Old archive bytes retain their
original interpretation; old decoders reject the new mapping explicitly.

Both NumPy and native encoders support the new mapping, including percentile
checks, all packing widths, strict actual-reconstruction checks, exceptions and
actual-byte selection against exact float32. A stale native extension falls back
to NumPy for midrise. Rebuild with `python setup.py build_ext --inplace` when
running from a checkout. The independent standard-library block decoder supports
both mappings. Detailed operations are in [the format contract](../docs/format.md).

The 73-test suite passes, including literal endpoints/all levels at 2/4/6/8/12/16
bits, ties, zero-bias rejection, negative zeros, NaN payloads, flags, RFI, overflow,
truncation, unknown mapping, native/reference agreement, old-codec regression,
stochastic reproducibility and independent MS export/source preservation.

## Matched full-MS measurement

Use the existing inspected 62,128 × 192 × 4 source. Both variants have independent
amplitude order per product, 128 × 192 × 1 quantization blocks, 24 workers,
4/8/12/16 candidates, max/.995 scale candidates, gzip level 1, and unchanged
engineering error ceilings: component MSE 1e-7, per-baseline MSE 2e-7, maximum
complex error .003 and complex bias .0003. Flags are preserved. No noise model
or science-approved tolerance is asserted.

| Measurement | Zero-preserving symmetric | Full-range midrise |
| --- | ---: | ---: |
| First-run whole archive bytes | 106,731,246 | 106,732,004 |
| Serialized visibility block bytes | 97,058,796 | 97,081,093 |
| Aggregate component MSE | 2.811694574e-8 | 2.805425666e-8 |
| 4 / 8 / 12 / 16 / exact block counts | 38 / 1556 / 344 / 6 / 0 | 40 / 1557 / 341 / 4 / 2 |
| Fresh-process encode seconds | 12.997 / 12.945 / 13.295 | 13.069 / 13.039 / 12.210 |
| Median encode seconds | 12.997 | 13.039 |

Both have 1,944 blocks and all six runs pass independent source verification.
The new codec reduces aggregate component MSE by about 0.223%, but **does not
reduce total bytes**. Its serialized visibility blocks are 22,297 bytes larger;
HDF5 chunk allocation makes physical payload size unchanged, while manifest and
quality storage account for the small archive difference. Repeat file lengths
differ slightly because timing/provenance JSON is stored; scientific dataset
hashes are checked for identical repetition.

All timing runs use sequential fresh processes and a previously read source on
the shared host, with no cache eviction or affinity. Imaging does not overlap
timed encoding. Timings include startup, inspection, exact columns, quantization,
error accounting, writing and integrity publication. Verification, export and
imaging are separate. These runs continue to meet the 15-second target; no
meaningful speed difference between codebooks is established.

The source has many exactly zero parallel-hand autocorrelation imaginary
components. With midrise, full-product XX/YY imaginary biases are 1.1083e-6 and
9.1734e-7, compared with -1.4032e-7 and 1.9476e-9 for the symmetric codec. All
configured block and whole-product limits still pass. Exact fallback is selected
for two midrise blocks. These results support an opt-in experiment, not replacing
the zero-preserving default.
Restricted to autocorrelation rows, the new XX/YY imaginary biases are
1.8288e-4 and 1.5567e-4; the whole-product means dilute this effect because only
352 of 62,128 rows are autocorrelations.

## Imaging and reconstruction

The new archive was exported and verified against the source. All 132 supported
non-DATA column entries and metadata pass exact checks, including UVW, flags,
weights, coordinates and original ordering. DATA errors reproduce encode-time
statistics. All 65 source file SHA-256 hashes are unchanged. The exported MS
reopens in python-casacore and is read successfully by WSClean; a separate new
CASA application run was not performed for this codec.

The pinned WSClean container uses 2048² pixels, 4 arcmin/pixel, Briggs 0,
multiscale CLEAN, Stokes I, one worker, and the reference restoring beam. CLEAN
reaches the 10,000-component limit; convergence is not claimed. Compared with
the original-MS single-thread reference, relative RMS differences over the valid
celestial footprint are 0.198534% dirty, 1.114193% restored, 10.726310% model and
2.477910% residual; the PSF difference is zero. Imaging took 258.65 s, separately
from the archive timing, and did not modify its input MS.

The prior zero-preserving reconstruction with the same measured visibility
errors/settings gives a 0.7137% restored-image difference. Thus the small decrease
in midrise visibility MSE does **not** establish an image-domain advantage. CLEAN
is nonlinear and the paths/auto-masks can change. WSClean uses the 61,776 cross
rows here; this image test does not validate autocorrelation imaginary offsets.
Weak polarization, calibration solutions, calibrated solar flux, narrow-feature
and temporal fidelity remain untested. No scientific acceptance threshold was
defined for these image metrics.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 visarchive benchmark INPUT.ms \
  --matrix configs/benchmark-midrise.json --output results/new-uint-comparison
visarchive encode INPUT.ms NEW.h5 --config configs/performance-midrise.json
visarchive decode NEW.h5 RESTORED.ms
visarchive verify NEW.h5 --source INPUT.ms --restored RESTORED.ms
```

The [paper](../docs/paper/varquantMS.md) contains equations and level/benchmark
figures. [The machine-readable report](../reports/uint-midrise-evaluation.json)
records full byte accounting, product errors, class-specific errors, timings,
source hashes, export verification and imaging evidence. Figure-only regeneration
requires no source MS: `python scripts/report_midrise.py --figures-only`.
