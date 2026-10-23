# Default precision ladder: 4/8/12/16 bits

On 2026-10-23 the default adaptive candidate set changed from `[2,4,6,8]` to
`[4,8,12,16]` **per real or imaginary component**. The engineering configuration
and current benchmark matrices use the new ladder. Fixed-bit smoke comparisons
now include 4, 8, 12 and 16. The archive's default mode remains lossless; adaptive
mode still requires explicit numeric error bounds. Exact complex64 blocks remain
available and compete on serialized size with quantized blocks and exceptions.

| Component bits | Complex sample payload | Ideal DATA-only reduction |
| --- | ---: | ---: |
| 4 | 1 byte | 8× |
| 8 | 2 bytes | 4× |
| 12 | 3 bytes | 2.67× |
| 16 | 4 bytes | 2× |

These ideal numbers exclude scales, headers, exceptions, maps, metadata and
secondary compression. Sixteen-bit codes are scaled integers, not float16;
float32 output does not restore quantization losses.

## Implementation and compatibility

Packing uses wide enough unsigned integers to retain all sixteen code bits;
decoding subtracts the codebook midpoint in signed int32. This avoids both
uint8 truncation and signed int16 wraparound. Twelve-bit codes are tightly
packed. The float32 scale and decoder's float32 multiplication are unchanged.
Codebook zero, reserved codes, padding, exceptions and exact escape rules remain.

The `uniform-symmetric-v1` family and block header are unchanged. Existing widths
retain their original bytes and arithmetic results. Readers still support 2/6
bits, and explicit legacy configurations can still reproduce earlier tests;
these widths are not defaults. Older decoders reject 12/16-bit blocks. The
standard-library reference decoder and literal conformance fixtures cover the
new widths independently of NumPy packing.

Tests cover high unsigned codes, negative and positive endpoints, exact byte
lengths, twelve-bit padding, reserved codes, truncated/trailing payloads,
protected special values, MS round trips and adaptive promotion to 12/16 bits
or exact storage under progressively tighter bounds.

## Grouping is unchanged

| Setting | Ordering | Maximum block shape (baselines × channels × correlations) |
| --- | --- | --- |
| Library defaults (if quantization enabled) | Original row order | 1 × 192 × 1 |
| `configs/engineering-adaptive.json` | Original row order | 32 × 64 × 1 |
| Earlier smallest measured archive | Independent amplitude order per product | 128 × 192 × 1 |

Amplitude order uses each baseline/product's median magnitude over valid channels
across the spectrum. Flagged/nonfinite samples are excluded from ordering
statistics. Sorting is stable with original row index as tie-breaker. Auto/cross
and no-valid-sample classes stay separate; blocks stop at those boundaries and
partial edges. Every block has one float32 scale shared by real/imaginary parts
and an independently chosen precision. Frequency segmentation controls which
channels share that scale; baseline spectra remain intact during sorting.

Other implemented choices are baseline length, projected UV orientation,
combined length/orientation bins, either ordered antenna endpoint, logarithmic
amplitude buckets, and class-preserving original order. None is silently selected
by changing precision defaults.

## Historical measurements

The previously reported 136,920,951-byte archive and associated WSClean image
differences used `[2,4,6,8]`. They are not measurements of the new ladder.
Historical reports retain their exact configurations; changing example files
does not retroactively change those results. Scientific tolerances remain
unapproved, and higher precision alone is not a science-fidelity guarantee.

## Full-MS regression with the new ladder

The 128-baseline × 192-channel, per-product amplitude configuration was rerun on
the original 59 MHz MS. Only `bits` changed relative to its earlier settings;
the numeric ceilings, grouping, scale candidates and flag preservation matched.
This is one regression profile, not a new sweep of all grouping methods.

| Ladder | Whole archive bytes | Source/archive ratio |
| --- | ---: | ---: |
| Earlier 2/4/6/8 | 136,920,951 | 2.991× |
| New 4/8/12/16 | 106,730,115 | 3.838× |

The archive is 30,190,836 bytes (22.05%) smaller for this configuration. Local
12/16-bit blocks replace the previous lossless blocks; this is not evidence that
higher-bit payloads are intrinsically smaller. Selection still counts actual
compressed block bytes, exceptions and headers.

| Selected component precision | Blocks | Share of complex samples by block precision |
| --- | ---: | ---: |
| 4 | 38 | 1.957% |
| 8 | 1,556 | 80.144% |
| 12 | 344 | 17.615% |
| 16 | 6 | 0.283% |
| Lossless components | 0 | 0% |

Exact exception values can still occur inside these quantized blocks. Overall
per-component MSE was `2.8116945739376e-8`; maximum complex error was
`0.0007737967887591865`, in input DATA units. All configured checks passed.
UVW/metadata/flags/weights and protected values passed their exact checks; the
exported visibilities matched archive reconstruction. All 65 original source
file hashes remained unchanged. All 55 tests and GitHub Actions passed.

WSClean used the same pinned container as the earlier comparisons: 2048×2048,
4-arcmin pixels, Stokes I, Briggs 0, multiscale, one thread, fixed restoring beam,
auto-mask 5, auto-threshold 1 and 10,000 iterations. Both reference and candidate
hit the iteration cap; convergence is not claimed. The imaging input MS hash
inventory was unchanged after the run, which used a private working copy.

| Image comparison against the original MS | Earlier ladder | New ladder |
| --- | ---: | ---: |
| Dirty-image relative RMS difference | 0.2045% | 0.1996% |
| Restored-image relative RMS difference | 0.6934% | 0.7137% |
| PSF relative RMS difference | 0% | 0% |

The new restored-image peak at the fixed reference position changed by
−0.03976%, with no peak-pixel shift. The fixed radius-four-pixel aperture sum
changed by −0.01177%. These are not identified solar peak/integrated flux
measurements. FITS units are software-labeled Jy/beam; absolute calibration is
not independently established. The comparison covers the finite celestial SIN
footprint and excludes invalid projection corners. Model and residual relative
RMS differences were 4.856% and 1.533%, respectively; nonlinear CLEAN responses
must not be summarized as only the dirty-image error. Q/U/V, calibration
solutions and temporal behavior were not newly tested.

The 22.05% smaller archive has slightly larger restored-image differences under
these controls. This is a rate/distortion observation, not proof that every
science product is preserved or that this grouping is optimal for the new ladder.

Machine-readable evidence:

- [Archive/configuration/quality and imaging summary](../reports/precision-ladder-4-8-12-16.json)
- [Full WSClean comparison](../reports/wsclean-precision-4-8-12-16.json)
- [Restored-image comparison plot](../reports/wsclean-precision-4-8-12-16.png)

Generated archives and exports are under ignored
`results/precision-4-8-12-16-20261023/`; imaging outputs and invocation logs are
under `results/wsclean-precision-4-8-12-16/`. The initial reference is
`results/wsclean-grouping/reference-single/`.
