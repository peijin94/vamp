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
