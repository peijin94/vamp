# Default component precision: 3, 4, 6, 8, 12, 16 bits

The default adaptive candidate ladder now includes **3/4/6/8/12/16 bits per real
component**, as requested on 24 October 2026. Both components of a complex sample
use the selected width. Exact complex64 remains a competing fallback; the default
encoding mode remains lossless until a lossy profile and error limits are supplied.

Three-bit codes are packed, not stored in eight-bit containers. Eight component
codes occupy three bytes; for example codes `[0,1,2,3,4,5,6,7]` produce hexadecimal
`88 c6 fa`. Partial final bytes have zero padding, checked by the decoder. A block
of N complex samples occupies `ceil(6*N/8)` code bytes before headers, exceptions
and secondary compression. This is an ideal 9.375% code payload relative to
complex64, before those costs and any lossless compression savings.

The Python encoder/decoder, native uniform and power encoders, and independent
reference decoder support three bits. Existing codec IDs and decoder arithmetic
are unchanged. Two-bit archives remain readable and can be explicitly requested;
two bits are not in the new default ladder. Older decoder implementations reject
unsupported widths. An older native extension triggers the portable encoder
when three-bit candidates are requested.

Validation covers literal packed bytes, every tail-byte alignment, invalid
padding, protected values, archive reconstruction, native/portable equivalence
and the independent decoder. **103 tests passed** locally after rebuilding the
native extension. Historical experiment configurations and their measured
results retain the ladders actually used. New 10%-payload experiments explicitly
use the new six-width ladder; their image quality must be measured independently.
