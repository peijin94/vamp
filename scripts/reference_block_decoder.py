"""Independent standard-library decoder for uniform and signed-power blocks.

This small conformance implementation returns little-endian complex64 bytes.
The complete HDF5-to-MS decoder is `python -m varquanms decode`.
"""

import argparse
import math
import os
from pathlib import Path
import struct
import zlib


def decode(blob, rows, channels):
    header = struct.Struct("<BBHIIf")
    if len(blob) < header.size:
        raise ValueError("Truncated header")
    bits, compression, mapping, count, exceptions, step = header.unpack_from(blob)
    if mapping not in (0, 1, 2) or count != rows * channels or not 0 < count <= 262144:
        raise ValueError("Invalid sample count/header")
    if bits == 0:
        if mapping or exceptions or step != 0:
            raise ValueError("Invalid exact header")
        expected = count * 8
    elif bits in (2, 4, 6, 8, 12, 16) and 0 <= exceptions <= count and 0 <= step < float("inf"):
        code_bytes = (2 * count * bits + 7) // 8
        expected = code_bytes + exceptions * 12
    else:
        raise ValueError("Unsupported block header")
    offset = header.size
    if mapping == 2:
        if len(blob) < offset+12:
            raise ValueError("Truncated power metadata")
        cr, ci, gamma = struct.unpack_from("<fff", blob, offset)
        if not math.isfinite(cr) or not math.isfinite(ci) or gamma not in (1., .5, .25):
            raise ValueError("Invalid power parameters")
        offset += 12
    body = blob[offset:]
    if compression == 1:
        decoder = zlib.decompressobj()
        body = decoder.decompress(body, expected + 1)
        if not decoder.eof or decoder.unconsumed_tail or decoder.unused_data:
            raise ValueError("Invalid compressed body")
    elif compression != 0:
        raise ValueError("Unknown compression ID")
    if len(body) != expected:
        raise ValueError("Body length mismatch")
    if bits == 0:
        return body
    tail = 2 * count * bits % 8
    if tail and body[code_bytes - 1] >> tail:
        raise ValueError("Nonzero bit padding")
    result = bytearray(count * 8)
    maximum = (2**bits - 1)/2 if mapping else 2**(bits - 1) - 1
    for component in range(count * 2):
        position = component * bits
        byte, shift = divmod(position, 8)
        span = (shift + bits + 7) // 8
        code = (int.from_bytes(body[byte:byte + span], "little") >> shift) & (2**bits - 1)
        if not mapping and code == 2**bits - 1:
            raise ValueError("Reserved code")
        if mapping == 2:
            u = (2*code-(2**bits-1))/(2**bits-1)
            magnitude = abs(u)
            if gamma != 1:
                magnitude = magnitude*magnitude
                if gamma == .25:
                    magnitude = magnitude*magnitude
            residual = math.copysign(magnitude, u)*step
            center = cr if component % 2 == 0 else ci
            value = center if step == 0 else residual+center
        else:
            value = (code-maximum)*step
        struct.pack_into("<f", result, component * 4, value)
    previous = -1
    for index in range(exceptions):
        position = struct.unpack_from("<I", body, code_bytes + index * 4)[0]
        if not previous < position < count:
            raise ValueError("Invalid exception position")
        previous = position
        source = code_bytes + exceptions * 4 + index * 8
        result[position * 8:position * 8 + 8] = body[source:source + 8]
    return bytes(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("block", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--channels", type=int, required=True)
    args = parser.parse_args()
    values = decode(args.block.read_bytes(), args.rows, args.channels)
    with args.output.open("xb") as stream:
        stream.write(values)


if __name__ == "__main__":
    main()
