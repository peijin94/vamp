import struct
from pathlib import Path
import tempfile
import unittest
import zlib

import numpy as np
import h5py

from varquanms.codec import error_metrics, within_bounds
from varquanms.experimental_mx import (FORMATS, HEADER, LEVELS, decode_block, decode_component,
                                      element_values, encode_batch, quantize_elements)


BOUNDS = {"max_component_mse": 1e-7, "max_baseline_component_mse": 2e-7,
          "max_complex_error": .003, "max_complex_bias": .0003, "max_p99_error": None}


class MicroScalingTests(unittest.TestCase):
    def test_ocp_literal_levels_and_midpoint_ties(self):
        np.testing.assert_array_equal(LEVELS[1], [0, .5, 1, 1.5, 2, 3, 4, 6])
        for fid, maximum, smallest in [(1, 6, .5), (2, 7.5, .125), (3, 28, .0625),
                                      (4, 448, 2**-9), (5, 57344, 2**-16), (6, 65504, 2**-24)]:
            self.assertEqual(float(LEVELS[fid][-1]), maximum)
            self.assertEqual(float(LEVELS[fid][1]), smallest)
            codes = np.arange(len(LEVELS[fid]), dtype=np.uint16)
            np.testing.assert_array_equal(quantize_elements(LEVELS[fid], fid), codes)
            midpoints = (LEVELS[fid][:-1].astype(np.float64)+LEVELS[fid][1:])/2
            expected = codes[:-1] + (codes[:-1] & 1)
            np.testing.assert_array_equal(quantize_elements(midpoints, fid), expected)
            negative = quantize_elements(-LEVELS[fid], fid)
            np.testing.assert_array_equal(element_values(negative, fid).view(np.uint32),
                                          (-LEVELS[fid]).view(np.uint32))

    def test_reserved_codes_rejected(self):
        for fid, code in [(4, 127), (5, 124), (6, 31744)]:
            with self.assertRaises(ValueError):
                element_values([code], fid)

    def test_mixed_scales_precisions_partial_edges_and_protected_bits(self):
        rng = np.random.default_rng(2048)
        source = (rng.normal(size=(5, 32, 7))+1j*rng.normal(size=(5, 32, 7))).astype(np.complex64)
        source *= np.logspace(-5, 1, 7)[None, None, :].astype(np.float32)
        flags = np.zeros(source.shape, dtype=bool)
        source.view(np.uint32).flat[:6] = [0x7fc12345, 0x80000000, 0x7f800000, 0xff800000, 0, 0x80000000]
        flags[1, 3, 2] = True
        source[2, 0, 0] = 1e30+1e-30j
        counts = np.array([32, 17, 1, 32, 31])
        for outer in ["unit", "max"]:
            for granularity in ["mega", "micro"]:
                records, _ = encode_batch(source, flags, counts, BOUNDS, outer_mode=outer, granularity=granularity)
                for i, record in enumerate(records):
                    decoded = decode_block(record.blob)
                    actual = source[i, :counts[i]]
                    protected = flags[i, :counts[i]] | ~np.isfinite(actual) | ((actual.real == 0)&np.signbit(actual.real)) | ((actual.imag == 0)&np.signbit(actual.imag))
                    self.assertEqual(actual[protected].tobytes(), decoded[protected].tobytes())
                    self.assertTrue(within_bounds(error_metrics(actual, decoded, ~flags[i, :counts[i]]), BOUNDS))
                    self.assertEqual(sum(record.formats.values()), actual.size*2)
                    self.assertEqual(record.accounting['stored_body']+HEADER.size, len(record.blob))

    def test_zero_budget_forces_exact_and_all_flagged_bits(self):
        rng = np.random.default_rng(17)
        source = (rng.normal(size=(1, 32, 24))+1j*rng.normal(size=(1, 32, 24))).astype(np.complex64)
        for flags, bounds in [(np.zeros(source.shape, bool), dict(BOUNDS, max_component_mse=0)),
                              (np.ones(source.shape, bool), BOUNDS)]:
            encoded, _ = encode_batch(source, flags, [32], bounds)
            self.assertEqual(decode_block(encoded[0].blob).tobytes(), source[0].tobytes())

    def test_independent_records_corruption_and_truncation(self):
        source = np.zeros((2, 32, 3), np.complex64)
        records, _ = encode_batch(source, np.zeros(source.shape, bool), [32, 9], BOUNDS)
        self.assertEqual(decode_block(records[1].blob).shape, (9, 3))
        blob = records[0].blob
        for bad in [blob[:4], blob[:-1], blob+b'\0', bytes([blob[0]^1])+blob[1:], blob[:-1]+bytes([blob[-1]^1])]:
            with self.assertRaises(ValueError):
                decode_block(bad)
        # Malicious length with a fresh CRC must still be rejected before decompression.
        header = list(HEADER.unpack_from(blob))
        header[10] = 1000000000
        header[-1] = 0
        head = HEADER.pack(*header)
        crc = zlib.crc32(blob[HEADER.size:], zlib.crc32(head[:-4]))
        bad = head[:-4]+struct.pack('<I',crc)+blob[HEADER.size:]
        with self.assertRaises(ValueError):
            decode_block(bad)

    def test_file_decoder_restores_product_and_baseline_order_and_rejects_missing_records(self):
        rng=np.random.default_rng(882)
        source=(rng.normal(size=(3,5,2))+1j*rng.normal(size=(3,5,2))).astype(np.complex64)
        permutations=[np.array([2,0,1],dtype='<u4'),np.array([1,2,0],dtype='<u4')]
        blobs=[]
        expected=np.empty_like(source)
        for p in range(2):
            padded=np.zeros((1,32,5),np.complex64)
            padded[0,:3]=source[permutations[p],:,p]
            records,_=encode_batch(padded,np.zeros_like(padded,dtype=bool),[3],BOUNDS)
            blobs.append(records[0].blob)
            expected[permutations[p],:,p]=decode_block(blobs[-1])
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'fixture.h5'
            with h5py.File(path,'w') as h5:
                h5.attrs.update(format='varquanms.experimental-mx-component/1',complete=True,
                                full_observation=True,shape=source.shape,variant='max-micro')
                for p in range(2):h5.create_dataset(f'permutations/{p}',data=permutations[p])
                h5.create_dataset('payload',data=np.frombuffer(b''.join(blobs),np.uint8))
                h5.create_dataset('offsets',data=np.asarray([0,len(blobs[0]),sum(map(len,blobs))],dtype='<u8'))
                h5.create_dataset('records',data=np.array([[0,0,3,0,5],[1,0,3,0,5]],dtype='<u4'))
            self.assertEqual(decode_component(path).tobytes(),expected.tobytes())
            with h5py.File(path,'r+') as h5:h5['records'][1,0]=0
            with self.assertRaisesRegex(ValueError,'Overlapping'):
                decode_component(path)


if __name__ == '__main__':
    unittest.main()
