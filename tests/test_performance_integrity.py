"""Row views, batched reads and cached power levels retain exact semantics."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import unittest

import numpy as np

from vamp.archive import _check_column
from vamp.config import Config, SUPPORTED_BITS
from vamp.container import Reader, Writer
from vamp.inspection import _open, row_window
from vamp.power import reconstruct
from vamp import codec
from test_inspection import make_ms, snapshot


class PerformanceIntegrityTests(unittest.TestCase):
    def test_readonly_window_bounds_context_and_exception_reset(self):
        with tempfile.TemporaryDirectory() as d:
            source=make_ms(Path(d)/'input.ms',pairs=((0,0),(0,1),(1,0),(1,1)))
            before=snapshot(source)
            with _open(source) as t:expected=t.getcol('ANTENNA1')[1:3]
            with self.assertRaisesRegex(RuntimeError,'deliberate'):
                with row_window(source,1,3,('DATA',)):
                    with _open(source) as t:
                        self.assertEqual(t.nrows(),2)
                        self.assertNotIn('DATA',t.colnames())
                        np.testing.assert_array_equal(t.getcol('ANTENNA1'),expected)
                        self.assertFalse(t.iswritable())
                    def rows():
                        with _open(source) as t:return t.nrows()
                    with ThreadPoolExecutor(max_workers=1) as pool:
                        self.assertEqual(pool.submit(rows).result(),4)
                    with self.assertRaisesRegex(ValueError,'Nested'):
                        with row_window(source,0,1):pass
                    raise RuntimeError('deliberate')
            with _open(source) as t:self.assertEqual(t.nrows(),4)
            for lo,hi in ((-1,2),(2,2),(1,5)):
                with self.assertRaisesRegex(ValueError,'window'):
                    with row_window(source,lo,hi):pass
            self.assertEqual(snapshot(source),before)

    def test_chunk_sized_exact_checks_detect_corrupt_tail(self):
        values=np.arange(1100*1024,dtype='<f4').reshape(1100,1024)
        class Source:
            def nrows(self):return len(values)
            def getcol(self,column,startrow,nrow):return values[startrow:startrow+nrow]
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'a.h5';w=Writer(path,Config())
            w.put('x/values',values,'auxiliary_columns')
            w.finish(dict(config=Config().to_dict(),tables=[]))
            spec=dict(encoding='dense',values='x/values',row_count=len(values),cell_shape=[1024],dtype='<f4')
            with Reader(path) as r:_check_column(r,spec,Source(),'x')
            values[-1,-1]+=1
            with Reader(path) as r:
                with self.assertRaisesRegex(ValueError,'reconstruction mismatch'):
                    _check_column(r,spec,Source(),'x')

    def test_every_cached_level_matches_literal_decoder_operations(self):
        for bits in SUPPORTED_BITS:
            q=np.arange(2**bits,dtype=np.uint16)
            maximum=2**bits-1
            for gamma in (1.,.5,.25):
                u=(2*q.astype(np.float64)-maximum)/maximum
                magnitude=np.abs(u)
                if gamma!=1:
                    magnitude=magnitude*magnitude
                    if gamma==.25:magnitude=magnitude*magnitude
                for radius,center in ((0.,-0.),(1.1234567,-2.34567),
                                      (np.finfo(np.float32).tiny,np.nextafter(np.float32(0),np.float32(1)))):
                    radius=float(np.float32(radius));center=float(np.float32(center))
                    expected=(np.copysign(magnitude,u)*radius+center).astype(np.float32)
                    if radius==0:expected.fill(center)
                    self.assertEqual(reconstruct(q,bits,radius,center,gamma).tobytes(),expected.tobytes())

    @unittest.skipIf(codec._native is None,'Optional native encoder unavailable')
    def test_native_search_preserves_every_power_level(self):
        for bits in SUPPORTED_BITS:
            for gamma in (1.,.5,.25):
                q=np.arange(2**bits,dtype=np.uint16)
                levels=reconstruct(q,bits,1.,0.,gamma)
                original=(levels+1j*levels).astype('c8').reshape(1,-1)
                cfg=Config(mode='fixed',codebook='power-midrise-v1',bits=[bits],
                           power_gammas=[gamma],scale_quantiles=[1.],compression='none',
                           max_component_mse=0.,max_complex_error=0.)
                encoded=codec.encode_block(original,np.zeros(original.shape,bool),cfg,cfg.bounds())
                self.assertEqual(encoded.bits,bits)
                self.assertEqual(encoded.decoded.tobytes(),original.tobytes())


if __name__=='__main__':unittest.main()
