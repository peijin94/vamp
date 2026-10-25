"""The standard experiment keeps independent I/V references and output paths."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
try:
    import stokes_validation as routine
finally:
    sys.path.pop(0)


class StokesRoutineTests(unittest.TestCase):
    def test_standard_paths_keep_both_polarizations_and_repeats(self):
        self.assertEqual(routine.STANDARD_STOKES,('I','V'))
        outputs={routine.image_paths(Path('runs'),'case',pol,repeat)
                 for pol in routine.STANDARD_STOKES for repeat in (False,True)}
        self.assertEqual(len(outputs),4)
        self.assertEqual(routine.reference_path('runs','I'),Path('runs/reference'))
        self.assertEqual(routine.reference_path('runs','V'),Path('runs/reference-v'))
        with self.assertRaises(ValueError):routine.image_paths('runs','case','Q')

    def test_existing_reference_rejects_wrong_stokes_or_changed_source(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);reference=root/'reference-v';reference.mkdir()
            (reference/'image-V-image.fits').touch()
            record=dict(status='complete',source_unchanged=True,source_sha256={'file':'hash'},
                parameters=dict(size=2048,scale='4amin',polarizations='V',data_column='CORRECTED_DATA',
                    clean_profile='calibrated',weighting='briggs 0',multiscale=True,threads=8,
                    niter=150000,mgain=.95,auto_mask=3,auto_threshold=.5,multiscale_scale_bias=.8,
                    local_rms=True,horizon_mask='10deg',reorder=True,memory_percent=40,update_model_required=False))
            path=reference/'complete.json';path.write_text(json.dumps(record))
            with patch.object(routine,'snapshot',return_value={'file':'hash'}):
                result=routine.ensure_references('source.ms',root,('V',),{})
                self.assertEqual(result,{'V':reference})
                record['parameters']['polarizations']='I';path.write_text(json.dumps(record))
                with self.assertRaisesRegex(ValueError,'source or calibrated imaging settings differ'):
                    routine.ensure_references('source.ms',root,('V',),{})
                record['parameters']['polarizations']='V';path.write_text(json.dumps(record))
            with patch.object(routine,'snapshot',return_value={'file':'changed'}):
                with self.assertRaisesRegex(ValueError,'source or calibrated imaging settings differ'):
                    routine.ensure_references('source.ms',root,('V',),{})


if __name__=='__main__':unittest.main()
