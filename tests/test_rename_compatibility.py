"""The software rename must not break archived data or the older entry points."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import h5py
from vamp.archive import encode, verify
from vamp.config import Config
from vamp.container import Reader
from varquanms.container import Reader as LegacyReader
from test_inspection import make_ms


class RenameCompatibilityTests(unittest.TestCase):
    def test_both_module_names_and_commands_remain_available(self):
        import vamp
        import varquanms
        self.assertEqual(vamp.__version__, varquanms.__version__)
        for name in ('vamp', 'varquanms'):
            run=subprocess.run([sys.executable,'-m',name,'--help'],capture_output=True,text=True,check=True)
            self.assertIn('usage: vamp',run.stdout)

    def test_archival_identifier_and_old_reader_remain_compatible(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            source=make_ms(root/'source.ms')
            archive=root/'archive.h5'
            encode(source,archive,Config())
            self.assertEqual(verify(archive,source=source)['status'],'passed')
            with h5py.File(archive,'r') as h5:
                self.assertEqual(h5.attrs['format'],'varquanms.hdf5')
            with Reader(archive) as current, LegacyReader(archive) as legacy:
                self.assertEqual(current.manifest['schema'],'varquanms.archive/1')
                self.assertEqual(legacy.manifest,current.manifest)

    def test_paper_title_and_repository_name(self):
        text=(Path(__file__).resolve().parents[1]/'docs/paper/vamp.md').read_text()
        title='Visibility-Adaptive Mixed Precision for Radio Interferometric Measurement Sets'
        self.assertIn(f"title: '{title}'",text)
        self.assertIn(f'# {title}',text)
        self.assertIn('https://github.com/peijin94/vamp',text)


if __name__=='__main__':
    unittest.main()
