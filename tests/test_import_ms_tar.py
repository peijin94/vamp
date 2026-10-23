import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location('import_ms_tar', Path(__file__).parents[1] / 'scripts/import_ms_tar.py')
IMPORTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPORTER)


class ImportTests(unittest.TestCase):
    def fixture(self, root, name='fixture.ms/table.dat', wrong_hash=False):
        source = root / 'supplied.tar'
        data = b'exact fixture payload\x00\xff'
        with tarfile.open(source, 'w') as archive:
            item = tarfile.TarInfo(name); item.size = len(data)
            archive.addfile(item, io.BytesIO(data))
        sha = IMPORTER.digest(source)
        manifest = {'root_directory': 'fixture.ms', 'archive_sha256': sha, 'archive_bytes': source.stat().st_size,
                    'file_count': 1, 'files': {'table.dat': {'bytes': len(data),
                    'sha256': '0'*64 if wrong_hash else hashlib.sha256(data).hexdigest()}}}
        source.with_name(source.name+'.manifest.json').write_text(json.dumps(manifest))
        source.with_name(source.name+'.sha256').write_text(sha+' supplied.tar\n')
        return source, data

    def run_import(self, root, source):
        with patch.object(sys, 'argv', ['import', str(source), '--raw', str(root/'raw'),
                          '--working', str(root/'working'), '--report', str(root/'report.json')]):
            IMPORTER.main()

    def test_copy_and_extraction_match_supplied_bytes(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); source, data = self.fixture(root)
            self.run_import(root, source)
            self.assertEqual((root/'working/fixture.ms/table.dat').read_bytes(), data)
            self.assertEqual(IMPORTER.digest(source), IMPORTER.digest(root/'raw/supplied.tar'))
            with self.assertRaisesRegex(ValueError, 'new'):
                self.run_import(root, source)

    def test_unsafe_member_and_bad_payload_do_not_publish_ms(self):
        for name, bad in [('fixture.ms/../../escape', False), ('fixture.ms/table.dat', True)]:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as d:
                root = Path(d); source, _ = self.fixture(root, name, bad)
                with self.assertRaises(ValueError):
                    self.run_import(root, source)
                self.assertFalse((root/'working/fixture.ms').exists())
                self.assertFalse((root/'report.json').exists())


if __name__ == '__main__':
    unittest.main()
