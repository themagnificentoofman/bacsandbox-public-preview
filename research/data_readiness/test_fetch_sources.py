"""Synthetic transfer fixtures only; no measurement data or live HTTP calls."""
import hashlib
import io
from pathlib import Path
import tempfile
import unittest

from research.data_readiness.fetch_sources import materialize, verify_bytes


class SourceTransferTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'source.xlsx'
        self.data = b'synthetic transfer fixture'
        self.entry = {'key': 'fixture/source.xlsx', 'size': len(self.data),
                      'sha256': hashlib.sha256(self.data).hexdigest(),
                      'md5': hashlib.md5(self.data).hexdigest(), 'url': 'https://example.invalid/fixture'}

    def test_all_integrity_fields_are_required(self):
        verify_bytes(self.data, self.entry)
        for field, value in [('size', 0), ('sha256', 'bad'), ('md5', 'bad')]:
            with self.subTest(field=field), self.assertRaises(ValueError):
                verify_bytes(self.data, self.entry | {field: value})

    def test_existing_source_needs_no_network(self):
        self.path.write_bytes(self.data)
        def no_network(*args, **kwargs):
            self.fail('existing valid source must not be downloaded again')
        self.assertEqual(materialize(self.entry, self.path, opener=no_network), 'verified_existing')

    def test_differing_existing_file_is_preserved(self):
        self.path.write_bytes(b'preserve unrelated bytes')
        with self.assertRaises(ValueError):
            materialize(self.entry, self.path)
        self.assertEqual(self.path.read_bytes(), b'preserve unrelated bytes')

    def test_invalid_download_is_not_installed(self):
        with self.assertRaises(ValueError):
            materialize(self.entry, self.path, opener=lambda *a, **k: io.BytesIO(b'invalid'))
        self.assertEqual(list(self.path.parent.iterdir()), [])
        with self.assertRaises(FileNotFoundError):
            materialize(self.entry, self.path, verify_only=True)

    def test_verified_download_is_installed_without_temporary_leftovers(self):
        self.assertEqual(materialize(self.entry, self.path, opener=lambda *a, **k: io.BytesIO(self.data)), 'downloaded_verified')
        self.assertEqual(self.path.read_bytes(), self.data)
        self.assertEqual(list(self.path.parent.iterdir()), [self.path])


if __name__ == '__main__':
    unittest.main()
