"""Tests for optional region-sector cache archiving and recovery."""
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import region_cache_archive as archive


class ArchiveTests(unittest.TestCase):
    def test_manifest_records_exact_git_objects(self):
        import json
        data = json.loads(archive.MANIFEST.read_text())
        self.assertEqual(data["baseline"], "d362f551357f80125649e583b22b6dbf8aeccb4e")
        self.assertEqual(len(data["entries"]), 6)
        self.assertEqual(sum(x["size"] for x in data["entries"]), 4_074_665)
        self.assertEqual(len({x["path"] for x in data["entries"]}), 6)
        for item in data["entries"]:
            self.assertTrue(item["path"].startswith("paper2/r0/theory/secdec_cache_"))
            self.assertEqual(len(item["git_sha"]), 40)
            self.assertEqual(len(item["sha256"]), 64)

    def test_wrong_blob_digest_or_size_rejected(self):
        payload = b"sector data"
        entry = {"path": "paper2/r0/theory/secdec_cache_test.pkl",
                 "git_sha": "a" * 40, "size": len(payload),
                 "sha256": hashlib.sha256(payload).hexdigest()}
        def fake_git(*args):
            if args[0] == "rev-parse":
                return ("a" * 40 + "\n").encode()
            return payload
        with patch.object(archive, "git", side_effect=fake_git):
            self.assertEqual(archive.inspect_record(entry, "some-commit"), payload)
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                archive.inspect_record({**entry, "sha256": "0" * 64}, "some-commit")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                archive.inspect_record({**entry, "size": 2}, "some-commit")

    def test_missing_cache_is_ignored_by_sector_cache(self):
        # The runtime cache loader must regenerate entries when the cache is absent.
        import sys
        import importlib
        sys.path.insert(0, str(archive.ROOT / "paper2/r0/theory"))
        try:
            cache = importlib.import_module("cache_utils")
            with tempfile.TemporaryDirectory() as d:
                _, entries = cache.sector_cache(Path(d) / "not-yet-created.pkl")
                self.assertEqual(entries, {})
        finally:
            sys.path.remove(str(archive.ROOT / "paper2/r0/theory"))


if __name__ == "__main__":
    unittest.main()
