"""Verify that cached replay cannot silently accept changed or unavailable evidence."""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest

from avi_replay import check_config


class CacheConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ["TEST_WORKDIR"])
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        inputs = {}
        for name in ("responses", "query", "metadata"):
            data = (name + " synthetic public fixture\n").encode()
            (self.root / name).write_bytes(data)
            inputs[name] = {"path": name, "sha256": hashlib.sha256(data).hexdigest()}
        self.cfg = {"schema_version": 1, "cache_only": True, "inputs": inputs,
                    "cutoff": 0.99379, "operator": ">=", "sdk_version": "0.9.0",
                    "score_scale": "direct AVI_SCORE quantiles layer"}
        self.path = self.root / "config.json"
        self.save()

    def save(self):
        self.path.write_text(json.dumps(self.cfg))

    def test_paths_relative_to_config_not_working_directory(self):
        resolved = check_config(self.path)
        self.assertEqual(resolved["responses"], self.root / "responses")

    def test_changed_response_rejected(self):
        (self.root / "responses").write_text("changed synthetic score\n")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            check_config(self.path)

    def test_missing_query_rejected(self):
        (self.root / "query").unlink()
        with self.assertRaisesRegex(ValueError, "Missing AVI input"):
            check_config(self.path)

    def test_threshold_scale_and_online_mode_changes_rejected(self):
        for key, value in [("cutoff", 0.8), ("score_scale", "raw score"), ("cache_only", False)]:
            with self.subTest(key=key):
                original = self.cfg[key]
                self.cfg[key] = value
                self.save()
                with self.assertRaises(ValueError):
                    check_config(self.path)
                self.cfg[key] = original


if __name__ == "__main__":
    unittest.main()
