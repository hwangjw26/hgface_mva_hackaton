"""Portability checks for the reorganized command-line entry points."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class LayoutTests(unittest.TestCase):
    def test_entry_points_outside_checkout(self):
        with tempfile.TemporaryDirectory(dir=os.environ["TEST_WORKDIR"]) as tmp:
            for name in ("reproduce", "verify_reproduction", "build_submission",
                         "prepare_native", "prepare_metadata"):
                result = subprocess.run(
                    [sys.executable, "-B", str(ROOT / "scripts" / (name + ".py")), "--help"],
                    cwd=tmp, capture_output=True, text=True, check=True)
                self.assertIn("usage:", result.stdout)

    def test_generated_commands_resolve_bundled_resources(self):
        with tempfile.TemporaryDirectory(dir=os.environ["TEST_WORKDIR"]) as tmp:
            native = Path(tmp) / "native"
            subprocess.run([sys.executable, "-B", str(ROOT / "scripts/prepare_native.py"),
                            "--config", str(ROOT / "examples/native_paths.json"),
                            "--output", str(native)], cwd=tmp, check=True, capture_output=True)
            commands = json.loads((native / "commands.json").read_text())
            for cmd in commands[:2]:
                for flag in ("--analysis", "--output"):
                    path = Path(cmd[cmd.index(flag) + 1])
                    self.assertTrue(path.is_file(), str(path))
                    self.assertEqual(path.parent, ROOT / "config")
            pairs = Path(tmp) / "pairs.tsv"
            pairs.write_text("gene_id\tallele_1_key\tallele_2_key\nNCBIGene:1\t1:1:A:G\t1:2:C:T\n")
            metadata = Path(tmp) / "metadata"
            subprocess.run([sys.executable, "-B", str(ROOT / "scripts/prepare_metadata.py"),
                            "--native-commands", str(native / "commands.json"),
                            "--pairs", str(pairs), "--sample-id", "SYNTHETIC",
                            "--java", "/tools/java", "--javac", "/tools/javac",
                            "--jar", "/tools/LIRICAL.jar", "--output", str(metadata)],
                           cwd=tmp, check=True, capture_output=True)
            compile_call = json.loads((metadata / "commands.json").read_text())[0]
            self.assertEqual(Path(compile_call[-1]), ROOT / "native/FullMetadata.java")
            self.assertTrue(Path(compile_call[-1]).is_file())
