"""Run all public synthetic tests; requires TEST_WORKDIR and BCFTOOLS."""
from pathlib import Path
import os
import shutil
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__":
    for key in ("TEST_WORKDIR", "BCFTOOLS"):
        if not os.environ.get(key):
            raise SystemExit(f"Set {key} before running tests")
    if not shutil.which(os.environ["BCFTOOLS"]):
        raise SystemExit("BCFTOOLS must identify an existing executable")
    sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts/avi")]
    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)
