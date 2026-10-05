"""Locate shared modules independently of the caller's working directory."""
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parents[1]

def run(module):
    sys.path.insert(0, str(ROOT / "src"))
    runpy.run_module(module, run_name="__main__")
