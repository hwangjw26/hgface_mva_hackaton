#!/usr/bin/env python3
"""Read public scorer metadata with an existing key; no variant query."""
import argparse
import json
from pathlib import Path

from alphagenome.atlas import atlas

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--key-file", required=True, type=Path)
p.add_argument("--output", required=True, type=Path)
a = p.parse_args()
client = atlas.create(a.key_file.read_text().strip(), timeout=30)
metadata = client.scorer_metadata()
result = {name: {"name": value.name, "is_signed": value.is_signed,
                 "tracks": value.track_metadata.reset_index().to_dict(orient="records")}
          for name, value in metadata.items()}
a.output.write_text(json.dumps(result, indent=2, default=str) + "\n")
print(json.dumps({"scorers": list(result), "track_counts": {k:len(v['tracks']) for k,v in result.items()}}))
