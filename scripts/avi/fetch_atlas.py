#!/usr/bin/env python3
"""Fetch only authorized minimal SNV requests; keep full score responses private."""
import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
from datetime import datetime, timezone

import grpc
from alphagenome.atlas import atlas
from alphagenome.data import genome

SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"]
FIELDS = ["chromosome", "position_1based", "reference_bases", "alternate_bases"]


def variant_key(variant):
    return f"{variant.chromosome.removeprefix('chr')}:{variant.position}:{variant.reference_bases}:{variant.alternate_bases}"


def serialize_response(result, expected):
    if set(result) != set(SCORERS):
        raise ValueError("Returned scorer set mismatch")
    serialized = {}
    for name, data in result.items():
        if list(data.obs.columns) != ["variant"] or len(data.obs) != 1:
            raise ValueError("Unexpected scorer observation unit")
        returned = variant_key(data.obs.iloc[0]["variant"])
        if returned != expected:
            raise ValueError("Returned exact allele mismatch")
        names = list(data.var["name"])
        expected_columns = 1 if name == "AVI_SCORE" else 18
        if data.shape != (1, expected_columns) or len(set(names)) != expected_columns:
            raise ValueError("Unexpected score matrix/feature identities")
        serialized[name] = {"returned_allele_key": returned, "track_names": names,
                            "values": data.X.tolist(),
                            "layers": {key: value.tolist() for key, value in data.layers.items()}}
    json.dumps(serialized, allow_nan=False)
    return serialized


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--key-file", required=True, type=Path)
    parser.add_argument("--query", type=Path)
    parser.add_argument("--query-sha256")
    parser.add_argument("--public-control", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    os.umask(0o077)
    if importlib.metadata.version("alphagenome") != "0.9.0":
        raise ValueError("This adapter is pinned to SDK0.9.0")
    if args.public_control:
        if args.query or args.query_sha256:
            parser.error("Public-control mode does not read patient requests")
        # Public UCSC import control, documented before this analysis.
        variants = [genome.Variant("chr21", 5010001, "G", "A")]
        digest = None
    else:
        if not args.query or not args.query_sha256:
            parser.error("Patient mode requires the approved minimal query file and its pinned SHA256")
        query_bytes = args.query.read_bytes()
        digest = hashlib.sha256(query_bytes).hexdigest()
        if digest != args.query_sha256:
            raise ValueError("Approved query file changed")
        with io.StringIO(query_bytes.decode("utf-8")) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FIELDS:
                raise ValueError("Only the four authorized request fields are allowed")
            rows = list(reader)
        variants = [genome.Variant(x["chromosome"], int(x["position_1based"]), x["reference_bases"], x["alternate_bases"]) for x in rows]
        if len(variants) != 201 or len({variant_key(v) for v in variants}) != 201:
            raise ValueError("Authorized request must contain201 distinct SNVs")
        if any(len(v.reference_bases) != 1 or len(v.alternate_bases) != 1 or set(v.reference_bases + v.alternate_bases) - set("ACGT") or v.reference_bases == v.alternate_bases for v in variants):
            raise ValueError("Non-SNV request is outside authorization")
    client = atlas.create(args.key_file.read_text().strip(), timeout=30)
    failures = 0
    with args.output.open("x") as handle:
        for index, variant in enumerate(variants, 1):
            row = {"requested_allele_key": variant_key(variant), "utc": datetime.now(timezone.utc).isoformat(),
                   "sdk_version": "0.9.0", "endpoint": "gdmscience.googleapis.com:443",
                   "scorers": SCORERS, "query_sha256": digest}
            try:
                result = client.query_variant(variant, requested_scorers=SCORERS)
            except (ValueError, IndexError, TimeoutError, grpc.RpcError) as error:
                # SDK already applies its pinned bounded retry policy to transient RPCs.
                # Service exception text is deliberately not persisted (may contain credentials).
                code = error.code().name if isinstance(error, grpc.RpcError) else "NON_RPC"
                row.update(status="API_ERROR", error_type=type(error).__name__, rpc_status=code)
                failures += 1
            except PermissionError:
                raise RuntimeError("API authentication/authorization failed; service details suppressed") from None
            else:
                row.update(status="RETURNED", response=serialize_response(result, variant_key(variant)))
            handle.write(json.dumps(row, allow_nan=False) + "\n")
            handle.flush()
            if index % 20 == 0 or index == len(variants):
                print(json.dumps({"processed": index, "total": len(variants), "API_errors": failures}), flush=True)
    if failures:
        raise RuntimeError(f"API returned explicit failures for{failures} requests; inspect protected responses")


if __name__ == "__main__":
    main()
