#!/usr/bin/env python3
"""Run deterministic bcftools normalization inside an allowed local root."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_within(path: Path, root: Path, label: str) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as error:
        raise RuntimeError(f"{label} must be under the declared root") from error
    return resolved


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    path.chmod(0o600)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--expected-fasta-sha256", required=True)
    parser.add_argument("--bcftools", required=True, type=Path)
    parser.add_argument("--regions")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--provenance-json", required=True, type=Path)
    parser.add_argument("--stdout-log", required=True, type=Path)
    parser.add_argument("--stderr-log", required=True, type=Path)
    args = parser.parse_args()

    if not 1 <= args.threads <= 4:
        raise RuntimeError("threads must be between 1 and 4")
    input_vcf = require_within(args.input_vcf, args.input_root, "input VCF")
    output_vcf = require_within(args.output_vcf, args.output_root, "output VCF")
    provenance = require_within(args.provenance_json, args.output_root, "provenance")
    stdout_log = require_within(args.stdout_log, args.output_root, "stdout log")
    stderr_log = require_within(args.stderr_log, args.output_root, "stderr log")
    args.output_root.mkdir(parents=True, exist_ok=True)
    for path in (output_vcf, provenance, stdout_log, stderr_log):
        path.parent.mkdir(parents=True, exist_ok=True)

    observed_fasta_sha = sha256(args.fasta)
    if observed_fasta_sha != args.expected_fasta_sha256:
        raise RuntimeError("approved FASTA SHA-256 mismatch")
    version_line = subprocess.run(
        [str(args.bcftools), "--version"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[0]
    if version_line != "bcftools 1.22":
        raise RuntimeError(f"bcftools version mismatch: {version_line}")

    command = [
        str(args.bcftools),
        "norm",
        "--fasta-ref", str(args.fasta),
        "--check-ref", "e",
        "--multiallelics", "-any",
        "--keep-sum", "AD",
        "--old-rec-tag", "HGFACE_ORIG",
        "--site-win", "1000000",
        "--sort", "lex",
        "--threads", str(args.threads),
        "--no-version",
        "--output-type", "z",
        "--write-index=tbi",
        "--output", str(output_vcf),
    ]
    if args.regions:
        command.extend(["--regions", args.regions])
    command.append(str(input_vcf))

    run_env = os.environ.copy()
    run_env["TMPDIR"] = os.environ["TMPDIR"]
    result = subprocess.run(command, capture_output=True, text=True, env=run_env)
    stdout_log.write_text(result.stdout)
    stderr_log.write_text(result.stderr)
    stdout_log.chmod(0o600)
    stderr_log.chmod(0o600)
    payload = {
        "bcftools_version": version_line.removeprefix("bcftools "),
        "command": command,
        "fasta_sha256": observed_fasta_sha,
        "input_sha256": sha256(input_vcf),
        "regions": args.regions or "ALL",
        "remote_query": False,
        "returncode": result.returncode,
        "status": "PASS" if result.returncode == 0 else "FAIL",
    }
    if result.returncode != 0:
        write_json(provenance, payload)
        print("bcftools normalization failed; inspect the protected stderr log", file=sys.stderr)
        raise SystemExit(result.returncode)

    index_path = Path(str(output_vcf) + ".tbi")
    if not output_vcf.is_file() or not index_path.is_file():
        raise RuntimeError("normalized VCF or tabix index missing")
    count_result = subprocess.run(
        [str(args.bcftools), "index", "--nrecords", str(output_vcf)],
        check=True,
        capture_output=True,
        text=True,
    )
    payload.update(
        {
            "normalized_record_count": int(count_result.stdout.strip()),
            "output_sha256": sha256(output_vcf),
            "output_tbi_sha256": sha256(index_path),
        }
    )
    output_vcf.chmod(0o600)
    index_path.chmod(0o600)
    write_json(provenance, payload)


if __name__ == "__main__":
    main()
