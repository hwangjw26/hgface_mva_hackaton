#!/usr/bin/env python3
"""Run deterministic cache-only VEP 115 annotation inside an approved root."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import resource
import subprocess
import time
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
        raise RuntimeError(f"{label} is outside its approved root") from error
    return resolved


def stream_record_count(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return sum(not line.startswith("#") for line in handle)


def indexed_record_count(bcftools: Path, path: Path) -> int:
    result = subprocess.run(
        [str(bcftools), "index", "--nrecords", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout.strip())


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    path.chmod(0o600)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--allowed-input-root", required=True, type=Path)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--expected-records", required=True, type=int)
    parser.add_argument("--require-input-index", action="store_true")
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--allowed-output-root", required=True, type=Path)
    parser.add_argument("--vep-prefix", required=True, type=Path)
    parser.add_argument("--expected-vep-sha256", required=True)
    parser.add_argument("--cache-root", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--expected-fasta-sha256", required=True)
    parser.add_argument("--clinvar", required=True, type=Path)
    parser.add_argument("--expected-clinvar-sha256", required=True)
    parser.add_argument("--expected-clinvar-tbi-sha256", required=True)
    parser.add_argument("--bcftools", required=True, type=Path)
    parser.add_argument("--mode", choices=["PUBLIC_SYNTHETIC", "PROTECTED_SHARD", "FULL_PATIENT"], required=True)
    parser.add_argument("--forks", type=int, default=1)
    parser.add_argument("--buffer-size", type=int, default=5000)
    parser.add_argument("--provenance-json", required=True, type=Path)
    parser.add_argument("--stdout-log", required=True, type=Path)
    parser.add_argument("--stderr-log", required=True, type=Path)
    args = parser.parse_args()

    if not 1 <= args.forks <= 8:
        raise RuntimeError("forks must be between 1 and 8")
    if not 500 <= args.buffer_size <= 50000:
        raise RuntimeError("buffer size must be between 500 and 50000")
    input_vcf = require_within(args.input_vcf, args.allowed_input_root, "input VCF")
    output_vcf = require_within(args.output_vcf, args.allowed_output_root, "output VCF")
    provenance = require_within(args.provenance_json, args.allowed_output_root, "provenance")
    stdout_log = require_within(args.stdout_log, args.allowed_output_root, "stdout log")
    stderr_log = require_within(args.stderr_log, args.allowed_output_root, "stderr log")
    if output_vcf.suffix != ".gz" or not output_vcf.name.endswith(".vcf.gz"):
        raise RuntimeError("output must end in .vcf.gz")
    for path in (output_vcf, provenance, stdout_log, stderr_log):
        path.parent.mkdir(parents=True, exist_ok=True)

    observed_input_sha = sha256(input_vcf)
    if observed_input_sha != args.expected_input_sha256:
        raise RuntimeError("input VCF SHA-256 mismatch")
    input_index = Path(str(input_vcf) + ".tbi")
    if args.require_input_index and not input_index.is_file():
        raise RuntimeError("required input TBI is missing")
    input_count = (
        indexed_record_count(args.bcftools, input_vcf)
        if args.require_input_index
        else stream_record_count(input_vcf)
    )
    if input_count != args.expected_records:
        raise RuntimeError("input record count mismatch")

    cache_leaf = args.cache_root / "homo_sapiens" / "115_GRCh38"
    cache_info = cache_leaf / "info.txt"
    if not cache_info.is_file() or args.fasta.parent.resolve() != cache_leaf.resolve():
        raise RuntimeError("VEP cache or cache-matched FASTA is unavailable")
    info_text = cache_info.read_text()
    required_cache_lines = {
        "assembly\tGRCh38",
        "source_assembly\tGRCh38.p14",
        "source_gencode\tGENCODE 49",
        "source_gnomADe\tv4.1",
        "source_gnomADg\tv4.1",
        "source_ClinVar\t202502",
        "source_sift\t6.2.1",
        "source_polyphen\t2.2.3",
    }
    if not required_cache_lines <= set(info_text.splitlines()):
        raise RuntimeError("VEP cache metadata contract mismatch")
    if sha256(args.fasta) != args.expected_fasta_sha256:
        raise RuntimeError("FASTA SHA-256 mismatch")
    clinvar_tbi = Path(str(args.clinvar) + ".tbi")
    if not args.clinvar.is_file() or not clinvar_tbi.is_file():
        raise RuntimeError("current ClinVar VCF or TBI is unavailable")
    if sha256(args.clinvar) != args.expected_clinvar_sha256:
        raise RuntimeError("current ClinVar SHA-256 mismatch")
    if sha256(clinvar_tbi) != args.expected_clinvar_tbi_sha256:
        raise RuntimeError("current ClinVar TBI SHA-256 mismatch")

    vep = args.vep_prefix / "bin" / "vep"
    if sha256(vep.resolve()) != args.expected_vep_sha256:
        raise RuntimeError("VEP executable SHA-256 mismatch")
    run_env = os.environ.copy()
    run_env["PATH"] = f"{args.vep_prefix / 'bin'}:/usr/bin:/bin"
    run_env["TMPDIR"] = os.environ["TMPDIR"]
    help_result = subprocess.run(
        [str(vep), "--help"], check=True, capture_output=True, text=True, env=run_env
    )
    version = re.search(r"ensembl-vep\s*:\s*([0-9.]+)", help_result.stdout)
    if not version or version.group(1) != "115.0":
        raise RuntimeError("VEP runtime is not exactly 115.0")
    bcftools_version = subprocess.run(
        [str(args.bcftools), "--version"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[0]
    if bcftools_version != "bcftools 1.22":
        raise RuntimeError("bcftools runtime is not exactly 1.22")

    ref_command = [
        str(args.bcftools), "norm", "--fasta-ref", str(args.fasta),
        "--check-ref", "e", "--do-not-normalize", "--output-type", "u",
        "--output", "/dev/null", str(input_vcf),
    ]
    ref_result = subprocess.run(
        ref_command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, env=run_env
    )
    if ref_result.returncode != 0:
        stderr_log.write_text(ref_result.stderr)
        stderr_log.chmod(0o600)
        raise RuntimeError("input REF agreement preflight failed")

    custom_spec = (
        f"{args.clinvar},ClinVarCurrent,vcf,exact,0,"
        "CLNSIG,CLNREVSTAT,CLNDN,CLNDISDB,CLNSIGCONF"
    )
    command = [
        str(vep),
        "--input_file", str(input_vcf),
        "--output_file", str(output_vcf),
        "--vcf", "--force_overwrite", "--no_stats",
        "--offline", "--cache", "--cache_version", "115",
        "--dir_cache", str(args.cache_root),
        "--species", "homo_sapiens", "--assembly", "GRCh38",
        "--fasta", str(args.fasta),
        "--symbol", "--numbers", "--canonical", "--mane", "--biotype",
        "--protein", "--uniprot", "--domains", "--variant_class",
        "--allele_number", "--hgvs", "--regulatory",
        "--af", "--af_1kg", "--af_gnomade", "--af_gnomadg",
        "--max_af", "--sift", "b", "--polyphen", "b",
        "--custom", custom_spec,
        "--compress_output", "bgzip",
        "--buffer_size", str(args.buffer_size),
    ]
    if args.forks > 1:
        command.extend(["--fork", str(args.forks)])

    start = time.monotonic()
    with stdout_log.open("w") as stdout_handle, stderr_log.open("w") as stderr_handle:
        if ref_result.stderr:
            stderr_handle.write(ref_result.stderr)
        result = subprocess.run(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            env=run_env,
        )
    elapsed = time.monotonic() - start
    stdout_log.chmod(0o600)
    stderr_log.chmod(0o600)
    if result.returncode != 0:
        write_json(
            provenance,
            {
                "status": "FAIL",
                "mode": args.mode,
                "returncode": result.returncode,
                "remote_coordinate_query": False,
                "vep_command": command,
            },
        )
        raise RuntimeError("VEP annotation failed; inspect the protected stderr log")
    if not output_vcf.is_file() or output_vcf.stat().st_size == 0:
        raise RuntimeError("VEP output is missing")

    subprocess.run(
        [str(args.bcftools), "index", "--tbi", "--force", str(output_vcf)],
        check=True,
        capture_output=True,
        text=True,
    )
    output_tbi = Path(str(output_vcf) + ".tbi")
    output_count = indexed_record_count(args.bcftools, output_vcf)
    if output_count != input_count:
        raise RuntimeError("VEP input/output record count mismatch")
    header = subprocess.run(
        [str(args.bcftools), "view", "--header-only", str(output_vcf)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if "##INFO=<ID=CSQ" not in header or '##VEP="v115.0"' not in header:
        raise RuntimeError("annotated VCF lacks the expected VEP/CSQ header")

    for path in (output_vcf, output_tbi, stdout_log, stderr_log):
        path.chmod(0o600)
    payload = {
        "status": "PASS",
        "mode": args.mode,
        "input_record_count": input_count,
        "output_record_count": output_count,
        "input_sha256": observed_input_sha,
        "output_sha256": sha256(output_vcf),
        "output_tbi_sha256": sha256(output_tbi),
        "vep_runtime_version": version.group(1),
        "vep_cache_version": "115",
        "vep_cache_assembly": "GRCh38.p14",
        "clinvar_current_release": "2026-08-22",
        "execution_mode": "offline_cache_only",
        "remote_coordinate_query": False,
        "filtering_options_present": False,
        "forks": args.forks,
        "buffer_size": args.buffer_size,
        "elapsed_seconds": round(elapsed, 3),
        "records_per_second": round(output_count / elapsed, 3),
        "maximum_resident_set_kb": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        "reference_preflight_command": ref_command,
        "reference_preflight_returncode": ref_result.returncode,
        "vep_command": command,
    }
    write_json(provenance, payload)


if __name__ == "__main__":
    main()
