#!/usr/bin/env python3
"""Tag primary source records and write a protected raw-allele source map."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path


PRIMARY_NUCLEAR = {str(value) for value in range(1, 23)} | {"X", "Y"}
MODULUS = 1 << 256


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_key(ordinal: int, chrom: str, pos: str, ref: str, alts: str) -> str:
    value = f"{ordinal}\t{chrom}\t{pos}\t{ref}\t{alts}"
    return hashlib.sha256(value.encode()).hexdigest()


def allele_token(record_key: str, alt_index: int) -> int:
    return int.from_bytes(hashlib.sha256(f"{record_key}:{alt_index}".encode()).digest(), "big")


def variant_class(ref: str, alt: str) -> str:
    if alt == "*":
        return "STAR"
    if alt.startswith("<") and alt.endswith(">"):
        return "SYMBOLIC"
    if "[" in alt or "]" in alt:
        return "BREAKEND"
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(ref) > len(alt):
        return "DELETION"
    if len(ref) < len(alt):
        return "INSERTION"
    return "COMPLEX"


def contig_class(chrom: str) -> str:
    if chrom in PRIMARY_NUCLEAR:
        return "PRIMARY_NUCLEAR"
    if chrom == "M":
        return "MITOCHONDRIAL"
    return "NON_PRIMARY"


def length_scope(ref: str, alt: str, kind: str) -> str:
    if kind in {"SYMBOLIC", "BREAKEND", "STAR"}:
        return f"OUT_OF_CORE_{kind}"
    return "CORE_SHORT" if max(len(ref), len(alt)) <= 50 else "OUT_OF_CORE_LONG"


def genotype_class(gt: str, alt_index: int) -> str:
    if not gt or gt == "." or "." in gt:
        return "MISSING"
    alleles = gt.replace("|", "/").split("/")
    copies = sum(value == str(alt_index) for value in alleles)
    if copies == 0:
        return "NOT_CARRIED"
    if copies == len(alleles):
        return "HOM_ALT_CARRIED"
    if copies == 1:
        return "HET_CARRIED"
    return "CARRIED_OTHER"


def parse_region(value: str | None):
    if value is None:
        return None
    chrom, interval = value.split(":", 1)
    start, end = (int(item) for item in interval.split("-", 1))
    return chrom, start, end


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-vcf", required=True, type=Path)
    parser.add_argument("--expected-raw-sha256", required=True)
    parser.add_argument("--tagged-primary-vcf", required=True, type=Path)
    parser.add_argument("--source-map", required=True, type=Path)
    parser.add_argument("--sensitive-root", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--bcftools", required=True, type=Path)
    parser.add_argument("--region")
    parser.add_argument("--expected-records", type=int)
    parser.add_argument("--expected-multiallelic-records", type=int)
    parser.add_argument("--provenance-json", required=True, type=Path)
    parser.add_argument("--stderr-log", required=True, type=Path)
    args = parser.parse_args()

    sensitive_root = args.sensitive_root.resolve()
    sensitive_root.mkdir(parents=True, exist_ok=True)
    sensitive_root.chmod(0o700)
    sensitive_paths = [
        args.tagged_primary_vcf, args.source_map, args.provenance_json, args.stderr_log
    ]
    for path in sensitive_paths:
        try:
            path.resolve().relative_to(sensitive_root)
        except ValueError as error:
            raise RuntimeError("sensitive output escaped the protected root") from error
        path.parent.mkdir(parents=True, exist_ok=True)

    observed_hash = sha256(args.raw_vcf)
    if observed_hash != args.expected_raw_sha256:
        raise RuntimeError("raw VCF SHA-256 mismatch")
    version = subprocess.run(
        [str(args.bcftools), "--version"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[0]
    if version != "bcftools 1.22":
        raise RuntimeError(f"bcftools version mismatch: {version}")

    view_command = [
        str(args.bcftools), "view", "--no-version", "--threads", "2",
        "--output-type", "z", "--output", str(args.tagged_primary_vcf), "-",
    ]
    process = subprocess.Popen(
        view_command,
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        env={**os.environ, "TMPDIR": os.environ["TMPDIR"]},
    )
    if process.stdin is None:
        raise RuntimeError("failed to open bcftools stdin")

    region = parse_region(args.region)
    counts: Counter[str] = Counter()
    grouped: Counter[tuple[str, str]] = Counter()
    source_xor = 0
    source_sum = 0
    record_ordinal = 0
    sample_columns: list[str] = []
    fieldnames = [
        "source_record_ordinal", "source_record_key", "source_chrom", "source_pos",
        "source_id", "source_ref", "source_alts", "source_alt_index", "source_alt",
        "source_qual", "source_filter", "source_gt", "source_dp", "source_gq",
        "source_ad_ref", "source_ad_alt", "contig_class", "source_variant_class",
        "length_scope", "genotype_class", "expected_normalization_disposition",
    ]
    with gzip.GzipFile(filename=str(args.source_map), mode="wb", mtime=0) as binary:
        with io.TextIOWrapper(binary, encoding="utf-8", newline="") as text_handle:
            writer = csv.DictWriter(text_handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            with open_text(args.raw_vcf) as input_handle:
                for line in input_handle:
                    if line.startswith("##"):
                        process.stdin.write(line)
                        continue
                    if line.startswith("#CHROM"):
                        process.stdin.write("##INFO=<ID=HGFACE_SRC_ORD,Number=1,Type=Integer,Description=Stable raw record ordinal>\n")
                        process.stdin.write("##INFO=<ID=HGFACE_SRC_KEY,Number=1,Type=String,Description=Stable raw record SHA256 key>\n")
                        process.stdin.write("##INFO=<ID=HGFACE_SRC_ALT_COUNT,Number=1,Type=Integer,Description=Original ALT count>\n")
                        process.stdin.write(line)
                        sample_columns = line.rstrip().split("\t")[9:]
                        if len(sample_columns) != 1:
                            raise RuntimeError("exactly one VCF sample is required")
                        continue
                    record_ordinal += 1
                    fields = line.rstrip("\n").split("\t")
                    chrom, pos_text, record_id, ref, alt_text, qual, filt, info = fields[:8]
                    pos = int(pos_text)
                    if region and not (chrom == region[0] and region[1] <= pos <= region[2]):
                        continue
                    counts["selected_records"] += 1
                    alts = alt_text.split(",")
                    if len(alts) > 1:
                        counts["multiallelic_records"] += 1
                    record_key = source_key(record_ordinal, chrom, pos_text, ref, alt_text)
                    format_keys = fields[8].split(":") if len(fields) > 8 else []
                    sample_values = fields[9].split(":") if len(fields) > 9 else []
                    sample = dict(zip(format_keys, sample_values))
                    gt = sample.get("GT", ".")
                    ad = sample.get("AD", ".").split(",")
                    for alt_index, alt in enumerate(alts, start=1):
                        kind = variant_class(ref, alt)
                        cclass = contig_class(chrom)
                        scope = length_scope(ref, alt, kind)
                        gclass = genotype_class(gt, alt_index)
                        if cclass == "PRIMARY_NUCLEAR":
                            expected = f"EXPECTED_NORMALIZED_{scope}"
                            counts["expected_normalized_alt_count"] += 1
                            token = allele_token(record_key, alt_index)
                            source_xor ^= token
                            source_sum = (source_sum + token) % MODULUS
                        else:
                            expected = f"NOT_NORMALIZED_{cclass}"
                        counts["selected_alt_count"] += 1
                        grouped[("contig_class", cclass)] += 1
                        grouped[("variant_class", kind)] += 1
                        grouped[("length_scope", scope)] += 1
                        grouped[("genotype_class", gclass)] += 1
                        grouped[("expected_disposition", expected)] += 1
                        writer.writerow(
                            {
                                "source_record_ordinal": record_ordinal,
                                "source_record_key": record_key,
                                "source_chrom": chrom,
                                "source_pos": pos_text,
                                "source_id": record_id,
                                "source_ref": ref,
                                "source_alts": alt_text,
                                "source_alt_index": alt_index,
                                "source_alt": alt,
                                "source_qual": qual,
                                "source_filter": filt,
                                "source_gt": gt,
                                "source_dp": sample.get("DP", "."),
                                "source_gq": sample.get("GQ", "."),
                                "source_ad_ref": ad[0] if ad else ".",
                                "source_ad_alt": ad[alt_index] if len(ad) > alt_index else ".",
                                "contig_class": cclass,
                                "source_variant_class": kind,
                                "length_scope": scope,
                                "genotype_class": gclass,
                                "expected_normalization_disposition": expected,
                            }
                        )
                    if chrom in PRIMARY_NUCLEAR:
                        tags = f"HGFACE_SRC_ORD={record_ordinal};HGFACE_SRC_KEY={record_key};HGFACE_SRC_ALT_COUNT={len(alts)}"
                        fields[7] = tags if info == "." else f"{info};{tags}"
                        process.stdin.write("\t".join(fields) + "\n")
                        counts["tagged_primary_records"] += 1

    process.stdin.close()
    stderr = process.stderr.read() if process.stderr else ""
    returncode = process.wait()
    args.stderr_log.write_text(stderr)
    if returncode != 0:
        raise RuntimeError("bcftools view failed while writing tagged primary VCF")
    subprocess.run(
        [str(args.bcftools), "index", "--tbi", str(args.tagged_primary_vcf)],
        check=True,
        capture_output=True,
        text=True,
    )

    if args.expected_records is not None and counts["selected_records"] != args.expected_records:
        raise RuntimeError("selected source-record count mismatch")
    if (
        args.expected_multiallelic_records is not None
        and counts["multiallelic_records"] != args.expected_multiallelic_records
    ):
        raise RuntimeError("multiallelic source-record count mismatch")

    summary = {
        "status": "PASS",
        "mode": "SHARD" if region else "FULL",
        "selected_record_count": counts["selected_records"],
        "selected_alt_count": counts["selected_alt_count"],
        "multiallelic_record_count": counts["multiallelic_records"],
        "tagged_primary_record_count": counts["tagged_primary_records"],
        "expected_normalized_alt_count": counts["expected_normalized_alt_count"],
        "expected_source_allele_xor": f"{source_xor:064x}",
        "expected_source_allele_sum": f"{source_sum:064x}",
        "grouped_counts": [
            {"dimension": dimension, "category": category, "count": count}
            for (dimension, category), count in sorted(grouped.items())
        ],
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    provenance = {
        "status": "PASS",
        "bcftools_version": version.removeprefix("bcftools "),
        "input_sha256": observed_hash,
        "mode": summary["mode"],
        "region": args.region or "ALL",
        "sample_count": len(sample_columns),
        "view_command": view_command,
    }
    args.provenance_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    for path in sensitive_paths + [Path(str(args.tagged_primary_vcf) + ".tbi")]:
        path.chmod(0o600)


if __name__ == "__main__":
    main()
