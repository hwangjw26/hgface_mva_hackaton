#!/usr/bin/env python3
"""Build a normalized-record-to-source-allele link and aggregate accounting."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import subprocess
from collections import Counter
from pathlib import Path


MODULUS = 1 << 256


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def allele_token(record_key: str, alt_index: int) -> int:
    return int.from_bytes(hashlib.sha256(f"{record_key}:{alt_index}".encode()).digest(), "big")


def parse_info(value: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    if value == ".":
        return parsed
    for item in value.split(";"):
        key, separator, field = item.partition("=")
        parsed[key] = field if separator else "true"
    return parsed


def source_alt_index(original: str | None) -> int:
    if not original or original == ".":
        return 1
    parts = original.split("|")
    if len(parts) == 4:
        return 1
    if len(parts) == 5:
        return int(parts[-1])
    raise RuntimeError("unexpected HGFACE_ORIG schema")


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


def length_scope(ref: str, alt: str, kind: str) -> str:
    if kind in {"SYMBOLIC", "BREAKEND", "STAR"}:
        return f"OUT_OF_CORE_{kind}"
    return "CORE_SHORT" if max(len(ref), len(alt)) <= 50 else "OUT_OF_CORE_LONG"


def genotype_class(gt: str) -> str:
    if not gt or gt == "." or "." in gt:
        return "MISSING"
    alleles = gt.replace("|", "/").split("/")
    copies = sum(value == "1" for value in alleles)
    if copies == 0:
        return "NOT_CARRIED"
    if copies == len(alleles):
        return "HOM_ALT_CARRIED"
    if copies == 1:
        return "HET_CARRIED"
    return "CARRIED_OTHER"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--normalized-vcf", required=True, type=Path)
    parser.add_argument("--prepare-summary", required=True, type=Path)
    parser.add_argument("--source-link", required=True, type=Path)
    parser.add_argument("--duplicate-keys", required=True, type=Path)
    parser.add_argument("--sensitive-root", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--provenance-json", required=True, type=Path)
    parser.add_argument("--bcftools", required=True, type=Path)
    args = parser.parse_args()

    sensitive_root = args.sensitive_root.resolve()
    for path in (args.source_link, args.duplicate_keys, args.provenance_json):
        try:
            path.resolve().relative_to(sensitive_root)
        except ValueError as error:
            raise RuntimeError("sensitive output escaped the protected root") from error
        path.parent.mkdir(parents=True, exist_ok=True)

    prepare = json.loads(args.prepare_summary.read_text())
    expected_count = int(prepare["expected_normalized_alt_count"])
    expected_xor = prepare["expected_source_allele_xor"]
    expected_sum = prepare["expected_source_allele_sum"]
    index_count = int(
        subprocess.run(
            [str(args.bcftools), "index", "--nrecords", str(args.normalized_vcf)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    )

    link_fields = [
        "normalized_record_ordinal", "normalized_chrom", "normalized_pos",
        "normalized_ref", "normalized_alt", "normalized_allele_key", "normalized_id",
        "normalized_filter", "normalized_gt", "normalized_ad", "normalized_dp",
        "normalized_gq", "normalized_variant_class", "length_scope",
        "normalization_disposition", "ref_agreement", "source_record_ordinal",
        "source_record_key", "source_alt_index", "source_alt_count", "original_record_tag",
    ]
    duplicate_fields = ["normalized_allele_key", "record_count"]
    counts: Counter[str] = Counter()
    grouped: Counter[tuple[str, str]] = Counter()
    observed_xor = 0
    observed_sum = 0
    previous_key: str | None = None
    previous_count = 0

    with gzip.GzipFile(filename=str(args.source_link), mode="wb", mtime=0) as link_binary:
        with io.TextIOWrapper(link_binary, encoding="utf-8", newline="") as link_text:
            link_writer = csv.DictWriter(link_text, delimiter="\t", fieldnames=link_fields, lineterminator="\n")
            link_writer.writeheader()
            with gzip.GzipFile(filename=str(args.duplicate_keys), mode="wb", mtime=0) as dup_binary:
                with io.TextIOWrapper(dup_binary, encoding="utf-8", newline="") as dup_text:
                    dup_writer = csv.DictWriter(dup_text, delimiter="\t", fieldnames=duplicate_fields, lineterminator="\n")
                    dup_writer.writeheader()
                    with gzip.open(args.normalized_vcf, "rt") as handle:
                        for line in handle:
                            if line.startswith("#"):
                                continue
                            fields = line.rstrip("\n").split("\t")
                            chrom, pos, record_id, ref, alt, _qual, filt, info_text = fields[:8]
                            if "," in alt:
                                raise RuntimeError("normalized VCF is not biallelic")
                            info = parse_info(info_text)
                            record_key = info.get("HGFACE_SRC_KEY")
                            source_ordinal = info.get("HGFACE_SRC_ORD")
                            source_alt_count = info.get("HGFACE_SRC_ALT_COUNT")
                            if not record_key or not source_ordinal or not source_alt_count:
                                raise RuntimeError("normalized record lacks stable source tags")
                            alt_index = source_alt_index(info.get("HGFACE_ORIG"))
                            if not 1 <= alt_index <= int(source_alt_count):
                                raise RuntimeError("source ALT index is out of range")

                            format_keys = fields[8].split(":") if len(fields) > 8 else []
                            sample_values = fields[9].split(":") if len(fields) > 9 else []
                            sample = dict(zip(format_keys, sample_values))
                            kind = variant_class(ref, alt)
                            scope = length_scope(ref, alt, kind)
                            disposition = f"NORMALIZED_{scope}"
                            normalized_key = f"{chrom}:{pos}:{ref}:{alt}"
                            counts["normalized_records"] += 1
                            token = allele_token(record_key, alt_index)
                            observed_xor ^= token
                            observed_sum = (observed_sum + token) % MODULUS
                            grouped[("variant_class", kind)] += 1
                            grouped[("length_scope", scope)] += 1
                            grouped[("genotype_class", genotype_class(sample.get("GT", ".")))] += 1
                            grouped[("normalization_disposition", disposition)] += 1

                            if previous_key is None:
                                previous_key, previous_count = normalized_key, 1
                            elif normalized_key == previous_key:
                                previous_count += 1
                            else:
                                if previous_count > 1:
                                    dup_writer.writerow({"normalized_allele_key": previous_key, "record_count": previous_count})
                                    counts["duplicate_key_groups"] += 1
                                    counts["duplicate_key_rows"] += previous_count
                                previous_key, previous_count = normalized_key, 1

                            link_writer.writerow(
                                {
                                    "normalized_record_ordinal": counts["normalized_records"],
                                    "normalized_chrom": chrom,
                                    "normalized_pos": pos,
                                    "normalized_ref": ref,
                                    "normalized_alt": alt,
                                    "normalized_allele_key": normalized_key,
                                    "normalized_id": record_id,
                                    "normalized_filter": filt,
                                    "normalized_gt": sample.get("GT", "."),
                                    "normalized_ad": sample.get("AD", "."),
                                    "normalized_dp": sample.get("DP", "."),
                                    "normalized_gq": sample.get("GQ", "."),
                                    "normalized_variant_class": kind,
                                    "length_scope": scope,
                                    "normalization_disposition": disposition,
                                    "ref_agreement": "PASS_BY_BCFTOOLS_CHECK_REF",
                                    "source_record_ordinal": source_ordinal,
                                    "source_record_key": record_key,
                                    "source_alt_index": alt_index,
                                    "source_alt_count": source_alt_count,
                                    "original_record_tag": info.get("HGFACE_ORIG", "."),
                                }
                            )
                    if previous_key is not None and previous_count > 1:
                        dup_writer.writerow({"normalized_allele_key": previous_key, "record_count": previous_count})
                        counts["duplicate_key_groups"] += 1
                        counts["duplicate_key_rows"] += previous_count

    observed_xor_hex = f"{observed_xor:064x}"
    observed_sum_hex = f"{observed_sum:064x}"
    reconciled = (
        counts["normalized_records"] == expected_count == index_count
        and observed_xor_hex == expected_xor
        and observed_sum_hex == expected_sum
    )
    if not reconciled:
        raise RuntimeError("normalized/source allele accounting mismatch")

    summary = {
        "status": "PASS",
        "normalized_record_count": counts["normalized_records"],
        "indexed_record_count": index_count,
        "expected_source_alt_count": expected_count,
        "source_reconciliation": "PASS",
        "duplicate_normalized_key_groups": counts["duplicate_key_groups"],
        "duplicate_normalized_key_rows": counts["duplicate_key_rows"],
        "grouped_counts": [
            {"dimension": dimension, "category": category, "count": count}
            for (dimension, category), count in sorted(grouped.items())
        ],
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    provenance = {
        "status": "PASS",
        "normalized_vcf_sha256": sha256(args.normalized_vcf),
        "source_link_sha256": sha256(args.source_link),
        "duplicate_keys_sha256": sha256(args.duplicate_keys),
        "source_allele_xor": observed_xor_hex,
        "source_allele_sum": observed_sum_hex,
    }
    args.provenance_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    for path in (args.source_link, args.duplicate_keys, args.provenance_json):
        path.chmod(0o600)


if __name__ == "__main__":
    main()
