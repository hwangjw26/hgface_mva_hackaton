#!/usr/bin/env python3
"""Materialize lossless VEP CSQ rows and deterministic allele-level views."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import os
import re
from collections import Counter
from pathlib import Path


SEVERITY_ORDER = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
]
SEVERITY_RANK = {term: rank for rank, term in enumerate(SEVERITY_ORDER)}

KEY_FIELDS = [
    "mode",
    "record_ordinal",
    "chrom",
    "pos",
    "id",
    "ref",
    "alt",
    "alt_index",
    "normalized_allele_key",
    "source_record_ordinal",
    "source_record_key",
    "source_alt_count",
    "source_alt_index",
    "source_allele_key",
]

ALLELE_FIELDS = KEY_FIELDS + [
    "qual",
    "filter",
    "sample_name",
    "format",
    "sample_value",
    "annotation_state",
    "csq_row_count",
    "gene_count",
    "feature_count",
    "consequence_terms_json",
    "highest_impact",
    "embedded_clinvar_state",
    "embedded_CLIN_SIG_json",
    "embedded_SOMATIC_json",
    "embedded_PHENO_json",
    "current_clinvar_state",
    "current_ClinVarCurrent_json",
    "current_CLNSIG_json",
    "current_CLNREVSTAT_json",
    "current_CLNDN_json",
    "current_CLNDISDB_json",
    "current_CLNSIGCONF_json",
    "frequency_state",
    "AF_json",
    "EAS_AF_json",
    "gnomADe_AF_json",
    "gnomADe_EAS_AF_json",
    "gnomADg_AF_json",
    "gnomADg_EAS_AF_json",
    "MAX_AF_json",
    "MAX_AF_POPS_json",
    "predictor_state",
    "SIFT_json",
    "PolyPhen_json",
    "vep_release",
    "vep_cache_assembly",
    "gnomad_release",
    "embedded_clinvar_release",
    "current_clinvar_release",
    "sift_release",
    "polyphen_release",
]

TRANSCRIPT_PREFIX_FIELDS = KEY_FIELDS + [
    "record_csq_row_index",
    "allele_csq_row_index",
]

REPRESENTATIVE_PREFIX_FIELDS = KEY_FIELDS + [
    "annotation_state",
    "representative_csq_row_index",
    "representative_allele_csq_row_index",
    "selection_basis",
    "severity_rank",
]

AGGREGATE_FIELDS = {
    "embedded_CLIN_SIG_json": "CLIN_SIG",
    "embedded_SOMATIC_json": "SOMATIC",
    "embedded_PHENO_json": "PHENO",
    "current_ClinVarCurrent_json": "ClinVarCurrent",
    "current_CLNSIG_json": "ClinVarCurrent_CLNSIG",
    "current_CLNREVSTAT_json": "ClinVarCurrent_CLNREVSTAT",
    "current_CLNDN_json": "ClinVarCurrent_CLNDN",
    "current_CLNDISDB_json": "ClinVarCurrent_CLNDISDB",
    "current_CLNSIGCONF_json": "ClinVarCurrent_CLNSIGCONF",
    "AF_json": "AF",
    "EAS_AF_json": "EAS_AF",
    "gnomADe_AF_json": "gnomADe_AF",
    "gnomADe_EAS_AF_json": "gnomADe_EAS_AF",
    "gnomADg_AF_json": "gnomADg_AF",
    "gnomADg_EAS_AF_json": "gnomADg_EAS_AF",
    "MAX_AF_json": "MAX_AF",
    "MAX_AF_POPS_json": "MAX_AF_POPS",
    "SIFT_json": "SIFT",
    "PolyPhen_json": "PolyPhen",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else path.open("r")


def parse_header(handle) -> tuple[list[str], str | None]:
    header = []
    for line in handle:
        if line.startswith("#"):
            header.append(line.rstrip("\n"))
        else:
            return header, line
    return header, None


def csq_schema_from_header(header: list[str]) -> list[str]:
    matches = [line for line in header if line.startswith("##INFO=<ID=CSQ,")]
    if len(matches) != 1 or "Format: " not in matches[0]:
        raise RuntimeError("annotated VCF must contain exactly one parseable CSQ header")
    schema_text = matches[0].split("Format: ", 1)[1].rsplit('">', 1)[0]
    fields = schema_text.split("|")
    required = {
        "Allele",
        "Consequence",
        "IMPACT",
        "Gene",
        "Feature",
        "ALLELE_NUM",
        "CANONICAL",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "gnomADe_AF",
        "gnomADg_AF",
        "CLIN_SIG",
        "ClinVarCurrent_CLNSIG",
    }
    if len(fields) != len(set(fields)) or not required <= set(fields):
        raise RuntimeError("CSQ header schema is duplicated or incomplete")
    return fields


def parse_info(text: str) -> dict[str, str]:
    result = {}
    if text == ".":
        return result
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item, ""
        result[key] = value
    return result


def parse_record(line: str) -> list[str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        raise RuntimeError("malformed VCF data row")
    return fields


def unique_values(annotations: list[dict[str, str]], field: str) -> list[str]:
    return sorted({row[field] for row in annotations if row.get(field) not in (None, "", ".")})


def compact_json(values: list[str]) -> str:
    return json.dumps(values, separators=(",", ":"), ensure_ascii=True)


def consequence_rank(annotation: dict[str, str]) -> int:
    terms = annotation.get("Consequence", "").split("&")
    return min((SEVERITY_RANK.get(term, len(SEVERITY_ORDER)) for term in terms), default=len(SEVERITY_ORDER))


def representative_key(item: tuple[int, int, dict[str, str]]) -> tuple:
    record_index, allele_index, annotation = item
    if annotation.get("MANE_SELECT"):
        transcript_class = 0
    elif annotation.get("MANE_PLUS_CLINICAL"):
        transcript_class = 1
    elif annotation.get("CANONICAL") == "YES" and annotation.get("BIOTYPE") == "protein_coding":
        transcript_class = 2
    else:
        transcript_class = 3
    return (
        transcript_class,
        consequence_rank(annotation),
        annotation.get("Feature", ""),
        record_index,
        allele_index,
    )


def selection_basis(annotation: dict[str, str]) -> str:
    if annotation.get("MANE_SELECT"):
        return "MANE_SELECT_THEN_SEVERITY"
    if annotation.get("MANE_PLUS_CLINICAL"):
        return "MANE_PLUS_CLINICAL_THEN_SEVERITY"
    if annotation.get("CANONICAL") == "YES" and annotation.get("BIOTYPE") == "protein_coding":
        return "CANONICAL_PROTEIN_CODING_THEN_SEVERITY"
    return "CONSEQUENCE_SEVERITY_THEN_FEATURE"


def source_identity(
    mode: str,
    record_ordinal: int,
    input_row: list[str],
    input_info: dict[str, str],
    alt_index: int,
) -> tuple[str, str, str, str]:
    if mode == "PUBLIC_SYNTHETIC":
        source_record_key = hashlib.sha256(
            (f"PUBLIC_SYNTHETIC\t{record_ordinal}\t" + "\t".join(input_row[:8])).encode()
        ).hexdigest()
        return str(record_ordinal), source_record_key, str(len(input_row[4].split(","))), str(alt_index)

    required = ("HGFACE_SRC_ORD", "HGFACE_SRC_KEY", "HGFACE_SRC_ALT_COUNT")
    if any(not input_info.get(field) for field in required):
        raise RuntimeError("protected input lacks stable HGFACE source tags")
    source_record_key = input_info["HGFACE_SRC_KEY"]
    if not re.fullmatch(r"[0-9a-f]{64}", source_record_key):
        raise RuntimeError("invalid HGFACE source key")
    source_alt_count = int(input_info["HGFACE_SRC_ALT_COUNT"])
    if len(input_row[4].split(",")) != 1:
        raise RuntimeError("protected normalized records must be biallelic")
    old_record = input_info.get("HGFACE_ORIG", "")
    if old_record:
        old_parts = old_record.split("|")
        if len(old_parts) != 5 or not old_parts[-1].isdigit():
            raise RuntimeError("invalid HGFACE_ORIG tag")
        source_alt_index = int(old_parts[-1])
    else:
        if source_alt_count != 1:
            raise RuntimeError("multiallelic source record lacks HGFACE_ORIG ALT index")
        source_alt_index = 1
    if not 1 <= source_alt_index <= source_alt_count:
        raise RuntimeError("source ALT index is outside source ALT count")
    return (
        input_info["HGFACE_SRC_ORD"],
        source_record_key,
        str(source_alt_count),
        str(source_alt_index),
    )


def deterministic_tsv(path: Path, fieldnames: list[str]):
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    writer = csv.DictWriter(text, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return raw, compressed, text, writer


def close_tsv(handles: tuple) -> None:
    raw, compressed, text, _writer = handles
    text.flush()
    text.close()
    if not compressed.closed:
        compressed.close()
    if not raw.closed:
        raw.close()


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--annotated-vcf", required=True, type=Path)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--expected-annotated-sha256", required=True)
    parser.add_argument("--expected-records", required=True, type=int)
    parser.add_argument("--expected-alleles", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument(
        "--mode",
        choices=["PUBLIC_SYNTHETIC", "PROTECTED_SHARD", "FULL_PATIENT"],
        required=True,
    )
    args = parser.parse_args()

    csv.field_size_limit(64 * 1024 * 1024)
    if not re.fullmatch(r"[a-z0-9_]+", args.prefix):
        raise RuntimeError("prefix must use lowercase letters, digits, or underscores")
    if sha256(args.input_vcf) != args.expected_input_sha256:
        raise RuntimeError("input VCF SHA-256 mismatch")
    if sha256(args.annotated_vcf) != args.expected_annotated_sha256:
        raise RuntimeError("annotated VCF SHA-256 mismatch")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.mode != "PUBLIC_SYNTHETIC":
        args.output_dir.chmod(0o700)

    output_paths = {
        "alleles": args.output_dir / f"{args.prefix}.allele_annotations.tsv.gz",
        "transcripts": args.output_dir / f"{args.prefix}.transcript_annotations.tsv.gz",
        "representatives": args.output_dir / f"{args.prefix}.representative_transcripts.tsv.gz",
        "summary": args.output_dir / f"{args.prefix}.annotation_summary.json",
        "schema_tsv": args.output_dir / f"{args.prefix}.csq_schema.tsv",
        "schema_json": args.output_dir / f"{args.prefix}.csq_schema.json",
        "provenance": args.output_dir / f"{args.prefix}.parser_provenance.json",
    }
    for path in output_paths.values():
        if path.exists():
            raise RuntimeError(f"refusing to overwrite existing output: {path.name}")
    temporary_paths = {name: Path(str(path) + ".tmp") for name, path in output_paths.items()}

    record_count = 0
    allele_count = 0
    transcript_count = 0
    representative_count = 0
    no_csq_count = 0
    consequence_counts: Counter[str] = Counter()
    feature_type_counts: Counter[str] = Counter()
    evidence_counts: Counter[str] = Counter()
    allele_digest = hashlib.sha256()
    transcript_digest = hashlib.sha256()

    allele_handles = transcript_handles = representative_handles = None
    try:
        with open_text(args.input_vcf) as input_handle, open_text(args.annotated_vcf) as annotated_handle:
            input_header, first_input = parse_header(input_handle)
            annotated_header, first_annotated = parse_header(annotated_handle)
            input_columns = next((line for line in input_header if line.startswith("#CHROM")), None)
            annotated_columns = next((line for line in annotated_header if line.startswith("#CHROM")), None)
            if input_columns is None or input_columns != annotated_columns:
                raise RuntimeError("input and annotated VCF column headers differ")
            column_fields = input_columns.lstrip("#").split("\t")
            sample_names = column_fields[9:]
            if args.mode == "PUBLIC_SYNTHETIC" and sample_names:
                raise RuntimeError("public synthetic fixture unexpectedly contains a sample")
            if args.mode != "PUBLIC_SYNTHETIC" and len(sample_names) != 1:
                raise RuntimeError("protected singleton input must contain exactly one sample")
            csq_schema = csq_schema_from_header(annotated_header)

            transcript_fields = TRANSCRIPT_PREFIX_FIELDS + csq_schema
            representative_fields = REPRESENTATIVE_PREFIX_FIELDS + csq_schema
            allele_handles = deterministic_tsv(temporary_paths["alleles"], ALLELE_FIELDS)
            transcript_handles = deterministic_tsv(temporary_paths["transcripts"], transcript_fields)
            representative_handles = deterministic_tsv(temporary_paths["representatives"], representative_fields)
            allele_writer = allele_handles[3]
            transcript_writer = transcript_handles[3]
            representative_writer = representative_handles[3]

            input_records = itertools.chain([first_input], input_handle) if first_input else iter(())
            annotated_records = itertools.chain([first_annotated], annotated_handle) if first_annotated else iter(())
            for record_ordinal, pair in enumerate(
                itertools.zip_longest(input_records, annotated_records), start=1
            ):
                input_line, annotated_line = pair
                if input_line is None or annotated_line is None:
                    raise RuntimeError("input and annotated VCF record counts differ")
                input_row = parse_record(input_line)
                annotated_row = parse_record(annotated_line)
                if input_row[:7] != annotated_row[:7] or input_row[8:] != annotated_row[8:]:
                    raise RuntimeError(f"VCF lockstep mismatch at record {record_ordinal}")
                record_count += 1
                input_info = parse_info(input_row[7])
                annotated_info = parse_info(annotated_row[7])
                for stable_tag in ("HGFACE_SRC_ORD", "HGFACE_SRC_KEY", "HGFACE_SRC_ALT_COUNT", "HGFACE_ORIG"):
                    if input_info.get(stable_tag) != annotated_info.get(stable_tag):
                        raise RuntimeError(f"source tag mismatch at record {record_ordinal}: {stable_tag}")

                annotations_by_alt: dict[int, list[tuple[int, dict[str, str]]]] = {
                    index: [] for index in range(1, len(input_row[4].split(",")) + 1)
                }
                csq_text = annotated_info.get("CSQ", "")
                if csq_text:
                    for record_csq_index, csq_value in enumerate(csq_text.split(","), start=1):
                        values = csq_value.split("|")
                        if len(values) > len(csq_schema):
                            raise RuntimeError(f"CSQ row wider than schema at record {record_ordinal}")
                        values.extend([""] * (len(csq_schema) - len(values)))
                        annotation = dict(zip(csq_schema, values))
                        allele_number = annotation.get("ALLELE_NUM", "")
                        if not allele_number.isdigit() or int(allele_number) not in annotations_by_alt:
                            raise RuntimeError(f"invalid ALLELE_NUM at record {record_ordinal}")
                        annotations_by_alt[int(allele_number)].append((record_csq_index, annotation))

                alts = input_row[4].split(",")
                for alt_index, alt in enumerate(alts, start=1):
                    allele_count += 1
                    source_ordinal, source_key, source_alt_count, source_alt_index = source_identity(
                        args.mode, record_ordinal, input_row, input_info, alt_index
                    )
                    key_values = {
                        "mode": args.mode,
                        "record_ordinal": str(record_ordinal),
                        "chrom": input_row[0],
                        "pos": input_row[1],
                        "id": input_row[2],
                        "ref": input_row[3],
                        "alt": alt,
                        "alt_index": str(alt_index),
                        "normalized_allele_key": f"{input_row[0]}:{input_row[1]}:{input_row[3]}:{alt}",
                        "source_record_ordinal": source_ordinal,
                        "source_record_key": source_key,
                        "source_alt_count": source_alt_count,
                        "source_alt_index": source_alt_index,
                        "source_allele_key": f"{source_key}:{source_alt_index}",
                    }
                    allele_digest.update(
                        ("\t".join(key_values[field] for field in KEY_FIELDS) + "\n").encode()
                    )
                    indexed_annotations = annotations_by_alt[alt_index]
                    annotations = [annotation for _index, annotation in indexed_annotations]
                    annotation_state = "ANNOTATED" if annotations else "NON_CALL_NO_CSQ"
                    if not annotations:
                        no_csq_count += 1

                    for allele_csq_index, (record_csq_index, annotation) in enumerate(
                        indexed_annotations, start=1
                    ):
                        transcript_count += 1
                        transcript_row = {
                            **key_values,
                            "record_csq_row_index": str(record_csq_index),
                            "allele_csq_row_index": str(allele_csq_index),
                            **annotation,
                        }
                        transcript_writer.writerow(transcript_row)
                        transcript_digest.update(
                            (
                                key_values["source_allele_key"]
                                + f"\t{allele_csq_index}\t"
                                + "|".join(annotation[field] for field in csq_schema)
                                + "\n"
                            ).encode()
                        )
                        for term in annotation.get("Consequence", "").split("&"):
                            if term:
                                consequence_counts[term] += 1
                        if annotation.get("Feature_type"):
                            feature_type_counts[annotation["Feature_type"]] += 1

                    aggregate_values = {
                        output_field: unique_values(annotations, csq_field)
                        for output_field, csq_field in AGGREGATE_FIELDS.items()
                    }
                    embedded_present = any(
                        aggregate_values[field]
                        for field in ("embedded_CLIN_SIG_json", "embedded_SOMATIC_json", "embedded_PHENO_json")
                    )
                    current_present = any(
                        aggregate_values[field]
                        for field in (
                            "current_ClinVarCurrent_json",
                            "current_CLNSIG_json",
                            "current_CLNREVSTAT_json",
                            "current_CLNDN_json",
                            "current_CLNDISDB_json",
                            "current_CLNSIGCONF_json",
                        )
                    )
                    frequency_present = any(
                        aggregate_values[field]
                        for field in (
                            "AF_json",
                            "EAS_AF_json",
                            "gnomADe_AF_json",
                            "gnomADe_EAS_AF_json",
                            "gnomADg_AF_json",
                            "gnomADg_EAS_AF_json",
                            "MAX_AF_json",
                        )
                    )
                    predictor_present = any(
                        aggregate_values[field] for field in ("SIFT_json", "PolyPhen_json")
                    )
                    for source, present in (
                        ("embedded_clinvar", embedded_present),
                        ("current_clinvar", current_present),
                        ("frequency", frequency_present),
                        ("predictor", predictor_present),
                    ):
                        evidence_counts[f"{source}_{'present' if present else 'missing'}"] += 1

                    consequences = sorted(
                        {
                            term
                            for annotation in annotations
                            for term in annotation.get("Consequence", "").split("&")
                            if term
                        },
                        key=lambda term: (SEVERITY_RANK.get(term, len(SEVERITY_ORDER)), term),
                    )
                    impacts = {annotation.get("IMPACT", "") for annotation in annotations}
                    highest_impact = next(
                        (impact for impact in ("HIGH", "MODERATE", "LOW", "MODIFIER") if impact in impacts),
                        "",
                    )
                    allele_row = {
                        **key_values,
                        "qual": input_row[5],
                        "filter": input_row[6],
                        "sample_name": sample_names[0] if sample_names else "",
                        "format": input_row[8] if len(input_row) > 8 else "",
                        "sample_value": input_row[9] if len(input_row) > 9 else "",
                        "annotation_state": annotation_state,
                        "csq_row_count": str(len(annotations)),
                        "gene_count": str(len({a["Gene"] for a in annotations if a.get("Gene")})),
                        "feature_count": str(len({a["Feature"] for a in annotations if a.get("Feature")})),
                        "consequence_terms_json": compact_json(consequences),
                        "highest_impact": highest_impact,
                        "embedded_clinvar_state": "PRESENT" if embedded_present else "MISSING",
                        "current_clinvar_state": "PRESENT" if current_present else "MISSING",
                        "frequency_state": "PRESENT" if frequency_present else "MISSING",
                        "predictor_state": "PRESENT" if predictor_present else "MISSING",
                        "vep_release": "115.0",
                        "vep_cache_assembly": "GRCh38.p14",
                        "gnomad_release": "v4.1",
                        "embedded_clinvar_release": "202502",
                        "current_clinvar_release": "2026-08-22",
                        "sift_release": "6.2.1",
                        "polyphen_release": "2.2.3",
                    }
                    allele_row.update(
                        {field: compact_json(values) for field, values in aggregate_values.items()}
                    )
                    allele_writer.writerow(allele_row)

                    representative_count += 1
                    if indexed_annotations:
                        candidates = [
                            (record_index, allele_index, annotation)
                            for allele_index, (record_index, annotation) in enumerate(indexed_annotations, start=1)
                        ]
                        record_index, allele_csq_index, selected = min(candidates, key=representative_key)
                        representative_row = {
                            **key_values,
                            "annotation_state": annotation_state,
                            "representative_csq_row_index": str(record_index),
                            "representative_allele_csq_row_index": str(allele_csq_index),
                            "selection_basis": selection_basis(selected),
                            "severity_rank": str(consequence_rank(selected)),
                            **selected,
                        }
                    else:
                        representative_row = {
                            **key_values,
                            "annotation_state": annotation_state,
                            "representative_csq_row_index": "",
                            "representative_allele_csq_row_index": "",
                            "selection_basis": "NON_CALL_NO_CSQ",
                            "severity_rank": "",
                            **{field: "" for field in csq_schema},
                        }
                    representative_writer.writerow(representative_row)

        if record_count != args.expected_records:
            raise RuntimeError(f"record count mismatch: {record_count} != {args.expected_records}")
        if allele_count != args.expected_alleles:
            raise RuntimeError(f"allele count mismatch: {allele_count} != {args.expected_alleles}")
        if representative_count != allele_count:
            raise RuntimeError("representative table is not one row per allele")

        for handles in (allele_handles, transcript_handles, representative_handles):
            close_tsv(handles)
        allele_handles = transcript_handles = representative_handles = None

        schema_payload = {
            "source": "annotated VCF ##INFO=<ID=CSQ header",
            "field_count": len(csq_schema),
            "fields": csq_schema,
        }
        write_json(temporary_paths["schema_json"], schema_payload)
        with temporary_paths["schema_tsv"].open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["field_index", "field_name", "schema_source"])
            for index, field in enumerate(csq_schema, start=1):
                writer.writerow([index, field, "VCF_CSQ_HEADER"])

        summary_payload = {
            "status": "PASS",
            "mode": args.mode,
            "input_record_count": record_count,
            "input_allele_count": allele_count,
            "allele_row_count": allele_count,
            "transcript_row_count": transcript_count,
            "representative_row_count": representative_count,
            "non_call_no_csq_count": no_csq_count,
            "csq_schema_field_count": len(csq_schema),
            "allele_identity_sha256": allele_digest.hexdigest(),
            "transcript_identity_sha256": transcript_digest.hexdigest(),
            "consequence_row_counts": dict(sorted(consequence_counts.items())),
            "feature_type_row_counts": dict(sorted(feature_type_counts.items())),
            "evidence_presence_counts": dict(sorted(evidence_counts.items())),
            "representative_view_is_filter": False,
        }
        write_json(temporary_paths["summary"], summary_payload)
        provenance_payload = {
            "status": "PASS",
            "mode": args.mode,
            "input_vcf_sha256": args.expected_input_sha256,
            "annotated_vcf_sha256": args.expected_annotated_sha256,
            "runtime_derived_csq_schema": True,
            "all_csq_rows_retained": True,
            "explicit_non_call_state": "NON_CALL_NO_CSQ",
            "representative_precedence": [
                "MANE_SELECT",
                "MANE_PLUS_CLINICAL",
                "CANONICAL_PROTEIN_CODING",
                "CONSEQUENCE_SEVERITY",
                "FEATURE_ID",
                "CSQ_ROW_INDEX",
            ],
            "source_lockstep_checked": True,
            "remote_coordinate_query": False,
        }
        write_json(temporary_paths["provenance"], provenance_payload)

        for name, final_path in output_paths.items():
            os.replace(temporary_paths[name], final_path)
            final_path.chmod(0o600)
    except Exception:
        for handles in (allele_handles, transcript_handles, representative_handles):
            if handles is not None:
                try:
                    close_tsv(handles)
                except Exception:
                    pass
        for path in temporary_paths.values():
            path.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    main()
