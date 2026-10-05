#!/usr/bin/env python3
"""Pure flag and transcript-selection logic for Task04 Increment 04."""

from __future__ import annotations

import json
import math
import re
from collections import defaultdict


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
SEVERITY_RANK = {term: index for index, term in enumerate(SEVERITY_ORDER)}
PROTEIN_OR_SPLICE_TERMS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}
MISSING = {None, "", ".", "NA", "None", "null"}


def compact_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def parse_int(value: str) -> int | None:
    if value in MISSING:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def numeric_tokens(value: str) -> list[float]:
    if value in MISSING:
        return []
    try:
        decoded = json.loads(value)
    except json.JSONDecodeError:
        decoded = [value]
    items = decoded if isinstance(decoded, list) else [decoded]
    numbers = []
    for item in items:
        for token in re.split(r"[&,]", str(item)):
            if token in MISSING:
                continue
            try:
                number = float(token)
            except ValueError:
                continue
            if math.isfinite(number) and number >= 0:
                numbers.append(number)
    return numbers


def gnomad_max_af(representative: dict[str, str]) -> float | None:
    values = []
    for field, value in representative.items():
        if field.endswith("_AF") and field.startswith(("gnomADe", "gnomADg")):
            values.extend(numeric_tokens(value))
    return max(values) if values else None


def consequence_terms(value: str) -> list[str]:
    try:
        decoded = json.loads(value)
    except json.JSONDecodeError:
        decoded = [value]
    terms = set()
    for item in decoded if isinstance(decoded, list) else [decoded]:
        terms.update(term for term in str(item).split("&") if term not in MISSING)
    return sorted(terms, key=lambda term: (SEVERITY_RANK.get(term, len(SEVERITY_ORDER)), term))


def format_sample(format_text: str, sample_text: str) -> dict[str, str]:
    if format_text in MISSING or sample_text in MISSING:
        return {}
    return dict(zip(format_text.split(":"), sample_text.split(":")))


def genotype_class(gt: str) -> str:
    if gt in MISSING or "." in gt:
        return "MISSING"
    alleles = gt.replace("|", "/").split("/")
    copies = sum(allele == "1" for allele in alleles)
    if copies == 0:
        return "NOT_CARRIED"
    if copies == len(alleles):
        return "HOM_ALT_CARRIED"
    if copies == 1:
        return "HET_CARRIED"
    return "CARRIED_OTHER"


def allele_depths(value: str) -> tuple[int | None, int | None]:
    if value in MISSING:
        return None, None
    parts = value.split(",")
    if len(parts) != 2:
        return None, None
    return parse_int(parts[0]), parse_int(parts[1])


def is_non_symbolic(ref: str, alt: str) -> bool:
    return bool(ref and alt) and not (
        alt == "*" or (alt.startswith("<") and alt.endswith(">")) or "[" in alt or "]" in alt
    )


def valid_identity(source: dict[str, str], allele: dict[str, str]) -> bool:
    expected = ":".join(
        [
            source.get("normalized_chrom", ""),
            source.get("normalized_pos", ""),
            source.get("normalized_ref", ""),
            source.get("normalized_alt", ""),
        ]
    )
    return bool(expected.strip(":")) and expected == source.get("normalized_allele_key") == allele.get(
        "normalized_allele_key"
    )


def derive_flags(
    source: dict[str, str], allele: dict[str, str], representative: dict[str, str]
) -> dict[str, str]:
    gt = source.get("normalized_gt", ".")
    gt_class = genotype_class(gt)
    ref_depth, alt_depth = allele_depths(source.get("normalized_ad", "."))
    dp = parse_int(source.get("normalized_dp", "."))
    gq = parse_int(source.get("normalized_gq", "."))
    allele_balance = None
    if ref_depth is not None and alt_depth is not None and ref_depth + alt_depth > 0:
        allele_balance = alt_depth / (ref_depth + alt_depth)

    missing_qc = []
    for name, value in (("GQ", gq), ("DP", dp), ("AD_REF", ref_depth), ("AD_ALT", alt_depth)):
        if value is None:
            missing_qc.append(name)
    high_confidence = (
        not missing_qc
        and source.get("normalized_filter") == "PASS"
        and gq >= 20
        and dp >= 10
        and alt_depth >= 3
        and (gt_class != "HET_CARRIED" or 0.20 <= allele_balance <= 0.80)
    )

    max_af = gnomad_max_af(representative)
    if max_af is None:
        af_state, rare = "AF_MISSING", False
    elif max_af <= 0.01:
        af_state, rare = "AF_PRESENT_RARE_RECESSIVE", True
    else:
        af_state, rare = "AF_PRESENT_ABOVE_RECESSIVE_THRESHOLD", False

    terms = consequence_terms(allele.get("consequence_terms_json", "[]"))
    protein_or_splice = bool(PROTEIN_OR_SPLICE_TERMS.intersection(terms))
    identity_valid = valid_identity(source, allele)
    non_symbolic = is_non_symbolic(
        source.get("normalized_ref", ""), source.get("normalized_alt", "")
    )
    pairable = (
        gt_class == "HET_CARRIED"
        and non_symbolic
        and source.get("length_scope") == "CORE_SHORT"
        and source.get("ref_agreement") == "PASS_BY_BCFTOOLS_CHECK_REF"
        and identity_valid
    )

    sample = format_sample(allele.get("format", "."), allele.get("sample_value", "."))
    if "|" in gt:
        phase_state = "PHASED_GT_NO_PARENTAL_TRANS_CONFIRMATION"
    elif sample.get("PGT") not in MISSING or sample.get("PID") not in MISSING:
        phase_state = "LOCAL_PHASE_TAG_NO_PARENTAL_TRANS_CONFIRMATION"
    else:
        phase_state = "UNPHASED_NO_PARENTS"

    return {
        "genotype_class": gt_class,
        "ad_ref": "" if ref_depth is None else str(ref_depth),
        "ad_alt": "" if alt_depth is None else str(alt_depth),
        "allele_balance": "" if allele_balance is None else f"{allele_balance:.8g}",
        "qc_missing_fields_json": compact_json(missing_qc),
        "high_confidence_call": str(high_confidence).lower(),
        "max_gnomad_af": "" if max_af is None else f"{max_af:.12g}",
        "gnomad_af_state": af_state,
        "rare_recessive": str(rare).lower(),
        "protein_or_splice_review": str(protein_or_splice).lower(),
        "task05_pairable": str(pairable).lower(),
        "normalized_identity_valid": str(identity_valid).lower(),
        "non_symbolic": str(non_symbolic).lower(),
        "phase_state": phase_state,
        "pgt": sample.get("PGT", ""),
        "pid": sample.get("PID", ""),
        "ps": sample.get("PS", ""),
    }


def row_severity(row: dict[str, str]) -> int:
    terms = row.get("Consequence", "").split("&")
    return min((SEVERITY_RANK.get(term, len(SEVERITY_ORDER)) for term in terms), default=len(SEVERITY_ORDER))


def transcript_key(row: dict[str, str]) -> tuple:
    if row.get("MANE_SELECT"):
        transcript_class = 0
    elif row.get("MANE_PLUS_CLINICAL"):
        transcript_class = 1
    elif row.get("CANONICAL") == "YES" and row.get("BIOTYPE") == "protein_coding":
        transcript_class = 2
    else:
        transcript_class = 3
    return transcript_class, row_severity(row), row.get("Feature", "")


def select_gene_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped = defaultdict(list)
    for row in rows:
        symbol = row.get("SYMBOL", "")
        if row.get("Feature_type") == "Transcript" and symbol not in MISSING:
            grouped[symbol].append(row)
    selected = []
    for symbol in sorted(grouped):
        gene_rows = grouped[symbol]
        best = min(gene_rows, key=transcript_key)
        terms = consequence_terms(compact_json([row.get("Consequence", "") for row in gene_rows]))
        selected.append(
            {
                **best,
                "gene_symbol": symbol,
                "gene_transcript_row_count": str(len(gene_rows)),
                "gene_consequence_terms_json": compact_json(terms),
                "gene_severity_rank": str(min(row_severity(row) for row in gene_rows)),
            }
        )
    return selected
