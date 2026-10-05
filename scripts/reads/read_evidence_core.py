#!/usr/bin/env python3
"""Frozen deterministic Task06 read-evidence rules."""

from __future__ import annotations

from collections.abc import Iterable, Mapping


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def classify_allele(values: Mapping[str, object]) -> str:
    required_true = ("reference_ok", "probe_ok", "representation_ok", "retrieval_ok")
    if not all(bool(values.get(key, False)) for key in required_true):
        return "READ_NONCALL"
    depth = int(values["qualified_depth"])
    alt = int(values["independent_alt"])
    if depth < 10:
        return "READ_NONCALL"
    fraction = alt / depth
    if (
        alt >= 5
        and 0.20 <= fraction <= 0.80
        and int(values["alt_forward"]) >= 1
        and int(values["alt_reverse"]) >= 1
        and int(values["supporting_lanes"]) >= 2
        and not bool(values.get("mapping_caveat", False))
    ):
        return "STRONG_READ_SUPPORT"
    if alt >= 3 and 0.10 <= fraction <= 0.90:
        return "READ_SUPPORT_WITH_CAVEAT"
    if alt <= 1 or fraction < 0.10:
        return "READ_CONTRADICTED"
    return "READ_NONCALL"


def classify_pair(allele_labels: Iterable[str]) -> str:
    labels = list(allele_labels)
    if len(labels) != 2:
        raise ValueError("a pair must contain exactly two allele labels")
    supported = {"STRONG_READ_SUPPORT", "READ_SUPPORT_WITH_CAVEAT"}
    if labels.count("STRONG_READ_SUPPORT") == 2:
        return "PAIR_BOTH_ALLELES_STRONG"
    if "READ_CONTRADICTED" in labels:
        return "PAIR_READ_CONTRADICTED"
    if all(label in supported for label in labels):
        return "PAIR_BOTH_ALLELES_SUPPORTED"
    if sum(label in supported for label in labels) == 1:
        return "PAIR_ONE_ALLELE_SUPPORTED"
    return "PAIR_READ_NONCALL"


def classify_phase(haplotypes: Mapping[str, int]) -> str:
    cis = int(haplotypes.get("ALT_ALT", 0)) >= 2
    trans = int(haplotypes.get("REF_ALT", 0)) >= 2 and int(haplotypes.get("ALT_REF", 0)) >= 2
    if cis and trans:
        return "PHASE_CONFLICT_UNRESOLVED"
    if cis:
        return "LOCAL_CIS_SUPPORTED"
    if trans:
        return "LOCAL_TRANS_SUPPORTED"
    return "PHASE_UNRESOLVED_NO_PARENTS"
