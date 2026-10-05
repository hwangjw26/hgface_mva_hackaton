"""Exact source excerpts; see canonical_provenance.json."""
import collections
import datetime
import hashlib
import json
import math
import re
from urllib.parse import unquote

PLP_ATOMS = {"pathogenic", "likely_pathogenic"}

BLB_ATOMS = {"benign", "likely_benign"}

STARRED_REVIEW_STATUS = {"criteria_provided,_single_submitter", "criteria_provided,_multiple_submitters,_no_conflicts",
                         "reviewed_by_expert_panel", "practice_guideline"}

def contig_rank(chrom: str) -> int:
    value = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if value.isdigit() and 1 <= int(value) <= 22:
        return int(value)
    ranks = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if value in ranks:
        return ranks[value]
    if value.startswith("SYN_"):
        return 100
    raise ValueError("chromosome is outside canonical contig order")

def allele_key(chrom: str, pos: int | str, ref: str, alt: str) -> str:
    chrom = chrom[3:] if chrom.startswith("chr") else chrom
    return f"{'MT' if chrom == 'M' else chrom}:{int(pos)}:{ref}:{alt}"

def order_tuple(key: str) -> tuple:
    chrom, pos, ref, alt = key.split(":")
    return contig_rank(chrom), int(pos), ref, alt, key

def canonical_pair(key_a: str, key_b: str) -> tuple[str, str]:
    if key_a == key_b:
        raise ValueError("self-pair")
    return tuple(sorted((key_a, key_b), key=order_tuple))

def pair_key(key_a: str, key_b: str) -> str:
    first, second = canonical_pair(key_a, key_b)
    return hashlib.sha256(f"{first}|{second}".encode()).hexdigest()

def _atoms(values: list[str]) -> set[str]:
    atoms: set[str] = set()
    for value in values:
        atoms.update(p for p in re.split(r"[&/,]", value.strip().lower().replace(" ", "_")) if p)
    return atoms

def clinvar_class(row: dict) -> str:
    current = json.loads(row.get("current_CLNSIG_json") or "[]")
    embedded = json.loads(row.get("embedded_CLIN_SIG_json") or "[]")
    revstat = {s.strip().lower().replace("&", ",") for s in json.loads(row.get("current_CLNREVSTAT_json") or "[]")}
    if row.get("current_clinvar_state") == "PRESENT":
        atoms = _atoms(current)
        if atoms and atoms <= PLP_ATOMS:
            return "PLP_NONCONFLICTING"
        starred = bool(revstat) and revstat <= STARRED_REVIEW_STATUS
        if atoms and atoms <= BLB_ATOMS and starred:
            return "BLB_NONCONFLICTING_STARRED"
        return "PRESENT_OTHER_OR_CONFLICTING"
    if row.get("embedded_clinvar_state") == "PRESENT":
        atoms = _atoms(embedded)
        return "PLP_NONCONFLICTING" if atoms and atoms <= PLP_ATOMS else "PRESENT_OTHER_OR_CONFLICTING"
    return "MISSING"
