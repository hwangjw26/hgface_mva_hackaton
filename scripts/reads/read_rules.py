import hashlib
import re

READ_SUPPORTED = {"STRONG_READ_SUPPORT", "READ_SUPPORT_WITH_CAVEAT"}

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

def pair_read_label(label_1: str | None, label_2: str | None) -> str:
    labels = [label_1 or "READ_NOT_ASSESSED", label_2 or "READ_NOT_ASSESSED"]
    if labels.count("STRONG_READ_SUPPORT") == 2:
        return "PAIR_BOTH_ALLELES_STRONG"
    if "READ_CONTRADICTED" in labels:
        return "PAIR_READ_CONTRADICTED"
    supported = sum(label in READ_SUPPORTED for label in labels)
    if supported == 2:
        return "PAIR_BOTH_ALLELES_SUPPORTED"
    if supported == 1:
        return "PAIR_ONE_ALLELE_SUPPORTED"
    if labels == ["READ_NOT_ASSESSED", "READ_NOT_ASSESSED"]:
        return "READ_NOT_ASSESSED"
    return "PAIR_READ_NONCALL"
