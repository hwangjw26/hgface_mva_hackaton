#!/usr/bin/env python3
"""Offline exact-allele AVI annotation; never uploads variants or changes ranking."""
import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter
from pathlib import Path

CUTOFF = 0.99379  # Qu & Hakonarson v1, Table 2; abstract rounds this to 0.9938.


def read_table(path, delimiter="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def key_parts(key):
    chrom, pos, ref, alt = key.split(":")
    chrom = chrom.removeprefix("chr")
    pos = int(pos)
    if pos < 1 or not ref or not alt or set(ref + alt) - set("ACGT") or ref == alt:
        raise ValueError("Invalid allele key")
    return chrom, pos, ref, alt


def canonical_key(key):
    return ":".join(map(str, key_parts(key)))


def snv(key):
    return len(key_parts(key)[2]) == len(key_parts(key)[3]) == 1


def check_reference(fasta, key):
    chrom, pos, ref, alt = key_parts(key)
    fasta_chrom = chrom if chrom in fasta.references else "chr" + chrom
    if fasta.fetch(fasta_chrom, pos - 1, pos - 1 + len(ref)).upper() != ref:
        raise ValueError("GRCh38 REF mismatch")


def score_annotation(raw_text, phred_text):
    raw, phred = float(raw_text), float(phred_text)
    if not math.isfinite(raw) or not math.isfinite(phred) or phred < 0:
        raise ValueError("Nonfinite/invalid AVI score")
    quantile = -math.expm1(-phred * math.log(10) / 10)
    # PHRED is a rounded rank transform, not a clinical probability.
    # Preserve uncertainty if last-digit rounding can cross the published cutoff.
    from decimal import Decimal
    half_unit = 0.5 * 10 ** Decimal(phred_text).as_tuple().exponent
    lo = -math.expm1(-max(0, phred - half_unit) * math.log(10) / 10)
    hi = -math.expm1(-(phred + half_unit) * math.log(10) / 10)
    flag = "ROUNDING_BOUNDARY" if lo <= CUTOFF <= hi else str(quantile >= CUTOFF).lower()
    return {"avi_status": "SCORED", "avi_raw": raw_text, "avi_phred": phred_text,
            "avi_quantile_derived": format(quantile, ".12g"), "avi_meets_0_99379": flag}


def lookup_score(tabix, key):
    chrom, pos, ref, alt = key_parts(key)
    source_chrom = "chr" + chrom
    if source_chrom not in tabix.contigs:
        raise ValueError("AVI resource chromosome unavailable")
    matches = []
    for line in tabix.fetch(source_chrom, pos - 1, pos):
        fields = line.split("\t")
        if len(fields) != 6:
            raise ValueError("AVI row does not match the pinned six-column schema")
        row_key = canonical_key(":".join(fields[:4]))
        if key_parts(row_key)[1] == pos and fields[2] != ref:
            raise ValueError("AVI resource REF disagrees with candidate/reference")
        if row_key == canonical_key(key):
            matches.append(score_annotation(fields[4], fields[5]))
    if len(matches) > 1:
        raise ValueError("Duplicate exact AVI allele")
    return matches[0] if matches else None


def summarize_pair(first, second):
    if any(x["avi_status"] != "SCORED" for x in (first, second)):
        return "INCOMPLETE"
    if any(x["avi_meets_0_99379"] == "ROUNDING_BOUNDARY" for x in (first, second)):
        return "ROUNDING_BOUNDARY"
    count = sum(x["avi_meets_0_99379"] == "true" for x in (first, second))
    return ("NEITHER_ABOVE", "ONE_ABOVE", "BOTH_ABOVE")[count]


def write_table(path, rows):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2 ** 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    import pysam  # Only reference/bulk lookup needs pysam; cache replay uses stdlib.

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--avi", type=Path, help="Local official bgzip TSV with existing .tbi")
    parser.add_argument("--allow-unavailable", action="store_true", help="Write pending inventory, NOT scored results")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    os.umask(0o077)
    if not args.avi and not args.allow_unavailable:
        parser.error("Local AVI resource required; use --allow-unavailable only for a pending inventory")
    if args.avi and args.allow_unavailable:
        parser.error("Choose either a real AVI resource or explicit pending mode")
    for path in [args.fasta, Path(str(args.fasta) + ".fai"), Path(str(args.fasta) + ".gzi")]:
        if not path.is_file():
            raise ValueError("Existing indexed GRCh38 bgzip FASTA required")
    inputs = [args.baseline / name for name in ("integrated_pairs.tsv", "allele_evidence.tsv", "predictions.csv")]
    input_hashes = {str(p): sha256(p) for p in inputs}
    pairs, contexts = read_table(inputs[0]), read_table(inputs[1])
    submission = read_table(inputs[2], ",")
    if (len(pairs), len(contexts), len(submission)) != (120, 238, 4):
        raise ValueError("Baseline counts changed")
    if len({p["review_id"] for p in pairs}) != 120:
        raise ValueError("Duplicate pair identifier")
    context_keys = {(x["gene_id"], x["allele_key"]) for x in contexts}
    if len(context_keys) != 238:
        raise ValueError("Duplicate gene–allele context")
    needed = {(p["gene_id"], p[f"allele_{n}_key"]) for p in pairs for n in (1, 2)}
    if needed != context_keys:
        raise ValueError("Pair/context universe mismatch")
    physical_keys = sorted({x["allele_key"] for x in contexts}, key=key_parts)
    if len({canonical_key(k) for k in physical_keys}) != len(physical_keys):
        raise ValueError("Canonical allele collision")
    tabix = None
    if args.avi:
        if not args.avi.is_file() or not Path(str(args.avi) + ".tbi").is_file():
            raise ValueError("Existing local AVI TSV/index required; remote queries are not supported")
        tabix = pysam.TabixFile(str(args.avi))
        if list(tabix.header) != ["#CHROM\tPOS\tREF\tALT\traw_score\tPHRED"]:
            raise ValueError("AVI header differs from pinned official schema")
    observations = {}
    with pysam.FastaFile(str(args.fasta)) as fasta:
        for key in physical_keys:
            check_reference(fasta, key)
            obs = {"avi_status": "", "avi_raw": "", "avi_phred": "", "avi_quantile_derived": "", "avi_meets_0_99379": ""}
            if not snv(key):
                obs["avi_status"] = "NON_SNV_OUTSIDE_PAPER_SCOPE"
            elif tabix is None:
                obs["avi_status"] = "PENDING_AVI_RESOURCE"
            else:
                score = lookup_score(tabix, key)
                if score is None:
                    obs["avi_status"] = "SNV_NOT_IN_RESOURCE"
                else:
                    obs.update(score)
            observations[key] = obs
    if tabix is not None:
        tabix.close()
    annotated_alleles = [{**x, **observations[x["allele_key"]]} for x in contexts]
    annotated_pairs = []
    for pair in pairs:
        row = dict(pair)
        a, b = (observations[pair[f"allele_{n}_key"]] for n in (1, 2))
        for n, allele in enumerate((a, b), 1):
            row.update({f"avi_{n}_{k.removeprefix('avi_')}": v for k, v in allele.items()})
        row["avi_pair_status"] = summarize_pair(a, b)
        annotated_pairs.append(row)
    final = []
    for rank, prediction in enumerate(submission, 1):
        target = {canonical_key(":".join(prediction[f"{f}_{n}"] for f in ("chrom", "pos", "ref", "alt"))) for n in (1, 2)}
        matches = [p for p in annotated_pairs if {canonical_key(p[f"allele_{n}_key"]) for n in (1, 2)} == target]
        if len(target) != 2 or len(matches) != 1:
            raise ValueError("Submission-to-pair binding is not unique")
        final.append({"submission_rank": rank, **matches[0]})
    args.output.mkdir(parents=True, exist_ok=False, mode=0o700)
    write_table(args.output / "allele_avi.tsv", annotated_alleles)
    write_table(args.output / "pair_avi.tsv", annotated_pairs)
    write_table(args.output / "current_four_avi.tsv", final)
    resource = None if args.avi is None else {"path": str(args.avi.resolve()), "bytes": args.avi.stat().st_size, "mtime_ns": args.avi.stat().st_mtime_ns, "index_sha256": sha256(Path(str(args.avi) + ".tbi"))}
    manifest = {"inputs_sha256": input_hashes, "avi_resource": resource,
                "fasta": {"path": str(args.fasta), "bytes": args.fasta.stat().st_size, "fai_sha256": sha256(Path(str(args.fasta) + ".fai")), "gzi_sha256": sha256(Path(str(args.fasta) + ".gzi"))},
                "cutoff": CUTOFF, "operator": ">= (explicit operational convention; paper fixed-threshold tables use >=)", "quantile": "1 - 10**(-PHRED/10); derived from static PHRED, not directly retrieved API AVI_QUANTILE", "code_sha256": sha256(Path(__file__))}
    (args.output / "provenance.private.json").write_text(json.dumps(manifest, indent=2) + "\n")
    summary = {"state": "PENDING_AVI_RESOURCE" if args.avi is None else "ANNOTATED", "pairs": len(pairs), "allele_contexts": len(contexts), "unique_alleles": len(observations), "allele_status": dict(Counter(x["avi_status"] for x in observations.values())), "pair_status": dict(Counter(x["avi_pair_status"] for x in annotated_pairs)), "current_four_status": dict(Counter(x["avi_pair_status"] for x in final)), "submission_unchanged": sha256(inputs[2]) == input_hashes[str(inputs[2])]}
    if any(sha256(p) != input_hashes[str(p)] for p in inputs):
        raise ValueError("Baseline input changed during annotation")
    (args.output / "aggregate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
