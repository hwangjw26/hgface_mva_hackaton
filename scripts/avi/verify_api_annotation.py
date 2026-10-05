#!/usr/bin/env python3
"""Reconstruct API-to-allele/pair joins without importing the producer."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def table(path, sep="\t"):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter=sep))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "responses", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    a = parser.parse_args()
    checks = []

    def check(label, ok):
        if not ok:
            raise ValueError("Verification failed: " + label)
        checks.append(label)

    originals = [table(a.baseline / "allele_evidence.tsv"), table(a.baseline / "integrated_pairs.tsv")]
    alleles, pairs = table(a.output / "allele_avi.tsv"), table(a.output / "pair_avi.tsv")
    for name, old, new in zip(["allele", "pair"], originals, [alleles, pairs]):
        check(name + " count", len(old) == len(new))
        for i, (before, after) in enumerate(zip(old, new)):
            check(name + " source fields " + str(i), all(after[k] == v for k, v in before.items()))
    responses = [json.loads(line) for line in a.responses.read_text().splitlines()]
    cached = {x["requested_allele_key"]: x for x in responses}
    check("201 unique successful responses", len(responses) == len(cached) == 201 and all(x["status"] == "RETURNED" for x in responses))
    physical = {}
    for i, row in enumerate(alleles):
        physical[row["allele_key"]] = row
        if row["allele_key"] not in cached:
            ref, alt = row["allele_key"].split(":")[2:]
            check("nonSNV scope " + str(i), (len(ref) != 1 or len(alt) != 1) and row["avi_status"] == "NON_SNV_OUTSIDE_PAPER_SCOPE")
            check("no nonSNV scores " + str(i), all(v == "" for k, v in row.items() if k.startswith("avi_") and k != "avi_status"))
            continue
        returned = cached[row["allele_key"]]["response"]
        check("returned identity " + str(i), all(x["returned_allele_key"] == row["allele_key"] for x in returned.values()))
        score, feature = returned["AVI_SCORE"], returned["AVI_SCORE_FEATURE_IMPORTANCE"]
        q, raw = score["layers"]["quantiles"][0][0], score["values"][0][0]
        check("raw/q/cutoff " + str(i), row["avi_status"] == "SCORED_API" and float(row["avi_raw"]) == raw and float(row["avi_quantile"]) == q and row["avi_meets_global"] == str(q >= 0.99379).lower())
        expected = dict(zip(feature["track_names"], feature["values"][0]))
        check("18 exact feature values " + str(i), len(expected) == 18 and json.loads(row["avi_feature_importance_json"]) == expected)
        check("top feature " + str(i), abs(expected[row["avi_top_feature"]]) == max(abs(x) for x in expected.values()))
        check("feature residual " + str(i), float(row["avi_feature_sum_minus_raw"]) == sum(feature["values"][0]) - raw)
        check("derived PHRED " + str(i), (q == 1 and row["avi_phred_derived"] == "INF") or (q < 1 and math.isclose(float(row["avi_phred_derived"]), -10 * math.log10(1-q), abs_tol=1e-9)))
    for i, row in enumerate(pairs):
        x, y = [physical[row[f"allele_{n}_key"]] for n in (1, 2)]
        for n, allele in enumerate((x, y), 1):
            check("pair allele join " + str(i) + ":" + str(n), all(row[f"avi_{n}_{k[4:]}"] == v for k, v in allele.items() if k.startswith("avi_")))
        expected = "INCOMPLETE" if any(v["avi_status"] != "SCORED_API" for v in [x, y]) else ["NEITHER_MEETS", "ONE_MEETS", "BOTH_MEET"][sum(float(v["avi_quantile"]) >= 0.99379 for v in [x, y])]
        check("pair state " + str(i), row["avi_pair_status"] == expected)
    final, submission = table(a.output / "current_four_avi.tsv"), table(a.baseline / "predictions.csv", ",")
    check("four recommendations", len(final) == len(submission) == 4)
    for i, (row, pred) in enumerate(zip(final, submission), 1):
        expected = {":".join([pred[f"chrom_{n}"].removeprefix("chr"), pred[f"pos_{n}"], pred[f"ref_{n}"], pred[f"alt_{n}"]]) for n in (1, 2)}
        bound = [p for p in pairs if p["review_id"] == row["review_id"]]
        check("complete final row " + str(i), len(bound) == 1 and {k:v for k,v in row.items() if k != "submission_rank"} == bound[0] and row["submission_rank"] == str(i) and {row["allele_1_key"],row["allele_2_key"]} == expected)
    summary = json.loads((a.output / "aggregate_summary.json").read_text())
    check("summary count units", (summary["unique_alleles"],summary["allele_contexts"],summary["pairs"]) == (len(physical),len(alleles),len(pairs)))
    for label, rows, field in [("allele_status",physical.values(),"avi_status"),("pair_status",pairs,"avi_pair_status"),("current_four_status",final,"avi_pair_status")]:
        check("summary " + label, summary[label] == dict(Counter(row[field] for row in rows)))
    check("summary SNV flags", summary["snv_meets_global"] == dict(Counter(row["avi_meets_global"] for row in physical.values() if row["avi_status"] == "SCORED_API")))
    check("completed API state", summary["state"] == "ANNOTATED_API")
    receipt = json.loads((a.output / "provenance.private.json").read_text())
    for path, digest in receipt["inputs_sha256"].items():
        check("unchanged input " + Path(path).name, hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest)
    check("protected directory", a.output.stat().st_mode & 0o777 == 0o700)
    for path in a.output.iterdir():
        check("protected file " + path.name, path.stat().st_mode & 0o777 == 0o600)
    print(json.dumps({"verification":"PASS","checks":len(checks),"scope":"all201 responses,238 contexts,120 pairs and four recommendations"}))


if __name__ == "__main__":
    main()
