#!/usr/bin/env python3
"""Join frozen API responses locally; preserve all baseline evidence and ranks."""
import argparse
import json
import math
import os
from collections import Counter
from pathlib import Path

from annotate_avi import CUTOFF, canonical_key, read_table, sha256, snv, write_table


def api_observation(record, expected_features):
    blank = {"avi_status": "API_ERROR", "avi_raw": "", "avi_quantile": "", "avi_phred_derived": "",
             "avi_meets_global": "", "avi_feature_importance_json": "", "avi_top_feature": "",
             "avi_feature_sum_minus_raw": ""}
    if record["status"] == "API_ERROR":
        return blank
    if record["status"] != "RETURNED" or record["sdk_version"] != "0.9.0":
        raise ValueError("Unexpected API response state/version")
    response = record["response"]
    if set(response) != {"AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"}:
        raise ValueError("Unexpected response scorers")
    score, importance = response["AVI_SCORE"], response["AVI_SCORE_FEATURE_IMPORTANCE"]
    if any(x["returned_allele_key"] != record["requested_allele_key"] for x in (score, importance)):
        raise ValueError("Returned/requested allele mismatch")
    if score["track_names"] != ["AVI_SCORE"] or importance["track_names"] != expected_features:
        raise ValueError("Response feature identities changed")
    if len(score["values"]) != 1 or len(score["values"][0]) != 1 or len(importance["values"]) != 1 or len(importance["values"][0]) != 18:
        raise ValueError("Unexpected response shape")
    quantiles = score["layers"].get("quantiles")
    if quantiles is None or len(quantiles) != 1 or len(quantiles[0]) != 1:
        raise ValueError("Missing/unexpected API AVI quantile")
    raw, quantile = score["values"][0][0], quantiles[0][0]
    values = importance["values"][0]
    if not all(math.isfinite(x) for x in [raw, quantile, *values]) or not 0 <= quantile <= 1:
        raise ValueError("Invalid/nonfinite API score")
    feature_values = dict(zip(expected_features, values, strict=True))
    top = max(feature_values, key=lambda name: abs(feature_values[name]))
    phred = "INF" if quantile == 1 else format(-10 * math.log10(1 - quantile), ".12g")
    return {"avi_status": "SCORED_API", "avi_raw": repr(raw), "avi_quantile": repr(quantile),
            "avi_phred_derived": phred, "avi_meets_global": str(quantile >= CUTOFF).lower(),
            "avi_feature_importance_json": json.dumps(feature_values, sort_keys=True),
            "avi_top_feature": top, "avi_feature_sum_minus_raw": repr(sum(values) - raw)}


def pair_state(a, b):
    if any(x["avi_status"] != "SCORED_API" for x in (a, b)):
        return "INCOMPLETE"
    count = sum(x["avi_meets_global"] == "true" for x in (a, b))
    return ("NEITHER_MEETS", "ONE_MEETS", "BOTH_MEET")[count]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "responses", "query", "metadata", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    os.umask(0o077)
    paths = [args.baseline / name for name in ("allele_evidence.tsv", "integrated_pairs.tsv", "predictions.csv")]
    inputs = paths + [args.responses, args.query, args.metadata]
    before = {str(x): sha256(x) for x in inputs}
    alleles, pairs, submission = read_table(paths[0]), read_table(paths[1]), read_table(paths[2], ",")
    if (len(alleles), len(pairs), len(submission)) != (238, 120, 4):
        raise ValueError("Baseline counts changed")
    contexts = {(x["gene_id"], x["allele_key"]) for x in alleles}
    if len(contexts) != 238 or contexts != {(p["gene_id"], p[f"allele_{n}_key"]) for p in pairs for n in (1, 2)}:
        raise ValueError("Gene/allele pair universe mismatch")
    expected_snv = {x["allele_key"] for x in alleles if snv(x["allele_key"])}
    query = read_table(args.query)
    query_keys = {canonical_key(":".join([x["chromosome"], x["position_1based"], x["reference_bases"], x["alternate_bases"]])) for x in query}
    if len(query) != 201 or query_keys != expected_snv:
        raise ValueError("API request differs from the baseline SNV set")
    records = [json.loads(line) for line in args.responses.read_text().splitlines()]
    if len(records) != 201 or {x["requested_allele_key"] for x in records} != expected_snv:
        raise ValueError("Missing/duplicated/unexpected API response")
    if any(x["query_sha256"] != before[str(args.query)] or x["endpoint"] != "gdmscience.googleapis.com:443" for x in records):
        raise ValueError("Response provenance mismatch")
    features = [x["name"] for x in json.loads(args.metadata.read_text())["AVI_SCORE_FEATURE_IMPORTANCE"]["tracks"]]
    obs = {x["requested_allele_key"]: api_observation(x, features) for x in records}
    for row in alleles:
        if row["allele_key"] not in obs:
            if snv(row["allele_key"]):
                raise ValueError("Unaccounted SNV")
            obs[row["allele_key"]] = {key: ("NON_SNV_OUTSIDE_PAPER_SCOPE" if key == "avi_status" else "") for key in next(iter(obs.values()))}
    annotated_alleles = [{**row, **obs[row["allele_key"]]} for row in alleles]
    annotated_pairs = []
    for row in pairs:
        a, b = (obs[row[f"allele_{n}_key"]] for n in (1, 2))
        result = dict(row)
        for n, item in enumerate((a, b), 1):
            result.update({f"avi_{n}_{key[4:]}": value for key, value in item.items()})
        result["avi_pair_status"] = pair_state(a, b)
        annotated_pairs.append(result)
    final = []
    for rank, row in enumerate(submission, 1):
        target = {canonical_key(":".join(row[f"{key}_{n}"] for key in ("chrom", "pos", "ref", "alt"))) for n in (1, 2)}
        found = [x for x in annotated_pairs if {x["allele_1_key"], x["allele_2_key"]} == target]
        if len(target) != 2 or len(found) != 1:
            raise ValueError("Submission binding is ambiguous")
        final.append({"submission_rank": rank, **found[0]})
    args.output.mkdir(mode=0o700, parents=True, exist_ok=False)
    for name, rows in [("allele_avi.tsv", annotated_alleles), ("pair_avi.tsv", annotated_pairs), ("current_four_avi.tsv", final)]:
        write_table(args.output / name, rows)
    summary = {"state": "ANNOTATED_API" if all(x["status"] == "RETURNED" for x in records) else "INCOMPLETE_API_ERRORS",
               "unique_alleles": len(obs), "allele_contexts": len(alleles), "pairs": len(pairs),
               "allele_status": dict(Counter(x["avi_status"] for x in obs.values())),
               "snv_meets_global": dict(Counter(x["avi_meets_global"] for x in obs.values() if x["avi_status"] == "SCORED_API")),
               "pair_status": dict(Counter(x["avi_pair_status"] for x in annotated_pairs)),
               "current_four_status": dict(Counter(x["avi_pair_status"] for x in final))}
    (args.output / "aggregate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    receipt = {"inputs_sha256": before, "cutoff": CUTOFF, "operator": ">=", "score_scale": "API quantiles layer of AVI_SCORE",
               "sdk": "0.9.0", "server_dataset_version": "not exposed by metadata; frozen response hashes pin this retrieval",
               "features": "Returned feature-importance values preserved; not assumed to sum to raw score or establish mechanism",
               "script_sha256": {x.name: sha256(x) for x in (Path(__file__), Path(__file__).with_name("annotate_avi.py"))}}
    (args.output / "provenance.private.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if any(sha256(path) != before[str(path)] for path in inputs):
        raise ValueError("Input changed during annotation")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
