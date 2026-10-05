"""Independent field/byte comparisons; reference decisions are verification-only."""
import argparse
import json
import os
from pathlib import Path
from io_utils import index, require, rows, sha, write_json


def verify(result, reference):
    checked = 0
    mismatches = []
    manifest = json.loads((reference.parent / "verification_sources.private.json").read_text())
    pinned = {r["name"]: r["sha256"] for r in manifest}
    require(len(pinned) == len(manifest) and set(pinned) == {p.name for p in reference.iterdir()}, "Verification reference file set changed")
    require(all(sha(reference / name) == digest for name, digest in pinned.items()), "Verification reference checksum mismatch")

    def compare(actual_name, expected_name, keys, fields):
        nonlocal checked
        actual, expected = index(rows(result / actual_name), keys), index(rows(reference / expected_name), keys)
        require(set(actual) == set(expected), "Verification identity mismatch: " + actual_name)
        for k in sorted(actual):
            for f in fields:
                require(f in actual[k] and f in expected[k], "Verification field missing: " + f)
                checked += 1
                if actual[k][f] != expected[k][f]:
                    mismatches.append({"table": actual_name, "key": k, "field": f,
                                       "actual": actual[k][f], "expected": expected[k][f]})
    pair_fields = ["allele_1_key", "allele_2_key", "v2_tier", "v2_phase", "positive_ar_association",
                   "t20_gene_returned", "t20_ar_returned", "t20_both_expected_gene_mapped",
                   "t20_both_score_route_eligible", "t20_both_plp_route_eligible",
                   "t20_both_any_native_route_eligible", "t20_native_pair_evidence_state"]
    pair_fields += ["t20_best_" + kind + "_" + field for kind in ["any", "ar"]
                    for field in ["disease", "full_disease_rank", "native_posterior", "genotype_lr"]]
    pair_fields += [f"audit_{i}_" + field for i in [1, 2] for field in [
        "normalized_gt", "high_confidence_call", "frequency_flag", "clinvar_class", "evidence_conflicts",
        "read_label", "read_contradiction_evaluable", "gene_protein_splice", "selected_protein_splice",
        "annotation_rows", "selected_consequence", "transcript_consequence_discordance"]]
    compare("integrated_pairs.tsv", "pairs.tsv", ["gene_id", "pair_key"], pair_fields)
    allele_fields = ["selected_transcript", "selected_consequence", "selected_mane", "selected_hgvsc", "selected_hgvsp",
                    "clinvar_class", "high_confidence_call", "max_gnomad_af", "t17_condition_ids_json",
                    "t20_expected_gene_mapped", "t20_native_passed_filter", "t20_native_optional_score_present",
                    "t20_native_optional_score", "t20_native_clinvar_2602", "t20_native_alt_count",
                    "t20_score_route_allele_count", "t20_plp_route_allele_count", "t20_any_native_route_eligible",
                    "t20_best_ar_target", "t20_condition_id_context"]
    compare("allele_evidence.tsv", "alleles.tsv", ["gene_id", "allele_key"], allele_fields)
    fields = [prefix + field for prefix in ["best_any_", "best_ar_", "fixed_baseline_any_", "fixed_baseline_ar_"]
              for field in ["disease", "full_disease_rank", "native_posterior", "genotype_lr"]]
    compare("pair_sensitivity.tsv", "sensitivity.tsv", ["case", "gene_id", "pair_key"], fields)
    compare("recomputed_ledger.tsv", "selection.tsv", ["gene_id", "pair_key"],
            ["exome_native_ar_rank", "genome_native_ar_rank", "existing_research_tier", "phase"])
    checked += 1
    if (result / "predictions.csv").read_bytes() != (reference / "predictions.csv").read_bytes():
        mismatches.append({"table": "predictions.csv", "field": "exact_bytes"})
    require(all(sha(reference / name) == digest for name, digest in pinned.items()), "Verification reference changed during comparison")
    write_json(result / "verification.private.json", {"status": "PASS" if not mismatches else "FAIL",
        "comparisons": checked, "mismatches": mismatches,
        "reference_sha256": {p.name: sha(p) for p in sorted(reference.iterdir()) if p.is_file()},
        "csv_sha256": sha(result / "predictions.csv")})
    require(not mismatches, "Verification differs; see protected verification.private.json")
    return checked


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--result", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    args = p.parse_args()
    os.umask(0o077)
    print(json.dumps({"status": "PASS", "comparisons": verify(args.result, args.reference)}))
