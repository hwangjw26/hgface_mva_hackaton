"""Keep disease-level results and exact-allele native eligibility separate."""
import collections
import json
from pathlib import Path
from canonical.condition_ids import condition_ids
from canonical.lirical import complete_rows, HIDDEN_TYPES
from evidence import inheritance
from io_utils import index, require, rows, write_tsv


def eligible(context, route):
    if not context["mapped"]:
        return False
    if route == "score":
        return (context["score_eligible_allele_count"] or 0) > 0
    if route == "plp":
        return (context["plp_allele_count"] or 0) > 0
    return eligible(context, "score") or eligible(context, "plp")


def best_fields(prefix, row):
    return {prefix + suffix: row[source] if row else "" for suffix, source in [
        ("disease", "disease_id"), ("full_disease_rank", "full_disease_rank"),
        ("native_posterior", "native_posttest_probability"), ("genotype_lr", "native_genotype_lr")]}


def reconcile(cfg, pairs, alleles, out):
    m = json.loads(Path(cfg["metadata"]).read_text())
    require(m["threshold"] == 0.8 and m["gene_retention_filter_used"] is False,
            "Unexpected native metadata settings")
    contexts = index(m["contexts"], ["gene", "key"])
    observed = index(m["allele_observations"], ["key"])
    expected = {(a["gene_id"], a["allele_key"]) for a in alleles}
    require(set(contexts) == expected, "Native metadata context coverage mismatch")
    require({k[0] for k in observed} == {a["allele_key"] for a in alleles}, "Native metadata allele coverage mismatch")
    ar = inheritance(cfg["hpoa"])
    by_gene, by_disease = {}, {}
    require(set(cfg["cases"]) == {"baseline", *[f"loo_{i:02d}" for i in range(1, 8)]}, "Eight frozen cases required")
    counts = {}
    for case, paths in cfg["cases"].items():
        source = json.loads(Path(paths["json"]).read_text())["analysisResults"]
        native = rows(paths["tsv"], "!")
        displayed = {r["diseaseCurie"]: int(r["rank"]) for r in native}
        require(len(displayed) == len(native), "Duplicate displayed disease")
        complete = complete_rows(source, displayed)
        require([r["disease_id"] for r in complete if r["native_genotype_match_type"] not in HIDDEN_TYPES] ==
                [r["diseaseCurie"] for r in native], "Native report filter/order mismatch")
        write_tsv(out / (case + ".full_diseases.tsv"), complete)
        by_gene[case] = collections.defaultdict(list)
        for d in complete:
            by_gene[case][d["gene_id"]].append(d)
        by_disease[case] = {d["disease_id"]: d for d in complete}
        counts[case] = len(complete)
    best_any = {g: rr[0] for g, rr in by_gene["baseline"].items()}
    best_ar = {g: next((r for r in rr if r["disease_id"] in ar), None) for g, rr in by_gene["baseline"].items()}
    allele_output = []
    for a in alleles:
        ctx, v = contexts[(a["gene_id"], a["allele_key"])], observed[(a["allele_key"],)]
        target = best_ar.get(a["gene_id"])
        target_id = target["disease_id"] if target else ""
        ids = condition_ids(json.loads(a["indexed_CLNDISDB_json"]))
        condition = ("NO_RETURNED_AR_TARGET" if not target_id else
                     "EXACT_ID_OVERLAP_NOT_ASSERTION_ADJUDICATION" if target_id in ids else
                     "NO_COMPARABLE_CONDITION_ID" if not ids else "NO_EXACT_ID_OVERLAP_NOT_EXCLUSION")
        allele_output.append({**a, "t17_condition_ids_json": json.dumps(ids),
            "t20_expected_gene_mapped": ctx["mapped"], "t20_native_passed_filter": v["passed"],
            "t20_native_optional_score_present": v["score"] is not None, "t20_native_optional_score": v["score"],
            "t20_native_clinvar_2602": v["clinvar"], "t20_native_alt_count": v["alt_count"],
            "t20_score_route_allele_count": ctx["score_eligible_allele_count"], "t20_plp_route_allele_count": ctx["plp_allele_count"],
            "t20_any_native_route_eligible": eligible(ctx, "any"), "t20_best_ar_target": target_id,
            "t20_condition_id_context": condition})
    pair_output, sensitivity, disease_contexts = [], [], []
    for p in pairs:
        gene = p["gene_id"]
        cc = [contexts[(gene, p[f"allele_{n}_key"])] for n in [1, 2]]
        mapped, n = all(x["mapped"] for x in cc), sum(eligible(x, "any") for x in cc)
        state = ("EXPECTED_GENE_ANNOTATION_MISSING" if not mapped else
                 "BOTH_NATIVE_ROUTES_ELIGIBLE_NOT_CLINICAL_CONFIRMATION" if n == 2 else
                 "ONE_NATIVE_ROUTE_ELIGIBLE" if n == 1 else "NEITHER_NATIVE_ROUTE_ELIGIBLE")
        pair_output.append({**p, "t20_gene_returned": gene in best_any, "t20_ar_returned": bool(best_ar.get(gene)),
            **best_fields("t20_best_any_", best_any.get(gene)), **best_fields("t20_best_ar_", best_ar.get(gene)),
            "t20_both_expected_gene_mapped": mapped, "t20_both_score_route_eligible": all(eligible(x, "score") for x in cc),
            "t20_both_plp_route_eligible": all(eligible(x, "plp") for x in cc), "t20_both_any_native_route_eligible": n == 2,
            "t20_native_pair_evidence_state": state, "t20_native_pair_posterior": "NOT_COMPUTED"})
        for case in cfg["cases"]:
            rr = by_gene[case].get(gene, [])
            any_row, ar_row = rr[0] if rr else None, next((r for r in rr if r["disease_id"] in ar), None)
            fixed_any = by_disease[case].get(best_any[gene]["disease_id"]) if gene in best_any else None
            fixed_ar = by_disease[case].get(best_ar[gene]["disease_id"]) if best_ar.get(gene) else None
            sensitivity.append({"case": case, "review_id": p["review_id"], "pair_key": p["pair_key"], "gene_id": gene,
                **best_fields("best_any_", any_row), **best_fields("best_ar_", ar_row),
                **best_fields("fixed_baseline_any_", fixed_any), **best_fields("fixed_baseline_ar_", fixed_ar)})
            disease_contexts.extend({"case": case, "review_id": p["review_id"], "pair_key": p["pair_key"],
                "is_positive_ar_annotation": d["disease_id"] in ar, **d} for d in rr)
    write_tsv(out / "allele_evidence.tsv", allele_output)
    write_tsv(out / "pair_sensitivity.tsv", sensitivity)
    write_tsv(out / "all_pair_disease_contexts.tsv", disease_contexts,
              ["case", "review_id", "pair_key", "is_positive_ar_annotation", *list(complete[0])])
    return pair_output, counts
