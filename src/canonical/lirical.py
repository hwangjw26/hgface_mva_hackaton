"""Exact source excerpts; see canonical_provenance.json."""
import collections
import datetime
import hashlib
import json
import math
import re
from urllib.parse import unquote

HIDDEN_TYPES = {"NO_VARIANTS_DETECTED_AD", "NO_VARIANTS_DETECTED_AR", "UNKNOWN"}

KNOWN_TYPES = HIDDEN_TYPES | {"ONE_P_OR_LP_CLINVAR_ALLELE_IN_AD", "TWO_P_OR_LP_CLINVAR_ALLELES_IN_AR",
    "ONE_DELETERIOUS_VARIANT_IN_AR", "HIGH_NUMBER_OF_OBSERVED_PREDICTED_PATHOGENIC_VARIANTS",
    "LIRICAL_GT_MODEL", "ONE_DELETERIOUS_CLINVAR_VARIANT_IN_AD", "TWO_DELETERIOUS_CLINVAR_VARIANTS_IN_AR"}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def complete_rows(results, displayed):
    require(len({r["diseaseId"] for r in results}) == len(results), "Duplicate native disease identity")
    ordered = sorted(results, key=lambda r: r["posttestProbability"], reverse=True)
    bounds = {}
    for i, r in enumerate(ordered, 1):
        bounds.setdefault(r["posttestProbability"], [i, i])[1] = i
    out = []
    for i, r in enumerate(ordered, 1):
        g = r["genotypeLR"]
        require(g["matchType"] in KNOWN_TYPES, "Unknown native match type")
        out.append({"full_disease_rank": i, "tie_rank_min": bounds[r["posttestProbability"]][0],
            "tie_rank_max": bounds[r["posttestProbability"]][1], "disease_id": r["diseaseId"],
            "gene_id": g["geneId"]["id"], "gene_symbol": g["geneId"]["symbol"],
            "native_pretest_probability": r["pretestProbability"], "native_composite_lr": r["compositeLR"],
            "native_posttest_probability": r["posttestProbability"], "native_genotype_lr": g["lr"],
            "native_genotype_match_type": g["matchType"], "native_filtered_tsv_rank": displayed.get(r["diseaseId"], "")})
    return out
