"""Reassemble observations without importing a previous tier or selection ledger."""
import collections
import json
import math
import re
from canonical import handoff_core as core
from canonical.identity import allele_key, clinvar_class, pair_key
from canonical.observations import frequency_flag, physical_phase, read_consensus
from io_utils import annotated_records, index, require, rows, write_tsv


def inheritance(hpoa):
    return {r["database_id"] for r in rows(hpoa, "#")
            if r["aspect"] == "I" and r["hpo_id"] == "HP:0000007"
            and "NOT" not in re.split(r"[|;, ]+", r["qualifier"])}


def read_observations(cfg):
    source = index(rows(cfg["read_input_alleles"]), ["allele_id"])
    ids = {k[0]: allele_key(r["chrom"], r["pos1"], r["ref"], r["alt"]) for k, r in source.items()}
    require(len(set(ids.values())) == len(ids), "Duplicate read-input allele")
    observations, phases = collections.defaultdict(list), collections.defaultdict(list)
    for r in rows(cfg["read_alleles_initial"]):
        observations[ids[r["allele_id"]]].append({"label": r["read_disposition"], "source": "initial_reads", "record": r})
    for r in rows(cfg["read_alleles_extended"]):
        observations[r["allele_key"]].append({"label": r["label_of_record"], "source": "extended_reads", "record": r})
    for r in rows(cfg["read_pairs_initial"]):
        phases[pair_key(ids[r["allele_1_id"]], ids[r["allele_2_id"]])].append(
            {"phase": r["phase_label"], "source": "initial_reads", "record": r})
    for r in rows(cfg["read_pairs_extended"]):
        phases[r["pair_key"]].append({"phase": r["phase_label_of_record"], "source": "extended_reads", "record": r})
    return observations, phases


def allele_observation(key, v, lineage, reads):
    fmt, annotations = v["format"], v["annotations"]
    a = {"allele_key": key, **lineage, "normalized_filter": v["filter"],
         **{"normalized_" + f.lower(): fmt.get(f, "") for f in ["GT", "AD", "DP", "GQ"]},
         **{f.lower(): fmt.get(f, "") for f in ["PGT", "PID", "PS"]}}
    representative = min(annotations, key=core.transcript_key) if annotations else {}
    invalid_af = False
    for field, value in representative.items():
        if field.endswith("_AF") and field.startswith(("gnomADe", "gnomADg")):
            for token in re.split(r"[&,]", value):
                if token in core.MISSING:
                    continue
                try:
                    number = float(token)
                    invalid_af |= not math.isfinite(number) or not 0 <= number <= 1
                except ValueError:
                    invalid_af = True
    flags = core.derive_flags(a, {}, representative)
    a.update({f: flags[f] for f in ["genotype_class", "allele_balance", "high_confidence_call", "max_gnomad_af", "gnomad_af_state"]})
    vals = lambda field: sorted({r.get(field, "") for r in annotations} - {"", "."})
    current = ["ClinVarCurrent", *["ClinVarCurrent_" + f for f in ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNDISDB", "CLNSIGCONF"]]]
    a["current_clinvar_state"] = "PRESENT" if any(vals(f) for f in current) else "MISSING"
    a["embedded_clinvar_state"] = "PRESENT" if any(vals(f) for f in ["CLIN_SIG", "SOMATIC", "PHENO"]) else "MISSING"
    a["embedded_CLIN_SIG_json"] = json.dumps(vals("CLIN_SIG"))
    for f in ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNDISDB", "CLNSIGCONF"]:
        a["current_" + f + "_json"] = json.dumps(vals("ClinVarCurrent_" + f))
        a["indexed_" + f + "_json"] = a["current_" + f + "_json"]
    a["clinvar_class"] = clinvar_class(a)
    sig = " ".join(vals("ClinVarCurrent_CLNSIG") + vals("ClinVarCurrent_CLNSIGCONF")).lower().replace(" ", "_")
    strong = bool(re.search(r"(?<![a-z])(?:likely_)?pathogenic(?![a-z])", sig)) and bool(re.search(r"(?<![a-z])(?:likely_)?benign(?![a-z])", sig))
    problems = ["CLINVAR_OPPOSING_STRONG_ASSERTIONS"] if strong else []
    a["read_label"], evaluable, issue = read_consensus(reads)
    if issue:
        problems.append(issue)
    a["read_contradiction_evaluable"] = str(evaluable).lower()
    a["read_sources_json"] = json.dumps(reads, sort_keys=True)
    a["frequency_flag"], issue = frequency_flag(a["max_gnomad_af"], 0.01)
    if invalid_af:
        a["frequency_flag"], issue = "UNKNOWN", "INVALID_AF"
    if issue:
        problems.append(issue)
    a["af_plp_conflict"] = str(a["frequency_flag"] == "ABOVE_BROAD_SCREEN" and a["clinvar_class"] == "PLP_NONCONFLICTING").lower()
    if a["af_plp_conflict"] == "true":
        problems.append("AF_VS_PLP")
    a["evidence_conflicts"] = ";".join(sorted(problems))
    return a


def reconcile_phase(a, b, previous):
    computed = physical_phase(a, b)
    historical = [r for r in previous if r["record"].get("phase_source") == "VCF_PGT_PID_SHARED_BLOCK"]
    labels = {r["phase"] for r in previous if r["phase"] in {"LOCAL_CIS_SUPPORTED", "LOCAL_TRANS_SUPPORTED"} and r not in historical}
    if computed != "UNRESOLVED":
        labels.add(computed)
    conflict = len(labels) > 1 or any(r["phase"] != computed for r in historical)
    return "PHASE_CONFLICT" if conflict else next(iter(labels)) if labels else "UNRESOLVED"


def assemble(cfg, pairs, out):
    keys = {p[f"allele_{i}_key"] for p in pairs for i in [1, 2]}
    records = annotated_records(cfg["annotated_vcf"], keys)
    lineage = index(rows(cfg["source_lineage"]), ["allele_key"])
    require(keys <= {k[0] for k in lineage}, "Source normalization lineage incomplete")
    read_sources, phase_sources = read_observations(cfg)
    aa = {key: allele_observation(key, records[key], lineage[(key,)], read_sources[key]) for key in sorted(keys)}
    gene_ens, ens_ncbi = collections.defaultdict(set), collections.defaultdict(set)
    for h in rows(cfg["hgnc"]):
        if h["status"] == "Approved" and h["entrez_id"].isdecimal() and h["ensembl_gene_id"]:
            gid, ens = "NCBIGene:" + h["entrez_id"], h["ensembl_gene_id"].split(".")[0]
            gene_ens[gid].add(ens)
            ens_ncbi[ens].add(gid)
    ar = inheritance(cfg["hpoa"])
    diseases = collections.defaultdict(list)
    for d in rows(cfg["gene_diseases"], "#"):
        gid = d["ncbi_gene_id"] if d["ncbi_gene_id"].startswith("NCBIGene:") else "NCBIGene:" + d["ncbi_gene_id"]
        diseases[gid].append(d)
    biology = collections.defaultdict(list)
    for b in rows(cfg["biology_gene_records"]):
        biology[b[cfg["biology_gene_field"]]].append(b)
    contexts, tx = {}, []
    for gid, key in sorted({(p["gene_id"], p[f"allele_{i}_key"]) for p in pairs for i in [1, 2]}):
        ens = gene_ens[gid]
        conflict = "MULTIPLE_ENSEMBL_MAPPINGS" if len(ens) > 1 else "HGNC_ENSEMBL_CONFLICT" if any(ens_ncbi[e] != {gid} for e in ens) else ""
        candidates = [r for r in records[key]["annotations"] if r.get("Feature_type") == "Transcript" and r.get("Gene", "").split(".")[0] in ens] if not conflict else []
        best = min(candidates, key=core.transcript_key) if candidates else {}
        terms = sorted({t for r in candidates for t in r.get("Consequence", "").split("&") if t})
        contexts[(gid, key)] = {"gene_id": gid, **aa[key], "ensembl_ids": ";".join(sorted(ens)),
            "gene_bridge_conflicts": conflict, "annotation_rows": len(candidates),
            "gene_consequences_json": json.dumps(terms),
            "gene_protein_splice": str(bool(set(terms) & core.PROTEIN_OR_SPLICE_TERMS)).lower(),
            "selected_protein_splice": str(bool(set(best.get("Consequence", "").split("&")) & core.PROTEIN_OR_SPLICE_TERMS)).lower(),
            "selected_transcript": best.get("Feature", ""), "selected_consequence": best.get("Consequence", ""),
            "selected_mane": best.get("MANE_SELECT", ""), "selected_hgvsc": best.get("HGVSc", ""),
            "selected_hgvsp": best.get("HGVSp", ""),
            "transcript_consequence_discordance": str(any(r["Consequence"] != best["Consequence"] for r in candidates)).lower()}
        tx.extend({"gene_id": gid, "allele_key": key, **r} for r in candidates)
    enriched = []
    for p in pairs:
        gg = [contexts[(p["gene_id"], p[f"allele_{i}_key"])] for i in [1, 2]]
        previous = phase_sources[p["pair_key"]]
        associations = diseases[p["gene_id"]]
        r = {**p, "source_reconciled": "true",
             "gene_id_resolved": str(bool(re.fullmatch(r"NCBIGene:\d+", p["gene_id"]))).lower(),
             "gene_bridge_conflicts": ";".join(sorted({g["gene_bridge_conflicts"] for g in gg} - {""})),
             "positive_ar_association": str(any(d["disease_id"] in ar for d in associations)).lower(),
             "all_disease_associations_json": json.dumps(associations, sort_keys=True),
             "biological_gene_records_json": json.dumps(biology[p["gene_symbol"]], sort_keys=True),
             "biology_scope": "Gene-level literature; no automatic allele pathogenicity or tier promotion",
             "audit_phase": reconcile_phase(*gg, previous), "phase_sources_json": json.dumps(previous, sort_keys=True)}
        for i, g in enumerate(gg, 1):
            for f in ["normalized_gt", "high_confidence_call", "frequency_flag", "clinvar_class", "evidence_conflicts", "read_label", "read_contradiction_evaluable", "af_plp_conflict", "gene_protein_splice", "selected_protein_splice", "annotation_rows", "selected_consequence", "transcript_consequence_discordance"]:
                r[f"audit_{i}_" + f] = g[f]
        enriched.append(r)
    write_tsv(out / "source_alleles.tsv", list(aa.values()))
    write_tsv(out / "transcripts.tsv", tx)
    return enriched, list(contexts.values())
