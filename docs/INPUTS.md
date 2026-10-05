# Explicit checkpoint inputs

`reproduction.private.json` binds each scientific input to a path and SHA-256 digest. It also pins bcftools, its version output, GRCh38 FASTA and existing FAI/GZI indexes. Relative paths resolve against the configuration directory. Never replace a file without documenting the new source/version and reviewing any changed result.

| Input | Required origin and content |
| --- | --- |
| Native exome/genome gene and variant TSVs | Exomiser 15.1.0, data 2602; preserve every returned mode and contributing flag |
| Annotated VCF | Normalized singleton VCF with VEP 115 CSQ, gnomAD v4.1 and frozen current ClinVar fields; retain every requested allele |
| Source lineage TSV | `allele_key,source_record_key,source_alt_index,source_alt_count`; primitive normalization provenance only |
| HGNC, HPOA, gene–disease associations | Frozen public snapshots; positive AR inheritance only, excluding NOT |
| Biological gene records | Source-linked curated literature matrix; gene key supplied in `biology_gene_field`; no score multiplier |
| Initial read targets/alleles/pairs and extended read ledgers | Exact recorded technical observations with retrieval/evaluability, allele identifiers and phase source; unassessed is not reference-only |
| Completed reviews JSON | Array of explicitly completed allele/condition review records; the frozen run has no such completed records |
| LIRICAL JSON and native TSV | Baseline plus seven leave-one-term-out cases; full original VCF, unrestricted gene/disease analysis |
| Native metadata JSON | Output of the supplied Java adapter using LIRICAL's own parser and Gene2Genotype functions; exact expected gene–allele contexts |

The frozen LIRICAL run is non-global, threshold 0.8, default genotype penalty, ENSEMBL transcripts, strict VCF input validation, no target-disease/gene restriction and omitted `--sdwndv`. “Full” refers to unrestricted input and full JSON reporting; it does not mean the `--global` flag. The full JSON contains disease rows hidden by the native TSV's deleterious-variant display filter.

Research rules are fixed in the source: DP≥10, GQ≥20, alternate depth≥3, heterozygous allele balance 0.20–0.80; broad AF screen 0.01. These are triage rules, not disease-specific pathogenicity cutoffs. The ClinVar aggregate class is not an ACMG clinical classification.

Verification references are kept in a separate directory. No reference tier, rank, selected-row list or final CSV is read by `reproduce.py`. Native scores, disease probabilities, allele eligibility, technical support and reviewed biological evidence remain distinct.

Supplementary private evidence preserves the later condition/allele review and targeted phase investigation. Those reports qualify interpretation; they do not silently replace computed phase or promote a tier.

## Optional AVI cache

`avi.private.json` independently pins three files in `avi_inputs/`: the authorized four-column SNV query, cached response JSONL and scorer metadata. Schema version 1 fixes cache-only execution, SDK 0.9.0, direct API quantiles, operator `>=` and cutoff 0.99379. No secret key is an input. The response hash pins the unavailable server snapshot version operationally; it does not identify a vendor model release.

`avi_verification_only/` contains four reference annotation files for delivery comparison. They are never read by the annotation stage. The baseline's newly generated `predictions.csv` is read only to label the current-four annotation; it is not a prior decision input or modified by AVI. Candidate-specific supplementary prose is explicitly curated and stored in `supporting_review/avi/`, separate from computed annotations.
