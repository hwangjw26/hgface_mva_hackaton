# Optional AVI annotation

This supplement applies the previously reviewed AlphaGenome Variant Impact method to the frozen candidate set. It does not introduce a selection rule, refit a threshold, infer phase or estimate causal probability.

## Method and outputs

The configuration pins the query, response JSONL and scorer metadata by SHA256. Exact GRCh38 chromosome, 1-based position, REF and ALT join scores to alleles; no gene or coordinate-only join is used. The official API retrieval used SDK 0.9.0 on 2026-10-05 and requested `AVI_SCORE` plus `AVI_SCORE_FEATURE_IMPORTANCE`. All 201 authorized SNVs returned scores; 37 indels remain `NON_SNV_OUTSIDE_PAPER_SCOPE`. All 238 gene–allele contexts and 120 pairs are retained.

The direct API `AVI_SCORE` quantiles layer is compared with `>=0.99379`, Table 2 of [Qu and Hakonarson, Applying AlphaGenome Variant Impact for SNV Prioritization](https://doi.org/10.64898/2026.09.26.754650), posted 2026-10-01. The abstract rounds this to 0.9938. The pair flag is BOTH_MEET, ONE_MEETS, NEITHER_MEETS or INCOMPLETE. It is annotation only, not a disease probability or exclusion rule. A descriptive missense-specific comparison in the private review is not an additional filter.

Raw scores, direct quantiles, derived PHRED, all 18 signed feature contributions and the largest absolute contributor are preserved. Contributions are approximate model explanations relative to a model baseline. Their sum is not forced to equal raw AVI; a large absolute contribution may be negative. See the [AlphaGenome Atlas paper](https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/alphagenome-atlas.pdf).

`avi/annotate_api.py` reuses key and I/O helpers from `avi/annotate_avi.py`. The only helper import change from the source run is a lazy pysam import: cached replay uses Python's standard library; the original reference/bulk command still requires pysam. `avi/verify_api_annotation.py` independently reconstructs observations without importing the producer. The included public-control fixture predates case analysis and contains no case query.

## Reproduction and provenance

Use the optional command in [README.md](../README.md). Paths in `avi.private.json` resolve against that configuration file, not the shell working directory. Modified or missing caches fail; they are never silently downloaded or filled with zero. The stage verifies existing baseline fields and submission binding and writes a protected `avi_stage_receipt.json`. The private `avi_verification_only/` files are immutable comparison references, never annotation inputs.

The frozen response includes each request time, endpoint, SDK version and returned allele. The service did not expose a dataset snapshot version, so response hashes and timestamps identify this retrieval. A new API request could return different data. The public code includes separately invoked acquisition scripts, but normal reproduction never invokes them. Keep the API key outside all archives.

## Interpretation limits

AVI shares AlphaMissense, conservation and consequence information with other predictors; independence from Exomiser or ClinVar cannot be assumed. The preprint's ClinVar-based evaluation is not prospective independent clinical validation. A high score does not establish disease relevance, pathogenicity of both alleles or trans; a low score does not establish benignity. No whole-genome noncoding discovery, indel scoring or tissue-specific functional assay was added. The four submitted pairs, their order, EPCR values, research tiers and unresolved phase remain unchanged.
