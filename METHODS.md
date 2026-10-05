# Exomiser–LIRICAL prioritization with supplementary AVI annotation

Participant: **hideonbocchi** · [GitHub repository](https://github.com/hwangjw26/hgface_mva_hackaton)
Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1

## Approach and data

We analyzed the supplied singleton WGS VCF using seven proband HPO terms; family history was reviewed separately. The organizers disclosed a compound-heterozygous answer structure. Case data require controlled access; supporting databases and literature were public.

Exomiser 15.1.0, hg38/phenotype release 2602, supplied coding/splice and genome-preset candidates. We retained 120 native AR genotypes with two contributing heterozygous alleles, preserving preset-specific ranks without restricting discovery to known MVA genes.

LIRICAL 2.4.1 analyzed the complete VCF in a baseline and seven leave-one-HPO-out runs: non-global, pathogenicity threshold 0.8, default genotype penalties, no gene/disease restriction, and no `--sdwndv`. Complete JSON supplied disease context; 30 pairs had an AR context. Disease return and native allele eligibility did not establish pair causality.

Automated annotation and source-based review considered exact-allele ClinVar records (2026-08-22), transcripts, population frequency, biological literature, reads and phase. Resources included bcftools 1.22, VEP 115/gnomAD v4.1 annotations, HPO 2026-06-23 and HGNC 2026-04-07. Gene-level literature, including MVA evidence, was distinguished from allele-specific support. Missing evidence remained unknown.

## Supplementary AVI analysis

Google AlphaGenome Atlas SDK 0.9.0 returned AVI scores and 18 feature contributions for 201 SNVs on 2026-10-05. Authorized queries sent only GRCh38 chromosome, position, REF and ALT. Thirty-seven indels remained unscored.

Direct API quantiles were compared with **>=0.99379**, the [Qu and Hakonarson preprint](https://doi.org/10.64898/2026.09.26.754650) Table 2 threshold. Among 93 two-SNV pairs, four met it for both alleles, 18 for one and 71 for neither; 27 pairs were incomplete. AVI annotated all 120 pairs without filtering, reranking or changing submission values. Its shared AlphaMissense/conservation inputs and ClinVar-derived benchmark preclude assuming independent evidence. Quantiles are not pathogenicity probabilities; feature contributions are not functional experiments.

## Selection and limitations

Four pairs met the provisional mechanism-review tier; none met the stronger allele-evidence tier. Ordering used the better native Exomiser AR rank, then exome rank, genome rank and identifiers. All four are `primary` alternatives; no incidental findings are submitted. Hidden answers and leaderboard feedback were not used.

EPCR 1.00, 0.75, 0.50 and 0.25 encode ordinal rank, not calibrated causal probabilities. All phases remain unresolved; three pairs have strong read support for both alleles. These are research hypotheses. Structural, repeat and mosaic mechanisms were not comprehensively assessed.

## Reproducibility and resources

Release v1.1.1 reorganizes v1.1.0 code without changing decisions. It reconstructs the baseline from native outputs and declared annotation/review checkpoints, then optionally adds independently verified AVI annotations from hash-pinned responses. Cached replay requires no API key or network request. The service exposes no dataset snapshot version. Large upstream runs were not repeated for packaging.

Recorded wall times: Exomiser 5m47s/10m28s; eight LIRICAL runs 34m32s; targeted phase extension 86m38s. These exclude other work; total effort and cost were not recorded.

**AI assistance:** OpenAI Codex; user-reported plan/tier “GPT x30” (unverified); account training and retention settings unverified.

We acknowledge the family and challenge organizers.
