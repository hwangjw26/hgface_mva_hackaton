# Reviewed ledger contract

The recomputed ledger is a private output of the complete replay, not a database of clinically confirmed variants. It retains all 120 native two-heterozygote AR pair hypotheses. Previous final tiers are verification references only. The standalone CSV renderer can also read the ledger schema below, but the complete replay recomputes these fields from source observations.

| Column | Source and meaning |
| --- | --- |
| review_id | Stable private review identifier |
| gene_id | Same-gene hypothesis identifier |
| pair_key | Stable canonical physical pair identity; used only as a deterministic final tie break |
| allele_1_key, allele_2_key | Normalized GRCh38 chromosome:position:REF:ALT; source-VCF and REF checked locally |
| existing_research_tier | Frozen research review state; not ACMG/AMP classification |
| exome_native_ar_rank | Native Exomiser coding/splice AR rank; empty if absent |
| genome_native_ar_rank | Native Exomiser genome AR rank; empty if absent |
| phase | Recorded phase disposition; this submission selects only unresolved-phase pairs |

Tier1 requires completed disease/mechanism-appropriate review of each exact allele. ClinVar P/LP is one possible route, not a mandatory requirement. Tier2 means technically acceptable and compatible with an AR disease connection, with protein/splice/clinical clues but incomplete allele/disease support. Tier3 retains exploratory hypotheses. Hold and exact-pair contradiction states remain in the ledger rather than disappearing.

The original evidence review separated: target condition and inheritance, transcript/consequence, population frequency or missingness, ClinVar condition/review/conflict information, gene versus exact-allele literature, read/call support, phase, supporting evidence, counterevidence and unresolved questions. None of those are independent calibrated likelihood factors by default. The private bundle preserves the selected exact-allele source dossiers and biological-evidence records.

The complete replay recalculates research tiers with the original queue/tier functions. A gene association or P/LP label alone does not establish reviewed allele support. The renderer then applies the approved order: tier, best native AR rank across two presets, exome AR rank, genome AR rank, pair identity, gene identity. Missing rank is never zero. No LIRICAL posterior is used as pair EPCR.
