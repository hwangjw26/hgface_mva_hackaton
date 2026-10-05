# Upstream execution and review boundaries

The verified release starts from native Exomiser/LIRICAL outputs and source observations. The programs below regenerate those inputs. Their large WGS, annotation and FASTQ computations were not repeated for this release. The protected `upstream_inputs/` directory supplies the historical parameter contracts, recorded commands and read target lists; `analysis_sources.PRIVATE.zip` preserves the executed workflow source. Paths in historical contracts must be rebound to authorized local resources while preserving versions, hashes and analysis parameters.

Use fresh private output directories, existing environments and an existing `TMPDIR`. Measure total local CPU utilization before substantial runs; this project's local-first threshold is 70%, with at most 48 cores per task. Never write generated indexes beside shared/raw resources.

## Exomiser and LIRICAL

The supplied `prepare_native.py` generates commands and private HPO/sample configuration from `examples/native_paths.json` after replacing placeholders. It does not launch analysis.

```bash
python -B scripts/prepare_native.py --config /private/native_paths.json --output /private/native-run
bash /private/native-run/run_native.sh
```

Use Exomiser 15.1.0 with 2602 hg38/phenotype data and ReMM v0.4; both original preset YAMLs and output options are included. LIRICAL is 2.4.1, with its frozen resource directory. Seven patient HPO terms yield baseline plus seven leave-one-out cases. The eighth clinical item was family history and remains outside the patient HPO list. Review extraction from the clinical document remains explicit, source-linked input.

LIRICAL uses the full supplied VCF, no target gene/disease restriction, no `--global`, threshold 0.8, default genotype penalty and omitted `--sdwndv`. Input validation `STRICT` is distinct from strict genotype handling. Native output commands request JSON/TSV/HTML and all variant display; the replay reads all JSON disease rows.

The native pair extractor can be run before downstream evidence to create metadata targets:

```python
import json, os
from pathlib import Path
import sys
sys.path.insert(0, "src")
import native_pairs
os.umask(0o077)
# Use absolute paths to the two gene/variant TSVs, bcftools and GRCh38 FASTA.
config = json.loads(Path("/private/native_pair_paths.json").read_text())
output = Path("/private/native-pairs")
output.mkdir(mode=0o700)  # must be new
native_pairs.run(config, output)
```

`native_pair_paths.json` has six keys: `exome_genes,exome_variants,genome_genes,genome_variants,bcftools,reference_fasta`. Require pre-existing reference indexes and verify source/tool hashes before this upstream call; the complete replay enforces those checks itself.

Generate exact LIRICAL metadata using the included Java adapter, which invokes native allele eligibility functions:

```bash
python -B scripts/prepare_metadata.py --native-commands /private/native-run/commands.json --pairs /private/native-pairs/native_pairs.tsv --sample-id SAMPLE --java /tools/java --javac /tools/javac --jar /tools/LIRICAL.jar --output /private/metadata-run
bash /private/metadata-run/run_metadata.sh
```

The Java runtime/compiler must match the installed LIRICAL distribution, including its `lib/` directory. This full-VCF parser pass is separate from disease scoring.

## Normalization and annotation

The original parameterized CLI programs are included under `scripts/normalization/` and `scripts/annotation/`. Invoke their `--help` for exact arguments. The private provenance JSONs retain original command arrays, resource hashes and record counts.

Run in order:

1. `prepare_normalization_input.py`: singleton VCF to tagged primary records and original allele map.
2. `run_bcftools_normalization.py`: GRCh38 reference checking, multiallelic splitting and normalization.
3. `build_normalized_source_link.py`: retain source-record identity, ALT index/count and normalized allele key.
4. `run_vep_annotation.py`: offline VEP 115/GRCh38 cache and frozen 2026-08-22 ClinVar VCF. It retains allele identity and CSQ observations.
5. `parse_vep_annotations.py`: full allele/transcript source tables when needed.

The replay can consume the full annotated VCF or a lossless exact-allele subset. It requires every selected allele exactly once. To form its lineage TSV, project `normalized_allele_key` (renamed `allele_key`), `source_record_key`, `source_alt_index` and `source_alt_count` from the normalized source-link table; keep one consistent source record per normalized allele and reject ambiguous mappings. No previous quality, tier or pathogenicity flag belongs in this input.

## Targeted read observations

Copy all files from `scripts/reads/` into a fresh private `<run-dir>/scripts/` and create `<run-dir>/data/sensitive/` and `<run-dir>/reports/`. The C prefilter source **must** be at `<run-dir>/scripts/lossless_seed_prefilter.c`; the prefilter stage compiles the binary there using the existing gcc/zlib toolchain.

The private upstream inputs include initial and extended target handoffs, the probe exclusion list, the FASTQ pair/integrity contracts and read evidence contract. Rebind their file/tool/reference paths to existing local resources; keep expected hashes. Copy the extended handoff to `<run-dir>/data/sensitive/read_target_pairs.tsv` and put `probe_excluded_alleles.tsv` next to it. It preserves previously declared low-complexity exclusions as noncalls.

For each fresh read run, call `run_targeted_fastq_validation.py` (initial targets) or `run_task09_read_validation.py` (extended targets) with `--run-dir`, `--candidate-handoff`, `--fastq-contract`, `--fastq-integrity`, `--contract` and `--reference`. Run `--stage` in this order: `prepare, smoke, prefilter, full, call, evidence`. Respect every smoke/coverage gate; do not bypass a failed gate. The initial historical run expects 21 target pairs/28 alleles. The extended wrapper derives target cardinality from its frozen handoff.

Afterward, `python -B EXTENDED_RUN/scripts/reconcile_read_phase.py --run-dir EXTENDED_RUN --source-alleles SOURCE_ALLELES --initial-read-dir INITIAL_RUN/data/sensitive` recreates the extended read/phase ledgers. The source allele table and exact original target list are private inputs. The main replay independently revalidates physical VCF phase and conflicting observations.

The later targeted interval phase extension is preserved in the private historical source and supporting review. Its result remains unresolved; this release does not infer trans from two heterozygous calls or rerun that investigation.

## Reviewed evidence

The biological literature matrix and completed-review ledger are curated source records. They cannot be regenerated as clinical conclusions by rerunning software. Preserve citations, evidence level, reviewer/date/rationale where present and explicit missing evidence. The frozen completed-review array is empty; no allele pair is silently promoted to Tier1.

Source copies and the limited portability edits are recorded in [upstream_provenance.json](../provenance/v1.1.0/upstream_provenance.json). The Java metadata adapter compiles and seven CLI entry points pass help checks. Read-ledger reconciliation was also replayed from recorded read observations: both compressed allele/phase ledgers match byte-for-byte and rerunning refuses existing outputs. Raw read retrieval/calling, full native scoring and full annotation were not repeated.

## AVI acquisition and frozen replay

The completed AVI supplement uses cached Google AlphaGenome Atlas SDK 0.9.0 responses. Acquisition code is included for traceability; it is not executed by `reproduce.py`. The response includes scorer identities, exact returned alleles and request times. A fresh upstream collection requires a separately authorized query, existing SDK environment and external key file:

```bash
python -B scripts/avi/inspect_atlas_metadata.py --key-file /secrets/key --output /private/new-metadata.json
python -B scripts/avi/fetch_atlas.py --key-file /secrets/key --query /private/authorized-query.tsv --query-sha256 AUTHORIZED_QUERY_SHA256 --output /private/new-responses.jsonl
```

The acquisition adapter enforces the frozen 201-SNV scope and four request columns. Do not rerun it just to reproduce this submission, and do not broaden its scope or overwrite the frozen cache. Provider data may change without an exposed snapshot version. A new cache would require new provenance and comparison rather than silently replacing existing evidence. The optional local stage is described in [AVI.md](AVI.md); its source-copy hashes and limited portability edits are recorded in [avi_provenance.json](../provenance/v1.1.0/avi_provenance.json).
