# Exomiser, LIRICAL and AVI analysis workflow

**Participant:** hideonbocchi

**Challenge:** Rare Disease, Real Kid — MVA Hackathon 2026, Track 1

**Code release:** 1.1.1 (organized from the validated v1.1.0 release)

This repository contains the scripts, decision rules and documentation used to
reconstruct an Exomiser-led variant-pair analysis, add LIRICAL disease context,
review recorded evidence and optionally annotate variants with AVI.
It produces research hypotheses, not clinical diagnoses.

[한국어 안내](README.ko.md) · [English methods report](METHODS.md) · [국문 보고서](METHODS.ko.md)

## Repository layout

```text
scripts/                 Command-line entry points
  normalization/         VCF normalization and source mapping
  annotation/            VEP execution and annotation parsing
  reads/                 Targeted read and phase workflows
  avi/                   AVI acquisition, cached annotation and verification
src/                     Shared reconstruction and evidence code
  canonical/             Preserved selection, queue and tier rules
config/                  Exomiser presets and output options
native/                  Java adapter for LIRICAL allele metadata
examples/                Public configuration templates with placeholders
tests/                   Synthetic tests and a public control fixture
docs/                    Input schemas, upstream steps and interpretation limits
provenance/              Original manifests, layout mapping and validation
```

Run commands from the repository root unless stated otherwise. The main Python
entry points also work when called by absolute path from another directory.

## What can be reproduced

The validated replay starts with **native Exomiser/LIRICAL outputs and declared
annotation, read and review checkpoints** in an authorized private bundle.
It recomputes evidence flags, research tiers, selection and the submission CSV.
Previous result tables are used only by the independent verifier.
The optional AVI stage reads hash-pinned cached responses; it makes no API calls.

Upstream VCF normalization, VEP, Exomiser/LIRICAL and targeted read scripts are
provided and documented in [docs/UPSTREAM.md](docs/UPSTREAM.md). Those large runs
were not repeated for this packaging release. Literature review and clinical
judgment are curated inputs, not conclusions generated automatically by software.
The later exploratory interval phase extension remains in the private historical
archive; this repository does not claim to publish every historical experiment.

## Requirements and private inputs

- Python 3.11 and bcftools 1.22 for frozen replay; no Python package install needed.
- Authorized v1.1.0 private input bundle, reference FASTA/indexes and pinned resources.
- A private writable output location and an existing scratch directory.
- AVI cached replay needs neither the AlphaGenome SDK nor an API key.

See [docs/INPUTS.md](docs/INPUTS.md). Adjust only machine-specific resource/tool
paths in the private configuration, retaining the expected hashes and versions.
Relative input paths resolve against the configuration file. Do not place the
private bundle inside this repository.

## Reproduce and verify

Replace `/private/...` with your authorized private paths. The output directory
must not already exist. Set `TMPDIR` to an existing approved scratch directory.

```bash
export TMPDIR=/private/scratch
python -B scripts/reproduce.py --config /private/reproduction.private.json --output /private/new-replay
python -B scripts/verify_reproduction.py --result /private/new-replay --reference /private/verification_only
```

To include AVI, use the following first command instead, then run the same verifier:

```bash
python -B scripts/reproduce.py --config /private/reproduction.private.json --avi-config /private/avi.private.json --output /private/new-replay
```

Both `verification.private.json` and, when enabled, `avi_stage_receipt.json` must
report PASS. Missing or changed caches fail rather than trigger a download.
Private outputs include `predictions.csv`, baseline tables and `avi/` annotations.
See [docs/AVI.md](docs/AVI.md) for the exact cutoff and interpretation boundaries.

## Analysis steps

| Step | Entry point | Purpose |
|---|---|---|
| 1 | `scripts/normalization/` | Normalize VCF and preserve source identity |
| 2 | `scripts/annotation/` | Run VEP and parse allele/transcript observations |
| 3 | `scripts/prepare_native.py` | Generate Exomiser and eight LIRICAL run commands |
| 4 | `scripts/prepare_metadata.py` | Generate exact LIRICAL allele metadata commands |
| 5 | `scripts/reads/` | Evaluate recorded targets and reconcile read/phase evidence |
| 6 | `scripts/reproduce.py` | Reconstruct candidates, evidence, tiers and CSV |
| 7 | `scripts/avi/` | Acquire responses separately or replay cached AVI annotations |
| 8 | `scripts/verify_reproduction.py` | Independently compare frozen outputs |

The two `prepare_*` scripts generate commands; they do not launch native analysis.
Full upstream arguments and resource requirements are in [docs/UPSTREAM.md](docs/UPSTREAM.md).

## Selection and interpretation

The frozen workflow retains all 120 native AR pairs. The submission rule includes
Tier1/Tier2 and orders them by native Exomiser AR rank; four Tier2 pairs satisfy it.
LIRICAL ranks and AVI scores do not alter this frozen selection. EPCR is an ordinal
encoding, not a calibrated probability. Phase remains unresolved. An AVI quantile
at least 0.99379 is a supplementary comparison, not a clinical pathogenicity call.
See [docs/EVIDENCE_REVIEW.md](docs/EVIDENCE_REVIEW.md).

## Run tests

```bash
export TEST_WORKDIR=/private/test-scratch
export BCFTOOLS=/tools/bcftools
python -B -W error::ResourceWarning scripts/run_tests.py
```

Tests use synthetic input and a public control, not patient data. Validation for
this layout is recorded in [provenance/VALIDATION.json](provenance/VALIDATION.json).
Original release records are retained in `provenance/v1.1.0/`.

## Submission and data boundaries

Use `hideonbocchi_exomiser_lirical_avi.csv` and `hideonbocchi_track1_report.md` for
the portal. The CSV stays private; [METHODS.md](METHODS.md) is the report source.
Anthropic Claude (Team plan for scientists) supported analysis; OpenAI Codex (ChatGPT Pro 20x) assisted with code and reporting. Account training and retention settings remain unverified.
Do not commit case inputs, exact candidate tables, raw reads, outputs, API caches
or keys. Repository publication does not submit predictions to the challenge.

See [delivery instructions](docs/DELIVERY.ko.md), [license](LICENSE.txt) and
[third-party notices](THIRD_PARTY.md). This repository retains the original
licenses; it does not redistribute Exomiser/LIRICAL binaries or reference databases.

The reviewed delivery CSV adds candidate-specific rationale and limitations only in `notes`. All other fields and row order match the frozen replay output. Replay and its byte-comparison verifier still produce/check the original CSV with `Phase unresolved.`; the reviewed notes are an editorial delivery step, not a selection rule. Keep the reviewed delivery CSV for portal upload.
