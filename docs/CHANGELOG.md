# Changes

## 1.1.1 — 2026-10-05

- Organize command-line entry points in `scripts/`, shared code in `src/`, tests in `tests/`, and supporting documents in `docs/`.
- Resolve config, Java adapter and AVI script paths relative to the repository, independent of the caller's working directory.
- Add English/Korean repository guides and participant name hideonbocchi.
- Preserve original release manifests under `provenance/v1.1.0/` and document old/new paths in `provenance/layout_map.json`.
- Pass 30 tests, 25,881 baseline comparisons and 2,024 AVI checks. Submission CSV and all four deterministic AVI tables match v1.1.0 byte-for-byte.

The scientific selection code, candidates, order, EPCR values, tiers and phase states are unchanged. This release replays existing checkpoints; it does not rerun large native analyses or acquire new AVI responses.
