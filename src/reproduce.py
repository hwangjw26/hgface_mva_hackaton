"""Reproduce the frozen Exomiser-led analysis from declared source checkpoints."""
import argparse
import collections
import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import traceback
import build_submission
import native_pairs
from canonical import queues, tiers
from evidence import assemble
from io_utils import require, sha, write_json, write_tsv
from reconcile_lirical import reconcile


def run(cfg, out):
    require(not out.exists(), "Output already exists; choose a fresh directory")
    out.mkdir(mode=0o700, parents=True)
    try:
        before = {}
        for name, item in cfg["inputs"].items():
            path = Path(item["path"])
            require(path.is_file(), "Declared input missing: " + name)
            before[name] = sha(path)
            require(before[name] == item["sha256"], "Declared input checksum mismatch: " + name)
        require(len(before) > 0, "An explicit input manifest is required")
        resource_before = {k: sha(v["path"]) for k, v in cfg["resources"].items()}
        require(all(resource_before[k] == v["sha256"] for k, v in cfg["resources"].items()), "Resource checksum mismatch")
        require(cfg["resources"]["reference_fasta"]["path"] == cfg["native"]["reference_fasta"] and
                cfg["resources"]["bcftools"]["path"] == cfg["native"]["bcftools"], "Normalization resource binding mismatch")
        fasta = cfg["native"]["reference_fasta"]
        for suffix in ([".fai", ".gzi"] if fasta.endswith(".gz") else [".fai"]):
            require(any(v["path"] == fasta + suffix for v in cfg["resources"].values()), "Reference index not pinned")
        version = subprocess.run([cfg["native"]["bcftools"], "--version"], check=True, capture_output=True, text=True).stdout
        require(version == cfg["bcftools_version"], "bcftools version mismatch")
        # Every scientific input must be in the manifest. Executables and FASTA are resources.
        bound = {str(Path(item["path"]).resolve()) for item in cfg["inputs"].values()}
        for section in ["native", "evidence"]:
            for name, value in cfg[section].items():
                if name not in {"bcftools", "reference_fasta", "biology_gene_field"}:
                    require(str(Path(value).resolve()) in bound, "Unmanifested analysis input: " + name)
        for value in [cfg["lirical"]["metadata"], cfg["lirical"]["hpoa"], cfg["completed_reviews"],
                      *[v for case in cfg["lirical"]["cases"].values() for v in case.values()]]:
            require(str(Path(value).resolve()) in bound, "Unmanifested LIRICAL/review input")
        pairs, presets = native_pairs.run(cfg["native"], out)
        enriched, alleles = assemble(cfg["evidence"], pairs, out)
        classified = tiers.classify_all(queues.classify_all(enriched), json.loads(Path(cfg["completed_reviews"]).read_text()))
        integrated, disease_counts = reconcile(cfg["lirical"], classified, alleles, out)
        write_tsv(out / "integrated_pairs.tsv", integrated)
        decisions = [{**r, "existing_research_tier": r["v2_tier"], "phase": r["v2_phase"]} for r in integrated]
        write_tsv(out / "recomputed_ledger.tsv", decisions)
        selected = build_submission.select(decisions)
        result = {"pairs": len(pairs), "preset_rows": len(presets), "allele_contexts": len(alleles),
                  "diseases_by_case": disease_counts, "returned_pairs": sum(r["t20_gene_returned"] for r in integrated),
                  "ar_pairs": sum(r["t20_ar_returned"] for r in integrated),
                  "tiers": dict(collections.Counter(r["v2_tier"] for r in integrated)), "selected_pairs": len(selected)}
        require(result == cfg["expected"], "Frozen analysis counts changed; inspect protected outputs")
        with (out / "predictions.csv").open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=build_submission.COLUMNS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(build_submission.render(selected, cfg["proband_id"]))
        require(all(sha(item["path"]) == before[name] for name, item in cfg["inputs"].items()), "Input changed during replay")
        require(all(sha(v["path"]) == resource_before[k] for k, v in cfg["resources"].items()), "Resource changed during replay")
        write_json(out / "run_receipt.json", {"status": "PASS", "summary": result, "input_sha256": before,
            "python": platform.python_version(), "config_sha256": cfg["_config_sha256"],
            "resource_sha256": resource_before, "bcftools_version": version,
            "code_sha256": {str(p.relative_to(Path(__file__).resolve().parents[1])): sha(p)
                            for p in Path(__file__).resolve().parents[1].rglob("*") if p.suffix in {".py", ".java", ".yml"}},
            "upstream_native_annotation_fastq_rerun": False})
    except Exception:
        with (out / "failure.log").open("x") as handle:
            traceback.print_exc(file=handle)
        raise SystemExit("Replay failed; details are in the protected output failure.log")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--avi-config", type=Path, help="Optional pinned local AVI cache; no API request")
    args = p.parse_args()
    os.umask(0o077)
    require(os.environ.get("TMPDIR") and Path(os.environ["TMPDIR"]).is_dir(), "Set TMPDIR to an existing scratch directory")
    cfg = json.loads(args.config.read_text())
    cfg["_config_sha256"] = sha(args.config)
    # Relative input paths resolve against the configuration file, never the shell cwd.
    base = args.config.resolve().parent
    for item in cfg["inputs"].values():
        item["path"] = str((base / item["path"]).resolve())
    for item in cfg["resources"].values():
        item["path"] = str((base / item["path"]).resolve())
    for section in ["native", "evidence"]:
        cfg[section] = {k: str((base / v).resolve()) if k != "biology_gene_field" else v for k, v in cfg[section].items()}
    cfg["completed_reviews"] = str((base / cfg["completed_reviews"]).resolve())
    for k in ["metadata", "hpoa"]:
        cfg["lirical"][k] = str((base / cfg["lirical"][k]).resolve())
    for case in cfg["lirical"]["cases"].values():
        for k in case:
            case[k] = str((base / case[k]).resolve())
    if args.avi_config:
        from avi_replay import check_config, run_avi
        check_config(args.avi_config)
    result = run(cfg, args.output.resolve())
    if args.avi_config:
        result = {**result, "avi": run_avi(args.avi_config.resolve(), args.output.resolve())}
    print(json.dumps(result, sort_keys=True))
