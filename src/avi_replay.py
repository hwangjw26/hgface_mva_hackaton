"""Run the optional AVI supplement from pinned local responses, without network access."""
import json
import subprocess
import sys
import traceback
from pathlib import Path

from io_utils import require, sha, write_json


def check_config(path):
    cfg = json.loads(path.read_text())
    require(cfg.get("schema_version") == 1 and cfg.get("cache_only") is True,
            "AVI requires the versioned cache-only configuration")
    require((cfg.get("cutoff"), cfg.get("operator"), cfg.get("sdk_version")) == (0.99379, ">=", "0.9.0"),
            "AVI method differs from the frozen supplement")
    require(cfg.get("score_scale") == "direct AVI_SCORE quantiles layer", "AVI score scale changed")
    require(set(cfg["inputs"]) == {"responses", "query", "metadata"}, "AVI requires three declared inputs")
    paths = {}
    for name, item in cfg["inputs"].items():
        source = (path.resolve().parent / item["path"]).resolve()
        require(source.is_file(), "Missing AVI input: " + name)
        require(sha(source) == item["sha256"], "AVI checksum mismatch: " + name)
        paths[name] = source
    return paths


def run_avi(config, baseline):
    paths = check_config(config)
    output = baseline / "avi"
    require(not output.exists(), "AVI output already exists")
    sources = {**paths, "config": config,
               **{name: baseline / name for name in ("allele_evidence.tsv", "integrated_pairs.tsv", "predictions.csv")}}
    before = {name: sha(path) for name, path in sources.items()}
    scripts = Path(__file__).resolve().parents[1] / "scripts/avi"
    commands = [
        ("avi_annotation.log", [sys.executable, "-B", str(scripts / "annotate_api.py"),
          "--baseline", str(baseline), "--responses", str(paths["responses"]),
          "--query", str(paths["query"]), "--metadata", str(paths["metadata"]), "--output", str(output)]),
        ("avi_verification.json", [sys.executable, "-B", str(scripts / "verify_api_annotation.py"),
          "--baseline", str(baseline), "--responses", str(paths["responses"]), "--output", str(output)]),
    ]
    try:
        for log_name, command in commands:
            with (baseline / log_name).open("x") as stdout, (baseline / (log_name + ".stderr")).open("x") as stderr:
                subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
        verification = json.loads((baseline / "avi_verification.json").read_text())
        require(verification["verification"] == "PASS", "Independent AVI verification did not pass")
        require(all(sha(path) == before[name] for name, path in sources.items()), "AVI source changed during replay")
        summary = json.loads((output / "aggregate_summary.json").read_text())
        require(summary["state"] == "ANNOTATED_API", "AVI results are incomplete")
        write_json(baseline / "avi_stage_receipt.json", {
            "status": "PASS", "cache_only": True, "network_requests": 0,
            "input_sha256": before, "config_sha256": before["config"],
            "verification": verification, "summary": summary,
            "script_sha256": {p.name: sha(p) for p in
                              [Path(__file__), scripts / "annotate_api.py", scripts / "annotate_avi.py", scripts / "verify_api_annotation.py"]},
            "annotation_sha256": {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()},
            "selection_rank_tier_phase_csv_changed": False})
    except Exception as error:
        with (baseline / "avi_stage_failure.log").open("x") as handle:
            traceback.print_exc(file=handle)
        write_json(baseline / "avi_stage_failure.json", {"status": "FAIL", "error_type": type(error).__name__})
        raise SystemExit("AVI supplement failed; inspect logs in the protected replay directory") from None
    return summary
