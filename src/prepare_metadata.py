"""Prepare exact native LIRICAL allele metadata commands; analysis is not launched."""
import argparse
import json
import os
from pathlib import Path
import shlex
from io_utils import require, rows, write_json


def prepare(commands, pairs, sample, java, javac, jar, out):
    require(not out.exists(), "Output already exists")
    out.mkdir(mode=0o700, parents=True)
    baseline = [c for c in commands if "yaml" in c and c[c.index("-x") + 1] == "baseline"]
    require(len(baseline) == 1, "Unique baseline native command required")
    command = baseline[0]
    cli = command[command.index("yaml") + 1:]
    config = json.loads(Path(cli[cli.index("-y") + 1]).read_text())
    require(config["sampleId"] == sample, "Sample binding mismatch")
    contexts = sorted({(r["gene_id"], r[f"allele_{n}_key"]) for r in pairs for n in [1, 2]})
    write_json(out / "input.json", {"cli_args": cli, "sample": sample, "vcf": config["vcf"],
        "contexts": [{"gene": g, "key": k} for g, k in contexts], "task19_scope_genes": []})
    classes = out.resolve() / "classes"
    source = Path(__file__).resolve().parents[1] / "native/FullMetadata.java"
    classpath = str(classes) + os.pathsep + str(Path(jar).resolve()) + os.pathsep + str(Path(jar).resolve().parent / "lib/*")
    calls = [[javac, "-cp", classpath, "-d", str(classes), str(source)],
             [java, "-Xmx8g", "-XX:ActiveProcessorCount=4", "-XX:-UsePerfData", "-cp", classpath,
              "org.monarchinitiative.lirical.cli.cmd.FullMetadata", str(out.resolve() / "input.json"),
              str(out.resolve() / "metadata.json")]]
    write_json(out / "commands.json", calls)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", "umask 077",
             ': "${TMPDIR:?Set TMPDIR to an existing scratch directory}"', 'test -d "$TMPDIR"',
             "test ! -e " + shlex.quote(str(classes)), "mkdir " + shlex.quote(str(classes)), "set -o noclobber"]
    for i, call in enumerate(calls):
        line = shlex.join(call) if i == 0 else shlex.quote(call[0]) + ' "-Djava.io.tmpdir=$TMPDIR" ' + shlex.join(call[1:])
        lines.append(line + " >" + shlex.quote(str(out.resolve() / f"metadata_{i}.log")) + " 2>&1")
    (out / "run_metadata.sh").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["native-commands", "pairs", "sample-id", "java", "javac", "jar", "output"]:
        p.add_argument("--" + name, required=True)
    a = p.parse_args()
    os.umask(0o077)
    prepare(json.loads(Path(a.native_commands).read_text()), rows(a.pairs), a.sample_id, a.java, a.javac, a.jar, Path(a.output))
    print("Prepared metadata commands; no native analysis executed.")
