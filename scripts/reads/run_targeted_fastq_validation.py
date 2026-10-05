#!/usr/bin/env python3
"""Protected Task06 targeted FASTQ validation orchestrator."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import threading
import time
import traceback
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from read_evidence_core import classify_allele, classify_pair, classify_phase, reverse_complement


EXPECTED_PAIR_COUNT = 21
EXPECTED_ALLELE_COUNT = 28
SMOKE_PAIRS_PER_LANE = 100_000
SMOKE_MATCH_FRACTION_MAX = 0.0005
MAX_LOCAL_THREADS = 16
PRIMARY_EXCLUDED_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


class GateFailure(RuntimeError):
    """A deterministic workflow gate failed."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    path.chmod(0o600 if "sensitive" in path.parts else 0o644)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    path.chmod(0o600 if "sensitive" in path.parts else 0o644)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_gzip_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            text = compressed.write
            header = "\t".join(fieldnames) + "\n"
            text(header.encode())
            for row in rows:
                line = "\t".join(str(row.get(field, "")) for field in fieldnames) + "\n"
                text(line.encode())
    path.chmod(0o600)


def natural_contig_key(contig: str) -> tuple[int, str]:
    value = contig.removeprefix("chr")
    if value.isdigit():
        return int(value), ""
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return order.get(value, 100), value


class CommandRunner:
    def __init__(self, log_dir: Path) -> None:
        self.log_dir = log_dir
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_dir.chmod(0o700)
        self.command_log = self.log_dir / "commands.log"
        self.lock = threading.Lock()
        self.counter = 0

    def _record(self, args: list[str]) -> tuple[int, Path]:
        with self.lock:
            self.counter += 1
            number = self.counter
            with self.command_log.open("a", encoding="utf-8") as handle:
                handle.write(f"{number:04d}\t{shlex.join(args)}\n")
            self.command_log.chmod(0o600)
        stderr_path = self.log_dir / f"command_{number:04d}.stderr.log"
        return number, stderr_path

    def run(self, args: list[str], stdout_path: Path | None = None) -> float:
        number, stderr_path = self._record(args)
        start = time.monotonic()
        stdout_handle = stdout_path.open("wb") if stdout_path else subprocess.DEVNULL
        try:
            with stderr_path.open("wb") as stderr_handle:
                result = subprocess.run(
                    args,
                    check=False,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    env={**os.environ, "TMPDIR": os.environ["TMPDIR"]},
                )
        finally:
            if stdout_path:
                stdout_handle.close()
                stdout_path.chmod(0o600)
            stderr_path.chmod(0o600)
        if result.returncode:
            raise GateFailure(f"command_{number:04d}_failed")
        return time.monotonic() - start

    def capture(self, args: list[str]) -> str:
        with self.lock:
            capture_number = self.counter + 1
        path = self.log_dir / f"capture_{capture_number:04d}.stdout.log"
        self.run(args, stdout_path=path)
        return path.read_text(encoding="utf-8", errors="strict")

    def pipeline(self, left: list[str], right: list[str]) -> float:
        left_number, left_stderr = self._record(left)
        right_number, right_stderr = self._record(right)
        start = time.monotonic()
        with left_stderr.open("wb") as left_error, right_stderr.open("wb") as right_error:
            first = subprocess.Popen(
                left,
                stdout=subprocess.PIPE,
                stderr=left_error,
                env={**os.environ, "TMPDIR": os.environ["TMPDIR"]},
            )
            assert first.stdout is not None
            second = subprocess.run(
                right,
                stdin=first.stdout,
                stdout=subprocess.DEVNULL,
                stderr=right_error,
                check=False,
                env={**os.environ, "TMPDIR": os.environ["TMPDIR"]},
            )
            first.stdout.close()
            first_status = first.wait()
        left_stderr.chmod(0o600)
        right_stderr.chmod(0o600)
        if first_status or second.returncode:
            raise GateFailure(f"pipeline_{left_number:04d}_{right_number:04d}_failed")
        return time.monotonic() - start


def tool_paths(contract: dict[str, object]) -> dict[str, Path]:
    return {name: Path(value["path"]) for name, value in contract["executables"].items()}


def cpu_utilization(runner: CommandRunner) -> float:
    output = runner.capture(["/usr/bin/mpstat", "1", "1"])
    average = next(line for line in reversed(output.splitlines()) if "Average:" in line and "all" in line)
    return 100.0 - float(average.split()[-1])


def validate_fastq_names(r1_header: str, r2_header: str) -> bool:
    def parts(header: str) -> tuple[str, str | None]:
        if not header.startswith("@"):
            raise GateFailure("malformed_fastq_header")
        fields = header[1:].rstrip().split()
        core = fields[0]
        mate: str | None = None
        if core.endswith("/1") or core.endswith("/2"):
            mate = core[-1]
            core = core[:-2]
        elif len(fields) > 1 and fields[1][:2] in {"1:", "2:"}:
            mate = fields[1][0]
        return core, mate

    core1, mate1 = parts(r1_header)
    core2, mate2 = parts(r2_header)
    return core1 == core2 and (mate1 in {None, "1"}) and (mate2 in {None, "2"})


def copy_fastq_prefix(r1: Path, r2: Path, out1: Path, out2: Path, maximum: int) -> int:
    count = 0
    with gzip.open(r1, "rt", encoding="ascii", newline="") as source1, gzip.open(r2, "rt", encoding="ascii", newline="") as source2:
        with gzip.open(out1, "wt", encoding="ascii", newline="", compresslevel=6) as target1, gzip.open(out2, "wt", encoding="ascii", newline="", compresslevel=6) as target2:
            for _ in range(maximum):
                records1 = [source1.readline() for _ in range(4)]
                records2 = [source2.readline() for _ in range(4)]
                if not records1[0] and not records2[0]:
                    break
                if not all(records1) or not all(records2) or not validate_fastq_names(records1[0], records2[0]):
                    raise GateFailure("fastq_pair_validation_failed")
                if not records1[0].startswith("@") or records1[2].rstrip() != "+" or records2[2].rstrip() != "+":
                    raise GateFailure("fastq_record_validation_failed")
                if len(records1[1].rstrip()) != len(records1[3].rstrip()) or len(records2[1].rstrip()) != len(records2[3].rstrip()):
                    raise GateFailure("fastq_sequence_quality_length_failed")
                target1.writelines(records1)
                target2.writelines(records2)
                count += 1
    out1.chmod(0o600)
    out2.chmod(0o600)
    return count


def count_fastq_pairs(path: Path) -> int:
    lines = 0
    with gzip.open(path, "rb") as handle:
        for _ in handle:
            lines += 1
    if lines % 4:
        raise GateFailure("extracted_fastq_record_count_failed")
    return lines // 4


def gzip_content_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_contract(path: Path, reference: Path) -> tuple[dict[str, object], dict[str, Path]]:
    contract = read_json(path)
    if contract.get("schema_version") != "task06.read_evidence_contract.v1":
        raise GateFailure("contract_schema_failed")
    if Path(contract["reference"]["path"]).resolve() != reference.resolve():
        raise GateFailure("reference_contract_path_failed")
    tools = tool_paths(contract)
    for name, tool in tools.items():
        if not tool.is_file() or sha256(tool) != contract["executables"][name]["sha256"]:
            raise GateFailure("executable_contract_failed")
    if sha256(reference) != contract["reference"]["sha256"]:
        raise GateFailure("reference_hash_failed")
    for item in contract["reference"]["indices"]:
        index = Path(item["path"])
        if not index.is_file() or sha256(index) != item["sha256"]:
            raise GateFailure("reference_index_hash_failed")
    return contract, tools


def resolve_fastqs(pair_contract: Path, integrity_contract: Path) -> list[dict[str, object]]:
    pairs = read_tsv(pair_contract)
    integrity = read_tsv(integrity_contract)
    if len(pairs) != 4 or any(row["pair_status"] != "PASS" for row in pairs):
        raise GateFailure("fastq_pair_contract_failed")
    if len(integrity) != 8 or any(row["status"] != "PASS" or row["sha256_match"].lower() != "true" for row in integrity):
        raise GateFailure("fastq_integrity_contract_failed")
    by_name = {row["file_name"]: row for row in integrity}
    rows: list[dict[str, object]] = []
    for index, pair in enumerate(sorted(pairs, key=lambda row: row["lane"]), start=1):
        values: dict[str, object] = {"lane_id": f"L{index:02d}", "source_lane": pair["lane"]}
        for read in ("r1", "r2"):
            entry = by_name.get(Path(pair[f"{read}_file"]).name)
            if entry is None:
                raise GateFailure("fastq_integrity_join_failed")
            path = Path(entry["absolute_path"])
            if not path.is_file():
                raise GateFailure("fastq_input_missing")
            values[f"{read}_path"] = str(path)
            values[f"{read}_expected_sha256"] = entry["observed_sha256"]
            values[f"{read}_expected_size"] = int(entry["size_bytes"])
        rows.append(values)
    return rows


def hash_fastq_row(row: dict[str, object]) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for read in ("r1", "r2"):
        path = Path(str(row[f"{read}_path"]))
        stat = path.stat()
        observed = sha256(path)
        results.append(
            {
                "lane_id": row["lane_id"],
                "read": read.upper(),
                "absolute_path": str(path),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "expected_sha256": row[f"{read}_expected_sha256"],
                "observed_sha256": observed,
                "hash_match": observed == row[f"{read}_expected_sha256"],
                "size_match": stat.st_size == row[f"{read}_expected_size"],
            }
        )
    return results


def load_candidates(path: Path) -> tuple[list[dict[str, str]], list[dict[str, object]], list[dict[str, object]]]:
    original = read_tsv(path)
    if len(original) != EXPECTED_PAIR_COUNT:
        raise GateFailure("candidate_pair_count_failed")
    allele_by_key: dict[tuple[str, int, str, str], dict[str, object]] = {}
    pairs: list[dict[str, object]] = []
    for pair_index, row in enumerate(original, start=1):
        allele_ids: list[str] = []
        for slot in (1, 2):
            key = (
                row[f"allele_{slot}_normalized_chrom"],
                int(row[f"allele_{slot}_normalized_pos"]),
                row[f"allele_{slot}_normalized_ref"].upper(),
                row[f"allele_{slot}_normalized_alt"].upper(),
            )
            if any(not value or re.search(r"[^ACGT]", value) for value in key[2:]) or key[2] == key[3]:
                raise GateFailure("candidate_representation_failed")
            if key not in allele_by_key:
                allele_by_key[key] = {
                    "chrom": key[0],
                    "pos1": key[1],
                    "ref": key[2],
                    "alt": key[3],
                    "variant_class": row[f"allele_{slot}_normalized_variant_class"],
                    "original_dp": row[f"allele_{slot}_normalized_dp"],
                    "original_ad": row[f"allele_{slot}_normalized_ad"],
                    "representation_ok": True,
                }
            allele_ids.append("")
        pairs.append({"pair_id": f"P{pair_index:03d}", "allele_keys": [
            (row["allele_1_normalized_chrom"], int(row["allele_1_normalized_pos"]), row["allele_1_normalized_ref"].upper(), row["allele_1_normalized_alt"].upper()),
            (row["allele_2_normalized_chrom"], int(row["allele_2_normalized_pos"]), row["allele_2_normalized_ref"].upper(), row["allele_2_normalized_alt"].upper()),
        ]})
    if len(allele_by_key) != EXPECTED_ALLELE_COUNT:
        raise GateFailure("candidate_allele_count_failed")
    ordered_keys = sorted(allele_by_key, key=lambda key: (natural_contig_key(key[0]), key[1], key[2], key[3]))
    id_by_key = {key: f"A{index:03d}" for index, key in enumerate(ordered_keys, start=1)}
    alleles: list[dict[str, object]] = []
    for key in ordered_keys:
        record = allele_by_key[key]
        record["allele_id"] = id_by_key[key]
        alleles.append(record)
    for pair in pairs:
        pair["allele_1_id"] = id_by_key[pair["allele_keys"][0]]
        pair["allele_2_id"] = id_by_key[pair["allele_keys"][1]]
        del pair["allele_keys"]
    return original, alleles, pairs


def read_fai(reference: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line in Path(str(reference) + ".fai").read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        lengths[fields[0]] = int(fields[1])
    return lengths


def resolve_contig(chrom: str, lengths: dict[str, int]) -> str:
    alternatives = [chrom]
    if chrom == "M":
        alternatives.append("MT")
    if chrom.startswith("chr"):
        alternatives.append(chrom[3:])
    else:
        alternatives.append("chr" + chrom)
    for value in alternatives:
        if value in lengths:
            return value
    raise GateFailure("candidate_contig_resolution_failed")


def fetch_reference_window(runner: CommandRunner, samtools: Path, reference: Path, contig: str, start0: int, end0: int) -> str:
    output = runner.capture([str(samtools), "faidx", str(reference), f"{contig}:{start0 + 1}-{end0}"])
    lines = [line.strip() for line in output.splitlines() if line and not line.startswith(">")]
    sequence = "".join(lines).upper()
    if len(sequence) != end0 - start0:
        raise GateFailure("reference_window_length_failed")
    return sequence


def tiled_probe_starts(start0: int, end0: int, contig_length: int, probe_length: int, tile_step: int, span: int) -> list[tuple[str, int]]:
    result: list[tuple[str, int]] = []
    left_first = max(0, start0 - span)
    left_last = start0 - probe_length
    if left_last >= left_first:
        starts = list(range(left_first, left_last + 1, tile_step))
        if not starts or starts[-1] != left_last:
            starts.append(left_last)
        result.extend(("LEFT", value) for value in starts)
    right_first = end0
    right_last = min(contig_length - probe_length, end0 + span - probe_length)
    if right_last >= right_first:
        starts = list(range(right_first, right_last + 1, tile_step))
        if not starts or starts[-1] != right_last:
            starts.append(right_last)
        result.extend(("RIGHT", value) for value in starts)
    return result


def build_probes(
    alleles: list[dict[str, object]],
    reference: Path,
    tools: dict[str, Path],
    sensitive: Path,
    runner: CommandRunner,
    probe_contract: dict[str, object],
) -> list[dict[str, object]]:
    lengths = read_fai(reference)
    probe_length = int(probe_contract["probe_length"])
    tile_step = int(probe_contract["tile_step"])
    span = int(probe_contract["flank_span"])
    probes: list[dict[str, object]] = []
    probe_index = 0
    for allele in alleles:
        contig = resolve_contig(str(allele["chrom"]), lengths)
        allele["reference_chrom"] = contig
        start0 = int(allele["pos1"]) - 1
        end0 = start0 + len(str(allele["ref"]))
        window_start = max(0, start0 - span)
        window_end = min(lengths[contig], end0 + span)
        window = fetch_reference_window(runner, tools["samtools"], reference, contig, window_start, window_end)
        observed_ref = window[start0 - window_start : end0 - window_start]
        reference_ok = observed_ref == allele["ref"]
        allele["reference_ok"] = reference_ok
        if not reference_ok:
            raise GateFailure("candidate_reference_match_failed")
        for side, probe_start in tiled_probe_starts(start0, end0, lengths[contig], probe_length, tile_step, span):
            local_start = probe_start - window_start
            sequence = window[local_start : local_start + probe_length]
            if len(sequence) != probe_length or "N" in sequence:
                continue
            if not (probe_start + probe_length <= start0 or probe_start >= end0):
                raise GateFailure("probe_overlap_failed")
            for orientation, oriented in (("FORWARD", sequence), ("REVERSE", reverse_complement(sequence))):
                probe_index += 1
                probes.append(
                    {
                        "probe_id": f"Q{probe_index:05d}",
                        "allele_id": allele["allele_id"],
                        "reference_chrom": contig,
                        "side": side,
                        "start0": probe_start,
                        "end0": probe_start + probe_length,
                        "orientation": orientation,
                        "sequence": oriented,
                    }
                )
    all_fasta = sensitive / "probes" / "all_probe_instances.fa"
    with all_fasta.open("w", encoding="ascii") as handle:
        for row in probes:
            handle.write(f">{row['probe_id']}\n{row['sequence']}\n")
    all_fasta.chmod(0o600)
    audit_sam = sensitive / "probes" / "full_reference_probe_audit.sam"
    runner.run(
        [str(tools["bwa"]), "mem", "-a", "-T", "0", "-k", "19", "-t", "4", str(reference), str(all_fasta)],
        stdout_path=audit_sam,
    )
    exact_hits: dict[str, set[tuple[str, int, str]]] = defaultdict(set)
    maximum_mapq: dict[str, int] = defaultdict(int)
    with audit_sam.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip().split("\t")
            flag = int(fields[1])
            if flag & 0x4:
                continue
            tags = {item.split(":", 2)[0]: item.split(":", 2)[2] for item in fields[11:] if item.count(":") >= 2}
            if fields[5] not in {f"{probe_length}M", f"{probe_length}="} or tags.get("NM") != "0":
                continue
            strand = "-" if flag & 0x10 else "+"
            exact_hits[fields[0]].add((fields[2], int(fields[3]), strand))
            maximum_mapq[fields[0]] = max(maximum_mapq[fields[0]], int(fields[4]))
    for row in probes:
        row["exact_hit_count"] = len(exact_hits[row["probe_id"]])
        row["maximum_mapq"] = maximum_mapq[row["probe_id"]]
        row["is_unique"] = row["exact_hit_count"] == 1 and row["maximum_mapq"] >= 20
    by_allele: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in probes:
        by_allele[str(row["allele_id"])].append(row)
    for allele in alleles:
        rows = by_allele[str(allele["allele_id"])]
        allele["left_probe_count"] = sum(row["side"] == "LEFT" for row in rows)
        allele["right_probe_count"] = sum(row["side"] == "RIGHT" for row in rows)
        allele["unique_probe_count"] = sum(bool(row["is_unique"]) for row in rows)
        allele["probe_ok"] = allele["left_probe_count"] > 0 and allele["right_probe_count"] > 0
        allele["mapping_caveat"] = allele["unique_probe_count"] == 0
    write_tsv(
        sensitive / "probes" / "probe_manifest.tsv",
        ["probe_id", "allele_id", "reference_chrom", "side", "start0", "end0", "orientation", "sequence", "exact_hit_count", "maximum_mapq", "is_unique"],
        probes,
    )
    adapters = sorted({str(row["sequence"]) for row in probes})
    adapter_fasta = sensitive / "probes" / "allele_neutral_probes.fa"
    with adapter_fasta.open("w", encoding="ascii") as handle:
        for index, sequence in enumerate(adapters, start=1):
            handle.write(f">T{index:05d}\n{sequence}\n")
    adapter_fasta.chmod(0o600)
    return probes


def write_targets(alleles: list[dict[str, object]], sensitive: Path, tools: dict[str, Path], runner: CommandRunner) -> None:
    bed = sensitive / "probes" / "targets.bed"
    target_tsv = sensitive / "probes" / "known_alleles.tsv"
    with bed.open("w", encoding="ascii") as bed_handle, target_tsv.open("w", encoding="ascii") as target_handle:
        for allele in alleles:
            start0 = int(allele["pos1"]) - 1
            end0 = start0 + len(str(allele["ref"]))
            bed_handle.write(f"{allele['reference_chrom']}\t{start0}\t{end0}\n")
            target_handle.write(f"{allele['reference_chrom']}\t{allele['pos1']}\t{allele['ref']},{allele['alt']}\n")
    bed.chmod(0o600)
    target_tsv.chmod(0o600)
    runner.run([str(tools["bgzip"]), "-f", str(target_tsv)])
    target_gz = Path(str(target_tsv) + ".gz")
    target_gz.chmod(0o600)
    runner.run([str(tools["tabix"]), "-f", "-s", "1", "-b", "2", "-e", "2", str(target_gz)])
    Path(str(target_gz) + ".tbi").chmod(0o600)


def prepare(args: argparse.Namespace, run_dir: Path, sensitive: Path, runner: CommandRunner) -> dict[str, object]:
    cpu_used = cpu_utilization(runner)
    disk = shutil.disk_usage(run_dir)
    if cpu_used >= 70.0:
        raise GateFailure("local_cpu_gate_failed")
    contract, tools = load_contract(args.contract.resolve(), args.reference.resolve())
    fastqs = resolve_fastqs(args.fastq_contract.resolve(), args.fastq_integrity.resolve())
    with ThreadPoolExecutor(max_workers=4) as pool:
        integrity_nested = list(pool.map(hash_fastq_row, fastqs))
    integrity = [row for group in integrity_nested for row in group]
    if not all(bool(row["hash_match"]) and bool(row["size_match"]) for row in integrity):
        raise GateFailure("pre_scan_fastq_integrity_failed")
    write_tsv(
        sensitive / "inputs" / "pre_scan_integrity.tsv",
        ["lane_id", "read", "absolute_path", "size_bytes", "mtime_ns", "expected_sha256", "observed_sha256", "hash_match", "size_match"],
        integrity,
    )
    write_tsv(
        sensitive / "inputs" / "fastq_manifest.tsv",
        ["lane_id", "source_lane", "r1_path", "r2_path", "r1_expected_sha256", "r2_expected_sha256", "r1_expected_size", "r2_expected_size"],
        fastqs,
    )
    original, alleles, pairs = load_candidates(args.candidate_handoff.resolve())
    probes = build_probes(alleles, args.reference.resolve(), tools, sensitive, runner, contract["probes"])
    if not probes:
        raise GateFailure("probe_generation_failed")
    write_targets(alleles, sensitive, tools, runner)
    write_tsv(
        sensitive / "inputs" / "candidate_alleles.tsv",
        ["allele_id", "chrom", "reference_chrom", "pos1", "ref", "alt", "variant_class", "original_dp", "original_ad", "representation_ok", "reference_ok", "left_probe_count", "right_probe_count", "unique_probe_count", "probe_ok", "mapping_caveat"],
        alleles,
    )
    write_tsv(sensitive / "inputs" / "candidate_pairs.tsv", ["pair_id", "allele_1_id", "allele_2_id"], pairs)
    source_fields = list(original[0])
    write_gzip_tsv(sensitive / "inputs" / "task05_handoff_snapshot.tsv.gz", source_fields, original)
    metadata = {
        "status": "PASS",
        "pair_count": len(pairs),
        "allele_count": len(alleles),
        "lane_count": len(fastqs),
        "fastq_file_count": len(integrity),
        "source_handoff_sha256": sha256(args.candidate_handoff.resolve()),
        "source_handoff_snapshot_sha256": sha256(sensitive / "inputs" / "task05_handoff_snapshot.tsv.gz"),
        "probe_instance_count": len(probes),
        "probe_sequence_count": len({str(row["sequence"]) for row in probes}),
        "unique_probe_instance_count": sum(bool(row["is_unique"]) for row in probes),
        "alleles_with_probe_caveat": sum(bool(row["mapping_caveat"]) for row in alleles),
        "alleles_with_missing_flank": sum(not bool(row["probe_ok"]) for row in alleles),
        "reference_sha256": contract["reference"]["sha256"],
        "contract_sha256": sha256(args.contract.resolve()),
    }
    write_json(sensitive / "inputs" / "prepare_metadata.json", metadata)
    ordinary = {
        "status": "PASS",
        "candidate_pair_count": len(pairs),
        "unique_allele_count": len(alleles),
        "lane_count": len(fastqs),
        "fastq_file_count": len(integrity),
        "fastq_hash_matches": sum(bool(row["hash_match"]) for row in integrity),
        "probe_instance_count": len(probes),
        "distinct_probe_sequence_count": len({str(row["sequence"]) for row in probes}),
        "unique_probe_instance_count": sum(bool(row["is_unique"]) for row in probes),
        "alleles_with_mapping_caveat": sum(bool(row["mapping_caveat"]) for row in alleles),
        "alleles_with_missing_flank": sum(not bool(row["probe_ok"]) for row in alleles),
        "variant_class_counts": dict(sorted(Counter(str(row["variant_class"]) for row in alleles).items())),
        "cpu_total_utilization_percent": round(cpu_used, 2),
        "local_cpu_gate": "PASS",
        "disk_free_bytes": disk.free,
        "file_descriptor_limit": os.sysconf("SC_OPEN_MAX"),
        "max_local_threads": MAX_LOCAL_THREADS,
        "patient_identity_in_ordinary_output": False,
    }
    write_json(run_dir / "reports" / "preflight_summary.json", ordinary)
    return ordinary


def cutadapt_command(cutadapt: Path, probes: Path, r1: Path, r2: Path, out1: Path, out2: Path, cores: int | None = None) -> list[str]:
    command = [
        str(cutadapt),
        "--quiet",
        "--action=none",
        "--no-indels",
        "-O",
        "31",
        "-e",
        "0.04",
        "-b",
        f"file:{probes}",
        "-B",
        f"file:{probes}",
        "--discard-untrimmed",
        "--pair-filter=both",
        "-o",
        str(out1),
        "-p",
        str(out2),
        str(r1),
        str(r2),
    ]
    if cores is not None:
        command.insert(1, str(cores))
        command.insert(1, "-j")
    return command


def smoke_lane(row: dict[str, str], sensitive: Path, cutadapt: Path, probes: Path, runner: CommandRunner) -> dict[str, object]:
    lane = row["lane_id"]
    lane_dir = sensitive / "smoke" / lane
    lane_dir.mkdir(parents=True, exist_ok=True)
    lane_dir.chmod(0o700)
    prefix1 = lane_dir / "prefix_R1.fastq.gz"
    prefix2 = lane_dir / "prefix_R2.fastq.gz"
    input_pairs = copy_fastq_prefix(Path(row["r1_path"]), Path(row["r2_path"]), prefix1, prefix2, SMOKE_PAIRS_PER_LANE)
    out1 = lane_dir / "extracted_R1.fastq.gz"
    out2 = lane_dir / "extracted_R2.fastq.gz"
    wall = runner.run(cutadapt_command(cutadapt, probes, prefix1, prefix2, out1, out2))
    extracted1 = count_fastq_pairs(out1)
    extracted2 = count_fastq_pairs(out2)
    if extracted1 != extracted2:
        raise GateFailure("smoke_mate_count_failed")
    return {
        "lane_id": lane,
        "input_pairs": input_pairs,
        "prefix_input_bytes": prefix1.stat().st_size + prefix2.stat().st_size,
        "extracted_pairs": extracted1,
        "match_fraction": extracted1 / input_pairs if input_pairs else 0.0,
        "output_bytes": out1.stat().st_size + out2.stat().st_size,
        "wall_seconds": round(wall, 3),
    }


def smoke(args: argparse.Namespace, run_dir: Path, sensitive: Path, runner: CommandRunner) -> dict[str, object]:
    contract = read_json(args.contract.resolve())
    tools = tool_paths(contract)
    fastqs = read_tsv(sensitive / "inputs" / "fastq_manifest.tsv")
    probes = sensitive / "probes" / "allele_neutral_probes.fa"
    if len(fastqs) != 4 or not probes.is_file():
        raise GateFailure("smoke_predecessor_failed")
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda row: smoke_lane(row, sensitive, tools["cutadapt"], probes, runner), fastqs))
    input_pairs = sum(int(row["input_pairs"]) for row in results)
    extracted_pairs = sum(int(row["extracted_pairs"]) for row in results)
    fraction = extracted_pairs / input_pairs if input_pairs else math.inf
    total_input_bytes = sum(int(row["r1_expected_size"]) + int(row["r2_expected_size"]) for row in fastqs)
    estimated_extracted_bytes = int(total_input_bytes * fraction * 2.0)
    full_bytes_by_lane = {row["lane_id"]: int(row["r1_expected_size"]) + int(row["r2_expected_size"]) for row in fastqs}
    estimated_lane_wall_seconds = {
        row["lane_id"]: float(row["wall_seconds"]) * full_bytes_by_lane[row["lane_id"]] / max(1, int(row["prefix_input_bytes"]))
        for row in results
    }
    disk_free = shutil.disk_usage(run_dir).free
    estimated_working_bytes = estimated_extracted_bytes * 10 + 20 * 1024**3
    gate = fraction <= SMOKE_MATCH_FRACTION_MAX and estimated_working_bytes < disk_free
    protected = {
        "status": "PASS" if gate else "FAIL",
        "smoke_pairs_per_lane": SMOKE_PAIRS_PER_LANE,
        "input_pairs": input_pairs,
        "extracted_pairs": extracted_pairs,
        "match_fraction": fraction,
        "maximum_match_fraction": SMOKE_MATCH_FRACTION_MAX,
        "estimated_extracted_bytes": estimated_extracted_bytes,
        "estimated_working_bytes": estimated_working_bytes,
        "estimated_parallel_full_scan_wall_seconds": max(estimated_lane_wall_seconds.values()),
        "disk_free_bytes": disk_free,
        "lanes": results,
    }
    write_json(sensitive / "smoke" / "smoke_summary.json", protected)
    ordinary = {
        "status": protected["status"],
        "lane_count": len(results),
        "input_pair_count": input_pairs,
        "extracted_pair_count": extracted_pairs,
        "probe_match_fraction": fraction,
        "maximum_probe_match_fraction": SMOKE_MATCH_FRACTION_MAX,
        "estimated_extracted_bytes": estimated_extracted_bytes,
        "estimated_working_bytes": estimated_working_bytes,
        "estimated_parallel_full_scan_wall_seconds": max(estimated_lane_wall_seconds.values()),
        "disk_free_bytes": disk_free,
        "lane_aggregates": results,
        "patient_identity_in_ordinary_output": False,
    }
    write_json(run_dir / "reports" / "protected_smoke_summary.json", ordinary)
    if not gate:
        raise GateFailure("bounded_smoke_gate_failed")
    return ordinary


def parse_metric_tsv(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    return {row["metric"]: row["value"] for row in rows}


def compile_seed_prefilter(run_dir: Path, runner: CommandRunner) -> tuple[Path, dict[str, str]]:
    source = run_dir / "scripts" / "lossless_seed_prefilter.c"
    binary = run_dir / "scripts" / "lossless_seed_prefilter"
    runner.run(["/usr/bin/gcc", "-O3", "-std=c11", "-Wall", "-Wextra", "-Werror", str(source), "-lz", "-o", str(binary)])
    binary.chmod(0o755)
    compiler = runner.capture(["/usr/bin/gcc", "--version"]).splitlines()[0]
    return binary, {
        "source_sha256": sha256(source),
        "binary_sha256": sha256(binary),
        "gcc_sha256": sha256(Path("/usr/bin/gcc")),
        "gcc_version": compiler,
    }


def seed_prefilter_lane(
    lane: str,
    binary: Path,
    probes: Path,
    r1: Path,
    r2: Path,
    out1: Path,
    out2: Path,
    stats: Path,
    runner: CommandRunner,
) -> dict[str, object]:
    wall = runner.run([str(binary), str(probes), str(r1), str(r2), str(out1), str(out2), str(stats)])
    out1.chmod(0o600)
    out2.chmod(0o600)
    stats.chmod(0o600)
    values = parse_metric_tsv(stats)
    if values.get("invalid_records") != "0":
        raise GateFailure("lossless_prefilter_fastq_validation_failed")
    return {
        "lane_id": lane,
        "input_pairs": int(values["input_pairs"]),
        "retained_pairs": int(values["retained_pairs"]),
        "retained_fraction": float(values["retained_fraction"]),
        "distinct_seeds": int(values["distinct_seeds"]),
        "wall_seconds": round(wall, 3),
        "output_bytes": out1.stat().st_size + out2.stat().st_size,
    }


def validate_prefilter_lane(
    row: dict[str, str],
    sensitive: Path,
    binary: Path,
    probes: Path,
    cutadapt: Path,
    runner: CommandRunner,
) -> dict[str, object]:
    lane = row["lane_id"]
    directory = sensitive / "smoke" / lane
    seed1 = directory / "seed_prefilter_R1.fastq.gz"
    seed2 = directory / "seed_prefilter_R2.fastq.gz"
    stats = directory / "seed_prefilter_stats.tsv"
    result = seed_prefilter_lane(
        lane,
        binary,
        probes,
        directory / "prefix_R1.fastq.gz",
        directory / "prefix_R2.fastq.gz",
        seed1,
        seed2,
        stats,
        runner,
    )
    accelerated1 = directory / "accelerated_extracted_R1.fastq.gz"
    accelerated2 = directory / "accelerated_extracted_R2.fastq.gz"
    cutadapt_wall = runner.run(cutadapt_command(cutadapt, probes, seed1, seed2, accelerated1, accelerated2, cores=4))
    direct1 = directory / "extracted_R1.fastq.gz"
    direct2 = directory / "extracted_R2.fastq.gz"
    content_identical = gzip_content_sha256(direct1) == gzip_content_sha256(accelerated1) and gzip_content_sha256(direct2) == gzip_content_sha256(accelerated2)
    direct_pairs = count_fastq_pairs(direct1)
    accelerated_pairs = count_fastq_pairs(accelerated1)
    result.update(
        {
            "cutadapt_wall_seconds": round(cutadapt_wall, 3),
            "direct_cutadapt_pairs": direct_pairs,
            "prefilter_then_cutadapt_pairs": accelerated_pairs,
            "content_identical": content_identical,
        }
    )
    return result


def prefilter_validation(args: argparse.Namespace, run_dir: Path, sensitive: Path, runner: CommandRunner) -> dict[str, object]:
    smoke_summary = read_json(sensitive / "smoke" / "smoke_summary.json")
    if smoke_summary.get("status") != "PASS":
        raise GateFailure("prefilter_smoke_predecessor_failed")
    contract = read_json(args.contract.resolve())
    tools = tool_paths(contract)
    binary, executable = compile_seed_prefilter(run_dir, runner)
    probes = sensitive / "probes" / "allele_neutral_probes.fa"
    fastqs = read_tsv(sensitive / "inputs" / "fastq_manifest.tsv")
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda row: validate_prefilter_lane(row, sensitive, binary, probes, tools["cutadapt"], runner), fastqs))
    proof_rows = []
    for overlap in range(31, 52):
        errors = int(0.04 * overlap)
        guaranteed_run = math.ceil((overlap - errors) / (errors + 1))
        proof_rows.append({"overlap_length": overlap, "maximum_mismatches": errors, "guaranteed_exact_run": guaranteed_run})
    guaranteed_minimum = min(int(row["guaranteed_exact_run"]) for row in proof_rows)
    equivalence = all(bool(row["content_identical"]) and row["direct_cutadapt_pairs"] == row["prefilter_then_cutadapt_pairs"] for row in results)
    proof_ok = guaranteed_minimum >= 15
    full_bytes = {row["lane_id"]: int(row["r1_expected_size"]) + int(row["r2_expected_size"]) for row in fastqs}
    prefix_bytes = {
        row["lane_id"]: (sensitive / "smoke" / row["lane_id"] / "prefix_R1.fastq.gz").stat().st_size
        + (sensitive / "smoke" / row["lane_id"] / "prefix_R2.fastq.gz").stat().st_size
        for row in fastqs
    }
    estimated_prefilter = max(float(row["wall_seconds"]) * full_bytes[row["lane_id"]] / prefix_bytes[row["lane_id"]] for row in results)
    maximum_retained_fraction = max(float(row["retained_fraction"]) for row in results)
    estimated_cutadapt = float(smoke_summary["estimated_parallel_full_scan_wall_seconds"]) * maximum_retained_fraction / 4.0
    estimated_total = estimated_prefilter + estimated_cutadapt
    status = "PASS" if equivalence and proof_ok else "FAIL"
    protected = {
        "status": status,
        "seed_length": 15,
        "minimum_cutadapt_overlap": 31,
        "maximum_error_rate": 0.04,
        "no_indels": True,
        "lossless_superset_proof": "At most e mismatches split matching positions into e+1 exact runs; for overlap 31-51 the shortest guaranteed run is 15.",
        "minimum_guaranteed_exact_run": guaranteed_minimum,
        "prefix_content_equivalence": equivalence,
        "lanes": results,
        "proof_rows": proof_rows,
        "executable": executable,
        "estimated_full_prefilter_seconds": estimated_prefilter,
        "estimated_full_cutadapt_seconds": estimated_cutadapt,
        "estimated_full_total_seconds": estimated_total,
    }
    write_json(sensitive / "smoke" / "lossless_prefilter_validation.json", protected)
    ordinary = {
        "status": status,
        "lane_count": len(results),
        "seed_length": 15,
        "minimum_guaranteed_exact_run": guaranteed_minimum,
        "prefix_content_equivalence": equivalence,
        "direct_extracted_pairs": sum(int(row["direct_cutadapt_pairs"]) for row in results),
        "accelerated_extracted_pairs": sum(int(row["prefilter_then_cutadapt_pairs"]) for row in results),
        "maximum_prefilter_retained_fraction": maximum_retained_fraction,
        "estimated_full_prefilter_seconds": estimated_prefilter,
        "estimated_full_cutadapt_seconds": estimated_cutadapt,
        "estimated_full_total_seconds": estimated_total,
        "source_sha256": executable["source_sha256"],
        "binary_sha256": executable["binary_sha256"],
        "patient_identity_in_ordinary_output": False,
    }
    write_json(run_dir / "reports" / "lossless_prefilter_validation.json", ordinary)
    if status != "PASS":
        raise GateFailure("lossless_prefilter_equivalence_failed")
    return ordinary


def extract_lane(row: dict[str, str], sensitive: Path, cutadapt: Path, probes: Path, runner: CommandRunner) -> dict[str, object]:
    lane = row["lane_id"]
    lane_dir = sensitive / "lanes" / lane
    lane_dir.mkdir(parents=True, exist_ok=True)
    lane_dir.chmod(0o700)
    out1 = lane_dir / "extracted_R1.fastq.gz"
    out2 = lane_dir / "extracted_R2.fastq.gz"
    input1 = Path(row.get("cutadapt_r1", row["r1_path"]))
    input2 = Path(row.get("cutadapt_r2", row["r2_path"]))
    cores = 4 if "cutadapt_r1" in row else None
    wall = runner.run(cutadapt_command(cutadapt, probes, input1, input2, out1, out2, cores=cores))
    extracted1 = count_fastq_pairs(out1)
    extracted2 = count_fastq_pairs(out2)
    if extracted1 != extracted2:
        raise GateFailure("full_extraction_mate_count_failed")
    return {
        "lane_id": lane,
        "input_bytes": int(row["r1_expected_size"]) + int(row["r2_expected_size"]),
        "extracted_pairs": extracted1,
        "output_bytes": out1.stat().st_size + out2.stat().st_size,
        "wall_seconds": round(wall, 3),
    }


def full_prefilter_lane(
    row: dict[str, str],
    sensitive: Path,
    binary: Path,
    probes: Path,
    runner: CommandRunner,
) -> dict[str, object]:
    lane = row["lane_id"]
    directory = sensitive / "lanes" / lane
    directory.mkdir(parents=True, exist_ok=True)
    directory.chmod(0o700)
    return seed_prefilter_lane(
        lane,
        binary,
        probes,
        Path(row["r1_path"]),
        Path(row["r2_path"]),
        directory / "prefilter_R1.fastq.gz",
        directory / "prefilter_R2.fastq.gz",
        directory / "prefilter_stats.tsv",
        runner,
    )


def align_lane(row: dict[str, str], reference: Path, tools: dict[str, Path], sensitive: Path, runner: CommandRunner) -> dict[str, object]:
    lane = row["lane_id"]
    lane_dir = sensitive / "lanes" / lane
    alignment_dir = sensitive / "alignments" / lane
    alignment_dir.mkdir(parents=True, exist_ok=True)
    alignment_dir.chmod(0o700)
    r1 = lane_dir / "extracted_R1.fastq.gz"
    r2 = lane_dir / "extracted_R2.fastq.gz"
    sam = alignment_dir / "full_reference.sam"
    name_bam = alignment_dir / "name_sorted.bam"
    fixmate_bam = alignment_dir / "fixmate.bam"
    coord_bam = alignment_dir / "coordinate_sorted.bam"
    rg = f"@RG\\tID:{lane}\\tSM:PROTECTED_SINGLETON\\tPL:ILLUMINA"
    wall = 0.0
    wall += runner.run([str(tools["bwa"]), "mem", "-t", "4", "-R", rg, str(reference), str(r1), str(r2)], stdout_path=sam)
    wall += runner.run([str(tools["samtools"]), "sort", "-n", "-@", "2", "-o", str(name_bam), str(sam)])
    wall += runner.run([str(tools["samtools"]), "fixmate", "-m", str(name_bam), str(fixmate_bam)])
    wall += runner.run([str(tools["samtools"]), "sort", "-@", "2", "-o", str(coord_bam), str(fixmate_bam)])
    return {
        "lane_id": lane,
        "alignment_wall_seconds": round(wall, 3),
        "sam_bytes": sam.stat().st_size,
        "coordinate_bam_bytes": coord_bam.stat().st_size,
    }


def parse_flagstat(text: str) -> dict[str, int]:
    result: dict[str, int] = {}
    labels = {
        "in total": "total_records",
        "primary mapped": "primary_mapped",
        "duplicates": "duplicate_records",
        "properly paired": "properly_paired",
    }
    for line in text.splitlines():
        for phrase, key in labels.items():
            if phrase in line and key not in result:
                result[key] = int(line.split("+", 1)[0].strip())
    return result


def full_scan(args: argparse.Namespace, run_dir: Path, sensitive: Path, runner: CommandRunner) -> dict[str, object]:
    cpu_used = cpu_utilization(runner)
    if cpu_used >= 70.0:
        raise GateFailure("full_scan_cpu_gate_failed")
    smoke_summary = read_json(sensitive / "smoke" / "smoke_summary.json")
    if smoke_summary.get("status") != "PASS":
        raise GateFailure("full_scan_smoke_predecessor_failed")
    prefilter_contract = read_json(sensitive / "smoke" / "lossless_prefilter_validation.json")
    if prefilter_contract.get("status") != "PASS" or not prefilter_contract.get("prefix_content_equivalence"):
        raise GateFailure("full_scan_prefilter_predecessor_failed")
    contract = read_json(args.contract.resolve())
    tools = tool_paths(contract)
    fastqs = read_tsv(sensitive / "inputs" / "fastq_manifest.tsv")
    probes = sensitive / "probes" / "allele_neutral_probes.fa"
    binary = run_dir / "scripts" / "lossless_seed_prefilter"
    if not binary.is_file() or sha256(binary) != prefilter_contract["executable"]["binary_sha256"]:
        raise GateFailure("full_scan_prefilter_binary_failed")
    with ThreadPoolExecutor(max_workers=4) as pool:
        prefilter = list(pool.map(lambda row: full_prefilter_lane(row, sensitive, binary, probes, runner), fastqs))
    cutadapt_rows = []
    for row in fastqs:
        updated = dict(row)
        updated["cutadapt_r1"] = str(sensitive / "lanes" / row["lane_id"] / "prefilter_R1.fastq.gz")
        updated["cutadapt_r2"] = str(sensitive / "lanes" / row["lane_id"] / "prefilter_R2.fastq.gz")
        cutadapt_rows.append(updated)
    with ThreadPoolExecutor(max_workers=4) as pool:
        extraction = list(pool.map(lambda row: extract_lane(row, sensitive, tools["cutadapt"], probes, runner), cutadapt_rows))
    with ThreadPoolExecutor(max_workers=4) as pool:
        alignment = list(pool.map(lambda row: align_lane(row, args.reference.resolve(), tools, sensitive, runner), fastqs))
    merged = sensitive / "alignments" / "merged_preduplicate.bam"
    marked = sensitive / "alignments" / "merged_marked.bam"
    lane_bams = [sensitive / "alignments" / row["lane_id"] / "coordinate_sorted.bam" for row in fastqs]
    merge_wall = runner.run([str(tools["samtools"]), "merge", "-f", "-@", "4", "-o", str(merged), *map(str, lane_bams)])
    runner.run([str(tools["samtools"]), "index", "-@", "4", str(merged)])
    markdup_wall = runner.run([str(tools["samtools"]), "markdup", "-@", "4", str(merged), str(marked)])
    index_wall = runner.run([str(tools["samtools"]), "index", "-@", "4", str(marked)])
    flagstat_text = runner.capture([str(tools["samtools"]), "flagstat", "-@", "4", str(marked)])
    flagstat = parse_flagstat(flagstat_text)
    protected = {
        "status": "PASS",
        "cpu_total_utilization_percent_before_full_scan": cpu_used,
        "lossless_prefilter": prefilter,
        "extraction": extraction,
        "alignment": alignment,
        "merge_wall_seconds": round(merge_wall, 3),
        "markdup_wall_seconds": round(markdup_wall, 3),
        "index_wall_seconds": round(index_wall, 3),
        "flagstat": flagstat,
        "marked_bam_bytes": marked.stat().st_size,
    }
    write_json(sensitive / "alignments" / "full_scan_summary.json", protected)
    ordinary = {
        "status": "PASS",
        "lane_count": len(fastqs),
        "cpu_total_utilization_percent_before_full_scan": round(cpu_used, 2),
        "lossless_prefilter": prefilter,
        "extraction": extraction,
        "alignment": alignment,
        "alignment_aggregate": flagstat,
        "marked_bam_bytes": marked.stat().st_size,
        "full_reference_mapping": True,
        "novel_calling_performed": False,
        "patient_identity_in_ordinary_output": False,
    }
    write_json(run_dir / "reports" / "full_scan_summary.json", ordinary)
    return ordinary


def call_known_alleles(
    bam: Path,
    output: Path,
    reference: Path,
    sensitive: Path,
    tools: dict[str, Path],
    runner: CommandRunner,
) -> float:
    mpileup = [
        str(tools["bcftools"]),
        "mpileup",
        "--no-version",
        "-A",
        "-Ou",
        "-f",
        str(reference),
        "-R",
        str(sensitive / "probes" / "targets.bed"),
        "-q",
        "20",
        "-Q",
        "20",
        "--ns",
        "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
        "-a",
        "FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR",
        str(bam),
    ]
    call = [
        str(tools["bcftools"]),
        "call",
        "--no-version",
        "-m",
        "-A",
        "-C",
        "alleles",
        "-T",
        str(sensitive / "probes" / "known_alleles.tsv.gz"),
        "-i",
        "-Ov",
        "-o",
        str(output),
    ]
    wall = runner.pipeline(mpileup, call)
    output.chmod(0o600)
    return wall


def minimal_allele_key(contig: str, pos1: int, ref: str, alt: str) -> tuple[str, int, str, str]:
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos1 += 1
    return contig, pos1, ref, alt


def parse_vcf(path: Path) -> dict[tuple[str, int, str, str], dict[str, object]]:
    records: dict[tuple[str, int, str, str], dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            alts = fields[4].split(",")
            format_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            sample = dict(zip(format_keys, sample_values))

            def numbers(tag: str) -> list[int]:
                value = sample.get(tag, ".")
                return [] if value in {"", "."} else [0 if item == "." else int(item) for item in value.split(",")]

            ad = numbers("AD")
            adf = numbers("ADF")
            adr = numbers("ADR")
            for index, alt in enumerate(alts, start=1):
                records[minimal_allele_key(fields[0], int(fields[1]), fields[3], alt)] = {
                    "dp": int(sample.get("DP", "0")) if sample.get("DP", ".") not in {"", "."} else 0,
                    "ref_count": ad[0] if ad else 0,
                    "alt_count": ad[index] if index < len(ad) else 0,
                    "ref_forward": adf[0] if adf else 0,
                    "alt_forward": adf[index] if index < len(adf) else 0,
                    "ref_reverse": adr[0] if adr else 0,
                    "alt_reverse": adr[index] if index < len(adr) else 0,
                }
    return records


def run_calls(args: argparse.Namespace, run_dir: Path, sensitive: Path, runner: CommandRunner) -> dict[str, object]:
    contract = read_json(args.contract.resolve())
    tools = tool_paths(contract)
    marked = sensitive / "alignments" / "merged_marked.bam"
    preduplicate = sensitive / "alignments" / "merged_preduplicate.bam"
    if not marked.is_file() or not preduplicate.is_file():
        raise GateFailure("calling_alignment_predecessor_failed")
    if not Path(str(preduplicate) + ".bai").is_file():
        runner.run([str(tools["samtools"]), "index", "-@", "4", str(preduplicate)])
    calls_dir = sensitive / "calls"
    primary = calls_dir / "primary_known_alleles.vcf"
    unfiltered = calls_dir / "preduplicate_known_alleles.vcf"
    timings = {
        "primary": call_known_alleles(marked, primary, args.reference.resolve(), sensitive, tools, runner),
        "preduplicate": call_known_alleles(preduplicate, unfiltered, args.reference.resolve(), sensitive, tools, runner),
    }
    lane_call_counts: dict[str, int] = {}
    for lane_row in read_tsv(sensitive / "inputs" / "fastq_manifest.tsv"):
        lane = lane_row["lane_id"]
        lane_bam = calls_dir / f"{lane}.marked.bam"
        lane_vcf = calls_dir / f"{lane}.known_alleles.vcf"
        runner.run([str(tools["samtools"]), "view", "-b", "-r", lane, "-o", str(lane_bam), str(marked)])
        runner.run([str(tools["samtools"]), "index", str(lane_bam)])
        timings[lane] = call_known_alleles(lane_bam, lane_vcf, args.reference.resolve(), sensitive, tools, runner)
        lane_call_counts[lane] = len(parse_vcf(lane_vcf))
    primary_records = parse_vcf(primary)
    allowed = {
        (row["reference_chrom"], int(row["pos1"]), row["ref"], row["alt"])
        for row in read_tsv(sensitive / "inputs" / "candidate_alleles.tsv")
    }
    novel = set(primary_records) - allowed
    if novel:
        raise GateFailure("novel_allele_emission_failed")
    protected = {
        "status": "PASS",
        "approved_allele_count": len(allowed),
        "primary_vcf_record_count": len(primary_records),
        "novel_allele_count": len(novel),
        "lane_vcf_record_counts": lane_call_counts,
        "wall_seconds": {key: round(value, 3) for key, value in timings.items()},
    }
    write_json(calls_dir / "call_summary.json", protected)
    ordinary = {
        "status": "PASS",
        "approved_allele_count": len(allowed),
        "constrained_call_record_count": len(primary_records),
        "novel_allele_count": 0,
        "lane_count": len(lane_call_counts),
        "known_allele_constraint": True,
        "patient_identity_in_ordinary_output": False,
    }
    write_json(run_dir / "reports" / "constrained_call_summary.json", ordinary)
    return ordinary


def truth(value: object) -> bool:
    return str(value).lower() in {"true", "1", "yes"}


def cigar_operations(cigar: str) -> list[tuple[int, str]]:
    operations = [(int(length), operation) for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
    if not operations or "".join(f"{length}{operation}" for length, operation in operations) != cigar:
        return []
    return operations


def observe_allele(fields: list[str], allele: dict[str, object]) -> tuple[str, int, int, bool] | None:
    flag = int(fields[1])
    if flag & PRIMARY_EXCLUDED_FLAGS or int(fields[4]) < 20 or fields[5] == "*" or fields[9] == "*":
        return None
    operations = cigar_operations(fields[5])
    if not operations:
        return None
    reference_start = int(fields[3]) - 1
    reference_end = reference_start + sum(length for length, operation in operations if operation in "MDN=X")
    target_start = int(allele["pos1"]) - 1
    target_end = target_start + len(str(allele["ref"]))
    if reference_start > target_start or reference_end < target_end:
        return None
    sequence = fields[9]
    qualities = [ord(value) - 33 for value in fields[10]] if fields[10] != "*" else []
    query_index = 0
    reference_index = reference_start
    observed: list[str] = []
    observed_qualities: list[int] = []
    for length, operation in operations:
        if operation in "M=X":
            for offset in range(length):
                position = reference_index + offset
                if target_start <= position < target_end:
                    observed.append(sequence[query_index + offset])
                    if qualities:
                        observed_qualities.append(qualities[query_index + offset])
            query_index += length
            reference_index += length
        elif operation == "I":
            if target_start < reference_index <= target_end:
                observed.extend(sequence[query_index : query_index + length])
                if qualities:
                    observed_qualities.extend(qualities[query_index : query_index + length])
            query_index += length
        elif operation in "DN":
            reference_index += length
        elif operation == "S":
            query_index += length
        elif operation in "HP":
            continue
    allele_string = "".join(observed).upper()
    if not observed_qualities or min(observed_qualities) < 20:
        return None
    if allele_string == str(allele["ref"]):
        state = "REF"
    elif allele_string == str(allele["alt"]):
        state = "ALT"
    else:
        return None
    clipped = operations[0][1] == "S" or operations[-1][1] == "S"
    return state, int(fields[4]), min(observed_qualities), clipped


def sam_records_for_region(
    marked_bam: Path,
    contig: str,
    start1: int,
    end1: int,
    samtools: Path,
    runner: CommandRunner,
) -> list[list[str]]:
    text = runner.capture([str(samtools), "view", str(marked_bam), f"{contig}:{start1}-{end1}"])
    return [line.split("\t") for line in text.splitlines() if line]


def fragment_observations(
    allele: dict[str, object],
    marked_bam: Path,
    samtools: Path,
    runner: CommandRunner,
) -> dict[str, object]:
    start1 = int(allele["pos1"])
    end1 = start1 + len(str(allele["ref"])) - 1
    records = sam_records_for_region(marked_bam, str(allele["reference_chrom"]), start1, end1, samtools, runner)
    by_fragment: dict[str, list[tuple[str, int, int, bool, str, str, bool]]] = defaultdict(list)
    for fields in records:
        observation = observe_allele(fields, allele)
        if observation:
            tags = {item.split(":", 2)[0]: item.split(":", 2)[2] for item in fields[11:] if item.count(":") >= 2}
            strand = "REVERSE" if int(fields[1]) & 0x10 else "FORWARD"
            by_fragment[fields[0]].append((*observation, strand, tags.get("RG", "UNKNOWN"), bool(int(fields[1]) & 0x40)))
    retained: list[tuple[str, int, int, bool, str, str, bool]] = []
    for observations in by_fragment.values():
        states = {value[0] for value in observations}
        if len(states) != 1:
            continue
        retained.append(max(observations, key=lambda value: (value[6], value[2], value[1])))
    alt = [value for value in retained if value[0] == "ALT"]
    ref = [value for value in retained if value[0] == "REF"]
    alt_lanes = Counter(value[5] for value in alt)
    return {
        "fragment_ref": len(ref),
        "fragment_alt": len(alt),
        "fragment_total": len(retained),
        "fragment_alt_forward": sum(value[4] == "FORWARD" for value in alt),
        "fragment_alt_reverse": sum(value[4] == "REVERSE" for value in alt),
        "fragment_supporting_lanes": len(alt_lanes),
        "fragment_lane_counts_json": json.dumps(dict(sorted(alt_lanes.items())), sort_keys=True, separators=(",", ":")),
        "alt_median_mapq": statistics.median(value[1] for value in alt) if alt else "",
        "alt_median_min_baseq": statistics.median(value[2] for value in alt) if alt else "",
        "alt_clipped_fragments": sum(value[3] for value in alt),
    }


def direct_phase(
    allele1: dict[str, object],
    allele2: dict[str, object],
    marked_bam: Path,
    samtools: Path,
    runner: CommandRunner,
) -> dict[str, object]:
    if allele1["reference_chrom"] != allele2["reference_chrom"]:
        return {"phase_label": "PHASE_UNRESOLVED_NO_PARENTS", "phase_reason": "DIFFERENT_CONTIGS", "haplotypes": {}}
    first = min(int(allele1["pos1"]), int(allele2["pos1"]))
    last = max(
        int(allele1["pos1"]) + len(str(allele1["ref"])) - 1,
        int(allele2["pos1"]) + len(str(allele2["ref"])) - 1,
    )
    if last - first > 1000:
        return {"phase_label": "PHASE_UNRESOLVED_NO_PARENTS", "phase_reason": "OUTSIDE_DIRECT_FRAGMENT_WINDOW", "haplotypes": {}}
    records = sam_records_for_region(marked_bam, str(allele1["reference_chrom"]), first, last, samtools, runner)
    fragments: dict[str, dict[str, list[str]]] = defaultdict(lambda: {"a1": [], "a2": []})
    for fields in records:
        flag = int(fields[1])
        if not flag & 0x2 or flag & PRIMARY_EXCLUDED_FLAGS:
            continue
        observation1 = observe_allele(fields, allele1)
        observation2 = observe_allele(fields, allele2)
        if observation1:
            fragments[fields[0]]["a1"].append(observation1[0])
        if observation2:
            fragments[fields[0]]["a2"].append(observation2[0])
    haplotypes: Counter[str] = Counter()
    for values in fragments.values():
        states1 = set(values["a1"])
        states2 = set(values["a2"])
        if len(states1) == 1 and len(states2) == 1:
            haplotypes[f"{next(iter(states1))}_{next(iter(states2))}"] += 1
    label = classify_phase(haplotypes)
    reason = "DIRECT_FRAGMENT_RULE_APPLIED" if haplotypes else "NO_QUALIFIED_SHARED_FRAGMENT"
    return {"phase_label": label, "phase_reason": reason, "haplotypes": dict(sorted(haplotypes.items()))}


def evidence_rationale(values: dict[str, object], label: str) -> str:
    details = [
        f"tier={label}",
        f"qualified_depth={values['qualified_depth']}",
        f"independent_alt={values['independent_alt']}",
        f"alt_forward={values['alt_forward']}",
        f"alt_reverse={values['alt_reverse']}",
        f"supporting_lanes={values['supporting_lanes']}",
        f"mapping_caveat={str(values['mapping_caveat']).lower()}",
    ]
    return ";".join(details)


def build_evidence(
    args: argparse.Namespace,
    sensitive: Path,
    output_dir: Path,
    runner: CommandRunner,
) -> dict[str, object]:
    contract = read_json(args.contract.resolve())
    tools = tool_paths(contract)
    allele_rows = read_tsv(sensitive / "inputs" / "candidate_alleles.tsv")
    pair_rows = read_tsv(sensitive / "inputs" / "candidate_pairs.tsv")
    source_rows = read_tsv(sensitive / "inputs" / "task05_handoff_snapshot.tsv.gz")
    if len(allele_rows) != EXPECTED_ALLELE_COUNT or len(pair_rows) != EXPECTED_PAIR_COUNT or len(source_rows) != EXPECTED_PAIR_COUNT:
        raise GateFailure("evidence_input_count_failed")
    primary = parse_vcf(sensitive / "calls" / "primary_known_alleles.vcf")
    preduplicate = parse_vcf(sensitive / "calls" / "preduplicate_known_alleles.vcf")
    lane_calls = {
        row["lane_id"]: parse_vcf(sensitive / "calls" / f"{row['lane_id']}.known_alleles.vcf")
        for row in read_tsv(sensitive / "inputs" / "fastq_manifest.tsv")
    }
    marked = sensitive / "alignments" / "merged_marked.bam"
    evidence: list[dict[str, object]] = []
    evidence_by_id: dict[str, dict[str, object]] = {}
    allele_by_id = {row["allele_id"]: row for row in allele_rows}
    for allele in allele_rows:
        key = (allele["reference_chrom"], int(allele["pos1"]), allele["ref"], allele["alt"])
        qualified = primary.get(key, {"dp": 0, "ref_count": 0, "alt_count": 0, "ref_forward": 0, "alt_forward": 0, "ref_reverse": 0, "alt_reverse": 0})
        raw = preduplicate.get(key, {"dp": 0, "ref_count": 0, "alt_count": 0, "ref_forward": 0, "alt_forward": 0, "ref_reverse": 0, "alt_reverse": 0})
        lane_values = {lane: records.get(key, {"dp": 0, "alt_count": 0}) for lane, records in lane_calls.items()}
        fragment = fragment_observations(allele, marked, tools["samtools"], runner)
        values: dict[str, object] = {
            "reference_ok": truth(allele["reference_ok"]),
            "probe_ok": truth(allele["probe_ok"]),
            "representation_ok": truth(allele["representation_ok"]),
            "retrieval_ok": int(fragment["fragment_total"]) > 0,
            "qualified_depth": int(fragment["fragment_total"]),
            "independent_alt": int(fragment["fragment_alt"]),
            "alt_forward": int(fragment["fragment_alt_forward"]),
            "alt_reverse": int(fragment["fragment_alt_reverse"]),
            "supporting_lanes": int(fragment["fragment_supporting_lanes"]),
            "mapping_caveat": truth(allele["mapping_caveat"]),
        }
        label = classify_allele(values)
        original_dp = int(float(allele["original_dp"])) if allele["original_dp"] else 0
        bcftools_discordant = int(fragment["fragment_alt"]) != int(qualified["alt_count"])
        reconciliation = "CONCORDANT"
        if bcftools_discordant:
            reconciliation = "INDEL_PADDED_AD_DISCORDANCE" if allele["variant_class"] in {"INSERTION", "DELETION"} else "COUNT_DISCORDANCE_REVIEW"
        row: dict[str, object] = {
            "allele_id": allele["allele_id"],
            **values,
            "qualified_ref": fragment["fragment_ref"],
            "qualified_alt_fraction": int(fragment["fragment_alt"]) / int(fragment["fragment_total"]) if int(fragment["fragment_total"]) else "",
            "bcftools_qualified_depth": qualified["dp"],
            "bcftools_qualified_ref": qualified["ref_count"],
            "bcftools_qualified_alt": qualified["alt_count"],
            "bcftools_alt_forward": qualified["alt_forward"],
            "bcftools_alt_reverse": qualified["alt_reverse"],
            "bcftools_supporting_lanes": sum(int(lane.get("alt_count", 0)) > 0 for lane in lane_values.values()),
            "bcftools_lane_counts_json": json.dumps(lane_values, sort_keys=True, separators=(",", ":")),
            "bcftools_preduplicate_depth": raw["dp"],
            "bcftools_preduplicate_ref": raw["ref_count"],
            "bcftools_preduplicate_alt": raw["alt_count"],
            "bcftools_fragment_alt_discordant": bcftools_discordant,
            "representation_reconciliation": reconciliation,
            "original_vcf_dp": original_dp,
            "original_vcf_ad": allele["original_ad"],
            "retrieval_dp_ratio": int(fragment["fragment_total"]) / original_dp if original_dp else "",
            **fragment,
            "read_disposition": label,
            "disposition_rationale": evidence_rationale(values, label),
        }
        evidence.append(row)
        evidence_by_id[str(allele["allele_id"])] = row
    evidence_fields = [
        "allele_id", "reference_ok", "probe_ok", "representation_ok", "retrieval_ok", "qualified_depth", "qualified_ref",
        "independent_alt", "qualified_alt_fraction", "alt_forward", "alt_reverse", "supporting_lanes", "mapping_caveat",
        "bcftools_qualified_depth", "bcftools_qualified_ref", "bcftools_qualified_alt", "bcftools_alt_forward", "bcftools_alt_reverse",
        "bcftools_supporting_lanes", "bcftools_lane_counts_json", "bcftools_preduplicate_depth", "bcftools_preduplicate_ref",
        "bcftools_preduplicate_alt", "bcftools_fragment_alt_discordant", "representation_reconciliation", "original_vcf_dp", "original_vcf_ad",
        "retrieval_dp_ratio", "fragment_ref", "fragment_alt", "fragment_total", "fragment_alt_forward", "fragment_alt_reverse",
        "fragment_supporting_lanes", "fragment_lane_counts_json", "alt_median_mapq", "alt_median_min_baseq",
        "alt_clipped_fragments", "read_disposition", "disposition_rationale",
    ]
    write_tsv(output_dir / "allele_read_evidence.tsv", evidence_fields, evidence)
    pair_evidence: list[dict[str, object]] = []
    phase_by_pair: dict[str, dict[str, object]] = {}
    for pair in pair_rows:
        first = evidence_by_id[pair["allele_1_id"]]
        second = evidence_by_id[pair["allele_2_id"]]
        pair_label = classify_pair([str(first["read_disposition"]), str(second["read_disposition"])])
        phase = direct_phase(
            allele_by_id[pair["allele_1_id"]],
            allele_by_id[pair["allele_2_id"]],
            marked,
            tools["samtools"],
            runner,
        )
        phase_by_pair[pair["pair_id"]] = phase
        pair_evidence.append(
            {
                "pair_id": pair["pair_id"],
                "allele_1_id": pair["allele_1_id"],
                "allele_2_id": pair["allele_2_id"],
                "allele_1_read_disposition": first["read_disposition"],
                "allele_2_read_disposition": second["read_disposition"],
                "pair_read_disposition": pair_label,
                "phase_label": phase["phase_label"],
                "phase_reason": phase["phase_reason"],
                "haplotype_counts_json": json.dumps(phase["haplotypes"], sort_keys=True, separators=(",", ":")),
                "parental_origin": "UNKNOWN_NO_PARENTS",
                "task05_rank_changed": False,
            }
        )
    pair_fields = [
        "pair_id", "allele_1_id", "allele_2_id", "allele_1_read_disposition", "allele_2_read_disposition",
        "pair_read_disposition", "phase_label", "phase_reason", "haplotype_counts_json", "parental_origin", "task05_rank_changed",
    ]
    write_tsv(output_dir / "pair_read_evidence.tsv", pair_fields, pair_evidence)
    appended = [
        "task06_allele_1_read_disposition", "task06_allele_2_read_disposition", "task06_pair_read_disposition",
        "task06_allele_1_qualified_depth", "task06_allele_2_qualified_depth", "task06_allele_1_independent_alt",
        "task06_allele_2_independent_alt", "task06_phase_label", "task06_phase_reason", "task06_haplotype_counts_json",
        "task06_parental_origin", "task06_fastq_read", "task06_rank_changed",
    ]
    handoff_rows: list[dict[str, object]] = []
    for source, pair, pair_result in zip(source_rows, pair_rows, pair_evidence, strict=True):
        first = evidence_by_id[pair["allele_1_id"]]
        second = evidence_by_id[pair["allele_2_id"]]
        updated: dict[str, object] = dict(source)
        updated.update(
            {
                "task06_allele_1_read_disposition": first["read_disposition"],
                "task06_allele_2_read_disposition": second["read_disposition"],
                "task06_pair_read_disposition": pair_result["pair_read_disposition"],
                "task06_allele_1_qualified_depth": first["qualified_depth"],
                "task06_allele_2_qualified_depth": second["qualified_depth"],
                "task06_allele_1_independent_alt": first["independent_alt"],
                "task06_allele_2_independent_alt": second["independent_alt"],
                "task06_phase_label": pair_result["phase_label"],
                "task06_phase_reason": pair_result["phase_reason"],
                "task06_haplotype_counts_json": pair_result["haplotype_counts_json"],
                "task06_parental_origin": "UNKNOWN_NO_PARENTS",
                "task06_fastq_read": True,
                "task06_rank_changed": False,
            }
        )
        handoff_rows.append(updated)
    handoff_fields = list(source_rows[0]) + appended
    write_gzip_tsv(output_dir / "task07_read_evidence_handoff.tsv.gz", handoff_fields, handoff_rows)
    summary = {
        "status": "PASS",
        "allele_count": len(evidence),
        "pair_count": len(pair_evidence),
        "allele_disposition_counts": dict(sorted(Counter(str(row["read_disposition"]) for row in evidence).items())),
        "pair_disposition_counts": dict(sorted(Counter(str(row["pair_read_disposition"]) for row in pair_evidence).items())),
        "phase_disposition_counts": dict(sorted(Counter(str(row["phase_label"]) for row in pair_evidence).items())),
        "technically_evaluable_alleles": sum(str(row["read_disposition"]) != "READ_NONCALL" for row in evidence),
        "supported_alleles": sum(str(row["read_disposition"]) in {"STRONG_READ_SUPPORT", "READ_SUPPORT_WITH_CAVEAT"} for row in evidence),
        "pairs_with_both_alleles_supported": sum(str(row["pair_read_disposition"]) in {"PAIR_BOTH_ALLELES_STRONG", "PAIR_BOTH_ALLELES_SUPPORTED"} for row in pair_evidence),
        "bcftools_fragment_alt_discordance_count": sum(bool(row["bcftools_fragment_alt_discordant"]) for row in evidence),
        "indel_representation_reconciliation_count": sum(row["representation_reconciliation"] == "INDEL_PADDED_AD_DISCORDANCE" for row in evidence),
        "parental_origin_inferred": False,
        "task05_rank_changed": False,
        "novel_alleles_emitted": 0,
    }
    write_json(output_dir / "evidence_summary.json", summary)
    return summary


def evidence_stage(args: argparse.Namespace, run_dir: Path, sensitive: Path, runner: CommandRunner) -> dict[str, object]:
    summary = build_evidence(args, sensitive, sensitive / "evidence", runner)
    ordinary = {
        **summary,
        "patient_identity_in_ordinary_output": False,
        "read_names_in_ordinary_output": False,
        "pathogenicity_claim": False,
        "causality_claim": False,
        "parental_origin_claim": False,
    }
    write_json(run_dir / "reports" / "read_evidence_summary.json", ordinary)
    return ordinary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("prepare", "smoke", "prefilter", "full", "call", "evidence", "replay"), required=True)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--candidate-handoff", required=True, type=Path)
    parser.add_argument("--fastq-contract", required=True, type=Path)
    parser.add_argument("--fastq-integrity", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    args = parser.parse_args()
    os.umask(0o077)
    run_dir = args.run_dir.resolve()
    sensitive = run_dir / "data" / "sensitive"
    if not run_dir.is_dir() or not sensitive.is_dir():
        print(json.dumps({"status": "FAIL", "error_type": "RUN_LAYOUT_MISSING"}, sort_keys=True))
        return 1
    for directory in [sensitive, *[path for path in sensitive.rglob("*") if path.is_dir()]]:
        directory.chmod(0o700)
    runner = CommandRunner(sensitive / "logs")
    try:
        if args.stage == "prepare":
            result = prepare(args, run_dir, sensitive, runner)
        elif args.stage == "smoke":
            result = smoke(args, run_dir, sensitive, runner)
        elif args.stage == "prefilter":
            result = prefilter_validation(args, run_dir, sensitive, runner)
        elif args.stage == "full":
            result = full_scan(args, run_dir, sensitive, runner)
        elif args.stage == "call":
            result = run_calls(args, run_dir, sensitive, runner)
        elif args.stage == "evidence":
            result = evidence_stage(args, run_dir, sensitive, runner)
        else:
            result = build_evidence(args, sensitive, sensitive / "replay", runner)
        print(json.dumps({"status": result["status"], "stage": args.stage}, sort_keys=True))
        return 0
    except Exception as error:
        failure = {
            "status": "FAIL",
            "stage": args.stage,
            "error_type": type(error).__name__,
            "error_code": str(error) if isinstance(error, GateFailure) else "unexpected_exception",
        }
        write_json(sensitive / "logs" / f"{args.stage}_failure.json", {**failure, "traceback": traceback.format_exc()})
        write_json(run_dir / "reports" / f"{args.stage}_failure.json", failure)
        print(json.dumps(failure, sort_keys=True))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
