"""Small explicit I/O helpers; patient records are never printed."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

csv.field_size_limit(10_000_000)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rows(path, comment=""):
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as handle:
        return list(csv.DictReader((line for line in handle if not comment or not line.startswith(comment)),
                                   delimiter="\t"))


def index(records, fields):
    result = {tuple(record[key] for key in fields): record for record in records}
    require(len(result) == len(records), "Duplicate input key")
    return result


def write_tsv(path, records, fields=None):
    fields = fields or list(dict.fromkeys(key for record in records for key in record))
    require(bool(fields), "Output schema missing")
    with Path(path).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def write_json(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, ensure_ascii=False)
        handle.write("\n")


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def annotated_records(path, wanted):
    records, schema = {}, None
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                schema = line.split("Format: ", 1)[1].rstrip().rstrip('\">').split("|")
            elif line.startswith("#CHROM"):
                require(len(line.rstrip().split("\t")) == 10, "Singleton annotated VCF required")
            elif not line.startswith("#"):
                c = line.rstrip("\n").split("\t")
                key = ":".join([c[0].removeprefix("chr"), c[1], c[3], c[4]])
                if key not in wanted:
                    continue
                require(key not in records and "," not in c[4], "Duplicate or unsplit annotation allele")
                info = dict(item.split("=", 1) for item in c[7].split(";") if "=" in item)
                require(schema is not None and "CSQ" in info, "VEP CSQ annotation missing")
                records[key] = {"filter": c[6], "format": dict(zip(c[8].split(":"), c[9].split(":"))),
                                "annotations": [dict(zip(schema, item.split("|"))) for item in info["CSQ"].split(",")]}
    require(set(records) == wanted, "Requested allele annotation coverage incomplete")
    return records
