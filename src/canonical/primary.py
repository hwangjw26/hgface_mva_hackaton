"""Exact source excerpts; see canonical_provenance.json."""
import collections
import datetime
import hashlib
import json
import math
import re
from urllib.parse import unquote

def require(ok, code):
    if not ok:
        raise ValueError(code)

def primary_genotypes(variants):
    grouped = collections.defaultdict(list)
    for row in variants:
        grouped[(row['preset'], 'NCBIGene:' + row['ENTREZ_GENE_ID'], row['MOI'])].append(row)
    result = []
    for (preset, gene, mode), records in grouped.items():
        keys = [(r['allele_key'], r['CONTRIBUTING_VARIANT'], r['GENOTYPE']) for r in records]
        require(len({x[0] for x in keys}) == len(keys), 'DUPLICATE_NATIVE_MODE_ALLELE')
        contributing = [r for r in records if r['CONTRIBUTING_VARIANT'] == '1']
        if mode == 'AR' and len(contributing) == 2 and all(r['GENOTYPE'] in {'0/1', '1/0', '0|1', '1|0'} for r in contributing):
            result.append((preset, gene, tuple(sorted(r['allele_key'] for r in contributing))))
    return result
