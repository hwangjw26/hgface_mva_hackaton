"""Exact source excerpts; see canonical_provenance.json."""
import collections
import datetime
import hashlib
import json
import math
import re
from urllib.parse import unquote

def require(ok, message):
    if not ok:
        raise ValueError(message)

def condition_ids(values):
    require(isinstance(values, list) and all(isinstance(v,str) for v in values), 'INVALID_CONDITION_ARRAY')
    ids = set()
    for value in values:
        for token in re.split(r'[&|,;\s]+', unquote(value)):
            m = re.fullmatch(r'(OMIM|MONDO|Orphanet|ORPHA):(?:(OMIM|MONDO|Orphanet|ORPHA):)?([0-9]+)', token)
            if m and (m[2] is None or m[1] == m[2]):
                ids.add(m[1]+':'+m[3])
    return sorted(ids)

def condition_state(ids, target):
    if not target:
        return 'TARGET_DISEASE_NOT_RETURNED'
    if target in ids:
        return 'EXACT_ID_OVERLAP_NOT_ASSERTION_ADJUDICATION'
    return 'NO_COMPARABLE_CONDITION_ID' if not ids else 'NO_EXACT_ID_OVERLAP_NOT_EXCLUSION'
