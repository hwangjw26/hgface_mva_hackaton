"""Exact source excerpts; see canonical_provenance.json."""
import collections
import datetime
import hashlib
import json
import math
import re
from urllib.parse import unquote

MISSING = {'', '.', 'NA', 'None', 'null'}

def frequency_flag(value, threshold):
    if value in MISSING:
        return 'UNKNOWN', ''
    try:
        af = float(value)
    except ValueError:
        return 'UNKNOWN', 'INVALID_AF'
    if not math.isfinite(af) or not 0 <= af <= 1:
        return 'UNKNOWN', 'INVALID_AF'
    return ('ABOVE_BROAD_SCREEN' if af > threshold else 'WITHIN_BROAD_SCREEN'), ''

def read_consensus(observations):
    labels = {r['label'] for r in observations}
    support = labels & {'STRONG_READ_SUPPORT', 'READ_SUPPORT_WITH_CAVEAT'}
    contra = [r for r in observations if r['label']=='READ_CONTRADICTED']
    evaluable = any(str(r['record'].get('retrieval_ok','')).lower()=='true' for r in contra)
    if support and evaluable:
        return 'READ_CONFLICT', False, 'READ_SUPPORT_VS_CONTRADICTION'
    if support:
        return ('STRONG_READ_SUPPORT' if 'STRONG_READ_SUPPORT' in support else 'READ_SUPPORT_WITH_CAVEAT'), False, ''
    if evaluable:
        return 'READ_CONTRADICTED', True, ''
    if contra:
        return 'READ_NONCALL', False, 'UNEVALUABLE_READ_CONTRADICTION'
    return ('READ_NONCALL' if 'READ_NONCALL' in labels else 'READ_NOT_ASSESSED'), False, ''

def physical_phase(a, b):
    """Conservative use of unchanged caller phase tags; never parental origin."""
    if a['allele_key'].split(':')[0] != b['allele_key'].split(':')[0]:
        return 'UNRESOLVED'
    if a.get('source_alt_count') != '1' or b.get('source_alt_count') != '1':
        return 'UNRESOLVED'
    if any(x.get('normalized_gt') not in {'0/1', '1/0', '0|1', '1|0'} for x in [a, b]):
        return 'UNRESOLVED'
    if a.get('pid', '') in MISSING or a.get('pid') != b.get('pid'):
        return 'UNRESOLVED'
    if a.get('pgt') not in {'0|1', '1|0'} or b.get('pgt') not in {'0|1', '1|0'}:
        return 'UNRESOLVED'
    return 'LOCAL_CIS_SUPPORTED' if a['pgt'] == b['pgt'] else 'LOCAL_TRANS_SUPPORTED'
