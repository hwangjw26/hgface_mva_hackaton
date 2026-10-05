"""Exact source excerpts; see canonical_provenance.json."""
import collections
import datetime
import hashlib
import json
import math
import re
from urllib.parse import unquote

VERSION = 'task15_research_tiers_v2'

T1 = 'TIER_1_REVIEWED_ALLELE_SUPPORT'

T2 = 'TIER_2_MECHANISM_REVIEW'

T3 = 'TIER_3_EXPLORATORY'

HOLD = 'HOLD_FOR_REVIEW'

CONTRA = 'CONTRADICTED_PAIR'

ORDER = [T1, T2, T3, HOLD, CONTRA]

LABELS = {T1: 'Tier1: allele 근거 검토 완료', T2: 'Tier2: 기전 검토 후보 (잠정)',
          T3: 'Tier3: 탐색 가설', HOLD: '보류: 추가 확인 필요', CONTRA: '정확한 쌍 설명에 반박 근거'}

ROUTES = {'MECHANISM_CHECKED_LOF', 'VALIDATED_FUNCTIONAL_WITH_CLINICAL_CONTEXT',
          'PRIMARY_CASE_SEGREGATION', 'DISEASE_SCOPED_REVIEWED_CLASSIFICATION'}

HOLD_QUEUES = {'IDENTITY_OR_CONFLICT_RESOLUTION', 'FREQUENCY_OR_CLINVAR_CONTEXT', 'CALL_QUALITY_REVIEW'}

def require(ok, code):
    if not ok:
        raise ValueError(code)

def nonempty(value):
    return isinstance(value, str) and bool(value.strip())

def valid_refs(value):
    return isinstance(value, list) and bool(value) and all(nonempty(x) for x in value)

def validate_review(row, review):
    """Validate a completed assessment's binding/provenance, not its scientific truth."""
    require(review.get('gene_id') == row['gene_id'] and review.get('pair_key') == row['pair_key'], 'REVIEW_PAIR_MISMATCH')
    for key in ['disease_id', 'rationale', 'reviewer', 'review_date']:
        require(nonempty(review.get(key)), 'REVIEW_METADATA_MISSING')
    datetime.date.fromisoformat(review['review_date'])
    require(review.get('status') == 'COMPLETE', 'REVIEW_NOT_COMPLETE')
    for key in ['ar_mechanism_compatible', 'phenotype_compatible', 'population_compatible', 'no_material_contradiction']:
        require(review.get(key) is True, 'REVIEW_COMPATIBILITY_NOT_ESTABLISHED')
    require(valid_refs(review.get('source_refs')), 'REVIEW_SOURCES_MISSING')
    alleles = review.get('alleles')
    require(isinstance(alleles, list) and len(alleles) == 2, 'REVIEW_REQUIRES_TWO_ALLELES')
    require({a.get('allele_key') for a in alleles} == {row['allele_1_key'], row['allele_2_key']}, 'REVIEW_ALLELE_MISMATCH')
    for a in alleles:
        require(a.get('disease_id') == review['disease_id'], 'REVIEW_DISEASE_MISMATCH')
        require(a.get('route') in ROUTES, 'UNSUPPORTED_REVIEW_ROUTE')
        require(a.get('mechanism_compatible') is True and a.get('transcript_relevant') is True, 'ALLELE_CONTEXT_NOT_REVIEWED')
        require(all(nonempty(a.get(k)) for k in ['transcript', 'rationale']), 'ALLELE_REVIEW_METADATA_MISSING')
        require(valid_refs(a.get('source_refs')) and valid_refs(a.get('evidence_units')), 'ALLELE_PRIMARY_EVIDENCE_MISSING')

def review_index(items, records):
    universe = {(r['gene_id'], r['pair_key']): r for r in items}
    require(len(universe) == len(items), 'DUPLICATE_HYPOTHESIS')
    result = {}
    for review in records:
        key = (review.get('gene_id'), review.get('pair_key'))
        require(key in universe and key not in result, 'UNKNOWN_OR_DUPLICATE_REVIEW')
        validate_review(universe[key], review)
        result[key] = review
    return result

def classify(row, review=None):
    if review is not None:
        validate_review(row, review)
    queue = row['workflow_queue']
    require(queue in HOLD_QUEUES | {'CONTRADICTED_TRANS_PAIR', 'AR_ALLELE_MECHANISM_REVIEW', 'GENE_DISEASE_OR_TARGET_LINK_REVIEW'}, 'UNKNOWN_SOURCE_QUEUE')
    if queue in HOLD_QUEUES:
        tier, reason = HOLD, queue
    elif queue == 'CONTRADICTED_TRANS_PAIR':
        tier, reason = CONTRA, 'SOURCE_QUALIFIED_EXACT_PAIR_CONTRADICTION'
    elif review is not None and row['positive_ar_association'] == 'true':
        tier, reason = T1, 'BOTH_ALLELES_AND_SAME_DISEASE_REVIEWED'
    elif queue == 'AR_ALLELE_MECHANISM_REVIEW':
        tier, reason = T2, 'ALLELE_AND_DISEASE_ADJUDICATION_PENDING'
    else:
        tier, reason = T3, 'GENE_DISEASE_TARGET_OR_ALLELE_LINK_PENDING'
    result = {'v2_rule_version': VERSION, 'v2_tier': tier, 'v2_label_ko': LABELS[tier], 'v2_reason': reason,
              'v2_phase': row['audit_phase'], 'v2_trans_confirmed': str(row['audit_phase'] == 'LOCAL_TRANS_SUPPORTED').lower(),
              'v2_review_state': 'COMPLETE_REVIEW_RECORD_APPLIED' if tier == T1 else 'COMPUTATIONAL_TRIAGE_EXPERT_REVIEW_PENDING',
              'v2_completed_review_present': str(review is not None).lower(),
              'v2_required_review': '' if tier == T1 else 'Same-disease AR mechanism, phenotype and population context; each exact allele/transcript and primary evidence; reviewer/date/rationale; material contradictions.',
              'v2_causal_probability': 'NOT_ESTIMATED', 'v2_clinical_classification': 'NOT_ASSIGNED'}
    for i in [1, 2]:
        if review is not None:
            a = next(a for a in review['alleles'] if a['allele_key'] == row[f'allele_{i}_key'])
            state = 'REVIEWED_' + a['route']
        else:
            state = 'ANNOTATION_ONLY' if row[f'audit_{i}_gene_protein_splice'] == 'true' else 'EFFECT_UNRESOLVED'
            if row[f'audit_{i}_clinvar_class'] == 'PLP_NONCONFLICTING':
                state += '_RECORDED_PLP_NOT_ADJUDICATED'
        result[f'v2_allele_{i}_evidence'] = state
    return result

def classify_all(items, reviews):
    indexed = review_index(items, reviews)
    return sorted(({**r, **classify(r, indexed.get((r['gene_id'], r['pair_key'])))} for r in items),
                  key=lambda r: (ORDER.index(r['v2_tier']), r['gene_id'], r['pair_key']))
