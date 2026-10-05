"""Exact source excerpts; see canonical_provenance.json."""
import collections
import datetime
import hashlib
import json
import math
import re
from urllib.parse import unquote

QUEUES = ['IDENTITY_OR_CONFLICT_RESOLUTION', 'CONTRADICTED_TRANS_PAIR', 'FREQUENCY_OR_CLINVAR_CONTEXT',
          'CALL_QUALITY_REVIEW', 'AR_ALLELE_MECHANISM_REVIEW', 'GENE_DISEASE_OR_TARGET_LINK_REVIEW']

DISPLAY = [QUEUES[i] for i in [0,4,3,5,2,1]]

LABELS = dict(zip(QUEUES,['ID·상충 근거 확인','해당 trans 쌍에 반대 근거','빈도·ClinVar 적용 맥락','호출 품질 확인','AR allele 기전 검토','유전자–질환·표적 연결 확인']))

CHANGE = dict(zip(QUEUES,[
    'Resolve the exact identifier/genotype/source conflict; then reapply every axis. A resolved bridge does not establish pathogenicity.',
    'Validated correction of the local phase or evaluable read/genotype contradiction could reopen this exact pair; an alternative allele or mechanism is a separate hypothesis.',
    'Disease-specific maximum credible AF, ancestry/coverage and documented exceptions, or condition/SCV-specific ClinVar reconciliation could change this workflow disposition.',
    'Independent alternate support or corrected original QC could move the pair into mechanism/target-link review; adequate validated reference-only evidence can contradict it.',
    'Disease-compatible effects for both alleles and reliable trans evidence strengthen the hypothesis; incompatible function, phenotype, or reliable cis weakens it. No automatic causal promotion.',
    'An exact gene-disease inheritance link and a validated variant-to-gene/protein/splice mechanism for the unresolved allele(s) could justify focused AR review.'
]))

def assess(r):
    flags=[]; missing=[]; questions=[]; identity=False; nonpair=False
    if r.get('source_reconciled')!='true' or r.get('gene_id_resolved')!='true':
        identity=True; flags.append('UNRESOLVED_IDENTITY')
    if r.get('gene_bridge_conflicts',''):
        identity=True; flags.append('GENE_BRIDGE_CONFLICT')
    if r.get('audit_phase')=='PHASE_CONFLICT': identity=True; flags.append('PHASE_CONFLICT')
    contra=r.get('audit_phase')=='LOCAL_CIS_SUPPORTED'
    if contra: flags.append('LOCAL_CIS_SUPPORTS_SAME_COPY')
    if r.get('audit_phase')=='UNRESOLVED': missing.append('TRANS_NOT_ESTABLISHED'); questions.append('Establish allele phase using reliable segregation or connected haplotypes; no parental origin is inferred.')
    frequency=False; poor=False; effects=[]
    for i in [1,2]:
        pre=f'audit_{i}_'; gt=r.get(pre+'normalized_gt','')
        if not re.fullmatch(r'[01](?:[/|][01])?',gt): identity=True; flags.append(f'A{i}_GENOTYPE_MISSING_OR_INVALID')
        elif gt not in {'0/1','1/0','0|1','1|0'}: nonpair=True; flags.append(f'A{i}_NOT_DIPLOID_HET')
        if r.get(pre+'evidence_conflicts',''):
            identity=True; flags.append(f'A{i}_'+r[pre+'evidence_conflicts'])
        read=r.get(pre+'read_label','READ_NOT_ASSESSED')
        if read=='READ_CONTRADICTED' and r.get(pre+'read_contradiction_evaluable')=='true': contra=True; flags.append(f'A{i}_EVALUABLE_READ_CONTRADICTION')
        if read in {'READ_NOT_ASSESSED','READ_NONCALL'}: missing.append(f'A{i}_'+read)
        if r.get(pre+'frequency_flag')=='UNKNOWN': missing.append(f'A{i}_AF_UNKNOWN'); questions.append(f'Allele {i}: verify population coverage/representation before inferring rarity.')
        if r.get(pre+'frequency_flag')=='ABOVE_BROAD_SCREEN': frequency=True; flags.append(f'A{i}_ABOVE_BROAD_AF_SCREEN')
        if r.get(pre+'clinvar_class')=='BLB_NONCONFLICTING_STARRED': frequency=True; flags.append(f'A{i}_RECORDED_CLINVAR_BLB')
        good=r.get(pre+'high_confidence_call')=='true' or read=='STRONG_READ_SUPPORT'
        if not good: poor=True; flags.append(f'A{i}_CALL_QC_REQUIRES_REVIEW'); questions.append(f'Allele {i}: inspect retrieval, depth and independent alternate support; an unassessed locus is not ref-only.')
        functional=r.get(pre+'gene_protein_splice')=='true'
        linked_plp=r.get(pre+'clinvar_class')=='PLP_NONCONFLICTING' and int(r.get(pre+'annotation_rows','0'))>0
        effects.append(functional or linked_plp)
        if not effects[-1]: missing.append(f'A{i}_GENE_SPECIFIC_EFFECT_UNRESOLVED'); questions.append(f'Allele {i}: establish the exact transcript/regulatory target and disease-compatible effect.')
        if r.get(pre+'transcript_consequence_discordance')=='true': flags.append(f'A{i}_TRANSCRIPT_CONSEQUENCES_DIFFER'); questions.append(f'Allele {i}: compare disease-relevant transcript with the recorded MANE/canonical effect; an alternative transcript consequence is not a functional assay.')
    ar=r.get('positive_ar_association')=='true'
    if not ar: missing.append('POSITIVE_AR_LINK_NOT_DOCUMENTED'); questions.append('Resolve disease/inheritance or target-gene association; a missing frozen association is not benign evidence.')
    if identity: queue=QUEUES[0]
    elif contra or nonpair: queue=QUEUES[1]
    elif frequency: queue=QUEUES[2]
    elif poor: queue=QUEUES[3]
    elif ar and all(effects): queue=QUEUES[4]
    else: queue=QUEUES[5]
    questions.insert(0,CHANGE[queue])
    if queue==QUEUES[4]: questions.append('Review exact phenotype/onset and disease mechanism, each allele/condition assertion and relevant functional evidence; native model scores remain separate.')
    return {'workflow_queue':queue,'workflow_label_ko':LABELS[queue],'evidence_flags':';'.join(flags),
            'missing_evidence':';'.join(missing),'next_questions':' | '.join(dict.fromkeys(questions)),
            'observation_that_changes_queue':CHANGE[queue], 'review_state':'ALL_RECORDED_AXES_AUDITED_EXPERT_ADJUDICATION_PENDING',
            'within_queue_order':'STABLE_IDENTIFIER_ORDER_NOT_BIOLOGICAL_RANK','causal_probability':'NOT_ESTIMATED',
            'clinical_pathogenicity_classification':'NOT_ASSIGNED'}

def classify_all(items):
    result=[{**r,**assess(r)} for r in items]
    return sorted(result,key=lambda r:(DISPLAY.index(r['workflow_queue']),r['gene_id'],r['pair_key']))
