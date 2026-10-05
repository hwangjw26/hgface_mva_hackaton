"""Task09 Increment04: reconcile reused Task06 and new targeted read/phase evidence for the frozen primary pairs.

Label of record: exact Task06 reuse where available, otherwise the Task09 targeted scan (same unchanged Task06 method).
Phase of record: Task06 label for exactly re-evaluated pairs, else new direct-fragment label; VCF PGT/PID physical phase is
used only when direct evidence is unresolved and both alleles share one PID block. Counts only to stdout.
"""
import collections
import csv
import gzip
import io
import json
import os
from pathlib import Path
import sys

import argparse
import read_rules as rules


def read(p):
    with (gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)) as s:
        return list(csv.DictReader(s, delimiter='\t'))


def write(p, rows, fields):
    with open(p, 'xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as gz:
        s = io.TextIOWrapper(gz, encoding='utf-8', newline='')
        w = csv.DictWriter(s, delimiter='\t', fieldnames=fields, lineterminator='\n'); w.writeheader(); w.writerows(rows); s.flush(); s.detach()
    os.chmod(p, 0o600)


def keyed(evidence_dir, inputs_dir):
    source = read(inputs_dir / 'candidate_alleles.tsv')
    ids = {r['allele_id']: rules.allele_key(r['chrom'], r['pos1'], r['ref'], r['alt']) for r in source}
    if len(ids) != len(source) or len(set(ids.values())) != len(ids):
        raise ValueError('Duplicate source allele identity')
    raw = read(evidence_dir / 'allele_read_evidence.tsv')
    alleles = {ids[r['allele_id']]: r for r in raw}
    if len(alleles) != len(raw):
        raise ValueError('Duplicate read evidence identity')
    raw = read(evidence_dir / 'pair_read_evidence.tsv')
    pairs = {rules.pair_key(ids[r['allele_1_id']], ids[r['allele_2_id']]): r for r in raw}
    if len(pairs) != len(raw):
        raise ValueError('Duplicate read evidence identity')
    return alleles, pairs


def vcf_phase(a, b):
    if not a['pid'] or not b['pid'] or a['pid'] != b['pid'] or not a['pgt'] or not b['pgt']:
        return 'NO_SHARED_VCF_PHASE_BLOCK'
    return 'VCF_PID_BLOCK_CIS' if a['pgt'] == b['pgt'] else 'VCF_PID_BLOCK_TRANS'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--source-alleles', required=True)
    parser.add_argument('--initial-read-dir', required=True)
    args = parser.parse_args()
    run = Path(args.run_dir).resolve(); sens = run / 'data/sensitive'; os.umask(0o077)
    outputs = [sens / 'allele_read_ledger.tsv.gz', sens / 'pair_read_phase_ledger.tsv.gz', run / 'data/read_phase_summary.json']
    if any(p.exists() for p in outputs):
        raise ValueError('Refusing existing read reconciliation output')
    targets = read(sens / 'read_target_pairs.tsv')
    source_rows = read(Path(args.source_alleles))
    source = {a['normalized_allele_key']: a for a in source_rows}
    if len(source) != len(source_rows):
        raise ValueError('Duplicate source annotation allele')
    new_a, new_p = keyed(sens / 'evidence', sens / 'inputs')
    old_a, old_p = keyed(Path(args.initial_read_dir) / 'evidence', Path(args.initial_read_dir) / 'inputs')
    keys = sorted({t[f'allele_{i}_normalized_allele_key'] for t in targets for i in (1, 2)}, key=rules.order_tuple)
    allele_rows = []
    for k in keys:
        n, o = new_a[k], old_a.get(k)
        record, src = (o['read_disposition'], 'TASK06_EXACT_REUSE') if o else (n['read_disposition'], 'TASK09_TARGETED_SCAN')
        allele_rows.append({'allele_key': k, 'label_of_record': record, 'label_source': src,
                            'task09_scan_label': n['read_disposition'], 'task06_label': o['read_disposition'] if o else '',
                            'reuse_concordant': ('' if not o else str(o['read_disposition'] == n['read_disposition']).lower()),
                            'retrieval_ok': n['retrieval_ok'], 'qualified_depth': n['qualified_depth'], 'independent_alt': n['independent_alt'],
                            'qualified_alt_fraction': n['qualified_alt_fraction'], 'alt_forward': n['alt_forward'], 'alt_reverse': n['alt_reverse'],
                            'supporting_lanes': n['supporting_lanes'], 'mapping_caveat': n['mapping_caveat'],
                            'original_vcf_dp': n['original_vcf_dp'], 'representation_reconciliation': n['representation_reconciliation'],
                            'disposition_rationale': n['disposition_rationale'],
                            'missing_coverage_semantics': 'NONCALL_NOT_CONTRADICTION' if n['read_disposition'] == 'READ_NONCALL' else ''})
    record = {r['allele_key']: r['label_of_record'] for r in allele_rows}
    pair_rows = []
    for t in targets:
        k1, k2 = t['allele_1_normalized_allele_key'], t['allele_2_normalized_allele_key']
        pk = t['pair_key']; n, o = new_p[pk], old_p.get(pk)
        direct, src = (o['phase_label'], 'TASK06_EXACT_REUSE') if o else (n['phase_label'], 'TASK09_DIRECT_FRAGMENT')
        vp = vcf_phase(source[k1], source[k2])
        phase = direct
        if direct == 'PHASE_UNRESOLVED_NO_PARENTS' and vp != 'NO_SHARED_VCF_PHASE_BLOCK':
            phase, src = ('LOCAL_CIS_SUPPORTED' if vp.endswith('CIS') else 'LOCAL_TRANS_SUPPORTED'), 'VCF_PGT_PID_SHARED_BLOCK'
        elif direct in ('LOCAL_CIS_SUPPORTED', 'LOCAL_TRANS_SUPPORTED') and vp != 'NO_SHARED_VCF_PHASE_BLOCK' and vp.endswith('CIS') != (direct == 'LOCAL_CIS_SUPPORTED'):
            phase, src = 'PHASE_CONFLICT_UNRESOLVED', 'DIRECT_VS_VCF_CONFLICT'
        c1, p1 = k1.split(':')[:2]; c2, p2 = k2.split(':')[:2]
        pair_rows.append({'hypothesis_id': t['hypothesis_id'], 'pair_key': pk, 'gene_id': t['gene_id'],
                          'allele_1_label_of_record': record[k1], 'allele_2_label_of_record': record[k2],
                          'pair_read_label': rules.pair_read_label(record[k1], record[k2]),
                          'task09_scan_pair_read_disposition': n['pair_read_disposition'],
                          'phase_label_of_record': phase, 'phase_source': src,
                          'task09_direct_phase': n['phase_label'], 'task09_direct_phase_reason': n['phase_reason'],
                          'task09_haplotype_counts_json': n['haplotype_counts_json'], 'task06_phase': o['phase_label'] if o else '',
                          'phase_reuse_concordant': ('' if not o else str(o['phase_label'] == n['phase_label']).lower()),
                          'vcf_physical_phase': vp,
                          'distance_bp': abs(int(p2) - int(p1)) if c1 == c2 else '', 'parental_origin': 'UNKNOWN_NO_PARENTS'})
    write(sens / 'allele_read_ledger.tsv.gz', allele_rows, list(allele_rows[0]))
    write(sens / 'pair_read_phase_ledger.tsv.gz', pair_rows, list(pair_rows[0]))
    cnt = lambda it: dict(sorted(collections.Counter(it).items()))
    dist = lambda d: '<=1kb' if d != '' and d <= 1000 else '1-10kb' if d != '' and d <= 10000 else '10-100kb' if d != '' and d <= 100000 else '>100kb'
    summary = {'alleles': len(allele_rows), 'pairs': len(pair_rows),
               'allele_label_source': cnt(r['label_source'] for r in allele_rows),
               'allele_labels_of_record': cnt(r['label_of_record'] for r in allele_rows),
               'task09_scan_labels': cnt(r['task09_scan_label'] for r in allele_rows),
               'reuse_concordance': cnt(r['reuse_concordant'] for r in allele_rows if r['reuse_concordant']),
               'pair_read_labels': cnt(r['pair_read_label'] for r in pair_rows),
               'phase_labels_of_record': cnt(r['phase_label_of_record'] for r in pair_rows),
               'phase_sources': cnt(r['phase_source'] for r in pair_rows),
               'task09_direct_phase_reasons': cnt(r['task09_direct_phase_reason'] for r in pair_rows),
               'phase_reuse_concordance': cnt(r['phase_reuse_concordant'] for r in pair_rows if r['phase_reuse_concordant']),
               'vcf_physical_phase': cnt(r['vcf_physical_phase'] for r in pair_rows),
               'distance_bins': cnt(dist(r['distance_bp']) for r in pair_rows),
               'parental_origin_inferred': False, 'novel_alleles_assessed': 0}
    with (run / 'data/read_phase_summary.json').open('x') as handle:
        handle.write(json.dumps(summary, indent=1, sort_keys=True) + '\n')
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
