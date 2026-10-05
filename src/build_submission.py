"""Reproduce a pair submission from a privately reviewed evidence ledger.

This program selects reviewed research tiers; it does not classify pathogenicity.
Only Python 3.11 standard-library modules are required. See README.md for scope.
"""
import argparse
import csv
import math
import os
from pathlib import Path
import re

COLUMNS = ['proband_id', 'chrom_1', 'pos_1', 'ref_1', 'alt_1',
           'chrom_2', 'pos_2', 'ref_2', 'alt_2', 'epcr', 'finding_type', 'notes']
TIERS = ['TIER_1_REVIEWED_ALLELE_SUPPORT', 'TIER_2_MECHANISM_REVIEW']
FIELDS = ['review_id', 'gene_id', 'pair_key', 'allele_1_key', 'allele_2_key',
          'existing_research_tier', 'exome_native_ar_rank', 'genome_native_ar_rank', 'phase']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def allele(key):
    match = re.fullmatch(r'(?:chr)?(1[0-9]|2[0-2]|[1-9]|X|Y|MT|M):([1-9][0-9]*):([ACGT]+):([ACGT]+)', key)
    require(match is not None, 'Invalid normalized allele key')
    chrom, pos, ref, alt = match.groups()
    require(ref != alt, 'REF must differ from ALT')
    return ('chr' + ('M' if chrom == 'MT' else chrom), int(pos), ref, alt)


def allele_order(v):
    contig = v[0][3:]
    return ((int(contig) if contig.isdigit() else {'X': 23, 'Y': 24, 'M': 25}[contig]), *v[1:])


def native_rank(value):
    if value == '':
        return math.inf
    require(value.isascii() and value.isdigit() and int(value) > 0, 'Invalid native rank')
    return int(value)


def select(rows):
    require(bool(rows), 'Empty reviewed ledger')
    seen_ids, seen_pairs = set(), set()
    eligible = []
    for row in rows:
        require(all(key in row for key in FIELDS), 'Missing reviewed ledger fields')
        require(row['review_id'] not in seen_ids, 'Duplicate review identifier')
        pair = frozenset(allele(row[k]) for k in ['allele_1_key', 'allele_2_key'])
        require(len(pair) == 2 and pair not in seen_pairs, 'Duplicate or degenerate physical pair')
        require(len({v[0] for v in pair}) == 1, 'A same-gene pair must share a chromosome')
        seen_ids.add(row['review_id'])
        seen_pairs.add(pair)
        e, g = native_rank(row['exome_native_ar_rank']), native_rank(row['genome_native_ar_rank'])
        require(math.isfinite(min(e, g)), 'At least one native rank required')
        if row['existing_research_tier'] in TIERS:
            require(row['phase'] == 'UNRESOLVED', 'This frozen submission contract requires unresolved phase')
            eligible.append(row)
    def order(row):
        e, g = native_rank(row['exome_native_ar_rank']), native_rank(row['genome_native_ar_rank'])
        return (TIERS.index(row['existing_research_tier']), min(e, g), e, g, row['pair_key'], row['gene_id'])
    return sorted(eligible, key=order)[:10]


def render(selected, proband_id):
    require(0 < len(selected) <= 10, 'Select between one and ten pairs')
    require(bool(proband_id.strip()), 'Missing proband identifier')
    result = []
    for rank, source in enumerate(selected, 1):
        row = {'proband_id': proband_id}
        variants = sorted([allele(source['allele_1_key']), allele(source['allele_2_key'])], key=allele_order)
        for n, variant in enumerate(variants, 1):
            row.update(dict(zip([f'chrom_{n}', f'pos_{n}', f'ref_{n}', f'alt_{n}'], variant)))
        row.update(epcr=f'{(len(selected) - rank + 1) / len(selected):.6f}', finding_type='primary',
                   notes='Phase unresolved.')
        result.append(row)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ledger', type=Path, required=True)
    parser.add_argument('--proband-id', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-rows', type=int, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    csv.field_size_limit(10000000)
    with args.ledger.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    selected = select(rows)
    require(len(selected) == args.expected_rows, 'Selected count differs from frozen submission decision')
    with args.output.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator='\n')
        writer.writeheader()
        writer.writerows(render(selected, args.proband_id))
    print(f'Wrote {len(selected)} candidate pairs; EPCR is ordinal, not calibrated.')


if __name__ == '__main__':
    main()
