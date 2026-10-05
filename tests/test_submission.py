"""Patient-free regression tests for the frozen submission contract."""
import copy
import unittest
from build_submission import allele, render, select


def fixture():
    rows = []
    for i, (tier, exome, genome) in enumerate([
        ('TIER_2_MECHANISM_REVIEW', '9', '2'),
        ('TIER_2_MECHANISM_REVIEW', '1', '1'),
        ('TIER_3_EXPLORATORY', '2', '3'),
        ('HOLD_FOR_REVIEW', '3', '4'),
        ('TIER_2_MECHANISM_REVIEW', '5', ''),
        ('TIER_2_MECHANISM_REVIEW', '6', '9'),
    ]):
        rows.append(dict(review_id=f'FAKE{i}', gene_id=f'GENE{i}', pair_key=f'PAIR{i}',
                         allele_1_key=f'1:{101 + i * 100}:A:G', allele_2_key=f'1:{121 + i * 100}:C:T',
                         existing_research_tier=tier, exome_native_ar_rank=exome,
                         genome_native_ar_rank=genome, phase='UNRESOLVED'))
    return rows


class SubmissionTests(unittest.TestCase):
    def test_frozen_order_and_ordinal_values(self):
        selected = select(list(reversed(fixture())))
        self.assertEqual([r['review_id'] for r in selected], ['FAKE1', 'FAKE0', 'FAKE4', 'FAKE5'])
        self.assertEqual([r['epcr'] for r in render(selected, 'SYNTHETIC')], ['1.000000', '0.750000', '0.500000', '0.250000'])

    def test_duplicate_swapped_pair_rejected(self):
        rows = fixture()
        row = copy.deepcopy(rows[0]); row['review_id'] = 'FAKE_DUP'
        row['allele_1_key'], row['allele_2_key'] = row['allele_2_key'], row['allele_1_key']
        with self.assertRaises(ValueError):
            select(rows + [row])

    def test_cis_cannot_be_rendered_as_unresolved(self):
        rows = fixture(); rows[0]['phase'] = 'LOCAL_CIS_SUPPORTED'
        with self.assertRaises(ValueError):
            select(rows)

    def test_invalid_and_missing_rank_rejected(self):
        for value in ['0', '-1', 'nan', '1.5', '']:
            rows = fixture(); rows[0]['exome_native_ar_rank'] = value; rows[0]['genome_native_ar_rank'] = ''
            with self.assertRaises(ValueError):
                select(rows)

    def test_allele_validation(self):
        for value in ['1:0:A:G', '1:1:A:A', '1:1:A:N', '23:10:A:T']:
            with self.assertRaises(ValueError):
                allele(value)
        self.assertEqual(allele('MT:2:A:G'), ('chrM', 2, 'A', 'G'))


if __name__ == '__main__':
    unittest.main()
