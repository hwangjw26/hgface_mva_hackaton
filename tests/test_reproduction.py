"""Synthetic boundary tests; no case inputs or remote access."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from canonical.condition_ids import condition_ids
from canonical.identity import pair_key
from canonical.observations import frequency_flag, physical_phase, read_consensus
from canonical.primary import primary_genotypes
from canonical import queues, tiers
from evidence import allele_observation, reconcile_phase
from io_utils import index
from reconcile_lirical import eligible, reconcile
from native_pairs import normalize
from reproduce import run


class Boundaries(unittest.TestCase):
    def workspace(self):
        root = Path(os.environ["TEST_WORKDIR"])
        root.mkdir(parents=True, exist_ok=True)
        temp = tempfile.TemporaryDirectory(dir=root)
        self.addCleanup(temp.cleanup)
        return Path(temp.name)

    def test_normalization_and_ref_failure(self):
        out = self.workspace()
        fasta = out / "ref.fa"
        fasta.write_text(">1\nCAAAAAG\n")
        Path(str(fasta) + ".fai").write_text("1\t7\t3\t7\t8\n")
        native = [{"CONTIG": "1", "START": "3", "REF": "AA", "ALT": "A"},
                  {"CONTIG": "1", "START": "1", "REF": "CA", "ALT": "C"}]
        result = normalize(native, out / "valid", os.environ["BCFTOOLS"], str(fasta))
        self.assertEqual(result, ["1:1:CA:C", "1:1:CA:C"])
        with self.assertRaises(subprocess.CalledProcessError):
            normalize([{"CONTIG": "1", "START": "2", "REF": "C", "ALT": "T"}],
                      out / "bad_ref", os.environ["BCFTOOLS"], str(fasta))

    def test_no_overwrite(self):
        out = self.workspace()
        sentinel = out / "sentinel"
        sentinel.write_text("unchanged")
        with self.assertRaises(ValueError):
            run({}, out)
        self.assertEqual(sentinel.read_text(), "unchanged")
        self.assertEqual(list(out.iterdir()), [sentinel])

    def test_unreturned_gene_retained(self):
        out = self.workspace()
        hpoa = out / "hpoa.tsv"
        hpoa.write_text("database_id\taspect\thpo_id\tqualifier\nOMIM:000001\tI\tHP:0000007\t\n")
        keys = ["1:11:A:C", "1:22:C:T"]
        m = {"threshold": .8, "gene_retention_filter_used": False,
             "contexts": [{"gene": "NCBIGene:1", "key": k, "mapped": False,
                           "score_eligible_allele_count": None, "plp_allele_count": None} for k in keys],
             "allele_observations": [{"key": k, "passed": True, "score": None, "clinvar": "UNKNOWN", "alt_count": 1} for k in keys]}
        meta = out / "metadata.json"
        meta.write_text(json.dumps(m))
        native = out / "native.json"
        native.write_text(json.dumps({"analysisResults": [{"diseaseId": "OMIM:000001", "posttestProbability": .2,
            "pretestProbability": .1, "compositeLR": 2, "genotypeLR": {"geneId": {"id": "NCBIGene:2", "symbol": "SYNTH"}, "lr": .5, "matchType": "UNKNOWN"}}]}))
        tsv = out / "native.tsv"
        tsv.write_text("diseaseCurie\trank\n")
        cfg = {"metadata": str(meta), "hpoa": str(hpoa), "cases": {case: {"json": str(native), "tsv": str(tsv)}
               for case in ["baseline"] + [f"loo_{i:02d}" for i in range(1, 8)]}}
        pairs = [{"gene_id": "NCBIGene:1", "pair_key": pair_key(*keys), "review_id": "SYNTH", "allele_1_key": keys[0], "allele_2_key": keys[1]}]
        alleles = [{"gene_id": "NCBIGene:1", "allele_key": k, "indexed_CLNDISDB_json": "[]"} for k in keys]
        result, counts = reconcile(cfg, pairs, alleles, out)
        self.assertEqual(len(result), 1)
        self.assertFalse(result[0]["t20_gene_returned"])
        self.assertFalse(result[0]["t20_ar_returned"])
        self.assertEqual(result[0]["t20_best_ar_disease"], "")
        self.assertEqual(set(counts.values()), {1})

    def test_malformed_af_is_flagged_before_reduction(self):
        for value in ["bad&0.0001", "nan&0.0001", "-0.1&0.0001"]:
            a = allele_observation("1:11:A:C", {"filter": "PASS", "format": {"GT": "0/1"},
                "annotations": [{"gnomADg_AF": value}]}, {"source_alt_count": "1"}, [])
            self.assertEqual(a["frequency_flag"], "UNKNOWN")
            self.assertIn("INVALID_AF", a["evidence_conflicts"])

    def test_pair_identity_order(self):
        self.assertEqual(pair_key("1:11:A:C", "1:22:C:T"), pair_key("1:22:C:T", "1:11:A:C"))

    def test_contributing_heterozygotes(self):
        rows = [{"preset": "exome", "ENTREZ_GENE_ID": "1", "MOI": "AR", "allele_key": k,
                 "CONTRIBUTING_VARIANT": "1", "GENOTYPE": "0/1"} for k in ["1:11:A:C", "1:22:C:T"]]
        self.assertEqual(len(primary_genotypes(rows)), 1)
        rows[1]["CONTRIBUTING_VARIANT"] = "0"
        self.assertEqual(primary_genotypes(rows), [])
        rows[1].update(CONTRIBUTING_VARIANT="1", GENOTYPE="1/1")
        self.assertEqual(primary_genotypes(rows), [])
        with self.assertRaises(ValueError):
            primary_genotypes([rows[0], rows[0]])

    def test_missing_is_not_benign(self):
        a = allele_observation("1:11:A:C", {"filter": "PASS", "format": {"GT": "0/1"}, "annotations": []},
                               {"source_alt_count": "1"}, [])
        self.assertEqual(a["frequency_flag"], "UNKNOWN")
        self.assertEqual(a["clinvar_class"], "MISSING")
        self.assertEqual(a["read_label"], "READ_NOT_ASSESSED")
        self.assertEqual(a["high_confidence_call"], "false")
        self.assertEqual(frequency_flag("nan", .01), ("UNKNOWN", "INVALID_AF"))

    def test_phase_requires_shared_valid_block(self):
        a = {"allele_key": "1:11:A:C", "source_alt_count": "1", "normalized_gt": "0/1", "pid": "BLOCK", "pgt": "0|1"}
        b = {**a, "allele_key": "1:22:C:T"}
        self.assertEqual(physical_phase(a, b), "LOCAL_CIS_SUPPORTED")
        self.assertEqual(physical_phase(a, {**b, "pgt": "1|0"}), "LOCAL_TRANS_SUPPORTED")
        self.assertEqual(physical_phase(a, {**b, "pid": "OTHER"}), "UNRESOLVED")
        self.assertEqual(physical_phase(a, {**b, "source_alt_count": "2"}), "UNRESOLVED")
        self.assertEqual(reconcile_phase(a, b, [{"phase": "LOCAL_TRANS_SUPPORTED", "record": {}}]), "PHASE_CONFLICT")

    def test_unevaluable_read_contradiction(self):
        self.assertEqual(read_consensus([{"label": "READ_CONTRADICTED", "record": {}}])[0], "READ_NONCALL")

    def test_condition_delimiters(self):
        self.assertEqual(condition_ids(["OMIM:000001&OMIM:000002|MONDO:0000003"]), ["MONDO:0000003", "OMIM:000001", "OMIM:000002"])

    def test_native_eligibility_missing_score(self):
        c = {"mapped": True, "score_eligible_allele_count": None, "plp_allele_count": 0}
        self.assertFalse(eligible(c, "any"))
        self.assertTrue(eligible({**c, "plp_allele_count": 1}, "any"))
        self.assertFalse(eligible({**c, "mapped": False, "plp_allele_count": 1}, "any"))

    def test_duplicate_join_rejected(self):
        with self.assertRaises(ValueError):
            index([{"key": "same"}, {"key": "same"}], ["key"])

    def test_cis_precludes_research_promotion(self):
        row = {"source_reconciled": "true", "gene_id_resolved": "true", "audit_phase": "LOCAL_CIS_SUPPORTED",
               "positive_ar_association": "true"}
        for n in [1, 2]:
            row.update({f"audit_{n}_" + k: v for k, v in {
                "normalized_gt": "0/1", "read_label": "STRONG_READ_SUPPORT", "high_confidence_call": "true",
                "frequency_flag": "WITHIN_BROAD_SCREEN", "gene_protein_splice": "true", "clinvar_class": "MISSING"}.items()})
        row.update(queues.assess(row))
        self.assertEqual(tiers.classify(row)["v2_tier"], "CONTRADICTED_PAIR")


if __name__ == "__main__":
    unittest.main()
