import copy
import json
import unittest
from pathlib import Path

from annotate_api import api_observation, pair_state


class ApiAnnotationTests(unittest.TestCase):
    def setUp(self):
        root = Path(__file__).resolve().parent
        self.record = json.loads((root / "fixtures/public_control_response.jsonl").read_text())
        self.features = self.record["response"]["AVI_SCORE_FEATURE_IMPORTANCE"]["track_names"]

    def test_public_control(self):
        result = api_observation(self.record, self.features)
        self.assertEqual(result["avi_status"], "SCORED_API")
        self.assertEqual(result["avi_meets_global"], "false")
        self.assertEqual(result["avi_top_feature"], "PHASTCONS_470_WAY")
        self.assertEqual(len(json.loads(result["avi_feature_importance_json"])), 18)

    def test_exact_cutoff(self):
        self.record["response"]["AVI_SCORE"]["layers"]["quantiles"] = [[0.99379]]
        self.assertEqual(api_observation(self.record, self.features)["avi_meets_global"], "true")
        self.record["response"]["AVI_SCORE"]["layers"]["quantiles"] = [[0.99378]]
        self.assertEqual(api_observation(self.record, self.features)["avi_meets_global"], "false")

    def test_wrong_returned_variant_rejected(self):
        self.record["response"]["AVI_SCORE"]["returned_allele_key"] = "1:1:A:C"
        with self.assertRaises(ValueError):
            api_observation(self.record, self.features)

    def test_missing_nan_or_invalid_quantile_rejected(self):
        for value in [None, [[float("nan")]], [[1.1]]]:
            row = copy.deepcopy(self.record)
            row["response"]["AVI_SCORE"]["layers"]["quantiles"] = value
            with self.assertRaises(ValueError):
                api_observation(row, self.features)

    def test_changed_features_rejected(self):
        with self.assertRaises(ValueError):
            api_observation(self.record, list(reversed(self.features)))

    def test_pair_error_not_false(self):
        low = api_observation(self.record, self.features)
        high = dict(low, avi_meets_global="true")
        error = api_observation({"status": "API_ERROR"}, self.features)
        self.assertEqual(error["avi_quantile"], "")
        self.assertEqual(pair_state(high, error), "INCOMPLETE")
        self.assertEqual(pair_state(high, low), "ONE_MEETS")
        self.assertEqual(pair_state(high, high), "BOTH_MEET")
        self.assertEqual(pair_state(low, low), "NEITHER_MEETS")


if __name__ == "__main__":
    unittest.main()
