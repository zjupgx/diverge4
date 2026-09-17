import unittest
from pathlib import Path

from diverge import Gu2001
from diverge.binding import Rvs, TypeOneAnalysis, get_colnames
from diverge.utils import CalPipe


ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = ROOT / "Tutorial" / "test_data"


class _Progress:
    def update(self, _amount):
        pass

    def set_description(self, _description):
        pass

    def refresh(self):
        pass


class _MethodAnalyzer:
    def __init__(self, *_args):
        pass

    def result(self):
        return "computed"


class RegressionTests(unittest.TestCase):
    def test_get_colnames_returns_a_list_copy(self):
        names = ("a", "b")
        self.assertEqual(get_colnames(names), ["a", "b"])

    def test_results_overrides_remain_properties(self):
        self.assertIsInstance(Rvs.results, property)
        self.assertIsInstance(TypeOneAnalysis.results, property)

    def test_calpipe_resolves_method_results(self):
        pipeline = object.__new__(CalPipe)
        pipeline.aln_file = "alignment"
        pipeline.tree_files = ()
        pipeline._results = {}
        config = {
            "class": _MethodAnalyzer,
            "result_attrs": [("result", "resolved_result")],
            "description": "dummy",
        }

        pipeline._execute_analysis("Dummy", config, _Progress())

        self.assertEqual(pipeline.resolved_result, "computed")
        self.assertEqual(pipeline._results["resolved_result"], "computed")

    def test_gu2001_failure_is_raised_in_python(self):
        with self.assertRaisesRegex(RuntimeError, "gu2001_compute"):
            Gu2001(
                str(TEST_DATA / "CASP.aln"),
                str(TEST_DATA / "cl1.tree"),
                str(TEST_DATA / "cl2.tree"),
                cluster_name=["A", "B"],
            )


if __name__ == "__main__":
    unittest.main()
