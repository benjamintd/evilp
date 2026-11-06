"""Engine-native tableaux to separately checked original-source GMI cuts."""
from fractions import Fraction as F
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

PRODUCER = Path(sys.argv.pop(1)).resolve()
CHECKER = Path(sys.argv.pop(1)).resolve()


def mps(a, lower, upper, lo, hi, objective, integer):
    rows = ["NAME ROOTCUT", "ROWS", " N OBJ"]
    for i, (l, u) in enumerate(zip(lower, upper)):
        rows.append(f" {'E' if l == u else 'L' if u is not None else 'G'} R{i}")
    rows.append("COLUMNS")
    inside = False
    marker = 0
    for j, is_integer in enumerate(integer):
        if bool(is_integer) != inside:
            rows.append(f" M{marker} 'MARKER' '{'INTORG' if is_integer else 'INTEND'}'")
            marker += 1
            inside = bool(is_integer)
        rows.append(f" X{j} OBJ {objective[j]}")
        for i, coefficients in enumerate(a):
            if F(str(coefficients[j])):
                rows.append(f" X{j} R{i} {coefficients[j]}")
    if inside:
        rows.append(f" M{marker} 'MARKER' 'INTEND'")
    rows.append("RHS")
    for i, (l, u) in enumerate(zip(lower, upper)):
        rows.append(f" RHS R{i} {u if u is not None else l}")
    ranges = [(i, F(str(u)) - F(str(l))) for i, (l, u) in enumerate(zip(lower, upper))
              if l is not None and u is not None and l != u]
    if ranges:
        rows.append("RANGES")
        for i, width in ranges:
            rows.append(f" RNG R{i} {float(width):.17g}")
    rows.append("BOUNDS")
    for j, (l, u) in enumerate(zip(lo, hi)):
        if l is not None and l == u:
            rows.append(f" FX B X{j} {l}")
            continue
        if l is None:
            rows.append(f" MI B X{j}")
        elif F(str(l)):
            rows.append(f" LO B X{j} {l}")
        if u is not None:
            rows.append(f" UP B X{j} {u}")
    return "\n".join(rows + ["ENDATA", ""])


class Pipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        retained = os.environ.get("MILPLAB_ROOT_CUT_ARTIFACTS")
        cls.root = None if retained else tempfile.TemporaryDirectory(prefix="root-cut-pipeline-")
        cls.directory = Path(retained) if retained else Path(cls.root.name)
        if retained:
            cls.directory.mkdir(parents=True, exist_ok=False)
        cls.proofs = cls.controls = cls.points = 0

    @classmethod
    def tearDownClass(cls):
        print(json.dumps({"checked_engine_proofs": cls.proofs,
                          "ownership_generation_controls": cls.controls,
                          "exact_feasible_points": cls.points}))
        if cls.root:
            cls.root.cleanup()

    def exercise(self, name, spec, *, count=1, overlay=None, profiles=True, no_exact=False):
        model = self.directory / (name + ".mps")
        model.write_text(mps(**spec))
        runs = itertools.product(range(2), range(2)) if profiles else [(0, 0)]
        outputs = []
        for dense, no_scaling in runs:
            prefix = self.directory / f"{name}-{dense}-{no_scaling}"
            argv = [str(PRODUCER), str(model), str(prefix), str(dense), str(no_scaling)]
            if overlay or no_exact:
                bounds = overlay or (spec["lo"][0], spec["hi"][0])
                argv += [str(bounds[0]), str(bounds[1]), str(int(no_exact))]
            process = subprocess.run(argv, text=True, capture_output=True, timeout=60)
            prefix.with_suffix(".producer.json").write_text(json.dumps({"argv": argv,
                "exit_code": process.returncode, "stdout": process.stdout, "stderr": process.stderr}, indent=2))
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            self.assertFalse(process.stderr)
            result = json.loads(process.stdout)
            self.assertEqual(result["cuts"], count, result)
            type(self).controls += result["controls"]
            for cut in result["published"]:
                checked = subprocess.run([str(CHECKER), "--model", str(model), "--gmi",
                                          cut["certificate"]], text=True, capture_output=True,
                                         timeout=30)
                self.assertEqual(checked.returncode, 0, checked.stdout + checked.stderr)
                self.assertFalse(checked.stderr)
                proof = json.loads(checked.stdout)
                Path(cut["certificate"] + ".checked.json").write_text(checked.stdout)
                self.assertTrue(proof["certifies_cut_validity"])
                self.assertTrue(proof["exact_tableau_identity"])
                for flag in ["certifies_original_lp_infeasibility", "certifies_original_mip_infeasibility",
                             "certifies_full_tree", "certifies_optimality",
                             "recorded_box_proven_implied_by_original_model"]:
                    self.assertFalse(proof[flag])
                text = Path(cut["certificate"]).read_text()
                fields = [line.split() for line in text.splitlines()]
                pub = {v[1]: F(v[2]) for v in fields if v[0] == "PUB"}
                self.assertEqual([pub[f"X{j}"] for j in range(len(spec["lo"]))],
                                 list(map(F, cut["coefficient"])))
                self.assertEqual([F(v[1]) for v in fields if v[0] == "PUB_RHS"],
                                 [F(cut["lower"])])
                self.check_points(spec, overlay, cut)
                outputs.append((model, cut, proof))
                type(self).proofs += 1
        return outputs

    def check_points(self, spec, overlay, cut):
        grids = []
        for j, integer in enumerate(spec["integer"]):
            values = map(F, range(-8, 19)) if integer else (F(k, 4) for k in range(-32, 73))
            lo, hi = spec["lo"][j], spec["hi"][j]
            if j == 0 and overlay:
                lo, hi = overlay
            grid = [x for x in values if (lo is None or x >= F(str(lo)))
                    and (hi is None or x <= F(str(hi)))]
            for row, rlo, rhi in zip(spec["a"], spec["lower"], spec["upper"]):
                if any(F(str(a)) for k, a in enumerate(row) if k != j):
                    continue
                grid = [x for x in grid if
                        (rlo is None or F(str(row[j]))*x >= F(str(rlo))) and
                        (rhi is None or F(str(row[j]))*x <= F(str(rhi)))]
            grids.append(grid)
        coefficients, rhs = list(map(F, cut["coefficient"])), F(cut["lower"])
        for point in itertools.product(*grids):
            feasible = True
            for row, lo, hi in zip(spec["a"], spec["lower"], spec["upper"]):
                activity = sum(F(str(a)) * x for a, x in zip(row, point))
                feasible &= ((lo is None or activity >= F(str(lo))) and
                             (hi is None or activity <= F(str(hi))))
            if feasible:
                self.assertGreaterEqual(sum(a*x for a, x in zip(coefficients, point)), rhs)
                type(self).points += 1

    def test_source_distances_and_scaling(self):
        cases = [
            ("lower", dict(a=[[2]], lower=[None], upper=[3], lo=[0], hi=[16], objective=[-1], integer=[1])),
            ("shifted", dict(a=[[2]], lower=[None], upper=[7], lo=[2], hi=[18], objective=[-1], integer=[1])),
            ("upper", dict(a=[[2]], lower=[-7], upper=[None], lo=[None], hi=[-2], objective=[1], integer=[1])),
            ("ranged", dict(a=[[2]], lower=[.5], upper=[3], lo=[0], hi=[16], objective=[-1], integer=[1])),
            ("nativebound", dict(a=[[2, 2]], lower=[None], upper=[3], lo=[0, 0], hi=[1, 1], objective=[-1, -1], integer=[1, 1])),
            ("mixed", dict(a=[[2, 1]], lower=[None], upper=[3], lo=[0, 0], hi=[16, 16], objective=[-1, 0], integer=[1, 0])),
            ("nonintegralanchor", dict(a=[[2, 2]], lower=[None], upper=[3], lo=[0, .25], hi=[16, 16.25], objective=[-1, 0], integer=[1, 1])),
            ("decimal", dict(a=[["1.000000000000000001"]], lower=[None], upper=[1.5], lo=[0], hi=[16], objective=[-1], integer=[1])),
            ("deletedrow", dict(a=[[0], [2]], lower=[0, None], upper=[0, 3], lo=[0], hi=[16], objective=[-1], integer=[1])),
        ]
        for name, spec in cases:
            with self.subTest(name=name):
                self.exercise(name, spec)

    def test_original_equalities_and_fixed_box(self):
        spec = dict(a=[[2, 3]], lower=[1.75], upper=[1.75], lo=[0, .25], hi=[16, .25], objective=[-1, 0], integer=[1, 0])
        for model, cut, proof in self.exercise("fixed", spec):
            self.assertTrue(proof["certifies_mip_in_box_infeasibility"])
            self.assertIn("FIX X1 ", Path(cut["certificate"]).read_text())

    def test_lp_feasible_integer_contradiction_is_scoped(self):
        spec = dict(a=[[2], [4]], lower=[1, 2], upper=[1, 2], lo=[0], hi=[16], objective=[-1], integer=[1])
        for model, cut, proof in self.exercise("parity", spec):
            self.assertTrue(proof["certifies_mip_in_box_infeasibility"])
            self.assertEqual(cut["coefficient"], ["0"])
            self.assertGreater(F(cut["lower"]), 0)

    def test_narrowed_session_box(self):
        spec = dict(a=[[2]], lower=[None], upper=[7], lo=[0], hi=[16], objective=[-1], integer=[1])
        self.exercise("overlay", spec, overlay=(2, 16))

    def test_multiple_cuts_and_deterministic_cap(self):
        spec = dict(a=[[2 if i == j else 0 for j in range(5)] for i in range(5)],
                    lower=[None]*5, upper=[3]*5, lo=[0]*5, hi=[16]*5,
                    objective=[-1]*5, integer=[1]*5)
        self.exercise("fourcuts", spec, count=4, profiles=False)

    def test_abstentions(self):
        base = dict(a=[[2]], lower=[None], upper=[4], lo=[0], hi=[16], objective=[-1], integer=[1])
        self.exercise("integral", base, count=0, profiles=False)
        self.exercise("continuous", dict(base, upper=[3], integer=[0]), count=0, profiles=False)
        self.exercise("noexact", dict(base, upper=[3]), count=0, profiles=False, no_exact=True)
        self.exercise("free", dict(base, lower=[1], upper=[1], lo=[None], hi=[None]), count=0, profiles=False)

    def test_actual_engine_proof_tampering(self):
        spec = dict(a=[[2, 1]], lower=[None], upper=[3], lo=[0, 0], hi=[16, 16], objective=[-1, 0], integer=[1, 0])
        model, cut, _ = self.exercise("tamper", spec, profiles=False)[0]
        original = Path(cut["certificate"]).read_text()
        for name, alter in [
            ("beta", lambda s: "\n".join("BETA 1000" if l.startswith("BETA ") else l for l in s.splitlines())),
            ("unsafe", lambda s: "\n".join("PUB_RHS 1000" if l.startswith("PUB_RHS ") else l for l in s.splitlines())),
            ("identity", lambda s: s.replace("TERM 0 ", "TERM 0 1000 # ", 1)),
            ("partial", lambda s: s.replace("END\n", "")),
        ]:
            path = self.directory / (name + ".cert")
            path.write_text(alter(original))
            result = subprocess.run([str(CHECKER), "--model", str(model), "--gmi", str(path)],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertFalse(json.loads(result.stdout)["certifies_cut_validity"])


if __name__ == "__main__":
    unittest.main()
