#!/usr/bin/env python3
"""Independent exact small-LP oracle and native solver acceptance contracts.

The oracle enumerates vertices using Fraction Gaussian elimination. It shares
neither a model reader nor a numerical algorithm with the C simplex engine.
All generated oracle problems have finite variable bounds, so a feasible
polytope has a vertex and a finite optimum.
"""
import fractions
import gzip
import itertools
import json
import math
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest

F = fractions.Fraction
BIN = str(pathlib.Path(sys.argv.pop(1)).resolve()) if __name__ == "__main__" else None


def decimal(value):
    """Print a terminating rational exactly, without a float conversion."""
    value = F(value)
    sign = "-" if value < 0 else ""
    num, den = abs(value.numerator), value.denominator
    places = 0
    while den != 1:
        common = math.gcd(den, 10)
        if common == 1:
            raise ValueError("nonterminating decimal in generated model")
        den //= common
        num *= 10 // common
        places += 1
    if not places:
        return sign + str(num)
    digits = str(num).rjust(places + 1, "0")
    return sign + digits[:-places] + "." + digits[-places:]


def linear_solve(matrix, rhs):
    """Solve a square rational system; a singular basis is not a vertex."""
    n = len(rhs)
    a = [[F(x) for x in row] + [F(b)] for row, b in zip(matrix, rhs)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        q = a[col][col]
        a[col] = [x / q for x in a[col]]
        for row in range(n):
            if row == col:
                continue
            q = a[row][col]
            if q:
                a[row] = [x - q * y for x, y in zip(a[row], a[col])]
    return [row[-1] for row in a]


def inequalities(model):
    """Return original row/bound inequalities as a*x <= b."""
    n = len(model["cost"])
    result = []
    for j, (lo, hi) in enumerate(model["bounds"]):
        unit = [F(int(k == j)) for k in range(n)]
        if lo is not None:
            result.append(([-x for x in unit], -F(lo)))
        if hi is not None:
            result.append((unit, F(hi)))
    for coeff, lo, hi in model["rows"]:
        coeff = list(map(F, coeff))
        if lo is not None:
            result.append(([-x for x in coeff], -F(lo)))
        if hi is not None:
            result.append((coeff, F(hi)))
    return result


def oracle(model):
    constraints = inequalities(model)
    n = len(model["cost"])
    best = None
    for chosen in itertools.combinations(constraints, n):
        point = linear_solve([a for a, _ in chosen], [b for _, b in chosen])
        if point is None or any(sum(a * x for a, x in zip(row, point)) > rhs
                                for row, rhs in constraints):
            continue
        value = F(model.get("offset", 0)) + sum(F(c) * x for c, x in zip(model["cost"], point))
        if best is None or (value > best if model.get("maximize") else value < best):
            best = value
    return best


def mps(model):
    n = len(model["cost"])
    text = ["NAME ORACLE", "OBJSENSE", " MAX" if model.get("maximize") else " MIN", "ROWS", " N OBJ"]
    rhs = []
    ranges = []
    for i, (_, lo, hi) in enumerate(model["rows"]):
        kind = "E" if lo is not None and lo == hi else "G" if hi is None else "L"
        text.append(f" {kind} r{i}")
        rhs.append((f"r{i}", lo if hi is None else hi))
        if lo is not None and hi is not None and lo != hi:
            ranges.append((f"r{i}", F(hi) - F(lo)))
    text.append("COLUMNS")
    for j in range(n):
        text.append(f" x{j} OBJ {decimal(model['cost'][j])}")
        for i, (coeff, _, _) in enumerate(model["rows"]):
            if coeff[j]:
                text.append(f" x{j} r{i} {decimal(coeff[j])}")
    text += ["RHS", " rhs OBJ " + decimal(-F(model.get("offset", 0)))]
    text.extend(f" rhs {name} {decimal(value)}" for name, value in rhs)
    if ranges:
        text.append("RANGES")
        text.extend(f" rng {name} {decimal(value)}" for name, value in ranges)
    text.append("BOUNDS")
    for j, (lo, hi) in enumerate(model["bounds"]):
        if lo is None and hi is None:
            text.append(f" FR b x{j}")
        elif lo is not None and lo == hi:
            text.append(f" FX b x{j} {decimal(lo)}")
        else:
            text.append(f" MI b x{j}" if lo is None else f" LO b x{j} {decimal(lo)}")
            if hi is not None:
                text.append(f" UP b x{j} {decimal(hi)}")
    return "\n".join(text + ["ENDATA", ""])


class Solver(unittest.TestCase):
    def direct(self, args):
        proc = subprocess.run([BIN] + list(map(str, args)), capture_output=True,
                              text=True, timeout=15)
        try:
            result = json.loads(proc.stdout)
        except ValueError:
            self.fail(f"non-JSON stdout: {proc.stdout!r}; stderr: {proc.stderr!r}")
        return proc.returncode, result

    def invoke(self, model, extra=(), *, compressed=False, relax=True):
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path = root / ("model.mps.gz" if compressed else "model.mps")
            text = model if isinstance(model, str) else mps(model)
            if compressed:
                with gzip.open(path, "wt") as f:
                    f.write(text)
            else:
                path.write_text(text)
            solution = root / "solution.sol"
            dual = root / "dual.sol"
            args = [BIN, str(path)] + (["--relax"] if relax else [])
            args += ["--solution", str(solution), "--dual", str(dual)] + list(extra)
            proc = subprocess.run(args, capture_output=True, text=True, timeout=15)
            self.assertIn(proc.returncode, (0, 2, 3, 4), proc.stderr)
            try:
                result = json.loads(proc.stdout)
            except ValueError:
                self.fail(f"non-JSON stdout: {proc.stdout!r}; stderr: {proc.stderr!r}")
            for key, value in result.items():
                if isinstance(value, float):
                    self.assertTrue(math.isfinite(value), (key, result))
            for key in ("iterations", "refactorizations"):
                if key in result:
                    self.assertGreaterEqual(result[key], 0, result)
            values = {}
            if solution.exists():
                for line in solution.read_text().splitlines():
                    if not line.strip() or line.lstrip().startswith("#"):
                        continue
                    name, value = line.split()
                    if name == "=obj=":
                        continue
                    self.assertNotIn(name, values)
                    values[name] = F(value)
            return proc.returncode, result, values, dual.read_text() if dual.exists() else None

    def compare(self, model):
        expected = oracle(model)
        code, result, values, _ = self.invoke(model)
        self.assertEqual(code, 0, (model, result))
        self.assertEqual(result["status"], "infeasible" if expected is None else "optimal", (model, result))
        if expected is None:
            return
        self.assertEqual(set(values), {f"x{i}" for i in range(len(model["cost"]))}, result)
        point = [values[f"x{i}"] for i in range(len(model["cost"]))]
        for a, b in inequalities(model):
            residual = sum(x * y for x, y in zip(a, point)) - b
            self.assertLessEqual(float(residual), 1e-7 * max(1.0, abs(float(b))), (model, point, a, b, result))
        value = F(model.get("offset", 0)) + sum(F(c) * x for c, x in zip(model["cost"], point))
        tolerance = 1e-7 * max(1.0, abs(float(expected)))
        self.assertAlmostEqual(float(value), float(expected), delta=tolerance, msg=str((model, result)))
        self.assertAlmostEqual(result["objective"], float(expected), delta=tolerance, msg=str((model, result)))
        return result

    def test_1000_generated_bounded_exact_oracle(self):
        rng = random.Random(20261106)
        for case in range(1000):
            n = 3 if case % 13 == 0 else rng.randint(1, 2)
            bounds = [(F(rng.randint(-4, -1)), F(rng.randint(1, 5))) for _ in range(n)]
            point = [F(rng.randint(int(lo * 2), int(hi * 2)), 2) for lo, hi in bounds]
            if case % 17 == 0:
                bounds[0] = (point[0], point[0])
            rows = []
            for _ in range(rng.randint(0, 3)):
                coeff = [F(rng.randint(-4, 4), rng.choice((1, 2))) for _ in range(n)]
                activity = sum(a * x for a, x in zip(coeff, point))
                mode = rng.choice(("L", "G", "E", "R"))
                slack = F(rng.randint(0, 5), 2)
                lo = activity - slack if mode in ("G", "R") else activity if mode == "E" else None
                hi = activity + slack if mode in ("L", "R") else activity if mode == "E" else None
                rows.append((coeff, lo, hi))
            if case % 11 == 0:
                rows.append(([F(int(j == 0)) for j in range(n)], bounds[0][1] + 1, None))
            model = {"cost": [F(rng.randint(-6, 6), rng.choice((1, 2))) for _ in range(n)],
                     "bounds": bounds, "rows": rows, "offset": F(rng.randint(-20, 20), 4),
                     "maximize": bool(case % 2)}
            with self.subTest(case=case):
                self.compare(model)

    def test_free_fixed_upper_only_and_redundant_rows(self):
        cases = [
            {"cost": [2, -1], "bounds": [(None, None), (3, 3)], "rows": [([1, 1], 5, 5)]},
            {"cost": [-1], "bounds": [(None, -2)], "rows": []},
            {"cost": [1], "bounds": [(None, None)], "rows": [([1], -3, None)]},
            {"cost": [1, 2], "bounds": [(0, 3), (0, 3)],
             "rows": [([1, 1], 2, 2), ([2, 2], 4, 4), ([0, 0], 0, 0), ([1, 1], 0, 3)]},
            {"cost": [1, -1], "bounds": [(2, 2), (-3, -3)], "rows": [([1, 1], -1, -1)], "offset": F(7, 10)},
        ]
        for model in cases:
            with self.subTest(model=model):
                self.compare(model)

    def test_infeasible_and_unbounded_statuses_are_claims(self):
        cases = [
            ({"cost": [1], "bounds": [(0, 1)], "rows": [([1], 2, None)]}, "infeasible"),
            ({"cost": [0], "bounds": [(0, 1)], "rows": [([0], 1, 1)]}, "infeasible"),
            ({"cost": [1], "bounds": [(None, None)], "rows": []}, "unbounded"),
            ({"cost": [-1], "bounds": [(0, None)], "rows": [([0], None, 1)]}, "unbounded"),
        ]
        for model, status in cases:
            with self.subTest(status=status, model=model):
                code, result, _, _ = self.invoke(model)
                self.assertEqual(code, 0, result)
                self.assertEqual(result["status"], status, result)
                self.assertFalse(result.get("certifies_optimality", False), result)

    def test_relaxation_does_not_enforce_integrality(self):
        model = """NAME RELAX
OBJSENSE
 MAX
ROWS
 N O
 L R
COLUMNS
 M0 'MARKER' 'INTORG'
 x O 1 R 2
 M1 'MARKER' 'INTEND'
RHS
 b R 1
BOUNDS
 BV b x
ENDATA
"""
        code, result, values, _ = self.invoke(model)
        self.assertEqual(code, 0, result)
        self.assertEqual(result["status"], "optimal", result)
        self.assertAlmostEqual(result["objective"], 0.5)
        self.assertEqual(values["x"], F(1, 2))
        code, result, _, _ = self.invoke(model, relax=False)
        self.assertIn(code, (2, 3), result)
        self.assertIn(result["status"], ("input_error", "unsupported"))

    def test_marker_defaults_and_explicit_bound_override(self):
        template = """NAME MARKERS
OBJSENSE
 MAX
ROWS
 N O
COLUMNS
 M0 'MARKER' 'INTORG'
 x O 1
 M1 'MARKER' 'INTEND'
BOUNDS
{bounds}ENDATA
"""
        for bounds, status, objective in [("", "optimal", 1),
                                           (" LO b x 2\n", "unbounded", None),
                                           (" UP b x 4\n", "optimal", 4)]:
            with self.subTest(bounds=bounds):
                code, result, _, _ = self.invoke(template.format(bounds=bounds))
                self.assertEqual(code, 0, result)
                self.assertEqual(result["status"], status, result)
                if objective is not None:
                    self.assertEqual(result["objective"], objective, result)

    def test_duplicate_coefficients_and_selected_objective(self):
        model = """NAME DUPLICATES
OBJSENSE
 MAX
OBJNAME
 COST
ROWS
 N DECOY
 N COST
 L CAP
COLUMNS
 x DECOY 77 COST 0.25
 x COST 0.75 CAP 1
 x CAP 1
RHS
 b CAP 1 COST 2
BOUNDS
 LO b x 0
 UP b x 3
ENDATA
"""
        code, result, values, _ = self.invoke(model)
        self.assertEqual(code, 0, result)
        self.assertEqual(result["status"], "optimal", result)
        self.assertEqual(values["x"], F(1, 2))
        self.assertAlmostEqual(result["objective"], -1.5)

    def test_gzip_d_exponent_and_max_offset(self):
        model = {"cost": [F(3, 10)], "bounds": [(-2, 4)], "rows": [],
                 "offset": F(-7, 2), "maximize": True}
        code, result, values, _ = self.invoke(mps(model).replace("0.3", "3D-1"), compressed=True)
        self.assertEqual(code, 0, result)
        self.assertEqual(result["status"], "optimal", result)
        self.assertAlmostEqual(result["objective"], -2.3)
        self.assertEqual(values["x0"], 4)

    def test_resource_limits_are_explicit(self):
        model = {"cost": [-1, -1], "bounds": [(0, 1), (0, 1)], "rows": []}
        for option, value, statuses in [
            ("--iteration-limit", "1", ("iteration_limit",)),
            ("--time-limit", "0.000000000001", ("time_limit",)),
            ("--max-dimension", "1", ("unsupported",)),
        ]:
            with self.subTest(option=option):
                code, result, _, _ = self.invoke(model, [option, value])
                self.assertIn(code, (3, 4), result)
                self.assertIn(result["status"], statuses, result)

    def test_parser_fails_closed(self):
        valid = mps({"cost": [1], "bounds": [(0, 2)], "rows": [([1], None, 1)]})
        malformed = ["", valid.replace("ENDATA", ""), valid.replace("x0 OBJ 1", "x0 MISSING 1"),
                     valid.replace("x0 OBJ 1", "x0 OBJ nan"), valid.replace("x0 OBJ 1", "x0 OBJ 1e9999")]
        for text in malformed:
            with self.subTest(model=text):
                code, result, _, _ = self.invoke(text)
                self.assertEqual(code, 2, result)
                self.assertEqual(result["status"], "input_error", result)
        code, result, _, _ = self.invoke(valid.replace("ENDATA", "QMATRIX\n x0 x0 1\nENDATA"))
        self.assertEqual(code, 3, result)
        self.assertEqual(result["status"], "unsupported", result)

    def test_equivalent_permutation_and_scaling(self):
        base = {"cost": [3, -2], "bounds": [(-2, 4), (-1, 3)],
                "rows": [([1, 2], -1, 4), ([-2, 1], None, 2)], "offset": F(13, 4)}
        original = self.compare(base)
        permuted = {**base, "cost": list(reversed(base["cost"])), "bounds": list(reversed(base["bounds"])),
                    "rows": [(list(reversed(a)), lo, hi) for a, lo, hi in reversed(base["rows"])]}
        self.assertAlmostEqual(self.compare(permuted)["objective"], original["objective"])
        scaled = {**base, "rows": [([-3 * x for x in a], None if hi is None else -3 * hi,
                                    None if lo is None else -3 * lo) for a, lo, hi in base["rows"]]}
        self.assertAlmostEqual(self.compare(scaled)["objective"], original["objective"])

    def test_fractional_vertices_and_scaled_rows(self):
        base = {"cost": [1, 2], "bounds": [(0, 3), (0, 3)],
                "rows": [([3, 0], 1, None), ([0, 7], 2, None)]}
        self.compare(base)
        for scale in (F(1, 10000), F(1000000)):
            scaled = {**base, "rows": [([scale * x for x in a],
                                        None if lo is None else scale * lo,
                                        None if hi is None else scale * hi)
                                       for a, lo, hi in base["rows"]]}
            with self.subTest(scale=scale):
                self.compare(scaled)

    def test_tiny_nonzero_equality_is_not_discarded(self):
        # A small coefficient is still a constraint: 1e-11*x = 0 pins x to
        # zero. Treating its artificial basis row as redundant fabricates an
        # improving ray for min -x. A safe explicit numerical failure is valid.
        model = {"cost": [-1], "bounds": [(0, None)],
                 "rows": [([F(1, 10**11)], 0, 0)]}
        code, result, values, _ = self.invoke(model)
        self.assertIn(result["status"], ("optimal", "numerical_failure"), result)
        if result["status"] == "optimal":
            self.assertEqual(code, 0, result)
            self.assertEqual(values["x0"], 0, result)
            self.assertEqual(result["objective"], 0, result)
        else:
            self.assertEqual(code, 4, result)

    def test_tiny_cost_with_large_bound_has_large_improvement(self):
        # An absolute reduced-cost cutoff loses a billion objective units here.
        model = {"cost": [F(-1, 10**11)], "bounds": [(0, 10**20)], "rows": []}
        code, result, values, _ = self.invoke(model)
        self.assertEqual(code, 0, result)
        self.assertEqual(result["status"], "optimal", result)
        self.assertEqual(values["x0"], 10**20, result)
        self.assertAlmostEqual(result["objective"], -10**9, delta=1e-6, msg=str(result))

    def test_small_explicit_cost_cannot_hide_an_improving_ray(self):
        # The x row creates a nonzero dual norm. It must not suppress a tiny
        # independent negative cost on the unbounded y direction.
        model = {"cost": [1, F('-1e-15')], "bounds": [(0, None), (0, None)],
                 "rows": [([1, 0], 1, None)]}
        for backend in ('sparse', 'dense'):
            code, result, _, _ = self.invoke(model, ['--backend', backend])
            self.assertEqual(code, 0, result)
            self.assertEqual(result['status'], 'unbounded', result)

    def test_constant_small_nonzero_equality_is_not_feasible(self):
        model = {"cost": [-1], "bounds": [(0, 1)],
                 "rows": [([0], F(1, 10**12), F(1, 10**12))]}
        code, result, _, _ = self.invoke(model)
        self.assertIn(result["status"], ("infeasible", "numerical_failure"), result)
        self.assertEqual(code, 0 if result["status"] == "infeasible" else 4, result)
        self.assertFalse(result.get("certifies_optimality", False), result)

    def test_cancelling_reduced_costs_do_not_cycle_or_choose_singular_basis(self):
        # Preserve generated cases 312 and 860 independently of generator order.
        # Removing every reduced-cost tolerance to handle tiny true costs once
        # promoted floating cancellation noise into entering-column decisions.
        cases = [
            {"cost": [-3, -3, 0], "bounds": [(-4, 2), (-3, 4), (-1, 3)],
             "rows": [([-3, 3, F(-3, 2)], F(9, 2), None)], "offset": F(3, 2)},
            {"cost": [3, F(-1, 2)], "bounds": [(-1, 2), (-4, 2)],
             "rows": [([0, -1], 3, None), ([0, F(3, 2)], F(-9, 2), None),
                      ([F(3, 2), F(3, 2)], None, F(-9, 4))], "offset": F(7, 2)},
        ]
        for model in cases:
            with self.subTest(model=model):
                self.compare(model)

    def test_degenerate_cycling_example(self):
        # Beale's example has a degenerate slack basis and cycles under the
        # classic largest-reduced-cost rule without an anti-cycling fallback.
        model = {"cost": [10, -57, -9, -24], "maximize": True,
                 "bounds": [(0, 10)] * 4,
                 "rows": [([F(1, 2), F(-11, 2), F(-5, 2), 9], None, 0),
                          ([F(1, 2), F(-3, 2), F(-1, 2), 1], None, 0),
                          ([1, 0, 0, 0], None, 1)]}
        self.compare(model)

    def test_basis_roundtrip_retains_redundant_row_policy(self):
        model = {"cost": [3, -2], "bounds": [(0, 3), (0, 3)],
                 "rows": [([1, 1], 2, 2), ([2, 2], 4, 4), ([0, 0], 0, 0)]}
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path = root / "model.mps"
            path.write_text(mps(model))
            basis = root / "first.bas"
            code, cold = self.direct([path, "--relax", "--basis", basis])
            self.assertEqual(code, 0, cold)
            self.assertEqual(cold["status"], "optimal", cold)
            self.assertGreater(cold["redundant_rows"], 0, cold)
            original = basis.read_bytes()
            code, warm = self.direct([path, "--relax", "--read-basis", basis,
                                      "--basis", root / "second.bas"])
            self.assertEqual(code, 0, warm)
            self.assertEqual(warm["status"], "optimal", warm)
            self.assertAlmostEqual(warm["objective"], float(oracle(model)))
            self.assertEqual(cold["objective"], warm["objective"])
            self.assertEqual(cold["redundant_rows"], warm["redundant_rows"])
            self.assertEqual(basis.read_bytes(), original)
            self.assertEqual((root / "second.bas").read_bytes(), original)

    def test_basis_rejects_wrong_model_and_tampered_ids(self):
        model = {"cost": [-1, -2], "bounds": [(0, 3), (0, 3)],
                 "rows": [([1, 1], None, 2)]}
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path = root / "model.mps"
            path.write_text(mps(model))
            basis = root / "basis.bas"
            code, result = self.direct([path, "--relax", "--basis", basis])
            self.assertEqual(code, 0, result)
            lines = basis.read_text().splitlines()
            header = lines[0].split()
            self.assertGreaterEqual(len(lines), 3)
            variants = {"wrong hash": " ".join([header[0], "0" * 64] + header[2:]) + "\n" + "\n".join(lines[1:]),
                        "wrong dimensions": " ".join(header[:2] + [str(int(header[2]) + 1)] + header[3:]) + "\n" + "\n".join(lines[1:]),
                        "truncated": "\n".join(lines[:-1]),
                        "trailing data": "\n".join(lines + ["unexpected"])}
            wrong_row = lines.copy()
            row, col = wrong_row[1].split()
            wrong_row[1] = f"{int(row) + 100} {col}"
            variants["wrong row ID"] = "\n".join(wrong_row)
            wrong_col = lines.copy()
            wrong_col[1] = f"{row} {int(header[3])}"
            variants["column out of range"] = "\n".join(wrong_col)
            duplicate_col = lines.copy()
            duplicate_col[2] = duplicate_col[2].split()[0] + " " + col
            variants["duplicate basis column"] = "\n".join(duplicate_col)
            for name, text in variants.items():
                with self.subTest(tampering=name):
                    bad = root / "tampered.bas"
                    bad.write_text(text + "\n")
                    code, result = self.direct([path, "--relax", "--read-basis", bad])
                    self.assertEqual(code, 2, result)
                    self.assertEqual(result["status"], "input_error", result)
            changed = root / "changed.mps"
            changed.write_text(mps({**model, "cost": [-2, -1]}))
            code, result = self.direct([changed, "--relax", "--read-basis", basis])
            self.assertEqual(code, 2, result)
            self.assertEqual(result["status"], "input_error", result)

    def test_output_files_are_created_fresh(self):
        model = {"cost": [-1], "bounds": [(0, 2)], "rows": []}
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path = root / "model.mps"
            path.write_text(mps(model))
            for option in ("--solution", "--dual", "--basis"):
                with self.subTest(option=option):
                    output = root / (option[2:] + ".existing")
                    output.write_text("preserve historical artifact\n")
                    code, result = self.direct([path, "--relax", option, output])
                    self.assertEqual(code, 2, result)
                    self.assertEqual(result["status"], "input_error", result)
                    self.assertEqual(output.read_text(), "preserve historical artifact\n")


if __name__ == "__main__":
    unittest.main()
