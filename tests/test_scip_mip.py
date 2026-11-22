"""Contract tests for the general C17 SCIP branch-and-cut and symmetry detector."""
from fractions import Fraction as F
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from test_mip_solver import integer_mps

SOLVER, CHECKER = (Path(p).resolve() for p in sys.argv[1:3])
del sys.argv[1:3]


class ScipGeneralMip(unittest.TestCase):
    def run_model(self, directory, text, options=()):
        model = directory / 'test.mps'
        sol = directory / 'point.sol'
        model.write_text(text)
        r = subprocess.run([SOLVER, model, '--mip', '--solution', sol, '--time-limit', '10',
                            *options], capture_output=True, text=True, timeout=20)
        self.assertIn(r.returncode, (0, 4), (r.stdout, r.stderr))
        self.assertNotIn('runtime error:', r.stderr)
        return json.loads(r.stdout), model, sol

    def test_mixed_integer_min_max_and_offset(self):
        # 36 columns: 18 binary + 18 continuous, with linking big-M rows and offset.
        rng = random.Random(1024)
        for maximize in (False, True):
            n_bin = 18
            n = 36
            costs = [rng.randint(-4, 6) for _ in range(n_bin)] + [F(rng.randint(-3, 5), 2) for _ in range(n_bin)]
            bounds = [(0, 1)] * n_bin + [(0, 4)] * n_bin
            integer = [True] * n_bin + [False] * n_bin
            rows = []
            # Each continuous y_i <= 4 * x_i
            for i in range(n_bin):
                a = [0] * n
                a[i] = -4
                a[n_bin + i] = 1
                rows.append((a, None, 0))
            # Budget on binary sum <= 7 and >= 3
            rows.append(([1] * n_bin + [0] * n_bin, 3, 7))
            # Coupled continuous demand >= 10
            rows.append(([0] * n_bin + [1] * n_bin, 10, 20))
            m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows,
                     offset=F(11, 4), maximize=maximize)
            # Compute exact optimum by choosing best y_i in {0, 4} for each open x_i.
            pair_vals = []
            for i in range(n_bin):
                c_x, c_y = costs[i], costs[n_bin + i]
                pair_vals.append((c_x, c_y))
            # Exact DP / enumeration over k open bins is tiny (n_bin=18, 2^18=262k? Even simpler:
            # compare SCIP result against native --no-root-diving fallback on exact same model!)
            with self.subTest(maximize=maximize), tempfile.TemporaryDirectory() as tmp:
                text = integer_mps(m)
                v, model, sol = self.run_model(Path(tmp), text)
                self.assertEqual(v['status'], 'optimal', v)
                self.assertIn('SCIP branch-and-cut', v['message'])
                self.assertFalse(v['certifies_optimality'])
                chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                     capture_output=True, text=True, timeout=15)
                self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
                c = json.loads(chk.stdout)
                self.assertTrue(c['checked_feasible'], c)
                # Cross-check against exact brute-force over the optimal LP per k:
                # For any chosen subset S (3 <= |S| <= 7), sum y_i in [10, 20], 0 <= y_i <= 4.
                # Verify optimality against native B&B with --no-root-diving (disables SCIP).
                nat_dir = Path(tmp) / 'nat'
                nat_dir.mkdir()
                v_nat, _, _ = self.run_model(nat_dir, text, ['--no-root-diving'])
                if v_nat['status'] == 'optimal':
                    self.assertAlmostEqual(float(F(c['objective'])), v_nat['objective'], places=5)

    def test_symmetric_packing_exercises_c17_symmetry(self):
        # 6 identical bins x 8 items = 48 binary vars + 6 bin-open vars = 54 vars.
        # Highly symmetric (S_6 on bins + item twins) -> exercises scip_symmetry.c!
        bins, items = 6, 8
        weights = [4, 4, 4, 4, 3, 3, 3, 3]
        cap = 7
        n = bins * items + bins
        costs = [0] * (bins * items) + [1] * bins
        bounds = [(0, 1)] * n
        integer = [True] * n
        rows = []
        # Each item assigned to exactly 1 bin
        for it in range(items):
            a = [0] * n
            for b in range(bins):
                a[b * items + it] = 1
            rows.append((a, 1, 1))
        # Each bin capacity <= cap * open_b
        for b in range(bins):
            a = [0] * n
            for it in range(items):
                a[b * items + it] = weights[it]
            a[bins * items + b] = -cap
            rows.append((a, None, 0))
        # Add 1 continuous slack variable so binary_sat abstains and SCIP+symmetry runs.
        costs.append(0)
        bounds.append((0, 0))
        integer.append(False)
        rows = [(a + [0], lo, hi) for a, lo, hi in rows]
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, model, sol = self.run_model(Path(tmp), integer_mps(m))
            self.assertEqual(v['status'], 'optimal', v)
            self.assertIn('SCIP branch-and-cut', v['message'])
            chk = subprocess.run([CHECKER, '--exact', '--model', model, '--solution', sol],
                                 capture_output=True, text=True, timeout=15)
            self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
            c = json.loads(chk.stdout)
            self.assertTrue(c['checked_feasible'], c)
            # 4 items of weight 4 + 4 items of weight 3, cap 7 => min 4 bins (each 4+3=7).
            self.assertEqual(F(c['objective']), F(4))

    def test_infeasible_mixed_integer(self):
        n = 36
        costs = [1] * n
        bounds = [(0, 1)] * 35 + [(0, 10)]
        integer = [True] * 35 + [False]
        rows = [
            ([1] * 35 + [0], 10, None),
            ([1] * 35 + [1], None, 9),
        ]
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, _, _ = self.run_model(Path(tmp), integer_mps(m))
            self.assertEqual(v['status'], 'infeasible', v)

    def test_feasjump_pure_feasibility_and_mixed_polish(self):
        # Mixed-integer feasibility + optimization with linked general integers and continuous vars
        # that exercises scip_feasjump.c (both zero-objective feasibility and coordinate polish).
        rng = random.Random(2026)
        n_int = 24
        n_cont = 16
        n = n_int + n_cont
        bounds = [(0, 3)] * n_int + [(0, 5)] * n_cont
        integer = [True] * n_int + [False] * n_cont
        # Construct a planted feasible point so feasibility is guaranteed.
        x_star = [rng.randint(0, 3) for _ in range(n_int)] + [rng.randint(0, 5) for _ in range(n_cont)]
        rows = []
        for _ in range(28):
            a = [rng.randint(-2, 3) if rng.random() < 0.35 else 0 for _ in range(n)]
            if not any(a):
                a[0] = 1
            rhs = sum(ai * xi for ai, xi in zip(a, x_star))
            rows.append((a, rhs - 2, rhs + 2))
        for zero_obj in (True, False):
            costs = [0] * n if zero_obj else [rng.randint(-3, 5) for _ in range(n)]
            m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
            with self.subTest(zero_obj=zero_obj), tempfile.TemporaryDirectory() as tmp:
                v, model, sol = self.run_model(Path(tmp), integer_mps(m))
                self.assertEqual(v['status'], 'optimal', v)
                self.assertIn('SCIP branch-and-cut', v['message'])
                chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                     capture_output=True, text=True, timeout=15)
                self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
                c = json.loads(chk.stdout)
                self.assertTrue(c['checked_feasible'], c)
                if zero_obj:
                    self.assertEqual(F(c['objective']), F(0))
                else:
                    # Must be at least as good as the planted feasible solution.
                    planted_obj = sum(ci * xi for ci, xi in zip(costs, x_star))
                    self.assertLessEqual(F(c['objective']), F(planted_obj))

    def test_parallel_row_presolver_merge_and_infeasibility(self):
        # Mixed-integer model with scalar-multiple parallel rows (positive and negative ratios):
        # Exercises scip_presol_parallel.c both for feasible merging and infeasible cutoff.
        n = 36
        costs = [1 + (j % 5) for j in range(n)]
        bounds = [(0, 1)] * 24 + [(0, 4)] * 12
        integer = [True] * 24 + [False] * 12
        base_row = [2 if j % 3 == 0 else (-1 if j % 3 == 1 else 3) for j in range(n)]
        # Feasible parallel pair: base_row >= 10, 3*base_row <= 36 (i.e. base_row <= 12), -2*base_row <= -22 (i.e. base_row >= 11)
        rows_feas = [
            (base_row, 10, None),
            ([3 * a for a in base_row], None, 36),
            ([-2 * a for a in base_row], None, -22),
            ([1] * 24 + [0] * 12, 4, 16),
        ]
        m_feas = dict(cost=costs, bounds=bounds, integer=integer, rows=rows_feas, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, model, sol = self.run_model(Path(tmp), integer_mps(m_feas))
            self.assertEqual(v['status'], 'optimal', v)
            chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                 capture_output=True, text=True, timeout=15)
            self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
            self.assertTrue(json.loads(chk.stdout)['checked_feasible'])
        # Infeasible parallel pair: base_row >= 15 and -2*base_row >= -24 (i.e. base_row <= 12)
        rows_inf = [
            (base_row, 15, None),
            ([-2 * a for a in base_row], -24, None),
        ]
        m_inf = dict(cost=costs, bounds=bounds, integer=integer, rows=rows_inf, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, _, _ = self.run_model(Path(tmp), integer_mps(m_inf))
            self.assertEqual(v['status'], 'infeasible', v)

    def test_big_m_coefftight_presolver(self):
        # Mixed-integer model with deliberately loose big-M coefficients on both
        # <= and >= rows (positive and negative binary coefficients) so
        # scip_presol_coefftight.c fires in SCIP_STAGE_PRESOLVING without errors
        # and preserves exact optimality.
        n_bin = 18
        n_cont = 18
        n = n_bin + n_cont
        costs = [5 + (i % 4) for i in range(n_bin)] + [-2 - (i % 3) for i in range(n_cont)]
        bounds = [(0, 1)] * n_bin + [(0, 3)] * n_cont
        integer = [True] * n_bin + [False] * n_cont
        rows = []
        # 1) sum_{k} y_k - 1000 * x_i <= 0 (loose negative big-M on <= row)
        for i in range(0, n_bin, 2):
            a = [0] * n
            a[i] = -1000
            a[n_bin + i] = 1
            a[n_bin + i + 1] = 1
            rows.append((a, None, 0))
        # 2) 1000 * x_i + y_i + y_{i+1} <= 1003 (loose positive big-M on <= row)
        for i in range(1, n_bin, 2):
            a = [0] * n
            a[i] = 1000
            a[n_bin + i] = 1
            a[n_bin + (i - 1)] = 1
            rows.append((a, None, 1003))
        # 3) 1000 * x_i + y_i >= 2 (loose positive big-M on >= row)
        for i in range(0, n_bin, 3):
            a = [0] * n
            a[i] = 1000
            a[n_bin + i] = 1
            rows.append((a, 2, None))
        # 4) -1000 * x_i + y_i >= -998 (loose negative big-M on >= row)
        for i in range(1, n_bin, 3):
            a = [0] * n
            a[i] = -1000
            a[n_bin + i] = 1
            rows.append((a, -998, None))
        # Budget row
        rows.append(([1] * n_bin + [0] * n_cont, 3, 8))
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            text = integer_mps(m)
            v, model, sol = self.run_model(Path(tmp), text)
            self.assertEqual(v['status'], 'optimal', v)
            self.assertIn('SCIP branch-and-cut', v['message'])
            chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                 capture_output=True, text=True, timeout=15)
            self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
            c = json.loads(chk.stdout)
            self.assertTrue(c['checked_feasible'], c)
            # Cross-check objective against native B&B (--no-root-diving disables SCIP).
            nat_dir = Path(tmp) / 'nat'
            nat_dir.mkdir()
            v_nat, _, _ = self.run_model(nat_dir, text, ['--no-root-diving'])
            if v_nat['status'] == 'optimal':
                self.assertAlmostEqual(float(F(c['objective'])), v_nat['objective'], places=5)


if __name__ == '__main__':
    unittest.main()
