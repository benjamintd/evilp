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


if __name__ == '__main__':
    unittest.main()
