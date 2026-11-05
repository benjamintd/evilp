#!/usr/bin/env python3
"""Cross the optimizer/checker boundary; neither output is its own oracle."""
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F
from test_lp_solver import mps, oracle

SOLVER = str(pathlib.Path(sys.argv.pop(1)).resolve())
CHECKER = str(pathlib.Path(sys.argv.pop(1)).resolve())


class Pipeline(unittest.TestCase):
    def solve(self, root, model, extra=()):
        path = root / 'model.mps'
        path.write_text(mps(model))
        solution, dual = root / 'solution.sol', root / 'dual.sol'
        p = subprocess.run([SOLVER, str(path), '--relax', '--solution', str(solution),
                            '--dual', str(dual)], text=True, capture_output=True, timeout=10)
        self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
        result = json.loads(p.stdout)
        self.assertEqual(result['status'], 'optimal', result)
        self.assertFalse(result['certifies_optimality'])
        q = subprocess.run([CHECKER, '--model', str(path), '--solution', str(solution),
                            '--dual', str(dual), '--relax', *extra], text=True,
                           capture_output=True, timeout=10)
        self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
        return json.loads(q.stdout)

    def test_exported_exact_and_fractional_witnesses(self):
        for maximum in (False, True):
            with tempfile.TemporaryDirectory() as td:
                model = {'cost': [3, 2, F('0.1')], 'offset': 5, 'maximize': maximum,
                         'bounds': [(0, 1), (0, 1), (None, None)],
                         'rows': [([2, 1, 0], None, 2), ([0, 0, 1], F('.3'), F('.3'))]}
                c = self.solve(pathlib.Path(td), model, ['--exact'])
                self.assertTrue(c['certifies_optimality'], c)
                self.assertEqual(F(c['objective']), oracle(model))
                self.assertEqual(F(c['certified_bound']), oracle(model))
        with tempfile.TemporaryDirectory() as td:
            model = {'cost': [1], 'bounds': [(0, 1)], 'rows': [([3], 1, None)]}
            c = self.solve(pathlib.Path(td), model)
            self.assertTrue(c['checked_feasible'], c)
            self.assertFalse(c['certifies_optimality'], c)
            self.assertEqual(c['bound_status'], 'valid', c)
            self.assertLessEqual(F(c['certified_bound']), F(1, 3))

    def test_60_generated_original_model_dual_bounds(self):
        rng = random.Random(841173)
        for case in range(60):
            maximum = bool(case % 2)
            model = {'cost': [rng.randint(-5, 5), rng.randint(-5, 5)],
                     'offset': F(rng.randint(-9, 9), 4), 'maximize': maximum,
                     'bounds': [(-2, 4), (-3, 5)],
                     'rows': [([rng.randint(1, 5), rng.randint(1, 5)],
                               F(rng.randint(-5, 0)), F(rng.randint(1, 7)))]}
            expected = oracle(model)
            self.assertIsNotNone(expected)
            with self.subTest(case=case), tempfile.TemporaryDirectory() as td:
                c = self.solve(pathlib.Path(td), model)
                self.assertTrue(c['checked_feasible'], c)
                self.assertEqual(c['bound_status'], 'valid', c)
                bound = F(c['certified_bound'])
                self.assertLessEqual(expected, bound) if maximum else self.assertLessEqual(bound, expected)
                self.assertLess(abs(float(bound - expected)), 1e-7)
                self.assertFalse(c['certifies_optimality'])

    def test_mixed_cost_scales_and_unbounded_direction(self):
        for upper in ('1e20', None):
            with tempfile.TemporaryDirectory() as td:
                root = pathlib.Path(td)
                model = {'cost': [1, F('-1e-11')], 'bounds': [(0, 1), (0, F(upper) if upper else None)],
                         'rows': []}
                path = root / 'model.mps'; path.write_text(mps(model))
                p = subprocess.run([SOLVER, str(path), '--relax'], text=True,
                                   capture_output=True, timeout=10)
                r = json.loads(p.stdout)
                self.assertEqual(p.returncode, 0, r)
                self.assertEqual(r['status'], 'optimal' if upper else 'unbounded', r)
                if upper:
                    self.assertAlmostEqual(r['objective'], -1e9, delta=1e-6)


if __name__ == '__main__':
    unittest.main()
