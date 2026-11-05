#!/usr/bin/env python3
"""Primal pricing contracts with exact fixtures and original-model checking.

--oracle-only validates the mathematical fixtures without executing a solver.
Full mode compares default and --bland-pricing policies using work counters.
"""
import argparse
import copy
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

ROOT = next(parent for parent in pathlib.Path(__file__).resolve().parents
            if (parent / 'CMakeLists.txt').is_file() and (parent / 'include/lp.h').is_file())
sys.path.insert(0, str(ROOT / 'tests'))
from test_lp_solver import inequalities, mps, oracle

SOLVER = CHECKER = None
TOLERANCE = F(1, 10**7)
COUNTERS = ('primal_accepted_pivots', 'primal_zero_steps', 'primal_short_steps',
            'pricing_bland_fallbacks')


def beale():
    return {'cost': [10, -57, -9, -24], 'maximize': True,
            'bounds': [(0, 10)] * 4,
            'rows': [([F(1, 2), F(-11, 2), F(-5, 2), 9], None, 0),
                     ([F(1, 2), F(-3, 2), F(-1, 2), 1], None, 0),
                     ([1, 0, 0, 0], None, 1)]}


def triangular(dimension):
    # Owned Klee-Minty-style triangular cube. The last row dominates every
    # positive objective coefficient, proving the upper bound 5**dimension.
    return {'cost': [2**(dimension-j-1) for j in range(dimension)], 'maximize': True,
            'bounds': [(0, 5**(j+1)) for j in range(dimension)],
            'rows': [([2**(i-j+1) if j < i else int(j == i)
                       for j in range(dimension)], None, 5**(i+1))
                     for i in range(dimension)]}


def weighted(columns, phase_one_zero=False):
    model = {'cost': [-j-1 for j in range(columns)], 'bounds': [(0, 1)] * columns,
             'rows': [([1] * columns, None, 1)]}
    if phase_one_zero:
        model['cost'].append(0)
        model['bounds'].append((0, 1))
        model['rows'][0][0].append(0)
        model['rows'].append(([0] * columns + [1], 0, 0))
    return model


def cancellation_models():
    # Retained algebraic cases, independent of generated test order or names.
    return [
        {'cost': [-3, -3, 0], 'bounds': [(-4, 2), (-3, 4), (-1, 3)],
         'rows': [([-3, 3, F(-3, 2)], F(9, 2), None)], 'offset': F(3, 2)},
        {'cost': [3, F(-1, 2)], 'bounds': [(-1, 2), (-4, 2)],
         'rows': [([0, -1], 3, None), ([0, F(3, 2)], F(-9, 2), None),
                  ([F(3, 2), F(3, 2)], None, F(-9, 4))], 'offset': F(7, 2)},
        {'cost': [1, F(-1, 10**11)], 'bounds': [(0, 1), (0, 10**20)], 'rows': []},
    ]


def random_models():
    rng = random.Random(202611091)
    for case in range(40):
        n = 3 if case % 10 == 0 else 2
        point = [F(rng.randint(-4, 6), 2) for _ in range(n)]
        rows = []
        for _ in range(case % 4):
            coefficients = [rng.randint(-4, 4) for _ in range(n)]
            activity = sum(a * x for a, x in zip(coefficients, point))
            slack = F(rng.randint(0, 4), 4)
            kind = rng.choice(('L', 'G', 'E', 'R'))
            rows.append((coefficients, activity if kind == 'E' else activity - slack
                         if kind in ('G', 'R') else None,
                         activity if kind == 'E' else activity + slack
                         if kind in ('L', 'R') else None))
        yield {'cost': [rng.randint(-5, 5) for _ in range(n)], 'bounds': [(-2, 3)] * n,
               'rows': rows, 'offset': F(rng.randint(-20, 20), 4), 'maximize': bool(case % 2)}


class FixtureOracles(unittest.TestCase):
    def feasible(self, model, point):
        for coefficients, rhs in inequalities(model):
            self.assertLessEqual(sum(a * x for a, x in zip(coefficients, point)), rhs)

    def test_beale_known_optimum(self):
        model = beale()
        self.assertEqual(oracle(model), F(1))
        self.feasible(model, list(map(F, [1, 0, 1, 0])))

    def test_triangular_known_optimum_and_bound(self):
        for n in range(2, 6):
            model = triangular(n)
            self.assertEqual(oracle(model), F(5**n))
            self.feasible(model, [F(0)] * (n-1) + [F(5**n)])
            last_row = model['rows'][-1][0]
            self.assertTrue(all(c <= a for c, a in zip(model['cost'], last_row)))
            self.assertTrue(all(c < a for c, a in zip(model['cost'][:-1], last_row[:-1])))

    def test_weighted_known_optimum(self):
        for n in (3, 4):
            self.assertEqual(oracle(weighted(n)), F(-n))
        for n in (8, 32):
            model = weighted(n)
            point = [F(0)] * (n-1) + [F(1)]
            self.feasible(model, point)
            self.assertTrue(all(c >= -n for c in model['cost']))
            self.assertEqual(sum(c * x for c, x in zip(model['cost'], point)), -n)

    def test_cancellation_and_generated_oracles_are_finite(self):
        for model in cancellation_models() + list(random_models()):
            self.assertIsNotNone(oracle(model))
        self.assertEqual(oracle(cancellation_models()[-1]), F(-10**9))


class PricingIntegration(unittest.TestCase):
    def solve(self, model, expected, bland=False, extra=()):
        with tempfile.TemporaryDirectory(prefix='pricing integration ') as td:
            root = pathlib.Path(td)
            path, solution = root / 'model.mps', root / 'solution.sol'
            path.write_text(mps(model))
            control = ['--bland-pricing'] if bland else []
            p = subprocess.run([SOLVER, str(path), '--relax', '--solution', str(solution),
                                '--iteration-limit', '1000', *control, *extra],
                               text=True, capture_output=True, timeout=20)
            self.assertEqual(p.returncode, 0, (model, p.stdout, p.stderr))
            result = json.loads(p.stdout)
            self.assertEqual(result['status'], 'optimal', result)
            self.assertTrue(result['has_primal'])
            self.assertFalse(result.get('certifies_optimality', False))
            for name in COUNTERS:
                self.assertIn(name, result, 'pricing integration must export this work counter')
                self.assertIsInstance(result[name], int)
                self.assertGreaterEqual(result[name], 0)
            self.assertLessEqual(result['primal_zero_steps'] + result['primal_short_steps'],
                                 result['primal_accepted_pivots'])
            if bland:
                self.assertEqual(result['pricing_bland_fallbacks'], 0)
            records = dict(line.split() for line in solution.read_text().splitlines())
            self.assertEqual(set(records), {'=obj='} | {f'x{i}' for i in range(len(model['cost']))})
            point = [F(records[f'x{i}']) for i in range(len(model['cost']))]
            for coefficients, rhs in inequalities(model):
                residual = sum(a * x for a, x in zip(coefficients, point)) - rhs
                self.assertLessEqual(residual, TOLERANCE, (model, point, result))
            value = F(model.get('offset', 0)) + sum(F(c) * x for c, x in zip(model['cost'], point))
            self.assertLessEqual(abs(value - expected), TOLERANCE, (model, point, value, expected, result))
            self.assertLessEqual(abs(F(str(result['objective'])) - expected), TOLERANCE)
            q = subprocess.run([CHECKER, '--model', str(path), '--solution', str(solution),
                                '--relax', '--lin-tol', '1e-7', '--int-tol', '0', '--obj-tol', '1e-7'],
                               text=True, capture_output=True, timeout=20)
            self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
            checked = json.loads(q.stdout)
            self.assertTrue(checked['checked_feasible'], checked)
            self.assertFalse(checked['certifies_optimality'], checked)
            self.assertLessEqual(abs(F(checked['objective']) - expected), TOLERANCE)
            return result

    def test_beale_fallback(self):
        for dense in (False, True):
            for unscaled in (False, True):
                extra = (['--backend', 'dense'] if dense else []) + (['--no-scaling'] if unscaled else [])
                for bland in (False, True):
                    with self.subTest(dense=dense, unscaled=unscaled, bland=bland):
                        result = self.solve(beale(), F(1), bland, extra)
                        self.assertGreaterEqual(result['primal_zero_steps'], 1)
                        if not bland:
                            self.assertGreaterEqual(result['pricing_bland_fallbacks'], 1)

    def test_triangular_family_correctness_without_universal_work_claim(self):
        for n in range(2, 6):
            for unscaled in (False, True):
                for bland in (False, True):
                    with self.subTest(n=n, unscaled=unscaled, bland=bland):
                        self.solve(triangular(n), F(5**n), bland, ['--no-scaling'] if unscaled else [])

    def test_weighted_work_counts(self):
        for n in (8, 32):
            for dense in (False, True):
                for deferred in (False, True):
                    extra = (['--backend', 'dense'] if dense else [])
                    if deferred:
                        extra += ['--bland-ratio', '--no-updates']
                    for bland in (False, True):
                        with self.subTest(n=n, dense=dense, deferred=deferred, bland=bland):
                            result = self.solve(weighted(n), F(-n), bland, extra)
                            self.assertEqual(result['primal_accepted_pivots'], n if bland else 1)
                            self.assertEqual(result['primal_zero_steps'], 0)
                            self.assertEqual(result['pricing_bland_fallbacks'], 0)

    def test_phase_reset_and_deferred_validation(self):
        n = 12
        for dense in (False, True):
            for deferred in (False, True):
                extra = ['--backend', 'dense'] if dense else []
                if deferred:
                    extra += ['--bland-ratio', '--no-updates']
                for bland in (False, True):
                    with self.subTest(dense=dense, deferred=deferred, bland=bland):
                        result = self.solve(weighted(n, phase_one_zero=True), F(-n), bland, extra)
                        self.assertEqual(result['primal_zero_steps'], 1)
                        self.assertEqual(result['primal_accepted_pivots'], n + 1 if bland else 2)
                        self.assertEqual(result['pricing_bland_fallbacks'], 0 if bland else 1)

    def test_unchanged_cancellation_and_true_tiny_cost_eligibility(self):
        for i, model in enumerate(cancellation_models()):
            expected = oracle(model)
            for unscaled in (False, True):
                for bland in (False, True):
                    with self.subTest(case=i, unscaled=unscaled, bland=bland):
                        self.solve(model, expected, bland, ['--no-scaling'] if unscaled else [])

    def test_generated_exact_oracle_comparison(self):
        for i, model in enumerate(random_models()):
            expected = oracle(model)
            for bland in (False, True):
                with self.subTest(case=i, bland=bland):
                    self.solve(copy.deepcopy(model), expected, bland)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oracle-only', action='store_true')
    parser.add_argument('--solver', type=pathlib.Path)
    parser.add_argument('--checker', type=pathlib.Path)
    args = parser.parse_args()
    global SOLVER, CHECKER
    classes = [FixtureOracles]
    if not args.oracle_only:
        if not args.solver or not args.checker:
            parser.error('full integration mode requires --solver and --checker')
        SOLVER, CHECKER = str(args.solver.resolve()), str(args.checker.resolve())
        classes.append(PricingIntegration)
    suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(cls) for cls in classes)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
