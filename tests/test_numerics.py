#!/usr/bin/env python3
"""Exact-oracle numerical regressions and equivalent-model scaling contracts.

These tests change mathematical units, not feasibility tolerances. Every finite
answer is evaluated with Fraction against both the submitted model and the
unscaled equivalent, and large objective constants are removed from the error
scale. The optimizer's numerical statuses are never treated as proof flags.
"""
import copy
import json
import math
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

from test_lp_solver import inequalities, mps, oracle
from test_mip_solver import generated_models as mip_models, integer_mps, mip_oracle

SOLVER = str(pathlib.Path(sys.argv.pop(1)).resolve()) if __name__ == '__main__' else None
CHECKER = (str(pathlib.Path(sys.argv.pop(1)).resolve())
           if __name__ == '__main__' and len(sys.argv) > 1 and not sys.argv[1].startswith('-')
           else None)
EXPONENTS = (-26, -20, -10, 0, 10, 20, 26)
ABS_TOL = F(1, 10**7)


def power2(exponent):
    return F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))


def scaled_rows(model, exponents):
    """Positive exact powers of two preserve every row's feasible set."""
    result = copy.deepcopy(model)
    result['rows'] = []
    for (coeff, lo, hi), exponent in zip(model['rows'], exponents):
        scale = power2(exponent)
        result['rows'].append(([F(a) * scale for a in coeff],
                               None if lo is None else F(lo) * scale,
                               None if hi is None else F(hi) * scale))
    return result


def lp_models(count=160):
    """Bounded, feasible-by-construction dyadic models plus exact contradictions."""
    rng = random.Random(20261108)
    for case in range(count):
        n = 3 if case % 9 == 0 else 2
        bounds = [(F(-2), F(3)) for _ in range(n)]
        point = [F(rng.randint(-4, 6), 2) for _ in range(n)]
        if case % 13 == 0:
            bounds[-1] = (point[-1], point[-1])
        rows = []
        for _ in range(1 + case % 3):
            coeff = [F(rng.randint(-4, 4)) for _ in range(n)]
            if not any(coeff):
                coeff[0] = 1
            activity = sum(a * x for a, x in zip(coeff, point))
            kind = rng.choice(('L', 'G', 'E', 'R'))
            slack = F(rng.randint(0, 4), 2)
            lo = activity - slack if kind in ('G', 'R') else activity if kind == 'E' else None
            hi = activity + slack if kind in ('L', 'R') else activity if kind == 'E' else None
            rows.append((coeff, lo, hi))
        if case % 19 == 0:
            rows.append(([F(int(j == 0)) for j in range(n)], 4, None))
        yield {'cost': [F(rng.randint(-5, 5), 2) for _ in range(n)],
               'bounds': bounds, 'rows': rows, 'offset': F(rng.randint(-8, 8), 4),
               'maximize': bool(case % 2)}


class Numerics(unittest.TestCase):
    def invoke(self, model, extra=(), *, integer=False):
        with tempfile.TemporaryDirectory(prefix='milp numerical oracle ') as td:
            root = pathlib.Path(td)
            path, solution = root / 'model.mps', root / 'solution.sol'
            path.write_text(integer_mps(model) if integer else mps(model))
            command = [SOLVER, str(path), '--mip' if integer else '--relax',
                       '--solution', str(solution), *map(str, extra)]
            p = subprocess.run(command, text=True, capture_output=True, timeout=30)
            try:
                result = json.loads(p.stdout)
            except ValueError:
                self.fail(f'non-JSON stdout: {p.stdout!r}, stderr: {p.stderr!r}')
            self.assertIn(p.returncode, (0, 2, 3, 4), (result, p.stderr))
            self.assertFalse(result.get('certifies_optimality', False), result)
            self.assertFalse(result.get('bound_certified', False), result)
            for key, value in result.items():
                if isinstance(value, float):
                    self.assertTrue(math.isfinite(value), (key, result))
            values = {}
            check = None
            if solution.exists():
                for line in solution.read_text().splitlines():
                    if not line.strip() or line.lstrip().startswith('#'):
                        continue
                    name, value = line.split()
                    if name != '=obj=':
                        self.assertNotIn(name, values)
                        values[name] = F(value)
                if CHECKER:
                    args = [CHECKER, '--model', str(path), '--solution', str(solution)]
                    if not integer:
                        args.append('--relax')
                    q = subprocess.run(args, text=True, capture_output=True, timeout=15)
                    self.assertEqual(q.returncode, 0, (model, result, q.stdout, q.stderr))
                    check = json.loads(q.stdout)
                    self.assertTrue(check['checked_feasible'], check)
            self.assertEqual(bool(result.get('has_primal')), bool(values), result)
            return p.returncode, result, values, check

    def point(self, model, values, result, *, integer=False):
        self.assertEqual(set(values), {f'x{j}' for j in range(len(model['cost']))}, result)
        point = [values[f'x{j}'] for j in range(len(model['cost']))]
        for coeff, rhs in inequalities(model):
            residual = sum(a * x for a, x in zip(coeff, point)) - rhs
            self.assertLessEqual(residual, ABS_TOL, (model, point, residual, result))
        if integer:
            for j, flag in enumerate(model['integer']):
                if flag:
                    self.assertEqual(point[j].denominator, 1, (j, point, result))
        value = F(model.get('offset', 0)) + sum(F(c) * x for c, x in zip(model['cost'], point))
        # A JSON double cannot preserve sub-ULP changes beside a large offset;
        # exact point cost is checked separately and remains the primary gate.
        self.assertAlmostEqual(result['objective'], float(value),
                               delta=max(1e-7, math.ulp(float(value))))
        return value

    def compare(self, base, transformed=None, extra=(), *, integer=False):
        submitted = transformed or base
        exact = mip_oracle if integer else oracle
        expected = exact(base)
        self.assertEqual(exact(submitted), expected, 'test transformation changed the optimum')
        code, result, values, _ = self.invoke(submitted, extra, integer=integer)
        self.assertEqual(code, 0, (base, submitted, result))
        self.assertEqual(result['status'], 'infeasible' if expected is None else 'optimal',
                         (base, submitted, result))
        if expected is None:
            self.assertFalse(result['has_primal'], result)
            return result
        actual = self.point(submitted, values, result, integer=integer)
        # Also check the original unit system. Tiny submitted row coefficients
        # must not make a materially infeasible unscaled point look feasible.
        self.point(base, values, result, integer=integer)
        constant = F(base.get('offset', 0))
        for cost, (lo, hi) in zip(base['cost'], base['bounds']):
            if lo is not None and lo == hi:
                constant += F(cost) * F(lo)
        tolerance = ABS_TOL * (1 + abs(expected - constant))
        self.assertLessEqual(abs(actual - expected), tolerance,
                             (base, submitted, actual, expected, result))
        return result

    def test_160_scaled_bounded_exact_lp_oracles(self):
        for case, model in enumerate(lp_models()):
            exponents = [EXPONENTS[(case + 3 * row) % len(EXPONENTS)]
                         for row in range(len(model['rows']))]
            with self.subTest(case=case, exponents=exponents):
                self.compare(model, scaled_rows(model, exponents))

    def test_60_scaled_exact_milp_oracles(self):
        for case, model in enumerate(mip_models(60)):
            exponents = [EXPONENTS[(case + 2 * row) % len(EXPONENTS)]
                         for row in range(len(model['rows']))]
            with self.subTest(case=case, exponents=exponents):
                self.compare(model, scaled_rows(model, exponents), integer=True)

    def test_row_scaling_equivalence_dense_sparse_and_cold_warm(self):
        model = {'cost': [3, 2, -1], 'integer': [True, True, False],
                 'bounds': [(0, 3), (-1, 2), (0, 2)],
                 'rows': [([2, 2, 1], None, F(9, 2)), ([1, -1, 2], -1, 3)],
                 'maximize': True, 'offset': F(3, 4)}
        transformed = scaled_rows(model, [-26, 26])
        for backend in ('sparse', 'dense'):
            for cold in (False, True):
                extra = ['--backend', backend] + (['--cold'] if cold else [])
                with self.subTest(backend=backend, cold=cold):
                    self.compare(model, transformed, extra, integer=True)

    def test_free_fixed_variables_and_ranged_equality_units(self):
        models = [
            {'cost': [2, -1], 'bounds': [(None, None), (3, 3)],
             'rows': [([1, 1], 5, 5), ([1, 0], -3, 4)]},
            {'cost': [1, -1], 'bounds': [(None, None), (0, 2)],
             'rows': [([1, 1], 1, 3), ([1, -1], -1, None)]},
            {'cost': [1, 2], 'bounds': [(0, 3), (0, 3)],
             'rows': [([1, 1], 2, 2), ([2, 2], 4, 4), ([0, 0], 0, 0)]},
        ]
        for i, model in enumerate(models):
            exponents = [-26, 26, -10][:len(model['rows'])]
            with self.subTest(case=i):
                self.compare(model, scaled_rows(model, exponents))

    def test_column_units_spanning_twelve_orders(self):
        for exponent in (20, 26):
            large, small = power2(exponent), power2(-exponent)
            for maximum in (False, True):
                model = {'cost': [3 * large, 2 * small],
                         'bounds': [(0, small), (0, large)], 'offset': 5,
                         'rows': [([large, small], F(5, 4), None),
                                  ([2 * large, small], None, 2)], 'maximize': maximum}
                with self.subTest(exponent=exponent, maximum=maximum):
                    self.compare(model)

    def test_large_constants_do_not_change_scaled_integer_choice(self):
        for fixed_constant in (False, True):
            model = {'cost': [0, 1], 'integer': [True, False],
                     'bounds': [(0, 1), (-2, 100)],
                     'rows': [([1, 1], 0, None),
                              ([-2**30, 1], -(2**30 - 32), None)],
                     'offset': 0 if fixed_constant else 10**12}
            if fixed_constant:
                model['cost'].append(10**12); model['integer'].append(False)
                model['bounds'].append((1, 1))
                model['rows'] = [(a + [0], lo, hi) for a, lo, hi in model['rows']]
            transformed = scaled_rows(model, [20, -20])
            for cold in (False, True):
                with self.subTest(fixed_constant=fixed_constant, cold=cold):
                    self.compare(model, transformed, ['--cold'] if cold else [], integer=True)

    def test_well_conditioned_ablation_modes_preserve_answer(self):
        model = {'cost': [-2, 1], 'integer': [True, False], 'bounds': [(0, 3), (0, 2)],
                 'rows': [([1, 3], 4, None)], 'offset': F(7, 4)}
        for integer in (False, True):
            for extra in (['--no-scaling'], ['--no-refinement'],
                          ['--no-scaling', '--no-refinement']):
                with self.subTest(integer=integer, extra=extra):
                    result = self.compare(model, extra=extra, integer=integer)
                    if '--no-scaling' in extra:
                        self.assertEqual(result['scaled_rows'], 0, result)
                        self.assertEqual(result['scaled_columns'], 0, result)
                    if '--no-refinement' in extra:
                        self.assertEqual(result['refinement_steps'], 0, result)

    def test_scaling_diagnostics_are_finite_and_explicit(self):
        model = {'cost': [1, 2], 'bounds': [(0, 4), (0, 4)],
                 'rows': [([1, 1], F(3, 2), None), ([2, -1], None, 3)]}
        result = self.compare(model, scaled_rows(model, [-26, 26]))
        for name in ('refinement_steps', 'scaled_rows', 'scaled_columns'):
            self.assertIn(name, result)
            self.assertIsInstance(result[name], int)
            self.assertGreaterEqual(result[name], 0)
        self.assertGreater(result['scaled_rows'], 0, result)
        self.assertIn('max_backward_error', result)
        self.assertGreaterEqual(result['max_backward_error'], 0, result)
        self.assertTrue(math.isfinite(result['max_backward_error']), result)

    def test_excluded_unsafe_entering_column_never_establishes_an_optimum(self):
        # Exact feasible set is x=0. Without scaling, the positive limiting
        # direction is too weak for a stable pivot. Rejecting its improving
        # column must remain an unresolved numerical outcome, not an optimum
        # or an unbounded-ray claim. Normal scaling solves the same model.
        model = {'cost': [-1], 'bounds': [(0, None)],
                 'rows': [([F(1, 10**13)], None, 0)]}
        self.compare(model)
        for flags in [[], ['--no-updates'], ['--bland-ratio']]:
            with self.subTest(flags=flags):
                code, result, values, check = self.invoke(model, ['--no-scaling', *flags])
                self.assertEqual(code, 4, result)
                self.assertEqual(result['status'], 'numerical_failure')
                self.assertIn('no validated primal pivot', result['message'])
                self.assertEqual(result['primal_candidate_rejections'], 1)
                self.assertFalse(values)
                self.assertIsNone(check)


if __name__ == '__main__':
    unittest.main()
