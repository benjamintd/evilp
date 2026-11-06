#!/usr/bin/env python3
"""Original-unit recovery regressions with independent exact point oracles.

An optimal status is accepted only with an exported point that satisfies both
the source decimal model and its equivalent normalized units. Large fixed
objective constants never enlarge the allowed error in the variable cost.
"""
import copy
import hashlib
import json
import math
import pathlib
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

from test_lp_solver import inequalities, mps

SOLVER = str(pathlib.Path(sys.argv.pop(1)).resolve())
CHECKER = str(pathlib.Path(sys.argv.pop(1)).resolve())
TOLERANCE = F(1, 10**7)
LARGE = 2**53


def row_units(model, factors):
    result = copy.deepcopy(model)
    result['rows'] = [([F(a) * factor for a in coefficients],
                       None if lower is None else F(lower) * factor,
                       None if upper is None else F(upper) * factor)
                      for (coefficients, lower, upper), factor in zip(model['rows'], factors)]
    return result


def column_units(model, factors):
    """Source coordinate x_j = factor_j * submitted coordinate u_j."""
    result = copy.deepcopy(model)
    result['cost'] = [F(a) * scale for a, scale in zip(model['cost'], factors)]
    result['bounds'] = [(None if lo is None else F(lo) / scale,
                         None if hi is None else F(hi) / scale)
                        for (lo, hi), scale in zip(model['bounds'], factors)]
    result['rows'] = [([F(a) * scale for a, scale in zip(coefficients, factors)], lower, upper)
                      for coefficients, lower, upper in model['rows']]
    return result


def fixed_shift_model(z=LARGE):
    # Rounding .1*z and subtracting in binary loses this genuine source gap.
    product = F(z, 10)
    rhs = F(-(-product.numerator * 4 // product.denominator), 4)
    if rhs == product:
        rhs += F(1, 4)
    gap = rhs - product
    return {'cost': [1, 0], 'bounds': [(0, 1), (z, z)],
            'rows': [([1, F(1, 10)], rhs, None)]}, gap


class PrimalRecovery(unittest.TestCase):
    def validate_point(self, model, point, result):
        for coefficients, rhs in inequalities(model):
            residual = sum(a * x for a, x in zip(coefficients, point)) - rhs
            self.assertLessEqual(residual, TOLERANCE, (model, point, residual, result))

    def solve(self, model, expected, *, normalized=None, coordinate_factors=None,
              extra=(), require_optimal=False, retained_slack_basis=False):
        with tempfile.TemporaryDirectory(prefix='primal recovery ') as td:
            root = pathlib.Path(td)
            path, solution = root / 'model.mps', root / 'solver.sol'
            path.write_text(mps(model))
            basis_arguments = []
            if retained_slack_basis:
                basis = root / 'retained.basis'
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                basis.write_text(f'MILPLAB_BASIS_1 {digest} 3 7 3\n0 1\n1 2\n2 3\n')
                # This minimized historical fixture names columns in the
                # explicit-bound representation. Retain that representation.
                basis_arguments = ['--keep-bound-rows', '--read-basis', str(basis)]
            p = subprocess.run([SOLVER, str(path), '--relax', '--solution', str(solution),
                                *basis_arguments, *extra],
                               text=True, capture_output=True, timeout=20)
            result = json.loads(p.stdout)
            self.assertFalse(result.get('certifies_optimality', False), result)
            self.assertFalse(result.get('certifies_infeasibility', False), result)
            self.assertIn(result['status'], ('optimal', 'numerical_failure'), result)
            if require_optimal:
                self.assertEqual(result['status'], 'optimal', (model, result))
            if result['status'] == 'numerical_failure':
                self.assertEqual(p.returncode, 4, (result, p.stderr))
                self.assertFalse(result['has_primal'], result)
                self.assertFalse(solution.exists(), result)
                return result
            self.assertEqual(p.returncode, 0, (result, p.stderr))
            self.assertTrue(result['has_primal'], result)
            self.assertLessEqual(F(str(result['primal_residual'])), TOLERANCE)
            values = {}
            for line in solution.read_text().splitlines():
                name, value = line.split()
                if name != '=obj=':
                    self.assertNotIn(name, values)
                    values[name] = F(value)
            self.assertEqual(set(values), {f'x{j}' for j in range(len(model['cost']))})
            point = [values[f'x{j}'] for j in range(len(model['cost']))]
            self.validate_point(model, point, result)
            original_point = [x * scale for x, scale in zip(point, coordinate_factors or [1] * len(point))]
            self.validate_point(normalized or model, original_point, result)
            objective = F(model.get('offset', 0)) + sum(F(a) * x for a, x in zip(model['cost'], point))
            self.assertLessEqual(abs(objective - expected), TOLERANCE,
                                 (model, point, objective, expected, result))
            # The JSON objective is a double; its spacing beside a large fixed
            # offset cannot replace the exact point-cost comparison above.
            self.assertLessEqual(abs(result['objective'] - float(objective)),
                                 max(1e-7, math.ulp(float(objective))))
            # Remove the rounded claimed objective so it cannot obscure an
            # independent original-point check beside a huge objective offset.
            witness = root / 'point.sol'
            witness.write_text(''.join(f'{name} {value}\n' for name, value in values.items()))
            q = subprocess.run([CHECKER, '--model', str(path), '--solution', str(witness),
                                '--format', 'rational', '--relax', '--lin-tol', '1e-7',
                                '--int-tol', '0', '--obj-tol', '1e-7'],
                               text=True, capture_output=True, timeout=20)
            self.assertEqual(q.returncode, 0, (model, result, q.stdout, q.stderr))
            checked = json.loads(q.stdout)
            self.assertTrue(checked['checked_feasible'], checked)
            self.assertFalse(checked['certifies_optimality'], checked)
            return result

    def test_well_conditioned_original_decimal_gap_is_recovered(self):
        # Formerly both backends returned x=0, objective 0, residual 0; the exact
        # source row is violated by 1/20. The remaining LP has a 1x1 unit basis.
        for backend in ('sparse', 'dense'):
            for z in (LARGE, LARGE + 2, LARGE + 4, LARGE + 6):
                model, expected = fixed_shift_model(z)
                with self.subTest(backend=backend, fixed=z, expected=expected):
                    self.solve(model, expected, extra=['--backend', backend], require_optimal=True)

    def test_tiny_row_units_do_not_hide_the_original_gap(self):
        model, expected = fixed_shift_model()
        for factor in (F(1, 2**80), F(1, 2**40), F(2**40)):
            for backend in ('sparse', 'dense'):
                with self.subTest(row_factor=factor, backend=backend):
                    self.solve(row_units(model, [factor]), expected, normalized=model,
                               extra=['--backend', backend])

    def test_retained_basis_negative_slack_is_restored_in_original_units(self):
        # This public retained-basis fixture uses only a unit row and dyadic
        # width. Its all-slack basis has exact canonical slack -2^-39 while
        # x=0 violates the source row by 1. Cost 0 makes the basis dual feasible,
        # and the feasible adjacent interval [1,2^39] has known optimum 0.
        model = {'cost': [0], 'bounds': [(0, 2**39)], 'rows': [([1], 1, None)]}
        for backend in ('sparse', 'dense'):
            with self.subTest(backend=backend):
                self.solve(model, F(0), extra=['--backend', backend], require_optimal=True,
                           retained_slack_basis=True)

    def test_column_units_and_scaling_controls_preserve_point_cost(self):
        model, expected = fixed_shift_model()
        controls = [(), ('--backend', 'dense'), ('--no-scaling',),
                    ('--no-refinement',), ('--no-exact',)]
        for factor in (F(1, 2**30), F(2**30)):
            submitted = column_units(model, [factor, 1])
            for extra in controls:
                with self.subTest(column_factor=factor, extra=extra):
                    self.solve(submitted, expected, normalized=model,
                               coordinate_factors=[factor, 1], extra=extra)

    def test_large_fixed_objective_offsets_do_not_mask_variable_cost_error(self):
        model, gap = fixed_shift_model()
        model['cost'][1], model['offset'] = F(7, 8), 10**18
        constant = F(10**18) + F(7 * LARGE, 8)
        for extra in ((), ('--backend', 'dense'), ('--no-scaling',)):
            with self.subTest(extra=extra):
                self.solve(model, constant + gap, extra=extra)

    def test_scaled_slack_cancellation_uses_original_rows_and_known_optimum(self):
        for scale in (10**6, 10**9, 10**12):
            for lower_side in (False, True):
                coefficients = [scale, 1 if lower_side else -1]
                rhs = scale + (1 if lower_side else -1)
                model = {'cost': [0, 1], 'bounds': [(0, 2), (0, 2)],
                         'rows': [([1, 0], 1, 1),
                                  (coefficients, rhs if lower_side else None,
                                   None if lower_side else rhs)]}
                for row_factor in (F(1, 2**40), F(1), F(2**40)):
                    with self.subTest(scale=scale, lower=lower_side, row_factor=row_factor):
                        self.solve(row_units(model, [1, row_factor]), F(1), normalized=model,
                                   require_optimal=scale == 10**6)

    def test_fixed_shift_summation_cancellation(self):
        model = {'cost': [1, 0, 0, 0],
                 'bounds': [(0, 2), (10**16, 10**16), (1, 1), (-10**16, -10**16)],
                 'rows': [([1, 1, 1, 1], 1, None)]}
        for extra in ((), ('--backend', 'dense'), ('--no-scaling',)):
            with self.subTest(extra=extra):
                self.solve(model, F(0), extra=extra)

    def test_shift_plus_structural_value_recovers_small_original_variable(self):
        for sign in (-1, 1):
            model = {'cost': [sign],
                     'bounds': [(-10**16, 2) if sign > 0 else (-2, 10**16)],
                     'rows': [([1], 1, None) if sign > 0 else ([1], None, -1)]}
            for extra in ((), ('--backend', 'dense'), ('--no-scaling',)):
                with self.subTest(sign=sign, extra=extra):
                    self.solve(model, F(1), extra=extra)

    def test_duplicate_source_objective_preserves_lp_and_milp_optimum(self):
        # The source sum is exactly -1, although summing the independently
        # rounded binary records gives0. Feasibility alone cannot catch the
        # wrong point x=0: the exact known optimum is x=1 with value -1.
        model = ('NAME COST\nROWS\n N O\nCOLUMNS\n'
                 ' x0 O 10000000000000000\n x0 O -10000000000000001\n'
                 'BOUNDS\n BV b x0\nENDATA\n')
        for mode in ('--relax', '--mip'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(prefix='source objective ') as td:
                root = pathlib.Path(td)
                path, solution = root / 'model.mps', root / 'solution.sol'
                path.write_text(model)
                p = subprocess.run([SOLVER, str(path), mode, '--solution', str(solution)],
                                   text=True, capture_output=True, timeout=20)
                self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
                result = json.loads(p.stdout)
                self.assertEqual(result['status'], 'optimal', result)
                self.assertEqual(F(str(result['objective'])), F(-1), result)
                self.assertFalse(result.get('certifies_optimality', False), result)
                values = dict(line.split() for line in solution.read_text().splitlines())
                self.assertEqual(F(values['x0']), F(1), (values, result))
                self.assertEqual(F(values['=obj=']), F(-1), (values, result))
                q = subprocess.run([CHECKER, '--model', str(path), '--solution', str(solution), '--exact'],
                                   text=True, capture_output=True, timeout=20)
                self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
                checked = json.loads(q.stdout)
                self.assertTrue(checked['checked_feasible'], checked)
                self.assertEqual(F(checked['objective']), F(-1))
                self.assertFalse(checked['certifies_optimality'], checked)


if __name__ == '__main__':
    unittest.main()
