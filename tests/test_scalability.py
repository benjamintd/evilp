#!/usr/bin/env python3
"""Sparse LP capacity and update contracts with exact generated optima.

Thousands of independent bound equations surround a three-variable coupled
core. Assertions use original points, an independent checker and algebraically
known objectives; refactor counts compare work, without timing thresholds.
"""
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

SOLVER = str(pathlib.Path(sys.argv.pop(1)).resolve())
CHECKER = str(pathlib.Path(sys.argv.pop(1)).resolve())
TOLERANCE = F(1, 10**7)


def instance(columns, active=0, integer=False):
    """Unique optimum (.5,.25,.75), then active ones and remaining zeros."""
    if columns < active + 3:
        raise ValueError('active variables must fit outside the coupled core')
    text = ['NAME SPARSE_CORE', 'ROWS', ' N OBJ', ' E A', ' E B', ' E C', 'COLUMNS']
    core = [(1, 2, 1), (1, -1, 2), (1, 1, -1)]
    for j in range(columns):
        name = f'x{j}'
        cost = (3, 2, 1)[j] if j < 3 else -1 if j < active + 3 else 1
        text.append(f' {name} OBJ {cost}')
        if j < 3:
            text.extend(f' {name} {row} {value}' for row, value in zip(('A', 'B', 'C'), core[j]))
    text += ['RHS', ' rhs A 1.5 B 1.5', ' rhs C .25', 'BOUNDS']
    for j in range(columns):
        if integer and j >= 3:
            text.append(f' BV b x{j}')
        else:
            text += [f' LO b x{j} 0', f' UP b x{j} 1']
    return '\n'.join(text + ['ENDATA', ''])


class Scalability(unittest.TestCase):
    def invoke(self, root, model, columns, active=0, extra=(), mode='--relax', expected=0):
        root.mkdir()
        path, solution = root / 'model.mps', root / 'solution.sol'
        path.write_text(model)
        p = subprocess.run([SOLVER, str(path), mode, '--solution', str(solution), *extra],
                           text=True, capture_output=True, timeout=120)
        self.assertEqual(p.returncode, expected, (p.stdout, p.stderr))
        result = json.loads(p.stdout)
        self.assertFalse(result.get('certifies_optimality', False), result)
        self.assertFalse(result.get('certifies_infeasibility', False), result)
        if expected:
            self.assertIn(result['status'], ('unsupported', 'input_error'), result)
            self.assertFalse(result.get('has_primal', False), result)
            self.assertFalse(result.get('has_incumbent', False), result)
            self.assertFalse(solution.exists(), result)
            return result
        self.assertEqual(result['status'], 'optimal', result)
        expected_objective = F(11, 4) - active
        self.assertLessEqual(abs(F(str(result['objective'])) - expected_objective), TOLERANCE)
        records = dict(line.split() for line in solution.read_text().splitlines())
        self.assertEqual(set(records), {'=obj='} | {f'x{j}' for j in range(columns)})
        values = [F(records[f'x{j}']) for j in range(columns)]
        truth = [F(1, 2), F(1, 4), F(3, 4)]
        objective = F(0)
        for j, value in enumerate(values):
            expected_value = truth[j] if j < 3 else F(int(j < active + 3))
            self.assertLessEqual(abs(value - expected_value), TOLERANCE, (j, value, result))
            self.assertGreaterEqual(value, -TOLERANCE)
            self.assertLessEqual(value, 1 + TOLERANCE)
            cost = (3, 2, 1)[j] if j < 3 else -1 if j < active + 3 else 1
            objective += cost * value
        for coefficients, rhs in [((1, 1, 1), F(3, 2)), ((2, -1, 1), F(3, 2)),
                                  ((1, 2, -1), F(1, 4))]:
            activity = sum(a * x for a, x in zip(coefficients, values))
            self.assertLessEqual(abs(activity - rhs), TOLERANCE, result)
        self.assertLessEqual(abs(objective - expected_objective), TOLERANCE)
        self.assertLessEqual(abs(F(records['=obj=']) - expected_objective), TOLERANCE)
        command = [CHECKER, '--model', str(path), '--solution', str(solution),
                   '--lin-tol', '1e-7', '--int-tol', '0', '--obj-tol', '1e-7']
        if mode == '--relax':
            command.append('--relax')
        q = subprocess.run(command, text=True, capture_output=True, timeout=120)
        self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
        checked = json.loads(q.stdout)
        self.assertTrue(checked['checked_feasible'], checked)
        self.assertFalse(checked['certifies_optimality'], checked)
        self.assertLessEqual(abs(F(checked['objective']) - expected_objective), TOLERANCE)
        return result

    def assert_sparse_core(self, result, columns):
        self.assertEqual(result['backend'], 'sparse_lu')
        self.assertEqual(result['columns'], columns)
        self.assertGreater(result['standard_columns'], 4096)
        self.assertGreater(result['standard_rows'], 4096)
        # Peak logical factor entries remain proportional to the sparse basis,
        # rather than its square. The coupled part has only three variables.
        self.assertLessEqual(result['factor_nnz_peak'], 8 * result['standard_rows'])
        self.assertLessEqual(result['factor_core_peak'], 3)
        self.assertGreater(result['factor_core_peak'], 0)
        self.assertGreaterEqual(result['singleton_eliminations_peak'], result['standard_rows'] - 3)

    def check_default_and_factor(self, root, model, columns, active=0, extra=()):
        root.mkdir(parents=True, exist_ok=True)
        result = self.invoke(root / 'default', model, columns, active, list(extra))
        if result.get('glpk_basis_status') == 'source_verified_point_and_bound':
            # The default point has already been checked on every original
            # variable. A separate basis-only solve must exercise all the
            # original sparse-factor structure assertions on the same input.
            self.assertEqual(result['refactorizations'], 0)
            self.assertGreater(result['source_support_scans'], 0)
            result = self.invoke(root / 'basis', model, columns, active, [*extra, '--glpk-basis'])
        self.assert_sparse_core(result, columns)

    def test_thousands_of_bounds_and_small_coupled_root(self):
        # Preserve a basis larger than the dense ceiling after implicit lower
        # rows disappear, and retain the original 2,500-column case separately.
        columns = 5000
        with tempfile.TemporaryDirectory(prefix='sparse capacity ') as td:
            self.check_default_and_factor(pathlib.Path(td) / 'implicit', instance(columns), columns)
            self.check_default_and_factor(pathlib.Path(td) / 'explicit', instance(2500), 2500,
                                          extra=['--keep-bound-rows'])

    def test_updates_reduce_refactors_without_changing_original_optimum(self):
        columns, active = 5000, 64
        model = instance(columns, active)
        with tempfile.TemporaryDirectory(prefix='basis update work ') as td:
            root = pathlib.Path(td)
            # Exercise the default starting policy on the same 64-active
            # problem. An imported optimal basis can already encode the
            # active upper bounds, leaving few simplex updates to measure.
            self.check_default_and_factor(root, model, columns, active)
            # Compare the two native update policies from the same cold basis.
            updated = self.invoke(root / 'updated', model, columns, active, ['--no-glpk-basis'])
            rebuilt = self.invoke(root / 'rebuilt', model, columns, active,
                                  ['--no-glpk-basis', '--no-updates'])
            self.assert_sparse_core(updated, columns)
            self.assert_sparse_core(rebuilt, columns)
            self.assertGreaterEqual(updated['iterations'], active)
            self.assertGreaterEqual(rebuilt['iterations'], active)
            self.assertGreaterEqual(updated['basis_updates'], active // 2)
            self.assertEqual(rebuilt['basis_updates'], 0)
            self.assertGreaterEqual(rebuilt['refactorizations'], updated['refactorizations'] + 16)
            self.assertGreaterEqual(rebuilt['refactorizations'], 2 * updated['refactorizations'])
            self.assertLessEqual(abs(F(str(updated['objective'])) - F(str(rebuilt['objective']))),
                                 TOLERANCE)

    def test_explicit_caps_dense_reference_and_cold_root(self):
        columns, active = 2500, 4
        model = instance(columns, active, integer=True)
        with tempfile.TemporaryDirectory(prefix='sparse caps ') as td:
            root = pathlib.Path(td)
            cold = self.invoke(root / 'cold', model, columns, active,
                               ['--cold', '--no-presolve'], mode='--mip')
            self.assertTrue(cold['cold'])
            self.assertTrue(cold['has_primal'])
            self.assertFalse(cold['presolve_enabled'])
            # A cap is an explicit unsupported/input outcome, never a claim of
            # mathematical infeasibility. Dense defaults must retain their cap.
            for i, extra in enumerate([['--max-dimension', '4096'], ['--backend', 'dense']]):
                result = self.invoke(root / f'cap-{i}', model, columns, active, extra, expected=3)
                self.assertEqual(result['status'], 'unsupported')
            result = self.invoke(root / 'cold-cap', model, columns, active,
                                 ['--cold', '--no-presolve', '--max-dimension', '4096'],
                                 mode='--mip', expected=3)
            self.assertEqual(result['status'], 'unsupported')
            refused = self.invoke(root / 'dense-million-cap', model, columns, active,
                                  ['--backend', 'dense', '--max-dimension', '1000000'], expected=2)
            self.assertIn('dense', refused['message'].lower())
            self.assertIn('4096', refused['message'])
            # Explicit dense reference still solves genuinely small models.
            self.invoke(root / 'small-dense', instance(12, 4), 12, 4, ['--backend', 'dense'])


if __name__ == '__main__':
    unittest.main()
