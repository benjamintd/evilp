#!/usr/bin/env python3
"""Genuine sub-ulp source negatives require a basis change, never clipping."""
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

from test_lp_solver import inequalities, mps
from test_mip_solver import integer_mps

SOLVER = str(pathlib.Path(sys.argv.pop(1)).resolve())
CHECKER = str(pathlib.Path(sys.argv.pop(1)).resolve())
TOL = F(1, 10**7)


def coefficient(power=18, numerator=1):
    # Each value rounds to binary1, while the original coefficient is >1.
    return F(1) + F(numerator, 10**power)


def loaded_model(blocks=1):
    n = 2 * blocks
    rows, costs = [], []
    for i in range(blocks):
        a = [F(0)] * n
        a[2*i], a[2*i+1] = 1, coefficient(18+i)
        rows.append((a, 1, 1))
        costs.extend([0, -1])
    return {'cost': costs, 'bounds': [(0, 1)] * n, 'rows': rows}


def loaded_columns(blocks):
    # Source equality rows first, then implicit lower/upper rows per variable.
    structural = 2 * blocks
    basis = list(range(0, structural, 2))
    for j in range(structural):
        lower_row = blocks + 2*j
        basis.extend([structural + lower_row,
                      j if j % 2 else structural + lower_row + 1])
    return basis


def natural_model(mip=False):
    # Phase I leaves x0=x1=1,y=0. Increasing y gives a binary ratio tie;
    # Bland leaves x0, while exact source x1 becomes -1e-18 at y=1.
    model = {'cost': [0, 0, -1], 'bounds': [(0, 1)] * 3,
             'rows': [([1, 0, 1], 1, 1), ([0, 1, coefficient()], 1, 1)]}
    if mip:
        model['cost'].append(0)
        model['bounds'].append((0, 1))
        model['rows'] = [(list(row) + [0], lo, hi) for row, lo, hi in model['rows']]
        model['rows'].append(([0, 0, 0, 1], 0, 0))
        model['integer'] = [False, False, False, True]
    return model


class SourceSignRepair(unittest.TestCase):
    def run_model(self, model, extra=(), blocks=None, require_optimal=True, retained_columns=None):
        with tempfile.TemporaryDirectory(prefix='source sign ') as td:
            root = pathlib.Path(td)
            path, solution = root / 'source.mps', root / 'point.sol'
            mip_mode = 'integer' in model
            path.write_text(integer_mps(model) if mip_mode else mps(model))
            source_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            arguments = []
            if blocks is not None or retained_columns is not None:
                basis = loaded_columns(blocks) if retained_columns is None else retained_columns
                size = len(basis)
                total_columns = len(model['cost']) + 2*size
                retained = root / 'retained.basis'
                retained.write_text(f'MILPLAB_BASIS_1 {source_hash} {size} {total_columns} {size}\n' +
                                    ''.join(f'{i} {j}\n' for i, j in enumerate(basis)))
                # Column IDs in this retained adversarial basis refer to the
                # explicit-bound representation, independent of the default.
                arguments = ['--keep-bound-rows', '--read-basis', str(retained)]
            p = subprocess.run([SOLVER, str(path), '--mip' if mip_mode else '--relax',
                                '--time-limit', '10', '--solution', str(solution),
                                *arguments, *map(str, extra)], text=True, capture_output=True, timeout=30)
            result = json.loads(p.stdout)
            self.assertFalse(result.get('certifies_optimality', False), result)
            self.assertFalse(result.get('certifies_infeasibility', False), result)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source_hash)
            if not require_optimal:
                self.assertEqual(p.returncode, 4, (result, p.stderr))
                self.assertIn(result['status'], ('numerical_failure', 'iteration_limit', 'time_limit'), result)
                self.assertFalse(result['has_primal'], result)
                self.assertFalse(solution.exists(), result)
                return result
            self.assertEqual(p.returncode, 0, (result, p.stderr))
            self.assertEqual(result['status'], 'optimal', result)
            self.assertTrue(result['has_primal'], result)
            values = {name: F(value) for name, value in
                      (line.split() for line in solution.read_text().splitlines()) if name != '=obj='}
            point = [values[f'x{j}'] for j in range(len(model['cost']))]
            for a, b in inequalities(model):
                self.assertLessEqual(sum(v*x for v, x in zip(a, point)) - b, TOL, (result, point))
            checked = subprocess.run([CHECKER, '--model', str(path), '--solution', str(solution),
                                      '--lin-tol', '1e-7', '--int-tol', '1e-7', '--obj-tol', '1e-7',
                                      *([] if mip_mode else ['--relax'])],
                                     capture_output=True, text=True, timeout=20)
            self.assertEqual(checked.returncode, 0, (result, checked.stdout, checked.stderr))
            self.assertTrue(json.loads(checked.stdout)['checked_feasible'])
            self.assertGreater(result['exact_primal_recoveries'], 0, result)
            return result

    def test_generated_source_negative_with_binary_zero_needs_real_pivot(self):
        for power in (18, 20, 24):
            model = loaded_model()
            model['rows'][0][0][1] = coefficient(power)
            self.assertEqual(float(coefficient(power)), 1)
            for flags in ((), ('--no-updates',), ('--backend', 'dense')):
                with self.subTest(power=power, flags=flags):
                    result = self.run_model(model, flags, blocks=1)
                    self.assertEqual(result['primal_repair_iterations'], 1, result)
                    self.assertGreaterEqual(result['exact_primal_attempts'], 2, result)

    def test_disabled_repair_still_rejects_every_source_negative(self):
        result = self.run_model(loaded_model(), ['--no-primal-repair'], blocks=1,
                                require_optimal=False)
        self.assertEqual(result['status'], 'numerical_failure')
        self.assertEqual(result['primal_repair_iterations'], 0)

    def test_upper_only_transform_and_maximization_use_canonical_source_sign(self):
        upper_only = {'cost': [0, -1], 'bounds': [(None, 0), (0, 1)],
                      'rows': [([-1, coefficient()], 1, 1)]}
        result = self.run_model(upper_only, retained_columns=[0, 3, 4, 1])
        self.assertEqual(result['primal_repair_iterations'], 1, result)
        maximized = loaded_model()
        maximized['cost'] = [0, 1]
        maximized['maximize'] = True
        result = self.run_model(maximized, blocks=1)
        self.assertEqual(result['primal_repair_iterations'], 1, result)

    def test_multiple_hints_are_reconstructed_after_each_pivot(self):
        result = self.run_model(loaded_model(2), blocks=2)
        self.assertEqual(result['primal_repair_iterations'], 2, result)
        self.assertEqual(result['primal_repairs'], 2, result)
        self.assertEqual(result['exact_primal_attempts'], 3, result)

    def test_source_repair_consumes_existing_pivot_budget(self):
        disabled = self.run_model(loaded_model(2), ['--no-primal-repair'], blocks=2,
                                  require_optimal=False)
        result = self.run_model(loaded_model(2), ['--iteration-limit', disabled['iterations'] + 1],
                                blocks=2, require_optimal=False)
        self.assertEqual(result['status'], 'iteration_limit', result)
        self.assertEqual(result['primal_repair_iterations'], 1, result)

    def test_natural_root_basis_with_presolve_warm_cold_and_update_controls(self):
        for mip_mode in (False, True):
            controls = [(), ('--no-updates',), ('--backend', 'dense')]
            if mip_mode:
                controls += [('--cold',), ('--no-presolve',)]
            for flags in controls:
                with self.subTest(mip=mip_mode, flags=flags):
                    result = self.run_model(natural_model(mip_mode), flags)
                    self.assertGreater(result['primal_repair_iterations'], 0, result)
                    if mip_mode and not flags:
                        self.assertGreater(result['presolve_rows'], 0, result)
                        self.assertGreater(result['warm_starts'], 0, result)

    def test_failed_source_sign_pivot_never_claims_infeasibility(self):
        # Binary x=1 appears feasible, but source (1-epsilon)*x=1 needs x>1.
        # There is no eligible dual column. A failed repair is not a proof.
        model = {'cost': [0], 'bounds': [(0, 1)],
                 'rows': [([F(1)-F(1, 10**18)], 1, 1)]}
        result = self.run_model(model, require_optimal=False)
        self.assertEqual(result['status'], 'numerical_failure', result)


if __name__ == '__main__':
    unittest.main()
