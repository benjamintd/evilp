"""Independent original-model checks for the optional bounded-variable basis lift."""
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from test_lp_solver import mps, oracle
from test_source_accuracy import case_text
from test_pricing_integration import random_models, cancellation_models

SOLVER, CHECKER = map(lambda p: str(Path(p).resolve()), sys.argv[1:3])
del sys.argv[1:3]


class BasisInitializer(unittest.TestCase):
    def test_original_lp_points_and_bounds(self):
        models = list(random_models()) + cancellation_models() + [
            # Reflected upper-only and split free columns, duplicate equalities,
            # a fixed column and nonbasic finite upper bounds.
            {'cost': [1, -2, 3, -4], 'bounds': [(None, 6), (None, None), (1, 1), (-2, 7)],
             'rows': [([1, 1, 0, 0], 2, 2), ([2, 2, 0, 0], 4, 4),
                      ([0, 1, 0, 0], -1, 3), ([0, 0, 1, 1], -1, 8)]},
            {'cost': [-3, 2], 'bounds': [(0, 5), (-2, 4)], 'rows': []},
            # GLPK's binary64 near-tie still needs native source pricing.
            {'cost': [1, Fraction('1.00000000000000001')], 'maximize': True,
             'bounds': [(0, 1), (0, 1)], 'rows': [([1, 1], 1, 1)]},
        ]
        accepted = 0
        for number, model in enumerate(models):
            with self.subTest(number=number), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                path, point, dual = [root / n for n in ['model.mps', 'point.sol', 'bound.sol']]
                path.write_text(mps(model))
                cmd = [SOLVER, str(path), '--relax', '--glpk-basis', '--solution', str(point),
                       '--dual', str(dual)]
                run = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
                self.assertEqual(run.returncode, 0, (run.stdout, run.stderr))
                result = json.loads(run.stdout)
                expected = oracle(model)
                if expected is None:
                    self.assertEqual(result['status'], 'infeasible')
                    self.assertFalse(point.exists())
                    continue
                self.assertEqual(result['status'], 'optimal', result)
                self.assertLessEqual(result['glpk_basis_iterations'], result['iterations'])
                accepted += result['glpk_bases_accepted']
                check = subprocess.run([CHECKER, '--model', str(path), '--solution', str(point),
                                        '--relax', '--dual', str(dual)], capture_output=True,
                                       text=True, timeout=10)
                self.assertEqual(check.returncode, 0, (check.stdout, check.stderr))
                checked = json.loads(check.stdout)
                self.assertTrue(checked['checked_feasible'])
                self.assertEqual(checked['bound_status'], 'valid')
                bound = Fraction(checked['certified_bound'])
                if model.get('maximize'):
                    self.assertGreaterEqual(bound, expected)
                else:
                    self.assertLessEqual(bound, expected)
                self.assertLessEqual(abs(Fraction(checked['objective']) - expected), Fraction('1e-7'))
                self.assertLessEqual(abs(bound - expected), Fraction('1e-7'))
        self.assertGreater(accepted, 15)

    def test_material_decimal_under_rounded_cost_tie(self):
        # The preferred integer point in original decimal arithmetic differs
        # by 1e-3 from the origin. GLPK sees a tie in binary64.
        text, _ = case_text('coupled')
        for mode in ('--relax', '--mip'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                path, point = [Path(directory) / n for n in ('model.mps', 'point.sol')]
                path.write_text(text)
                run = subprocess.run([SOLVER, str(path), mode, '--glpk-basis', '--solution',
                                      str(point)], capture_output=True, text=True, timeout=10)
                self.assertEqual(run.returncode, 0, (run.stdout, run.stderr))
                result = json.loads(run.stdout)
                self.assertEqual(result['status'], 'optimal', result)
                check = subprocess.run([CHECKER, '--model', str(path), '--solution', str(point),
                                        *(['--relax'] if mode == '--relax' else [])],
                                       capture_output=True, text=True, timeout=10)
                self.assertEqual(check.returncode, 0, (check.stdout, check.stderr))
                self.assertEqual(Fraction(json.loads(check.stdout)['objective']), Fraction('0.001'))

    def test_shared_pivot_limit(self):
        model = {'cost': [1, 1, 1], 'bounds': [(0, 2)] * 3,
                 'rows': [([int(i == j) for i in range(3)], 1, None) for j in range(3)]}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'model.mps'
            path.write_text(mps(model))
            run = subprocess.run([SOLVER, str(path), '--relax', '--glpk-basis',
                                  '--iteration-limit', '1'], capture_output=True,
                                 text=True, timeout=10)
            result = json.loads(run.stdout)
            self.assertEqual(result['status'], 'iteration_limit')
            self.assertEqual(result['iterations'], 1)
            self.assertEqual(result['glpk_basis_iterations'], 1)


if __name__ == '__main__':
    unittest.main()
