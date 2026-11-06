#!/usr/bin/env python3
"""Objective claims must describe the decimal coordinates actually exported."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from test_source_accuracy import case_text


class SourceSolution(unittest.TestCase):
    def check_solution(self, model, solution, objective, mode, claimed=True):
        records = dict(line.split() for line in solution.read_text().splitlines())
        if claimed:
            self.assertEqual(F(records['=obj=']), objective, records)
        else:
            self.assertNotIn('=obj=', records)
        checked = subprocess.run([str(CHECKER), '--model', str(model), '--solution',
                                  str(solution), *(['--relax'] if mode == '--relax' else [])],
                                 capture_output=True, text=True, timeout=30)
        self.assertEqual(checked.returncode, 0, (checked.stdout, checked.stderr))
        validation = json.loads(checked.stdout)
        self.assertTrue(validation['checked_feasible'])
        self.assertEqual(F(validation['objective']), objective)
        self.assertEqual(F(validation['objective_error']), 0)
        return records

    def solve(self, text, objective, mode):
        with tempfile.TemporaryDirectory(prefix='milplab-source-solution-') as temporary:
            root = Path(temporary)
            model, solution = root / 'model.mps', root / 'point.sol'
            model.write_text(text)
            solved = subprocess.run([str(SOLVER), str(model), mode, '--time-limit', '10',
                                     '--solution', str(solution)],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(solved.returncode, 0, (solved.stdout, solved.stderr))
            result = json.loads(solved.stdout)
            self.assertEqual(result['status'], 'optimal', result)
            return self.check_solution(model, solution, objective, mode)

    def test_small_improvement_beside_large_source_offset(self):
        magnitude = '1000000000000000.000000000000000001'
        for large in [False, True]:
            for sign in [-1, 1]:
                for mode in ['--relax', '--mip']:
                    with self.subTest(large=large, sign=sign, mode=mode):
                        text, _ = case_text('coupled', large=large)
                        rhs = ('-' if sign > 0 else '') + magnitude
                        text = text.replace('RHS\n', f'RHS\n rhs OBJ {rhs}\n')
                        objective = sign * F(magnitude) + F(1, 1000)
                        records = self.solve(text, objective, mode)
                        self.assertEqual(F(records['X']), F(10**15))
                        self.assertEqual(F(records['Y']), F(10**15))

    def test_score_exported_decimal_coordinates(self):
        # The binary representation of .3 differs from the decimal witness.
        # A large coefficient amplifies that difference beyond 1e-7.
        for value in ['.3', '-.3', '3e-200']:
            with self.subTest(value=value):
                with tempfile.TemporaryDirectory(prefix='milplab-point-export-') as temporary:
                    root = Path(temporary)
                    model, solution = root / 'model.mps', root / 'point.sol'
                    text = ('NAME exported_decimal\nROWS\n N OBJ\n L R\nCOLUMNS\n'
                            ' X OBJ 1000000000000000\nRHS\n rhs R 0\nBOUNDS\n'
                            f' FX b X {value}\nENDATA\n')
                    model.write_text(text)
                    # Exercise serialization directly. These are point checks,
                    # not claims that the optimizer completes every fixed LP.
                    exported = subprocess.run([str(DRIVER), str(model), str(solution), value],
                                              capture_output=True, text=True, timeout=30)
                    self.assertEqual(exported.returncode, 0, exported.stderr)
                    records = self.check_solution(model, solution, F(value) * 10**15, '--relax')
                    self.assertEqual(F(records['X']), F(value))

    def test_unavailable_source_omits_optional_objective_claim(self):
        with tempfile.TemporaryDirectory(prefix='milplab-point-no-claim-') as temporary:
            root = Path(temporary)
            model, solution = root / 'model.mps', root / 'point.sol'
            model.write_text('NAME no_objective_source\nROWS\n N OBJ\n L R\nCOLUMNS\n'
                             ' X OBJ 1000000000000000\nRHS\n rhs R 0\nBOUNDS\n'
                             ' FX b X .3\nENDATA\n')
            exported = subprocess.run([str(DRIVER), str(model), str(solution), '.3',
                                       'without-objective-source'],
                                      capture_output=True, text=True, timeout=30)
            self.assertEqual(exported.returncode, 0, exported.stderr)
            self.check_solution(model, solution, F(3 * 10**14), '--relax', claimed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--checker', type=Path, required=True)
    parser.add_argument('--driver', type=Path, required=True)
    args, remaining = parser.parse_known_args()
    SOLVER, CHECKER = args.solver.resolve(), args.checker.resolve()
    DRIVER = args.driver.resolve()
    unittest.main(argv=[__file__, *remaining])
