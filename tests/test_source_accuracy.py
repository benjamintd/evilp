#!/usr/bin/env python3
"""Material source-decimal errors must not become wrong optimal answers."""
from fractions import Fraction as F
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import unittest


def case_text(kind, large=False, fixed=False):
    columns = 300 if large else (2 if kind == 'coupled' else 1)
    rows = 200 if large else 1
    names = ['X', 'Y'][:2 if kind == 'coupled' else 1]
    text = ['NAME source_accuracy', 'OBJSENSE', ' MAX', 'ROWS', ' N OBJ']
    text += [' E LINK' if kind == 'coupled' else ' L LIMIT']
    text += [f' L Z{i}' for i in range(1, rows)]
    text += ['COLUMNS', " M0 'MARKER' 'INTORG'"]
    if kind == 'row':
        text += [' X OBJ 1', ' X LIMIT 1.000000000000000001', ' X LIMIT -1']
    elif kind == 'objective':
        text += [' X OBJ 1.000000000000000001', ' X OBJ -1']
    else:
        text += [' X OBJ 1.000000000000000001 LINK 1', ' Y OBJ -1 LINK -1']
    text += [" M1 'MARKER' 'INTEND'"]
    text += [f' P{j} OBJ 0' for j in range(columns - len(names))]
    text += ['RHS', ' rhs LINK 0' if kind == 'coupled' else ' rhs LIMIT 0', 'BOUNDS']
    for name in names:
        text += [f' FX b {name} 1000000000000000'] if fixed else [
            f' LO b {name} 0', f' UP b {name} 1000000000000000']
    text += [f' FX b P{j} 0' for j in range(columns - len(names))]
    text += ['ENDATA', '']
    return '\n'.join(text), names


class SourceAccuracy(unittest.TestCase):
    def exercise(self, kind, large=False, fixed=False, mode='--mip', flags=()):
        # The exact one-dimensional oracles use t=X, or X=Y=t. The source
        # coefficient difference is 10^-18, while every input binary double
        # in that difference cancels to zero. Row case forces t=0; objective
        # cases maximize 10^-18*t over [0,10^15], hence optimum 1/1000.
        text, names = case_text(kind, large, fixed)
        target = F(0) if kind == 'row' else F(10**15)
        optimum = target if kind == 'row' else F(1, 10**18) * target
        with tempfile.TemporaryDirectory(prefix='milplab-source-accuracy-') as temporary:
            root = Path(temporary)
            model, solution = root / 'model.mps', root / 'point.sol'
            model.write_text(text)
            command = [str(SOLVER), str(model), mode, '--time-limit', '10',
                       '--solution', str(solution), *flags]
            run = subprocess.run(command, text=True, capture_output=True, timeout=30)
            self.assertEqual(run.returncode, 0, (command, run.stdout, run.stderr))
            result = json.loads(run.stdout)
            self.assertEqual(result['status'], 'optimal', result)
            self.assertFalse(result['certifies_optimality'])
            point = dict(line.split() for line in solution.read_text().splitlines()
                         if line.strip() and not line.startswith('#'))
            for name in names:
                self.assertEqual(F(point[name]), target, (kind, result, point))
            checked = subprocess.run([str(CHECKER), '--model', str(model), '--solution',
                                      str(solution), *(['--relax'] if mode == '--relax' else [])],
                                     text=True, capture_output=True, timeout=30)
            self.assertEqual(checked.returncode, 0, (checked.stdout, checked.stderr))
            validated = json.loads(checked.stdout)
            self.assertTrue(validated['checked_feasible'])
            self.assertEqual(F(validated['objective']), optimum)
            self.assertAlmostEqual(result['objective'], float(optimum), delta=1e-7)

    def test_duplicate_row_coefficients(self):
        for large in [False, True]:
            for flags in [[], ['--no-presolve'], ['--cold'], ['--no-scaling']]:
                with self.subTest(large=large, flags=flags):
                    self.exercise('row', large, flags=flags)

    def test_duplicate_objective_coefficients(self):
        for large in [False, True]:
            for fixed in [False, True]:
                for mode in ['--relax', '--mip']:
                    with self.subTest(large=large, fixed=fixed, mode=mode):
                        self.exercise('objective', large, fixed, mode)

    def test_joint_objective_requires_real_optimization(self):
        for large in [False, True]:
            for fixed in [False, True]:
                for mode in ['--relax', '--mip']:
                    with self.subTest(large=large, fixed=fixed, mode=mode):
                        self.exercise('coupled', large, fixed, mode)

    def test_joint_objective_warm_cold_scaling_and_heuristic_controls(self):
        for flags in [['--cold'], ['--no-scaling'], ['--backend', 'dense'],
                      ['--no-root-diving', '--no-root-cuts'], ['--no-presolve']]:
            with self.subTest(flags=flags):
                self.exercise('coupled', large=True, flags=flags)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--checker', type=Path, required=True)
    args, remaining = parser.parse_known_args()
    SOLVER, CHECKER = args.solver.resolve(), args.checker.resolve()
    unittest.main(argv=[__file__, *remaining])
