#!/usr/bin/env python3
"""Independent exact LP-support bound and rational witness contracts."""
import fractions
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest

BIN = str(pathlib.Path(sys.argv.pop(1)).resolve())
F = fractions.Fraction


class Certificates(unittest.TestCase):
    def call(self, model, dual=None, solution=None, extra=(), expected=0):
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            (root / 'model.mps').write_text(model)
            args = [BIN, '--model', str(root / 'model.mps')]
            if dual is not None:
                (root / 'dual.txt').write_text(dual)
                args += ['--dual', str(root / 'dual.txt')]
            if solution is not None:
                (root / 'solution.sol').write_text(solution)
                args += ['--solution', str(root / 'solution.sol')]
            p = subprocess.run(args + list(extra), text=True, capture_output=True)
            self.assertEqual(p.returncode, expected, (p.stdout, p.stderr))
            return json.loads(p.stdout)

    def fractional(self, maximize=False, bound=' FR b x\n'):
        return ('NAME FRAC\nOBJSENSE\n ' + ('MAX' if maximize else 'MIN') +
                '\nROWS\n N O\n ' + ('L' if maximize else 'G') +
                ' R\nCOLUMNS\n x O 2 R 3\nRHS\n b O -5 R 1\nBOUNDS\n' +
                bound + 'ENDATA\n')

    def test_exact_fraction_min_max_and_offset(self):
        for maximize in (False, True):
            result = self.call(self.fractional(maximize),
                               'R ' + ('-2/3' if maximize else '2/3') + '\n',
                               '=obj= 17/3\nx 1/3\n',
                               ['--format', 'rational', '--exact', '--relax'])
            self.assertEqual(result['certified_bound'], '17/3')
            self.assertEqual(result['objective'], '17/3')
            self.assertEqual(result['bound_direction'], 'upper' if maximize else 'lower')
            self.assertTrue(result['certifies_optimality'])
            self.assertTrue(result['relaxation'])

    def test_bound_only_does_not_invent_primal(self):
        result = self.call(self.fractional(), 'R +2/3\n')
        self.assertEqual(result['bound_status'], 'valid')
        self.assertEqual(result['certified_bound'], '17/3')
        self.assertFalse(result['checked_feasible'])
        self.assertFalse(result['certifies_optimality'])
        self.assertIsNone(result['objective'])

    def test_selected_objective_and_ignored_nonconstraint_row(self):
        model = ('NAME SELECT\nOBJNAME COST\nROWS\n N IGNORE\n N COST\n G R\n'
                 'COLUMNS\n x IGNORE 100 COST 2\n x R 3\nRHS\n b IGNORE 999 COST -5\n'
                 ' b R 1\nBOUNDS\n FR b x\nENDATA\n')
        result = self.call(model, 'R 2/3\n', 'x +1/3\n', ['--format', 'rational', '--exact'])
        self.assertEqual(result['certified_bound'], '17/3')
        self.assertTrue(result['certifies_optimality'])
        self.call(model, 'IGNORE 1\n', expected=2)

    def test_fraction_input_requires_explicit_format(self):
        self.call(self.fractional(), solution='x 1/3\n', expected=2)
        self.call(self.fractional(), solution='x 1/0\n', extra=['--format', 'rational'], expected=2)

    def test_residual_box_correction_and_empty_dual(self):
        model = ('NAME BOX\nROWS\n N O\nCOLUMNS\n x O -3\n'
                 'BOUNDS\n LO b x -1\n UP b x 2\nENDATA\n')
        result = self.call(model, '# all row multipliers zero\n', 'x 2\n', ['--exact'])
        self.assertEqual(result['certified_bound'], '-6')
        self.assertTrue(result['certifies_optimality'])

    def test_signed_ranges(self):
        for sense, ran, lo, hi in [('L', -2, 3, 5), ('G', -2, 5, 7),
                                    ('E', 2, 5, 7), ('E', -2, 3, 5)]:
            model = (f'NAME RANGES\nROWS\n N O\n {sense} R\nCOLUMNS\n x O 1 R 1\n'
                     f'RHS\n b R 5\nRANGES\n g R {ran}\nBOUNDS\n FR b x\nENDATA\n')
            result = self.call(model, 'R 1\n', f'x {lo}\n', ['--exact'])
            self.assertEqual(F(result['certified_bound']), lo)
            self.assertTrue(result['certifies_optimality'])
            # A valid but looser multiplier must not be mistaken for a proof.
            model = model.replace('x O 1 R 1', 'x O -1 R 1')
            result = self.call(model, 'R -1\n', f'x {hi}\n', ['--exact'])
            self.assertEqual(F(result['certified_bound']), -hi)
            self.assertTrue(result['certifies_optimality'])

    def test_selected_infinite_sides_are_unavailable(self):
        result = self.call(self.fractional(), 'R -1\n', expected=4)
        self.assertEqual(result['bound_status'], 'unavailable')
        self.assertIsNone(result['certified_bound'])
        self.assertFalse(result['certifies_optimality'])
        result = self.call(self.fractional(), '', expected=4)
        self.assertIn('column x', result['bound_message'])
        self.assertIsNone(result['certified_bound'])

    def test_relaxation_does_not_certify_integer_model(self):
        model = ('NAME INTEGER\nROWS\n N O\n G R\nCOLUMNS\n x O 1 R 2\n'
                 'RHS\n b R 1\nBOUNDS\n BV b x\nENDATA\n')
        result = self.call(model, 'R 1/2\n', 'x 1/2\n',
                           ['--format', 'rational', '--exact'], expected=1)
        self.assertEqual(result['integrality_violations'], 1)
        self.assertFalse(result['certifies_optimality'])
        result = self.call(model, 'R 1/2\n', 'x 1/2\n',
                           ['--format', 'rational', '--exact', '--relax'])
        self.assertTrue(result['certifies_optimality'])
        self.assertTrue(result['relaxation'])

    def test_tolerance_feasibility_is_not_exact_proof(self):
        model = ('NAME TOL\nROWS\n N O\n E R\nCOLUMNS\n x O 0 R 1\n'
                 'RHS\n b R 1\nENDATA\n')
        result = self.call(model, '', 'x 1.000001\n')
        self.assertTrue(result['checked_feasible'])
        self.assertEqual(result['certified_bound'], '0')
        self.assertFalse(result['certifies_optimality'])
        self.call(model, '', 'x 1.000001\n', ['--exact'], expected=1)

    def test_bad_certificates_fail_closed(self):
        for dual in ['R 1\nR 1\n', 'missing 1\n', 'O 1\n', 'R nan\n', 'R 1/0\n',
                     'R 1/2/3\n', 'R -/2\n', 'R 2/-3\n', 'R\n', 'R 1 extra\n']:
            result = self.call(self.fractional(), dual, expected=2)
            self.assertIsNone(result['certified_bound'])
            self.assertFalse(result['certifies_optimality'])

    def test_claimed_objective_mismatch_cannot_certify(self):
        result = self.call(self.fractional(), 'R 2/3\n', '=obj= 6\nx 1/3\n',
                           ['--format', 'rational', '--exact'], expected=1)
        self.assertFalse(result['certifies_optimality'])

    def test_highs_lp_primal_dual_and_basis_blocks(self):
        model = ('NAME HIGHS\nROWS\n N O\n G R\nCOLUMNS\n x O 2 R 1\n y O 3 R 1\n'
                 'RHS\n b O -1.25 R 1.5\nBOUNDS\n UP b x 4\n UP b y 4\nENDATA\n')
        primal = ('Model status\nOptimal\n\n# Primal solution values\nFeasible\nObjective 4.25\n'
                  '# Columns 2\nx 1.5\ny 0\n# Rows 1\nR 1.5\n')
        dual = ('\n# Dual solution values\nFeasible\n# Columns 2\nx 0\ny 1\n'
                '# Rows 1\nR 2\n')
        basis = ('\n# Basis\nHiGHS_basis_file v2\nValid\n# Columns 2\nx 1\ny 0\n'
                 '# Rows 1\nR 0\n')
        for suffix in [dual, basis, dual + basis]:
            result = self.call(model, solution=primal + suffix,
                               extra=['--format', 'highs', '--exact'])
            self.assertEqual(result['objective'], '17/4')
            self.assertTrue(result['checked_feasible'])
        # Later Feasible/Columns records cannot replace absent primal data.
        self.call(model, solution='Model status\nOptimal\n' + dual + basis,
                  extra=['--format', 'highs'], expected=2)
        self.call(model, solution=primal + primal,
                  extra=['--format', 'highs'], expected=2)

    def test_randomized_exact_support_formula(self):
        rng = random.Random(20261106)
        for _ in range(120):
            n, m = 3, 2
            a = [[rng.randint(-4, 4) for _ in range(n)] for _ in range(m)]
            c = [rng.randint(-4, 4) for _ in range(n)]
            x = [F(rng.randint(-4, 4), 2) for _ in range(n)]
            lows = [v - rng.randint(0, 3) for v in x]
            ups = [v + rng.randint(0, 3) for v in x]
            activity = [sum(F(v) * z for v, z in zip(row, x)) for row in a]
            rl = [v - rng.randint(0, 3) for v in activity]
            ru = [v + rng.randint(0, 3) for v in activity]
            y = [F(rng.randint(-5, 5), rng.randint(1, 5)) for _ in range(m)]
            maximize = rng.choice([False, True])
            sign = -1 if maximize else 1
            offset = rng.randint(-5, 5)
            residual = [sign * c[j] - sum(y[i] * a[i][j] for i in range(m)) for j in range(n)]
            lower = sign * offset + sum(v * (rl[i] if v >= 0 else ru[i]) for i, v in enumerate(y))
            lower += sum(v * (lows[j] if v >= 0 else ups[j]) for j, v in enumerate(residual))
            expected = lower * sign
            objective = offset + sum(v * z for v, z in zip(c, x))
            dec = lambda v: str(float(v))
            model = ('NAME RANDOM\nOBJSENSE\n ' + ('MAX' if maximize else 'MIN') +
                     '\nROWS\n N O\n' + ''.join(f' L r{i}\n' for i in range(m)) + 'COLUMNS\n' +
                     ''.join(f' x{j} O {c[j]}\n' + ''.join(f' x{j} r{i} {a[i][j]}\n' for i in range(m)) for j in range(n)) +
                     f'RHS\n b O {-offset}\n' + ''.join(f' b r{i} {dec(ru[i])}\n' for i in range(m)) +
                     'RANGES\n' + ''.join(f' g r{i} {dec(ru[i]-rl[i])}\n' for i in range(m)) + 'BOUNDS\n' +
                     ''.join(f' LO b x{j} {dec(lows[j])}\n UP b x{j} {dec(ups[j])}\n' for j in range(n)) + 'ENDATA\n')
            dual = ''.join(f'r{i} {v}\n' for i, v in enumerate(y))
            primal = '=obj= ' + str(objective) + '\n' + ''.join(f'x{j} {v}\n' for j, v in enumerate(x))
            result = self.call(model, dual, primal, ['--format', 'rational', '--exact', '--relax'])
            self.assertEqual(F(result['certified_bound']), expected)
            self.assertLessEqual(lower, sign * objective)
            self.assertEqual(result['certifies_optimality'], expected == objective)


unittest.main()
