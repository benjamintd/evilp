#!/usr/bin/env python3
"""Independent exact original-model Farkas certificate contracts."""
import fractions
import gzip
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest

BIN = str(pathlib.Path(sys.argv.pop(1)).resolve())
F = fractions.Fraction
HEADER = 'MILPLAB_FARKAS_1\n'


class FarkasCertificates(unittest.TestCase):
    def call(self, model, certificate=HEADER, extra=(), expected=0, compressed=False):
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            model_path = root / ('model.mps.gz' if compressed else 'model.mps')
            if compressed:
                model_path.write_bytes(gzip.compress(model.encode()))
            else:
                model_path.write_text(model)
            cert_path = root / 'farkas.txt'
            cert_path.write_bytes(certificate if isinstance(certificate, bytes)
                                  else certificate.encode())
            args = [BIN, '--model', str(model_path), '--farkas', str(cert_path)]
            p = subprocess.run(args + list(extra), text=True, capture_output=True)
            self.assertEqual(p.returncode, expected, (p.stdout, p.stderr))
            if expected == 64:
                return None
            result = json.loads(p.stdout)
            self.assertFalse(result['checked_feasible'])
            self.assertFalse(result['certifies_optimality'])
            self.assertEqual(result['bound_status'], 'absent')
            self.assertIsNone(result['certified_bound'])
            self.assertIsNone(result['objective'])
            self.assertEqual(result['certifies_infeasibility'], expected == 0)
            self.assertEqual(result['certificate_status'],
                             {0: 'valid', 1: 'invalid', 2: 'error',
                              3: 'unsupported', 4: 'unavailable'}[expected])
            if expected in (2, 3, 4):
                self.assertIsNone(result['contradiction_margin'])
            return result

    @staticmethod
    def single(sense='G', rhs='2', bounds=' UP b x 1\n', coefficient='1'):
        return (f'NAME SINGLE\nROWS\n N O\n {sense} R\nCOLUMNS\n'
                f' x O 123 R {coefficient}\nRHS\n b R {rhs}\n'
                f'BOUNDS\n{bounds}ENDATA\n')

    def test_lower_row_upper_column_and_objective_independence(self):
        for maximize in (False, True):
            model = self.single().replace('ROWS\n', 'OBJSENSE\n ' +
                                           ('MAX' if maximize else 'MIN') + '\nROWS\n')
            model = model.replace(' b R 2', ' b O 999999999999999999999999 R 2')
            result = self.call(model, HEADER + 'R +1\n', compressed=maximize)
            self.assertEqual(result['certificate_kind'], 'row_combination')
            self.assertEqual(result['contradiction_margin'], '1')
            self.assertEqual(result['objective_sense'], 'max' if maximize else 'min')

    def test_upper_row_lower_column(self):
        result = self.call(self.single('L', '-1', ''), HEADER + 'R -2/3\n')
        self.assertEqual(result['contradiction_margin'], '2/3')

    def test_fixed_and_one_sided_original_columns(self):
        cases = [
            (self.single('G', '-3', ' FX b x -2\n', '2'), 'R 1\n'),
            (self.single('G', '0', ' MI b x\n UP b x -1\n'), 'R 1\n'),
            (self.single('L', '1', ' LO b x 2\n PL b x\n'), 'R -1\n'),
        ]
        for model, ray in cases:
            self.assertEqual(self.call(model, HEADER + ray)['contradiction_margin'], '1')

    def test_all_signed_range_sides(self):
        for sense, ran, lo, hi in [('L', -2, 3, 5), ('L', 2, 3, 5),
                                    ('G', -2, 5, 7), ('G', 2, 5, 7),
                                    ('E', 2, 5, 7), ('E', -2, 3, 5), ('E', 0, 5, 5)]:
            for sign in (1, -1):
                bounds = (f' UP b x {lo-1}\n' if sign > 0 else f' LO b x {hi+1}\n')
                model = self.single(sense, '5', bounds)
                model = model.replace('BOUNDS\n', f'RANGES\n g R {ran}\nBOUNDS\n')
                result = self.call(model, HEADER + f'R {sign}\n')
                self.assertEqual(result['contradiction_margin'], '1')

    @staticmethod
    def equalities():
        return ('NAME CANCEL\nOBJNAME COST\nROWS\n N IGNORED\n N COST\n E A\n E B\n'
                'COLUMNS\n x COST 999999 A 3\n x B 9 IGNORED 10000\n'
                ' z COST -1 IGNORED 99\nRHS\n b A 6 B 9\n'
                ' b COST -100 IGNORED 500\nBOUNDS\n FR b x\n FR b z\nENDATA\n')

    def test_equalities_exact_zero_on_free_columns(self):
        result = self.call(self.equalities(), HEADER + 'A 1/3\nB -1/9\n')
        self.assertEqual(result['contradiction_margin'], '1')
        # An opposite combination has a finite negative margin, not a proof.
        result = self.call(self.equalities(), HEADER + 'A -1/3\nB 1/9\n', expected=1)
        self.assertEqual(result['contradiction_margin'], '-1')

    def test_nonzero_free_column_residual_never_rounded(self):
        for perturbation in (F(1, 10**100), -F(1, 10**100)):
            ray = HEADER + f'A 1/3\nB {F(-1, 9) + perturbation}\n'
            result = self.call(self.equalities(), ray,
                               ['--lin-tol', '1e100', '--int-tol', '1e100',
                                '--obj-tol', '1e100'], expected=4)
            self.assertIn('column x', result['certificate_message'])

    def test_infinite_selected_row_or_column_side(self):
        for model, ray, message in [
            (self.single(), 'R -1\n', 'row R'),
            (self.single('L', '-1'), 'R 1\n', 'row R'),
            (self.single(bounds=' PL b x\n'), 'R 1\n', 'column x'),
            (self.single('L', '-1', ' MI b x\n UP b x 1\n'), 'R -1\n', 'column x'),
        ]:
            result = self.call(model, HEADER + ray, expected=4)
            self.assertIn(message, result['certificate_message'])

    def test_zero_and_negative_margins_ignore_primal_tolerances(self):
        for rhs, margin in [('1', '0'), ('0.999999999999999999999999999999',
                                         '-1/1000000000000000000000000000000')]:
            for extra in ((), ['--lin-tol', '1e100', '--int-tol', '1e100',
                               '--obj-tol', '1e100', '--relax'], ['--exact']):
                result = self.call(self.single(rhs=rhs), HEADER + 'R 1\n', extra, expected=1)
                self.assertEqual(result['contradiction_margin'], margin)
                self.assertEqual(result['linear_tolerance'], '0')
                self.assertEqual(result['integer_tolerance'], '0')
                self.assertEqual(result['objective_tolerance'], '0')

    def test_any_positive_exact_margin_is_a_proof(self):
        result = self.call(self.single(rhs='1.000000000000000000000000000001'),
                           HEADER + 'R 1\n')
        self.assertEqual(result['contradiction_margin'],
                         '1/1000000000000000000000000000000')

    def test_constant_rows_with_zero_combined_coefficients(self):
        for sense, rhs, multiplier in [('G', '1', '1'), ('L', '-1', '-1'),
                                        ('E', '1', '1'), ('E', '-1', '-1')]:
            model = self.single(sense, rhs, ' FR b x\n', '0')
            result = self.call(model, HEADER + f'R {multiplier}\n')
            self.assertEqual(result['contradiction_margin'], '1')

    def test_zero_ray_and_omitted_rows(self):
        for ray in (HEADER, HEADER + 'R 0\n'):
            self.assertEqual(self.call(self.single(), ray, expected=1)['contradiction_margin'], '0')
        model = self.single().replace(' G R\n', ' G R\n G UNUSED\n')
        self.assertEqual(self.call(model, HEADER + 'R 1\n')['contradiction_margin'], '1')

    def test_empty_original_column_box_is_explicit_proof(self):
        model = self.single(bounds=' LO b x 3\n UP b x 2\n')
        result = self.call(model)
        self.assertEqual(result['certificate_kind'], 'column_bounds')
        self.assertEqual(result['contradiction_margin'], '1')
        self.assertIn('column x', result['certificate_message'])
        # Contradictory model bounds cannot excuse an ill-formed certificate.
        self.call(model, HEADER + 'UNKNOWN 1\n', expected=2)

    def test_integer_rounding_is_not_an_lp_infeasibility_proof(self):
        model = ('NAME INTEGER\nROWS\n N O\nCOLUMNS\n x O 1\nBOUNDS\n'
                 ' LI b x 0.1\n UI b x 0.9\nENDATA\n')
        result = self.call(model, expected=1)
        self.assertEqual(result['contradiction_margin'], '0')

    def test_binary_domain_preserves_original_box(self):
        model = self.single(rhs='1.5', bounds=' BV b x\n')
        self.assertEqual(self.call(model, HEADER + 'R 1\n')['contradiction_margin'], '1/2')

    def test_duplicate_coefficients_sum_exactly(self):
        model = self.equalities().replace(' x B 9 IGNORED 10000',
                                          ' x B 4 IGNORED 10000\n x B 5')
        result = self.call(model, HEADER + 'A 1/3\nB -1/9\n')
        self.assertEqual(result['contradiction_margin'], '1')

    def test_header_comments_and_decimal_scientific_multipliers(self):
        for certificate in ['\n # proof\n$ comment\n' + HEADER + 'R +1.0D+0 $ row\n',
                            HEADER.rstrip() + ' $ version\nR 1', HEADER + 'R 2/2\r\n']:
            self.assertEqual(self.call(self.single(), certificate)['contradiction_margin'], '1')

    def test_malformed_certificates_fail_closed(self):
        for certificate in ['', '# no header\n', 'R 1\n', 'MILPLAB_FARKAS_2\nR 1\n',
                            'MILPLAB_FARKAS_1 extra\nR 1\n', HEADER + HEADER,
                            HEADER + 'R 1\nR 1\n', HEADER + 'UNKNOWN 1\n',
                            HEADER + 'O 1\n', HEADER + 'R\n', HEADER + 'R 1 extra\n',
                            HEADER + 'R nan\n', HEADER + 'R Inf\n', HEADER + 'R 1/0\n',
                            HEADER + 'R 1/2/3\n', HEADER + 'R -/2\n', HEADER + 'R 2/-3\n',
                            HEADER + 'R 2.1/3\n', HEADER + 'R ' + '1' * 4097 + '\n',
                            HEADER + '#' + 'x' * 16384 + '\n',
                            (HEADER + 'R 1\x00trailing\n').encode()]:
            self.call(self.single(), certificate, expected=2)

    def test_farkas_mode_is_exclusive(self):
        for extra in (['--inspect'], ['--dual', 'unused'], ['--solution', 'unused']):
            self.call(self.single(), extra=extra, expected=64)

    def test_model_parse_failures_remain_errors(self):
        self.call(self.single().replace('ENDATA\n', ''), HEADER + 'R 1\n', expected=2)
        self.call(self.single().replace(' UP b x 1', ' SC b x 1'),
                  HEADER + 'R 1\n', expected=3)

    def test_randomized_exact_support_oracle(self):
        rng = random.Random(20261107)
        for _ in range(180):
            n, m = 4, 3
            a = [[rng.randint(-5, 5) for _ in range(n)] for _ in range(m)]
            lower = [F(rng.randint(-8, 0), 4) for _ in range(n)]
            upper = [v + F(rng.randint(0, 12), 4) for v in lower]
            rl = [F(rng.randint(-12, 12), 4) for _ in range(m)]
            ru = [v + F(rng.randint(0, 12), 4) for v in rl]
            y = [F(rng.randint(-8, 8), rng.randint(1, 7)) for _ in range(m)]
            q = [sum(y[i] * a[i][j] for i in range(m)) for j in range(n)]
            row_support = sum(v * (rl[i] if v >= 0 else ru[i]) for i, v in enumerate(y))
            box_support = sum(v * (upper[j] if v >= 0 else lower[j]) for j, v in enumerate(q))
            expected = row_support - box_support
            decimal = lambda v: str(float(v))  # All model denominators divide four.
            model = ('NAME RANDOM\nOBJSENSE\n ' + rng.choice(['MIN', 'MAX']) +
                     '\nROWS\n N O\n' + ''.join(f' L r{i}\n' for i in range(m)) +
                     'COLUMNS\n' + ''.join(f' x{j} O {rng.randint(-5, 5)}\n' +
                         ''.join(f' x{j} r{i} {a[i][j]}\n' for i in range(m)) for j in range(n)) +
                     'RHS\n b O 42\n' + ''.join(f' b r{i} {decimal(ru[i])}\n' for i in range(m)) +
                     'RANGES\n' + ''.join(f' g r{i} {decimal(ru[i]-rl[i])}\n' for i in range(m)) +
                     'BOUNDS\n' + ''.join(f' LO b x{j} {decimal(lower[j])}\n'
                     f' UP b x{j} {decimal(upper[j])}\n' for j in range(n)) + 'ENDATA\n')
            ray = HEADER + ''.join(f'r{i} {v}\n' for i, v in enumerate(y))
            result = self.call(model, ray, expected=0 if expected > 0 else 1)
            self.assertEqual(F(result['contradiction_margin']), expected)


unittest.main()
