#!/usr/bin/env python3
"""Source pricing: actual optimization, source witnesses, state reuse and caps.

An expected numerical resource abstention is a control-flow test only. It is
never counted as a completed optimization problem or independent proof.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from test_source_accuracy import case_text


def is_prime(n):
    """Deterministic Miller-Rabin oracle for the generated 64-bit integers."""
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        x = pow(a % n, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def rational_growth_model(count, finite=True):
    # At the binary optimum X_i=1,Z=0, the exact row multipliers are 1/p_i.
    # Pricing Z sums pairwise-prime denominators. The 150-row control stays
    # below8192bits;190 rows exceed it. T=U is the material source-cost trigger.
    primes, p = [], 10**15 + 1
    while len(primes) < count:
        if is_prime(p):
            primes.append(p)
        p += 2
    text = ['NAME rational_growth', 'OBJSENSE', ' MAX', 'ROWS', ' N OBJ', ' E LINK']
    text += [f' L R{i}' for i in range(count)]
    text += ['COLUMNS', ' T OBJ 1.000000000000000001 LINK 1', ' U OBJ -1 LINK -1']
    text += [f' X{i} OBJ 1 R{i} {p}' for i, p in enumerate(primes)]
    text += [' Z OBJ 0'] + [f' Z R{i} 1' for i in range(count)]
    text += ['RHS', ' rhs LINK 0'] + [f' rhs R{i} {p}' for i, p in enumerate(primes)]
    text += ['BOUNDS', ' UP b T 1000000000000000', ' UP b U 1000000000000000',
             ' UP b Z 1']
    if finite:
        text += [f' UP b X{i} 2' for i in range(count)]
    text += ['ENDATA', '']
    denominator_bits = sum((F(1, p) for p in primes), F(0)).denominator.bit_length()
    return '\n'.join(text), denominator_bits


def signed_model(source, maximum, negative, offset):
    low, high = ('-1000000000000000', '0') if negative else ('0', '1000000000000000')
    cost = '1.000000000000000001' if source else '1.5'
    constant = F('1000000.000000000000000001') if offset else F(0)
    text = ['NAME signed_source', 'OBJSENSE', ' MAX' if maximum else ' MIN',
            'ROWS', ' N OBJ', ' E LINK']
    text += [f' L PAD{i}' for i in range(199)]
    text += ['COLUMNS', " M0 'MARKER' 'INTORG'", f' X OBJ {cost} LINK 1',
             ' Y OBJ -1 LINK -1', " M1 'MARKER' 'INTEND'"]
    text += [f' P{i} OBJ 0' for i in range(298)] + ['RHS', ' rhs LINK 0']
    if offset:
        text += [' rhs OBJ -1000000.000000000000000001']
    text += ['BOUNDS', f' LO b X {low}', f' UP b X {high}',
             f' LO b Y {low}', f' UP b Y {high}']
    text += [f' FX b P{i} 0' for i in range(298)] + ['ENDATA', '']
    target = F(high if maximum else low)
    return '\n'.join(text), target, (F(cost) - 1) * target + constant


class SourcePricing(unittest.TestCase):
    def checked(self, model, solution, relax):
        command = [str(CHECKER), '--model', str(model), '--solution', str(solution)]
        checked = subprocess.run(command + (['--relax'] if relax else []),
                                 text=True, capture_output=True, timeout=30)
        self.assertEqual(checked.returncode, 0, (checked.stdout, checked.stderr))
        result = json.loads(checked.stdout)
        self.assertTrue(result['checked_feasible'])
        self.assertFalse(result['certifies_optimality'])
        return F(result['objective'])

    def solve(self, root, text, mode='--relax', flags=()):
        root.mkdir()
        model, solution = root / 'model.mps', root / 'point.sol'
        model.write_text(text)
        command = [str(SOLVER), str(model), mode, '--time-limit', '10',
                   '--solution', str(solution), *flags]
        run = subprocess.run(command, text=True, capture_output=True, timeout=30)
        return run, json.loads(run.stdout), model, solution

    def test_session_bounds_restore_and_original_model_ownership(self):
        text, _ = case_text('coupled', large=True)
        with tempfile.TemporaryDirectory(prefix='source-pricing-session-') as temporary:
            root = Path(temporary)
            model = root / 'model.mps'
            model.write_text(text)
            run = subprocess.run([str(SESSION), str(model), str(root)],
                                 text=True, capture_output=True, timeout=60)
            self.assertEqual(run.returncode, 0, (run.stdout, run.stderr))
            records = [json.loads(line) for line in run.stdout.splitlines()]
            self.assertEqual(len(records), 40)
            self.assertGreaterEqual(sum(record['warm'] for record in records), 8)
            solutions = list(root.glob('*.sol'))
            self.assertEqual(len(solutions), 40)
            for solution in solutions:
                point = dict(line.split() for line in solution.read_text().splitlines()
                             if line.strip() and not line.startswith('#'))
                target = F(point['X'])
                self.assertEqual(F(point['Y']), target)
                self.assertIn(target, (F(10**15), F(5 * 10**14), F(25 * 10**13)))
                self.assertEqual(self.checked(model, solution, True), target / 10**18)

    def test_sign_sense_reference_constants_and_point_fast_path(self):
        with tempfile.TemporaryDirectory(prefix='source-pricing-signed-') as temporary:
            root = Path(temporary)
            for source in (False, True):
                for maximum in (False, True):
                    for negative in (False, True):
                        for offset in (False, True):
                            text, target, expected = signed_model(source, maximum, negative, offset)
                            for mode in ('--relax', '--mip'):
                                name = f'{source}-{maximum}-{negative}-{offset}-{mode[2:]}'
                                flags = ['--no-presolve'] if mode == '--mip' else []
                                with self.subTest(case=name):
                                    run, result, model, solution = self.solve(root / name, text,
                                                                             mode, flags)
                                    self.assertEqual(run.returncode, 0, (run.stdout, run.stderr))
                                    self.assertEqual(result['status'], 'optimal')
                                    point = dict(line.split() for line in solution.read_text().splitlines()
                                                 if line.strip() and not line.startswith('#'))
                                    self.assertEqual(F(point['X']), target)
                                    self.assertEqual(F(point['Y']), target)
                                    self.assertEqual(self.checked(model, solution, mode == '--relax'),
                                                     expected)
                                    if source:
                                        self.assertGreater(result['source_support_scans'], 0)
                                    elif not offset:
                                        if result.get('glpk_basis_status') == 'source_verified_point_and_bound':
                                            # The direct proposal has its own rigorous support scan.
                                            # Exercise the native binary pricing path on the same input.
                                            self.assertGreater(result['source_support_scans'], 0)
                                            run, result, model, solution = self.solve(
                                                root / (name + '-basis'), text, mode, ['--glpk-basis'])
                                            self.assertEqual(run.returncode, 0, (run.stdout, run.stderr))
                                            self.assertEqual(result['status'], 'optimal')
                                            self.assertEqual(self.checked(model, solution, True), expected)
                                        self.assertEqual(result['source_support_scans'], 0)
                                        self.assertGreater(result['source_support_binary_paths'], 0)

    def test_rational_growth_completes_below_cap_and_abstains_above(self):
        with tempfile.TemporaryDirectory(prefix='source-pricing-bits-') as temporary:
            for count, finite in ((150, False), (150, True), (190, False), (190, True)):
                text, bits = rational_growth_model(count, finite)
                root = Path(temporary) / f'{count}-{finite}'
                run, result, model, solution = self.solve(root, text)
                self.assertGreater(result['source_pricing_attempts'], 0)
                if count == 150:
                    self.assertLessEqual(bits, 8192)
                    self.assertEqual(run.returncode, 0, (run.stdout, run.stderr))
                    self.assertEqual(result['status'], 'optimal')
                    self.assertEqual(self.checked(model, solution, True), F(count) + F(1, 1000))
                    if not finite:
                        self.assertGreater(result['source_dual_closures'], 0)
                else:
                    self.assertGreater(bits, 8192)
                    self.assertEqual(run.returncode, 4, (run.stdout, run.stderr))
                    self.assertEqual(result['status'], 'numerical_failure')
                    self.assertEqual(result['message'], 'Bounded source objective pricing unavailable')
                    self.assertFalse(result['has_primal'])
                    self.assertFalse(result['certifies_optimality'])
                    self.assertFalse(solution.exists())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--checker', type=Path, required=True)
    parser.add_argument('--session-driver', type=Path, required=True)
    args, remaining = parser.parse_known_args()
    SOLVER, CHECKER, SESSION = (path.resolve() for path in
                                (args.solver, args.checker, args.session_driver))
    unittest.main(argv=[__file__, *remaining])
