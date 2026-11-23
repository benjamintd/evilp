"""Check the binary translation against independent enumeration and the original MPS."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import itertools
import json
import random
import shlex
import subprocess
import sys
import tempfile
import unittest

from test_lp_solver import decimal
from test_mip_solver import integer_mps, mip_oracle

SOLVER, CHECKER = map(lambda p: str(Path(p).resolve()), sys.argv[1:3])
del sys.argv[1:3]


def model(cost, rows, maximize=False):
    return dict(cost=cost, bounds=[(0, 1)] * len(cost), integer=[True] * len(cost),
                rows=rows, offset=F(15, 8), maximize=maximize)


class BinarySAT(unittest.TestCase):
    def solve(self, data, options=()):
        with tempfile.TemporaryDirectory() as tmp:
            path, point = (Path(tmp) / name for name in ('original.mps', 'point.sol'))
            path.write_text(integer_mps(data))
            cmd = [SOLVER, str(path), '--mip', '--binary-sat', '--solution', str(point), *options]
            run = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            self.assertIn(run.returncode, (0, 4), (run.stdout, run.stderr))
            self.assertNotIn('runtime error:', run.stderr)
            r = json.loads(run.stdout)
            value = None
            if point.exists():
                c = subprocess.run([CHECKER, '--model', str(path), '--solution', str(point)],
                                   capture_output=True, text=True, timeout=10)
                self.assertEqual(c.returncode, 0, (c.stdout, c.stderr))
                value = F(json.loads(c.stdout)['objective'])
            self.assertFalse(r['certifies_optimality'])
            return r, value

    def test_truth_tables_against_enumeration(self):
        rng = random.Random(8321)
        calls = 0
        for i in range(60):
            n = rng.randrange(4, 9)
            rows = []
            for j in range(rng.randrange(1, 8)):
                coeff = [F(rng.randint(-5, 5), (2, 4, 10)[i % 3]) for _ in range(n)]
                rhs = F(rng.randint(-6, 7), 4)
                lo, hi = [(None, rhs), (rhs, None), (rhs, rhs), (rhs - 2, rhs + 1)][j % 4]
                rows.append((coeff, lo, hi))
            m = model([F(rng.randrange(-7, 8), 4) for _ in range(n)], rows, bool(i % 2))
            expected = mip_oracle(m)
            with self.subTest(i=i):
                r, value = self.solve(m)
                self.assertEqual(r['status'], 'infeasible' if expected is None else 'optimal', r)
                self.assertEqual(value, expected)
                calls += r['sat_calls']
        self.assertGreater(calls, 15)

    def test_long_signed_at_most_one(self):
        n = 70
        # At most one selected, with a heavier literal forbidden outright.
        m = model([0] * (n-2) + [-2, -4], [([1] * (n-1) + [2], None, 1)])
        r, v = self.solve(m)
        self.assertGreater(r['sat_calls'], 0)
        self.assertEqual(v, F(15,8)-2)

    def test_nonrepresentable_cost_and_sub_double_boundary(self):
        delta = F(1,10**21)
        for m in [model([1 + delta, 1], [([1,1],1,None)]),
                  model([-1,-2], [([1+delta,1],None,1)]),
                  model([1,1], [([1,1],1+delta,None)])]:
            r, value = self.solve(m)
            self.assertEqual(value, mip_oracle(m))
            self.assertGreater(r['sat_calls'],0)

    def test_minimum_int64_coefficient_uses_exact_fallback(self):
        n = 24
        data = model([1] * n, [([1] * n, 5, None),
                               ([-(1 << 63), 1] + [0] * (n - 2), None, 0)])
        result, checked = self.solve(data)
        self.assertGreater(result['sat_calls'], 0)
        self.assertEqual((result['status'], checked), ('optimal', data['offset'] + 5))

    def test_large_relaxations_keep_the_full_witness(self):
        n = 540
        def row(entries, lo, hi):
            a = [0] * n
            for j, coeff in entries:
                a[j] = coeff
            return a, lo, hi
        rows = [row([(0,1),(1,1)],1,None), row([(0,1),(2,1)],1,None)]
        rows.extend(row([(j,1),(j+1,1)],1,1) for j in range(2,n-1))
        m = model([3,1] + [0]*(n-2),rows)
        r, value = self.solve(m)
        self.assertEqual((r['status'],value),('optimal',F(15,8)+1))
        self.assertGreater(r['sat_calls'],0)

    def test_many_prices_and_fixed_binary(self):
        n = 44
        costs = [(-1 if j % 2 else 2) for j in range(20)] + [0] * (n-20)
        # Signed and sub-double source costs must choose the correct face.
        costs[0] = F(1, 10**18)
        for maximize in (True, False):
            rows = [([0] * 20 + [1] * 24, 1, None)]
            m = model(costs, rows, maximize)
            m['bounds'][3] = (1, 1)
            r, value = self.solve(m)
            expected = m['offset'] + costs[3] + sum(
                c for j, c in enumerate(costs) if j != 3 and (c > 0 if maximize else c < 0))
            self.assertEqual((r['status'], value), ('optimal', expected))
            self.assertEqual(r['sat_calls'], 1)
            self.assertEqual(r['bound_origin'], 'objective_box')

    def test_relaxed_long_row_cannot_close_the_face(self):
        # The free variables are equal. Their sum cannot be 12 on the
        # cheapest face. Moving the first priced variable permits all zeros.
        n = 41
        rows = [([12] + [0] * 16 + [1] * 24, 12, 12)]
        for j in range(17, n-1):
            coeff = [0] * n
            coeff[j], coeff[j+1] = 1, -1
            rows.append((coeff, 0, 0))
        m = model([1] * 17 + [0] * 24, rows)
        r, value = self.solve(m)
        self.assertEqual((r['status'], value), ('optimal', m['offset'] + 1))

    def test_disjoint_and_overlapping_cover_bound(self):
        n = 64
        rows = []
        for start in range(0, n, 2):
            a = [0] * n
            for k in range(4):
                a[(start + k) % n] = 1
            rows.append((a, 1, None))
        r, value = self.solve(model([1] * n, rows))
        self.assertEqual((r['status'], value), ('optimal', F(15, 8) + n // 4))
        self.assertIn(r['sat_calls'], (1, 2))
        self.assertEqual(r['bound_origin'], 'row_penalty')

    def test_signed_row_charges_against_enumeration(self):
        # Three equal copies of six independent bits cross the PB objective
        # threshold while retaining a small exact oracle. The final two
        # variables have a fixed priced value and a free zero-cost value.
        rng = random.Random(34019)
        for case in range(16):
            cost = [F(rng.choice([-4, -2, 1, 3]), 2) for _ in range(18)] + [2, 0]
            planted = [rng.randrange(2) for _ in range(6)] * 3 + [1, 1]
            rows = []
            for i in range(8):
                a = [F(rng.randint(-3, 3), 10) for _ in cost]
                b = sum(x * c for x, c in zip(planted, a))
                lo = b - F(rng.randrange(4), 10) if i % 3 != 0 else None
                hi = b + F(rng.randrange(4), 10) if i % 3 != 1 else None
                rows.append((a, lo, hi))
            for j in range(6, 18):
                a = [0] * len(cost)
                a[j], a[j % 6] = 1, -1
                rows.append((a, 0, 0))
            m = model(cost, rows, maximize=bool(case % 2))
            m['bounds'][-2] = (1, 1)
            attainable = []
            for bits in itertools.product((0, 1), repeat=7):
                x = bits[:6] * 3 + (1, bits[6])
                acts = [sum(c * v for c, v in zip(a, x)) for a, _, _ in rows[:8]]
                if all((lo is None or lo <= y) and (hi is None or y <= hi)
                       for y, (_, lo, hi) in zip(acts, rows)):
                    attainable.append(m['offset'] + sum(c * v for c, v in zip(cost, x)))
            expected = (max if m['maximize'] else min)(attainable)
            with self.subTest(case=case):
                result, checked = self.solve(m)
                self.assertEqual((result['status'], checked), ('optimal', expected))

    def test_fractional_cost_per_unit_on_overlapping_weighted_rows(self):
        # Five units of penalty buy eight units of activity (5 + 3).
        # At the cheapest ratio 3/5, a block needs ceil(8 * 3/5) = 5.
        # The third row crosses both blocks and reuses their priced bits.
        gains = [5, 1, 3] * 10
        left = gains + [0] * 30
        right = [0] * 30 + gains
        overlap = gains[:15] + [0] * 30 + gains[15:]
        for maximize in (False, True):
            sign = -1 if maximize else 1
            cost = [sign * w for w in [3, 1, 2] * 20]
            data = model(cost, [(a, 8, None) for a in (left, right, overlap)], maximize)
            result, checked = self.solve(data)
            self.assertGreater(result['sat_calls'], 0)
            self.assertEqual((result['status'], checked), ('optimal', data['offset'] + sign * 10))

    def test_equality_bound_in_binary_and_bounded_integer_encodings(self):
        for upper in (1, 3):
            data = model([1] * 5, [([1] * 5, 2, 2)])
            data['bounds'] = [(0, upper)] * 5
            with self.subTest(upper=upper):
                r, value = self.solve(data)
                self.assertEqual(value, mip_oracle(data))
                self.assertEqual((r['search_route'], r['bound_origin'], r['best_bound']),
                                 ('binary_sat', 'source_equalities', value))

    def test_objective_cutoff_and_interrupted_cutoff_in_runner(self):
        # At most one of the five bits can be set. Two decisions leave an
        # incumbent and an open native frontier. Four decisions close via SAT.
        n = 5
        rows = [([int(k in (i, j)) for k in range(n)], None, 1)
                for i in range(n) for j in range(i)]
        data = model([-1] * n, rows)
        optimum = mip_oracle(data)
        for limit in (2, 4):
            with self.subTest(limit=limit):
                record, cli = self.run_in_runner(data, optimum, ['--node-limit', str(limit)])
                self.assertEqual(record['claimed_status'], cli['status'])
                for field in ('search_route', 'closure_reason', 'bound_origin'):
                    self.assertEqual(record[field], cli[field])
                self.assertTrue(record['checker']['checked_feasible'])
                self.assertEqual(F(record['checker']['objective']), optimum)
                self.assertFalse(record['checker']['certifies_optimality'])
                if limit == 2:
                    self.assertEqual((cli['status'], cli['search_route'], cli['closure_reason'],
                                      cli['bound_origin']), ('node_limit', 'native', 'none', 'open_nodes'))
                    self.assertLess(cli['best_bound'], optimum)
                else:
                    self.assertEqual((cli['status'], cli['search_route'], cli['closure_reason'],
                                      cli['bound_origin']), ('optimal', 'binary_sat', 'matching_bound', 'sat_cutoff_unsat'))
                    self.assertEqual(cli['best_bound'], optimum)

    def run_in_runner(self, data, optimum, options=()):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'model.mps'
            source.write_text(integer_mps(data))
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            suite = root / 'suite.tsv'
            suite.write_text('id\tpath\tsha256\treference_status\treference_objective\n'
                             f'case\t{source}\t{digest}\t=opt=\t{decimal(optimum)}\n')
            wrapper = root / 'solver'
            args = ' '.join(shlex.quote(arg) for arg in options)
            wrapper.write_text(f'#!/bin/sh\nexec {shlex.quote(SOLVER)} "$@" --binary-sat {args}\n')
            wrapper.chmod(0o755)
            sha = hashlib.sha256(wrapper.read_bytes()).hexdigest()
            solvers = root / 'solvers.tsv'
            solvers.write_text('solver\tbinary\tsha256\tversion\n'
                               f'lab\t{wrapper}\t{sha}\ttest\n')
            out = root / 'runs'
            subprocess.run([Path(SOLVER).with_name('milpbench'), '--suite', suite,
                            '--solvers', solvers, '--checker', CHECKER, '--mip',
                            '--out', out, '--limit', '10'], check=True, capture_output=True)
            record = json.loads((out / 'runs.jsonl').read_text())
            cli = json.loads((out / 'case/lab/seed-1/solver.log').read_text(), parse_float=F)
            return record, cli

    def test_closing_bound_encloses_exact_witness_after_serialization(self):
        for maximize in (False, True):
            for optimum in (F(21, 10), F(-21, 10), F(2**53 + 1), F(-2**53 - 1), F(25, 8)):
                data = model([1] * 5, [([1] * 5, 2, 2)], maximize=maximize)
                data['offset'] = optimum - 2
                with self.subTest(maximize=maximize, optimum=optimum):
                    record, cli = self.run_in_runner(data, optimum)
                    self.assertEqual((cli['status'], cli['bound_origin']), ('optimal', 'source_equalities'))
                    self.assertEqual(F(record['checker']['objective']), optimum)
                    self.assertTrue(record['checker']['checked_feasible'])
                    self.assertFalse(record['checker']['certifies_optimality'])
                    bound = F(record['reported_dual_bound'])
                    self.assertEqual(bound, cli['best_bound'])
                    for published in (bound, F.from_float(float(bound))):
                        if maximize:
                            self.assertGreaterEqual(published, optimum)
                        else:
                            self.assertLessEqual(published, optimum)
                    if optimum == F(25, 8):
                        self.assertEqual(bound, optimum)

    def test_short_budget(self):
        n = 7
        rows = []
        for i in range(n):
            for j in range(i):
                a = [0] * n
                a[i] = a[j] = 1
                rows.append((a,None,1))
        m = model([-1]*n,rows)
        r, value = self.solve(m, ['--node-limit','1'])
        if r['status'] == 'optimal':
            self.assertEqual(value,F(15,8)-1)
        else:
            self.assertIn(r['status'], ('node_limit','time_limit'))
        self.assertLessEqual(r['sat_calls'],1)

    def test_large_fixed_activities_cancel_in_each_encoding(self):
        # Ten positive products exceed INT64_MAX before the other ten cancel.
        k = 10**9
        fixed = [k] * 10 + [-k] * 10
        for upper, free in [(1, [1, 1, 0, 0, 0]), (2, [1, 1, 0, 0, 0]),
                            (2, [1] * 5), (35, [-1, 1, 0, 0, 0])]:
            for sign in (1, -1):
                coeff = [sign * c for c in fixed + free]
                data = model([0] * 20 + [1, 1, 0, 0, 0], [(coeff, sign, sign)])
                data['bounds'] = [(k, k)] * 20 + [(0, upper)] + [(0, min(upper, 2))] * 4
                with self.subTest(upper=upper, free=free, sign=sign):
                    result, point = self.solve(data)
                    self.assertEqual(point, mip_oracle(data))
                    self.assertEqual((result['status'], result['search_route']), ('optimal', 'binary_sat'))
                    self.assertEqual(result['best_bound'], point)

    def test_fixed_activity_beyond_int64(self):
        k = 10**9
        # RHS minus the fixed activity lies outside int64_t. The two row
        # directions are respectively redundant and impossible.
        for sign in (-1, 1):
            coeff = [sign * k] * 10 + [1, -1]
            data = model([0] * 10 + [1, 1], [(coeff, None, 0), ([0] * 10 + [1, 1], 1, None)])
            data['bounds'] = [(k, k)] * 10 + [(0, 2)] * 2
            with self.subTest(sign=sign):
                result, value = self.solve(data)
                optimum = mip_oracle(data)
                self.assertEqual(value, optimum)
                self.assertEqual(result['status'], 'infeasible' if optimum is None else 'optimal')

    def test_translated_ranged_row_against_enumeration(self):
        # A signed ranged row selected from the differential oracle. Both
        # forms have the same 510 possible points and identical row slacks.
        coeff = [-4, -4, 3]
        costs = [-4, 6, 6]
        origin = [-10**9, 10**9 - 40, 10**9 - 40]
        rhs = 2999999933
        for shift in (origin, [0, 0, 0]):
            translation = [a - b for a, b in zip(origin, shift)]
            side = rhs - sum(a * b for a, b in zip(coeff, translation))
            data = model(costs, [(coeff, side - 4, side)], maximize=True)
            data['bounds'] = [(lo, lo + width) for lo, width in zip(shift, [33, 4, 2])]
            data['offset'] += sum(a * b for a, b in zip(costs, translation))
            with self.subTest(shift=shift):
                r, value = self.solve(data)
                self.assertEqual(value, mip_oracle(data))
                self.assertEqual((r['search_route'], r['best_bound']), ('binary_sat', value))

    def test_source_rhs_beyond_machine_integer(self):
        # Both sides of this decimal equality round to the fixed activity in
        # double precision. The residual quarter-unit determines the optimum.
        k = 10**9
        for sign in (1, -1):
            row = [sign * (k // 4)] * 40 + [sign * F(1, 4), sign * F(1, 2)]
            rhs = sign * (10 * k * k + F(1, 4))
            data = model([0] * 40 + [1, 1], [(row, rhs, rhs)])
            data['bounds'] = [(k, k)] * 40 + [(0, 2)] * 2
            with self.subTest(sign=sign):
                r, value = self.solve(data)
                self.assertEqual(value, mip_oracle(data))
                self.assertEqual((r['search_route'], r['best_bound']), ('binary_sat', value))

    def test_no_backend_on_continuous_columns(self):
        m = model([1, 1], [([1, 1], F(1,2), None)])
        m['integer'][0] = False
        r, v = self.solve(m)
        self.assertEqual(r['sat_calls'], 0)
        self.assertEqual(v, F(15,8) + F(1,2))

    def test_per_row_lcm_rational_scaling_sat_and_unsat(self):
        # 20 binary vars with 0.333 (= 333/1000) coefficients; gcd(333)=333 => sum <= floor(1000/333)=3.
        n = 20
        w = [F(333, 1000)] * n
        m_sat = model([-1] * 10 + [0] * 10, [(w, None, 1)])
        r, v = self.solve(m_sat)
        self.assertEqual(r['status'], 'optimal', r)
        self.assertEqual(v, F(15, 8) - 3)
        self.assertGreater(r['sat_calls'], 0)
        # Infeasible: 0.333 * sum x <= 1 (max 3) AND sum x >= 4
        m_unsat = model([1] + [0] * (n - 1), [(w, None, 1), ([1] * n, 4, None)])
        r2, _ = self.solve(m_unsat)
        self.assertEqual(r2['status'], 'infeasible', r2)
        self.assertGreater(r2['sat_calls'], 0)

    def test_cheapest_face_general_integer_penalty_slacks(self):
        # 12 binary vars + 4 unbounded-above nonnegative integer slacks summed into 1 total var.
        # Setting total = 0 propagates all 4 slacks to 0 and solves binary face in Kissat.
        n = 17
        costs = [0] * 16 + [1]
        bounds = [(0, 1)] * 12 + [(0, None)] * 5
        integer = [True] * n
        rows = [
            ([0] * 12 + [1, 1, 1, 1, -1], 0, 0),
            ([1] * 6 + [0] * 6 + [1, 0, 0, 0, 0], 3, None),
            ([0] * 6 + [1] * 6 + [0, 1, 0, 0, 0], 3, None),
            ([1] * 12 + [0] * 5, None, 6),
        ]
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
        r, v = self.solve(m)
        self.assertEqual(r['status'], 'optimal', r)
        self.assertEqual(v, F(0))
        self.assertGreater(r['sat_calls'], 0)

    def test_pb_adder_large_weighted_constraint(self):
        # 36 binary variables with weighted constraint capacity > 512 (triggers pb_adder_at_most)
        # and > 16 objective terms (triggers cached-CNF pb_mode binary search).
        n = 36
        weights = [37 + 13 * j for j in range(n)]
        cap = 1500  # > 512, gcd(weights) = 1 -> forces pb_adder_at_most
        costs = [-(1 + (j % 3)) for j in range(n)]
        rows = [
            (weights, None, cap),
            ([1 if j < 18 else -2 for j in range(n)], None, 3),
            ([1 if j % 2 == 0 else -1 for j in range(n)], -2, 2),
        ]
        m = model(costs, rows, maximize=False)
        r, v = self.solve(m)
        self.assertEqual(r['status'], 'optimal', r)
        self.assertGreater(r['sat_calls'], 0)
        r_ref, v_ref = self.solve(m, ['--no-binary-sat'])
        self.assertEqual(r_ref['status'], 'optimal', r_ref)
        self.assertEqual(v, v_ref)

    def test_small_integer_pb_optimization_and_wide_binary_encoding(self):
        # Pure-integer model mixing small-domain (order-encoded) and wide-domain (binary-encoded,
        # width = 75 > 32) variables with a non-constant objective.
        costs = [-3, -5, -1, 2, -2, -4]
        bounds = [(0, 6), (0, 5), (0, 75), (0, 4), (0, 1), (0, 1)]
        integer = [True] * 6
        rows = [
            ([2, 3, 1, -2, 4, 5], None, 55),
            ([1, 1, 0, 1, 1, 1], 3, 10),
            ([0, 0, 1, -10, 0, 0], None, 40),
        ]
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=F(3, 2), maximize=False)
        r, v = self.solve(m)
        self.assertEqual(r['status'], 'optimal', r)
        self.assertGreater(r['sat_calls'], 0)
        # Verify exact optimality against native B&B (--no-binary-sat --no-scip)
        r_nat, v_nat = self.solve(m, ['--no-binary-sat', '--no-scip'])
        self.assertEqual(r_nat['status'], 'optimal', r_nat)
        self.assertEqual(v, v_nat)

    def test_pb_totalizer_unit_cardinality_and_descent(self):
        # 36 binary vars with unit-weight cardinality rows (n >= 16, cap >= 2 -> GAC Totalizer)
        # and unit objective (> 16 terms -> pb_mode adaptive linear-exponential descent).
        n = 36
        costs = [1] * n
        rows = [
            ([1] * n, 14, 24),
            ([1 if j < 24 else 0 for j in range(n)], None, 9),
            ([1 if j >= 12 else 0 for j in range(n)], None, 10),
            ([1 if j % 2 == 0 else -1 for j in range(n)], -1, 1),
            ([1 if j % 3 == 0 else -1 if j % 3 == 1 else 0 for j in range(n)], 0, 2),
        ]
        m = model(costs, rows, maximize=False)
        r, v = self.solve(m)
        self.assertEqual(r['status'], 'optimal', r)
        self.assertGreater(r['sat_calls'], 0)
        r_ref, v_ref = self.solve(m, ['--no-binary-sat'])
        self.assertEqual(r_ref['status'], 'optimal', r_ref)
        self.assertEqual(v, v_ref)


if __name__ == '__main__':
    unittest.main()
