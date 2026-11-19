"""Check the binary translation against independent enumeration and the original MPS."""
from fractions import Fraction as F
from pathlib import Path
import json
import random
import subprocess
import sys
import tempfile
import unittest

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
        # Verify exact optimality against native B&B (--no-binary-sat --no-root-diving)
        r_nat, v_nat = self.solve(m, ['--no-binary-sat', '--no-root-diving'])
        self.assertEqual(r_nat['status'], 'optimal', r_nat)
        self.assertEqual(v, v_nat)


if __name__ == '__main__':
    unittest.main()
