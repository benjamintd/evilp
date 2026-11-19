"""Enumerate small identical-bin assignments and check original integer witnesses."""
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction

from test_mip_solver import integer_mps

SOLVER, CHECKER = map(lambda s: str(pathlib.Path(s).resolve()), sys.argv[1:3])
del sys.argv[1:3]


def packing(demand, weights, capacity, bins, maximize=False):
    types = len(demand)
    n = bins * (types + 1)
    rows = [([0] * n, d, d) for d in demand] + [([0] * n, None, 0) for _ in range(bins)]
    bounds = [(0, 1)] * bins
    for b in range(bins):
        rows[types + b][0][b] = -capacity
        for t in range(types):
            j = bins + b * types + t
            rows[t][0][j] = 1
            rows[types + b][0][j] = weights[t]
            bounds.append((0, max(1, capacity // weights[t])))
    charge = Fraction(-3 if maximize else 3, 4)
    return {'rows': rows, 'cost': [charge] * bins + [0] * (n - bins),
            'bounds': bounds, 'integer': [True] * n, 'offset': Fraction(1, 2),
            'maximize': maximize}


def minimum_bins(demand, weights, capacity):
    # Independent finite-state enumeration over individual item placements.
    items = [w for d, w in zip(demand, weights) for _ in range(d)]
    loads = {()}
    for size in items:
        nxt = set()
        for state in loads:
            nxt.add(tuple(sorted((*state, size))))
            for j, used in enumerate(state):
                if used + size <= capacity:
                    new = list(state)
                    new[j] += size
                    nxt.add(tuple(sorted(new)))
        loads = nxt
    return min(map(len, loads))


class Patterns(unittest.TestCase):
    def solve(self, text, *flags):
        with tempfile.TemporaryDirectory() as td:
            path = pathlib.Path(td) / 'packing.mps'
            sol = pathlib.Path(td) / 'point.sol'
            path.write_text(text)
            r = subprocess.run([SOLVER, path, '--mip', '--solution', sol, *flags],
                               capture_output=True, text=True, timeout=15)
            data = json.loads(r.stdout)
            self.assertNotIn('runtime error:', r.stderr)
            value = None
            if sol.exists():
                check = subprocess.run([CHECKER, '--model', path, '--solution', sol, '--exact'],
                                       capture_output=True, text=True, timeout=10)
                self.assertEqual(check.returncode, 0, check.stdout)
                witness = json.loads(check.stdout)
                self.assertTrue(witness['checked_feasible'])
                self.assertFalse(witness['certifies_optimality'])
                value = Fraction(witness['objective'])
            return data, value

    def test_enumeration_and_both_objective_senses(self):
        rng = random.Random(67027)
        for case in range(24):
            demand = [rng.randint(1, 3) for _ in range(2 + case % 2)]
            weights = [rng.randint(1, 5) for _ in demand]
            cap = rng.randint(5, 9)
            expected = minimum_bins(demand, weights, cap)
            m = packing(demand, weights, cap, sum(demand), bool(case % 2))
            data, obj = self.solve(integer_mps(m))
            self.assertEqual(data['status'], 'optimal', data)
            self.assertIn('pattern lower bound', data['message'], data)
            self.assertEqual(obj, m['offset'] + expected * m['cost'][0])

    def test_nonuniform_and_budget_fall_back(self):
        m = packing([2, 3], [3, 2], 5, 4)
        original = integer_mps(m)
        for altered in [original.replace('x0 OBJ 0.75', 'x0 OBJ 0.750000000000000001'),
                        original.replace('x4 r2 3', 'x4 r2 3.000000000000000001')]:
            self.assertNotEqual(original, altered)
            result, _ = self.solve(altered, '--node-limit', '30')
            self.assertNotIn('pattern lower bound', result['message'])
        result, _ = self.solve(original, '--node-limit', '1')
        self.assertNotIn('pattern lower bound', result['message'])
        self.assertLessEqual(result['nodes'], 1)

    def test_insufficient_bin_supply(self):
        result, obj = self.solve(integer_mps(packing([3, 2], [4, 3], 5, 2)))
        self.assertEqual(result['status'], 'infeasible')
        self.assertIsNone(obj)


if __name__ == '__main__':
    unittest.main()
