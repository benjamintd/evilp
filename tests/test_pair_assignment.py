"""Compare recovered product costs with independent permutation enumeration."""
import itertools
import json
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from test_mip_solver import integer_mps, mip_oracle


SOLVER, CHECKER = [str(Path(arg).resolve()) for arg in sys.argv[1:3]]
del sys.argv[1:3]


def instance(n, rng, maximize=False, continuous=False):
    choices = list(itertools.product(range(n), repeat=2))
    pairs = [(u, v) for u, v in itertools.combinations(choices, 2)
             if u[0] != v[0] and u[1] != v[1]]
    columns = [(u,) for u in choices] + pairs
    rng.shuffle(columns)
    costs = [rng.randint(-15, 15) for _ in columns]
    rows = []
    for axis, group in itertools.product(range(2), range(n)):
        rows.append(([int(len(c) == 1 and c[0][axis] == group) for c in columns], 1, 1))
    for u in choices:
        for group in range(n):
            if group != u[0]:
                rows.append(([-1 if c == (u,) else
                              int(len(c) == 2 and u in c and any(v[0] == group for v in c))
                              for c in columns], 0, 0))
    rng.shuffle(rows)
    m = dict(rows=rows, integer=[len(c) == 1 or not continuous for c in columns],
             bounds=[(0, 1)] * len(columns), cost=costs, offset=Fraction(3, 8), maximize=maximize)
    values = []
    for permutation in itertools.permutations(range(n)):
        chosen = set(enumerate(permutation))
        values.append(m['offset'] + sum(v for c, v in zip(columns, costs) if chosen.issuperset(c)))
    return m, (max if maximize else min)(values)


class PairAssignment(unittest.TestCase):
    def solve(self, m, limit=100000):
        with tempfile.TemporaryDirectory() as temp:
            model, witness = Path(temp)/'random.mps', Path(temp)/'w.sol'
            model.write_text(integer_mps(m))
            run = subprocess.run([SOLVER, model, '--mip', '--node-limit', str(limit),
                                  '--solution', witness], capture_output=True, text=True, timeout=10)
            self.assertIn(run.returncode, (0, 4), (run.stdout, run.stderr))
            data = json.loads(run.stdout)
            value = None
            if witness.exists():
                check = subprocess.run([CHECKER, '--model', model, '--solution', witness, '--exact'],
                                       text=True, capture_output=True, timeout=5, check=True)
                checked = json.loads(check.stdout)
                self.assertTrue(checked['checked_feasible'])
                value = Fraction(checked['objective'])
            return data, value

    def test_permutations(self):
        for seed in range(24):
            model, expected = instance(2 + seed % 4, random.Random(seed),
                                       maximize=bool(seed & 1), continuous=bool(seed & 2))
            with self.subTest(seed=seed):
                result, value = self.solve(model)
                self.assertIn('Pair assignment search exhausted', result['message'])
                self.assertFalse(result['certifies_optimality'])
                self.assertEqual((result['status'], value), ('optimal', expected))

    def test_mixed_partner_groups(self):
        model, _ = instance(3, random.Random(70))
        links = [a for a, lo, _ in model['rows'] if lo == 0]
        a = links[0]
        owner = a.index(-1)
        b = next(other for other in links if other is not a and other[owner] == -1)
        # Both remaining groups supply partners. An owner with two selected
        # partners cannot satisfy this extra row if all products are set to one.
        model['rows'].append(([-1 if j == owner else a[j] + b[j] for j in range(len(a))], 0, 0))
        result, value = self.solve(model)
        self.assertNotIn('Pair assignment', result['message'])
        self.assertEqual(result['status'], 'optimal')
        self.assertIsNotNone(value)

    def test_full_width_mask(self):
        n = 16
        rows = []
        for j in range(n):
            a = [int(i == j) for i in range(n + 1)]
            rows.extend([(a, 1, 1), (a[:], 1, 1)])
        for j in [0, n - 1]:
            a = [0] * (n + 1)
            a[j], a[n] = -1, 1
            rows.append((a, 0, 0))
        model = dict(rows=rows, integer=[True]*n + [False], bounds=[(0, None)]*(n+1),
                     cost=[1000000000]*n + [-1000000000], maximize=True)
        result, value = self.solve(model)
        self.assertIn('Pair assignment search exhausted', result['message'])
        self.assertEqual(value, (n-1)*1000000000)

    def test_node_limit(self):
        model, _ = instance(4, random.Random(10))
        result, _ = self.solve(model, limit=3)
        self.assertNotIn('Pair assignment', result['message'])
        self.assertLessEqual(result['nodes'], 3)

    def test_near_matches(self):
        for change in ['merged_product', 'coefficient', 'cost', 'unlinked', 'missing_product',
                       'missing_owner', 'fractional_rhs', 'fractional_box']:
            model, _ = instance(2, random.Random(30))
            zero_rows = [r for r, (_, lo, _) in enumerate(model['rows']) if lo == 0]
            r = zero_rows[0]
            terms = model['rows'][r][0]
            aux = terms.index(1)
            if change == 'merged_product':
                # Two distinct columns for one pair can split the flow and cost.
                for a, _, _ in model['rows']:
                    a.append(a[aux])
                model['integer'].append(True)
                model['cost'].append(-100)
                model['bounds'].append((0, 1))
            elif change == 'coefficient':
                terms[aux] = 1 + Fraction(1, 10**20)
            elif change == 'cost':
                model['cost'][aux] += Fraction(1, 10**20)
            elif change == 'unlinked':
                terms[aux] = 0
            elif change == 'missing_product':
                # The cheapest relaxed matching cannot satisfy an empty link.
                owner = terms.index(-1)
                for a, _, _ in model['rows']:
                    a.pop(aux)
                model['cost'][owner] = -100
                for field in ['integer', 'bounds', 'cost']:
                    model[field].pop(aux)
            elif change == 'missing_owner':
                terms[terms.index(-1)] = 0
            elif change == 'fractional_rhs':
                # A genuine change large enough for the ordinary solver as well.
                model['rows'][r] = (terms, Fraction(1, 2), Fraction(1, 2))
            else:
                model['bounds'][aux] = (0, Fraction(1, 2))
            with self.subTest(change=change):
                # For the deliberately sub-double coefficient only check abstention.
                # The general numerical solver has its own feasibility tolerances.
                if change == 'coefficient':
                    model['cost'][aux] = 100
                result, value = self.solve(model)
                self.assertNotIn('Pair assignment', result['message'])
                if change != 'coefficient':
                    self.assertEqual(value, mip_oracle(model))


if __name__ == '__main__':
    unittest.main()
