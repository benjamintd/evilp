"""Exact-cover costs, incomplete searches, source-decimal eligibility and native MIP routing."""
import json
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from test_mip_solver import integer_mps

DRIVER, SOLVER, CHECKER = [str(Path(p).resolve()) for p in sys.argv[1:4]]
del sys.argv[1:4]


def model(rows, masks, costs, maximize=False):
    return {'cost': costs, 'integer': [True] * len(costs), 'bounds': [(0, 1)] * len(costs),
            'rows': [([int(bool(mask & (1 << i))) for mask in masks], 1, 1) for i in range(rows)],
            'maximize': maximize}


def oracle(rows, masks, costs, maximize):
    # Include/exclude each original column. Keep the best cost for each covered
    # set; no dual prices, ordering by reduced cost or branch pruning.
    best = {0: 0}
    pick = max if maximize else min
    for mask, cost in zip(masks, costs):
        for state, previous in list(best.items()):
            if state & mask:
                continue
            new = state | mask
            best[new] = pick(best.get(new, previous + cost), previous + cost)
    return best.get((1 << rows) - 1)


class Partition(unittest.TestCase):
    def run_case(self, data, limit=100000, mode=0):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'cover.mps'
            path.write_text(integer_mps(data))
            p = subprocess.run([DRIVER, str(path), str(limit), str(mode)], capture_output=True,
                               text=True, timeout=10)
            self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
            r = json.loads(p.stdout)
            self.assertLessEqual(r['nodes'], limit)
            if r['point']:
                for coeff, _, _ in data['rows']:
                    self.assertEqual(sum(coeff[j] for j in r['columns']), 1)
            else:
                self.assertFalse(r['optimal'])
            return r

    def test_subset_oracle(self):
        rng = random.Random(8842)
        for case in range(64):
            rows = rng.randint(2, 8)
            masks = [rng.randrange(1, 1 << rows) for _ in range(rng.randint(4, 22))]
            if case % 5:
                masks += [1 << i for i in range(rows)]
            masks += masks[:2]
            costs = [rng.randint(-35, 65) for _ in masks]
            maximize = bool(case % 2)
            data = model(rows, masks, costs, maximize)
            target = oracle(rows, masks, costs, maximize)
            with self.subTest(case=case):
                r = self.run_case(data, mode=case % 4)
                if target is None:
                    self.assertFalse(r['point'])
                else:
                    self.assertTrue(r['optimal'], r)
                    self.assertEqual(sum(costs[j] for j in r['columns']), target)

    def test_limits_last_row_and_extreme_cost(self):
        masks = [1 << j for j in range(64)] + [(1 << 64) - 1]
        costs = [1 << 30] * 64 + [-(1 << 30)]
        data = model(64, masks, costs)
        self.assertEqual(self.run_case(data, mode=2)['columns'], [64])
        self.assertFalse(self.run_case(data, limit=1)['optimal'])
        self.assertFalse(self.run_case(data, mode=4)['optimal'])

    def test_wide_decimal_costs(self):
        rng = random.Random(952)
        for case in range(20):
            masks = [1 << i for i in range(6)] + [rng.randrange(1, 64) for _ in range(14)]
            costs = [Fraction(rng.randint(-5000000000, 8000000000), 1000000)
                     for _ in masks]
            maximize = bool(case % 2)
            result = self.run_case(model(6, masks, costs, maximize), mode=case % 4)
            self.assertTrue(result['optimal'])
            self.assertEqual(sum(costs[j] for j in result['columns']),
                             oracle(6, masks, costs, maximize))

    def test_nonmatching_source_and_domains(self):
        for kind in ('coefficient', 'cost', 'rhs', 'bound', 'continuous', 'empty'):
            data = model(2, [1, 2, 3], [1, 1, 1])
            near_one = Fraction('1.0000000000000000001')
            if kind == 'coefficient':
                data['rows'][0][0][0] = near_one
            elif kind == 'rhs':
                data['rows'][0] = (data['rows'][0][0], near_one, near_one)
            elif kind == 'cost':
                data['cost'][0] = near_one
            elif kind == 'bound':
                data['bounds'][0] = (0, near_one)
            elif kind == 'continuous':
                data['integer'][1] = False
            else:
                for row, _, _ in data['rows']:
                    row[2] = 0
            with self.subTest(kind=kind):
                r = self.run_case(data)
                self.assertFalse(r['point'])
                self.assertEqual(r['nodes'], 0)

    def test_normal_mip_entry(self):
        masks = [3, 5, 6, 1, 2, 4] * 12
        costs = [2, 2, 2, 3, 3, 3] * 12
        data = model(3, masks, costs)
        with tempfile.TemporaryDirectory() as directory:
            path, point = [Path(directory) / s for s in ('model.mps', 'point.sol')]
            path.write_text(integer_mps(data))
            p = subprocess.run([SOLVER, str(path), '--mip', '--solution', str(point)],
                               capture_output=True, text=True, timeout=10)
            self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
            r = json.loads(p.stdout)
            self.assertEqual(r['status'], 'optimal')
            self.assertGreater(r['partition_nodes'], 0)
            self.assertEqual(r['nodes'], 1 + r['partition_nodes'])
            self.assertFalse(r['certifies_optimality'])
            c = subprocess.run([CHECKER, '--model', str(path), '--solution', str(point)],
                               capture_output=True, text=True, timeout=10)
            self.assertEqual(c.returncode, 0, (c.stdout, c.stderr))
            self.assertEqual(Fraction(json.loads(c.stdout)['objective']), 5)


if __name__ == '__main__':
    unittest.main()
