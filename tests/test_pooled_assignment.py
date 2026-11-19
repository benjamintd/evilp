"""Exact enumeration of opening patterns and capacity-constrained assignments."""
import itertools
import json
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from test_mip_solver import integer_mps

DRIVER, SOLVER, CHECKER = (str(Path(p).resolve()) for p in sys.argv[1:4])
del sys.argv[1:4]


def instance(seed, permute=True):
    rng = random.Random(seed)
    # Repeated types cross quota groups, some group/type pairs are missing.
    types = [0, 1, 2, 0, 2, 1, 1, 0]
    groups = [0, 0, 1, 1, 2, 2, 0, 1]
    quota = [rng.randint(1, 2), rng.randint(1, 2), 1]
    cap = rng.randint(2, 3)
    clients = min(6, sum(quota) * cap)
    costs = [[F(rng.randrange(-18000, 20000), 10000) for _ in range(3)]
             for _ in range(clients)]
    # Enumerate open sites, then assign clients directly to physical sites.
    # The oracle uses neither class flow nor dual multipliers.
    maximize = seed % 2 == 1
    pick = max if maximize else min
    optimum = None
    for bits in itertools.product((0, 1), repeat=len(types)):
        if any(sum(bits[f] for f, g in enumerate(groups) if g == i) != q
               for i, q in enumerate(quota)):
            continue
        sites = [f for f, b in enumerate(bits) if b]

        @lru_cache(None)
        def assign(d, left):
            if d == clients:
                return F(0)
            values = []
            for i, f in enumerate(sites):
                if left[i]:
                    rest = list(left)
                    rest[i] -= 1
                    values.append(costs[d][types[f]] + assign(d + 1, tuple(rest)))
            return pick(values)
        goal = assign(0, (cap,) * len(sites))
        optimum = goal if optimum is None else pick(goal, optimum)
    offset = F(7, 8)
    columns = [(d, f) for d in range(clients) for f in range(len(types))]
    columns += [(None, f) for f in range(len(types))]
    if permute:
        rng.shuffle(columns)
    rows = []
    for f in range(len(types)):
        rows.append(([(-cap if d is None else 1) if j == f else 0 for d, j in columns], None, 0))
    for i in range(clients):
        rows.append(([int(d == i) for d, _ in columns], 1, 1))
    for g, target in enumerate(quota):
        rows.append(([int(d is None and groups[f] == g) for d, f in columns], target, target))
    if permute:
        rng.shuffle(rows)
    return dict(cost=[0 if d is None else costs[d][types[f]] for d, f in columns],
                rows=rows, integer=[d is None or rng.choice((True, False)) for d, _ in columns],
                bounds=[(0, 1)]*len(columns), maximize=maximize, offset=offset), optimum + offset


class PooledAssignment(unittest.TestCase):
    def run_model(self, model, nodes=10000, iters=10000, normal=False):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            path, point = p/'case.mps', p/'source.sol'
            path.write_text(integer_mps(model))
            cmd = ([SOLVER, path, '--mip', '--solution', point] if normal else
                   [DRIVER, path, str(nodes), str(iters), point])
            run = subprocess.run(cmd, text=True, capture_output=True, timeout=10)
            self.assertEqual(run.returncode, 0, (run.stdout, run.stderr))
            record = json.loads(run.stdout)
            if not normal:
                self.assertLessEqual(record['nodes'], nodes)
            value = None
            if point.exists():
                check = subprocess.run([CHECKER, '--model', path, '--solution', point, '--exact'],
                                       text=True, capture_output=True, timeout=10)
                self.assertEqual(check.returncode, 0, check.stdout)
                c = json.loads(check.stdout)
                self.assertTrue(c['checked_feasible'])
                value = F(c['objective'])
            return record, value

    def test_exhaustive(self):
        branched = False
        for i in range(18):
            m, optimum = instance(i)
            with self.subTest(seed=i):
                r, cost = self.run_model(m)
                self.assertEqual(r['result'], 2, r)
                self.assertEqual(cost, optimum)
                branched |= r['nodes'] > 1
                if r['nodes'] > 1:
                    limited, _ = self.run_model(m, nodes=1)
                    self.assertLess(limited['result'], 2)
            if i == 0:
                run, cost = self.run_model(m, normal=True)
                self.assertEqual(run['status'], 'optimal')
                self.assertIn('Pooled assignment', run['message'])
                self.assertEqual(cost, optimum)
        self.assertTrue(branched)

    def test_source_controls(self):
        m, _ = instance(1, False)
        # A tiny original change to the cost or capacity cannot be rounded away.
        original = m['cost'][0]
        m['cost'][0] += F(1, 10**19)
        self.assertEqual(self.run_model(m)[0]['result'], 0)
        m['cost'][0] = original
        a, _, _ = m['rows'][0]
        a[0] += F(1, 10**19)
        self.assertEqual(self.run_model(m)[0]['result'], 0)
        a[0] = 1
        m['integer'][-1] = False
        self.assertEqual(self.run_model(m)[0]['result'], 0)
        m['integer'][-1] = True
        self.assertEqual(self.run_model(m, iters=1)[0]['result'], 0)

    def test_missing_edges_and_group_quotas(self):
        m, _ = instance(2, False)
        # Coalesce two physical assignments into the same facility/client
        # position, leaving another absent. Column count alone still matches.
        m['rows'][0][0][1] = 1
        m['rows'][1][0][1] = 0
        self.assertEqual(self.run_model(m)[0]['result'], 0)
        m, _ = instance(3, False)
        row, _, _ = m['rows'][-1]
        m['rows'][-1] = (row, sum(row) + 1, sum(row) + 1)
        r, value = self.run_model(m)
        self.assertLess(r['result'], 2)
        self.assertIsNone(value)
        m, _ = instance(4, False)
        m['rows'][0][0][-8] -= 1
        self.assertEqual(self.run_model(m)[0]['result'], 0)


if __name__ == '__main__':
    unittest.main()
