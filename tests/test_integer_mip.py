"""Original-model oracles for the C17 integer branch-and-cut integration."""
import itertools
import json
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path

from test_mip_solver import integer_mps

SOLVER, CHECKER = (Path(p).resolve() for p in sys.argv[1:3])
del sys.argv[1:3]


class IntegerSearch(unittest.TestCase):
    def model(self, seed, maximize=False):
        rng = random.Random(seed)
        n = 36
        groups = [j % 4 for j in range(n)]
        rng.shuffle(groups)
        bounds = [(-1, 2)] * n
        costs = [rng.randrange(-5, 7) for _ in range(n)]
        rows = []
        for g in range(4):
            members = [j for j in range(n) if groups[j] == g]
            for a, b in zip(members, members[1:]):
                row = [0] * n
                row[a], row[b] = -3, 3
                rows.append((row, 0, 0))
            bounds[members[-1]] = (0, 2)
        rows.append(rows[2])
        rng.shuffle(rows)
        weight = [F(rng.randrange(-8, 13), 8) for _ in range(n)]
        rows.append((weight, -1, 3))
        return dict(cost=costs, bounds=bounds, rows=rows, maximize=maximize,
                    integer=[True]*n, offset=7), groups

    def run_model(self, directory, text, options=()):
        model = directory/'test.mps'
        sol = directory/'point.sol'
        model.write_text(text)
        r = subprocess.run([SOLVER, model, '--mip', '--solution', sol, '--time-limit', '6',
                            *options], capture_output=True, text=True, timeout=12)
        return json.loads(r.stdout), model, sol

    def test_exact_integer_objectives_and_lift(self):
        for seed in range(12):
            maximize = seed % 2 == 1
            m, groups = self.model(seed, maximize)
            values = []
            weights, lo, hi = m['rows'][-1]
            for assignment in itertools.product(range(3), repeat=4):
                x = [assignment[g] for g in groups]
                if lo <= sum(a*b for a,b in zip(weights,x)) <= hi:
                    values.append(m['offset']+sum(a*b for a,b in zip(m['cost'],x)))
            expected = (max if maximize else min)(values)
            with self.subTest(seed=seed), tempfile.TemporaryDirectory() as tmp:
                v, model, sol = self.run_model(Path(tmp), integer_mps(m))
                self.assertEqual(v['status'], 'optimal', v)
                self.assertIn('SCIP reported optimality', v['message'])
                self.assertFalse(v['certifies_optimality'])
                check = subprocess.run([CHECKER,'--exact','--model',model,'--solution',sol],
                                       capture_output=True,text=True)
                c = json.loads(check.stdout)
                self.assertTrue(c['checked_feasible'],c)
                self.assertEqual(F(c['objective']),expected)

    def test_subdouble_coefficient_and_node_limit(self):
        m, _ = self.model(291)
        text = integer_mps(m)
        for tag, options in [('subdouble', []), ('node', ['--node-limit','1'])]:
            with self.subTest(tag=tag), tempfile.TemporaryDirectory() as tmp:
                problem = text
                if tag == 'subdouble':
                    lines = problem.splitlines()
                    start = lines.index('COLUMNS')
                    j = next(j for j in range(start,len(lines)) if ' OBJ ' in lines[j])
                    fields = lines[j].split()
                    lines[j] = ' '+ ' '.join(fields[:2])+' '+str(abs(int(fields[2])))+'.00000000000000001'
                    problem = '\n'.join(lines)+'\n'
                v, _, _ = self.run_model(Path(tmp), problem, options)
                self.assertNotIn('SCIP reported optimality', v['message'])
                if tag == 'node':
                    self.assertLessEqual(v['nodes'],1)


if __name__ == '__main__':
    unittest.main()
