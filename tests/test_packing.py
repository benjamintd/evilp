"""Check capacity lower bounds and original points by enumerating assignment choices."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import random
import subprocess
import sys
import tempfile
import unittest

from test_mip_solver import integer_mps

DRIVER, SOLVER, CHECKER = (str(Path(x).resolve()) for x in sys.argv[1:4])
del sys.argv[1:4]


def model(items, caps, choices, maximize=False):
    gates = len(caps)
    n = gates + len(choices)
    rows = [([0]*n, 1, 1) for _ in range(items)]
    for b, cap in enumerate(caps):
        row = [0]*n
        row[b] = -cap
        rows.append((row, None, 0))
    for j, (covered, load) in enumerate(choices, gates):
        for item in covered:
            rows[item][0][j] = 1
        for b, weight in load:
            rows[items+b][0][j] = weight
    charge = F(-3 if maximize else 3, 10)
    return dict(rows=rows, bounds=[(0,1)]*n, integer=[True]*n, cost=[charge]*gates+[0]*len(choices),
                maximize=maximize, offset=F(7, 8))


def oracle(n, caps, options):
    best = None
    for count in range(1, min(n,len(options))+1):
        for selected in combinations(options, count):
            items = [i for covered, _ in selected for i in covered]
            if sorted(items) != list(range(n)):
                continue
            used = [0]*len(caps)
            for _, entries in selected:
                for b, w in entries:
                    used[b] += w
            if all(x <= cap for x, cap in zip(used, caps)):
                price = sum(x > 0 for x in used)
                best = price if best is None else min(best, price)
    return best


class ResourcePacking(unittest.TestCase):
    def solve(self, m, limit=200, normal=False):
        with tempfile.TemporaryDirectory() as td:
            path, point = Path(td)/'original.mps', Path(td)/'point.sol'
            path.write_text(integer_mps(m))
            args = [SOLVER, path, '--mip', '--solution', point] if normal else [DRIVER,path,str(limit)]
            p = subprocess.run(args, text=True, capture_output=True, timeout=10)
            self.assertEqual(p.returncode, 0, (p.stdout,p.stderr))
            self.assertNotIn('runtime error:', p.stderr)
            data = json.loads(p.stdout)
            if normal:
                check = subprocess.run([CHECKER,'--model',path,'--solution',point,'--exact'],
                                       text=True,capture_output=True,timeout=10)
                self.assertEqual(check.returncode,0,check.stdout)
                return data,F(json.loads(check.stdout)['objective'])
            self.assertLessEqual(data['attempts'],limit)
            if data['solved']:
                for coeff,lo,hi in m['rows']:
                    value = sum(coeff[j] for j in data['columns'])
                    self.assertTrue(lo is None or value >= lo)
                    self.assertTrue(hi is None or value <= hi)
            return data

    def test_assignment_oracle(self):
        rng = random.Random(9408)
        solved = 0
        for case in range(50):
            n = 3 + case%2
            caps = [rng.randint(4, 8) for _ in range(3)]
            options = []
            for i in range(n):
                w = rng.randrange(1, 5)
                for b in range(3):
                    if b == i%3 or rng.randrange(3) != 0:
                        options.append(([i],[(b,w+rng.randrange(2))]))
            if case%2:
                options.append(([0,1],[(0,rng.randrange(1,8))]))
            if case%3==0:
                options.append(([0],[(0,1),(1,1)]))
            m = model(n,caps,options,case%2)
            expected = oracle(n,caps,options)
            r = self.solve(m)
            if r['solved']:
                solved += 1
                self.assertEqual(sum(j<len(caps) for j in r['columns']),expected)
        self.assertGreater(solved,15)

    def test_source_and_bundled_columns(self):
        options = [([i],[(b,4)]) for i in range(2) for b in range(2)]
        # The bundled column permits sharing a bin despite a singleton sum of eight.
        options.append(([0,1],[(0,4)]))
        self.assertEqual(self.solve(model(2,[4,4],options))['attempts'],0)
        for kind in ('cost','unit','demand','capacity','bound','side','charge','continuous'):
            m = model(2,[4,4],options[:-1])
            delta = F(1,10**22)
            if kind=='cost': m['cost'][0]+=delta
            if kind=='unit': m['rows'][0][0][2]+=delta
            if kind=='demand': m['rows'][0]=(m['rows'][0][0],1+delta,1+delta)
            if kind=='capacity': m['rows'][2][0][0]+=delta
            if kind=='bound': m['bounds'][0]=(0,1-delta)
            if kind=='side': m['rows'][2]=(m['rows'][2][0],None,delta)
            if kind=='charge': m['cost'][0]=-m['cost'][0]
            if kind=='continuous': m['integer'][2]=False
            with self.subTest(kind=kind):
                self.assertEqual(self.solve(m)['attempts'],0)

    def test_cli_and_budget(self):
        options = [([i],[(b,2)]) for i in range(4) for b in range(3)]
        m = model(4,[4,4,3],options,True)
        self.assertEqual(self.solve(m,0)['solved'],0)
        r,value = self.solve(m,normal=True)
        self.assertEqual((r['status'],value),('optimal',F(7,8)-F(6,10)))
        self.assertGreater(r['packing_attempts'],0)
        self.assertFalse(r['certifies_optimality'])


if __name__ == '__main__':
    unittest.main()
