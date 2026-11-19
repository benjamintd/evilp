"""Check projection, integer objective rounding and original witnesses."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import random
import subprocess
import sys
import tempfile
import unittest

from test_mip_solver import integer_mps

SOLVER, CHECKER = (str(Path(s).resolve()) for s in sys.argv[1:3])
del sys.argv[1:3]


def formulation(prices, members, demand, maximize=False):
    # members[j] lists the coverage rows for each auxiliary of gate j.
    aux = sum(map(len, members))
    n = len(prices) + aux
    rows = [[0]*n for _ in range(len(prices)+len(demand))]
    j = len(prices)
    for g, cols in enumerate(members):
        rows[g][g] = -len(cols)
        for covered in cols:
            rows[g][j] = 1
            for i in covered:
                rows[len(prices)+i][j] = 1
            j += 1
    sides = [(None, 0)]*len(prices) + [(d, None) for d in demand]
    return dict(cost=[-x if maximize else x for x in prices] + [0]*aux,
                bounds=[(0,1)]*n, integer=[True]*n, offset=F(3,8), maximize=maximize,
                rows=[(r, lo, hi) for r, (lo,hi) in zip(rows,sides)])


def enumerate_gates(prices, members, demand):
    costs = []
    for x in product((0, 1), repeat=len(prices)):
        out = [0]*len(demand)
        for val, cols in zip(x,members):
            if val:
                for rows in cols:
                    for i in rows:
                        out[i] += 1
        if all(a >= b for a,b in zip(out,demand)):
            costs.append(sum(v*c for v,c in zip(x,prices)))
    return min(costs) if costs else None


class CoverProjection(unittest.TestCase):
    def run_model(self, m):
        with tempfile.TemporaryDirectory() as tmp:
            path, sol = (Path(tmp)/s for s in ('model.mps','solution.sol'))
            path.write_text(integer_mps(m))
            r = subprocess.run([SOLVER,path,'--mip','--solution',sol],capture_output=True,
                               text=True,timeout=15)
            self.assertIn(r.returncode, (0, 4), r.stderr)
            self.assertNotIn('runtime error:', r.stderr)
            rec = json.loads(r.stdout)
            val = None
            if sol.is_file():
                c = subprocess.run([CHECKER,'--exact','--model',path,'--solution',sol],
                                   capture_output=True,text=True,timeout=10)
                self.assertEqual(c.returncode,0,c.stdout)
                val = F(json.loads(c.stdout)['objective'])
            return rec, val

    def test_generated_against_enumeration(self):
        rng = random.Random(839261)
        closures = 0
        for i in range(35):
            n, r = 5+i%3, 2+i%4
            prices = [rng.randint(1,8) for _ in range(n)]
            members = [[rng.sample(range(r),rng.randrange(r+1)) for _ in range(2)]
                       for _ in range(n)]
            demand = [1]*r
            optimum = enumerate_gates(prices,members,demand)
            if optimum is None:
                continue
            m = formulation(prices,members,demand,maximize=bool(i%2))
            result,value = self.run_model(m)
            self.assertEqual(result['status'],'optimal',result)
            self.assertEqual(value,m['offset'] + (-optimum if i%2 else optimum))
            closures += 'Activated covering' in result['message']
        self.assertGreaterEqual(closures, 10)

    def test_fractional_relaxation_keeps_bound_open(self):
        # Triangle covering has LP cost 3 and integer optimum 4.
        m = formulation([2]*3,[[[0,1]],[[1,2]],[[0,2]]],[1]*3)
        r, v = self.run_model(m)
        self.assertNotIn('Activated covering', r['message'])
        self.assertEqual(v,m['offset']+4)

    def test_insufficient_capacity_and_decimal_coefficient(self):
        for tiny in (False,True):
            m = formulation([1, 7],[[[0],[1]],[[0,1]]],[1,1])
            m['rows'][0][0][0] = -(F(19999999999999999,10**16) if tiny else 1)
            r, _ = self.run_model(m)
            self.assertNotIn('Activated covering', r['message'])

if __name__ == '__main__':
    unittest.main()
