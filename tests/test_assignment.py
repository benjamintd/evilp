"""Exact enumeration, independent decimal witness checks and structural refusals."""
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from test_mip_solver import integer_mps

SOLVER,CHECKER=map(lambda p:str(Path(p).resolve()),sys.argv[1:3])
del sys.argv[1:3]


def model(costs, maximize=False):
    n=len(costs)
    m=n*n
    rows=[]
    for i in range(n):
        rows.append(([int(k//n==i) for k in range(m)],1,1))
    for j in range(n):
        rows.append(([int(k%n==j) for k in range(m)],1,1))
    return dict(cost=sum(costs,[]), integer=[True]*m,bounds=[(0,1)]*m,rows=rows,
                maximize=maximize,offset=F(1,8))


class Assignment(unittest.TestCase):
    def solve(self,m,*flags):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'original.mps'
            sol=Path(t)/'point.sol'
            path.write_text(integer_mps(m) if isinstance(m,dict) else m)
            r=subprocess.run([SOLVER,path,'--mip','--solution',sol,*flags],
                             capture_output=True,text=True,timeout=15)
            data=json.loads(r.stdout)
            self.assertNotIn('runtime error:',r.stderr)
            value=None
            if sol.exists():
                c=subprocess.run([CHECKER,'--exact','--model',path,'--solution',sol],
                                 capture_output=True,text=True,timeout=10)
                self.assertEqual(c.returncode,0,(data,c.stdout))
                check=json.loads(c.stdout)
                self.assertTrue(check['checked_feasible'])
                self.assertFalse(check['certifies_optimality'])
                value=F(check['objective'])
            return data,value

    def test_permutations_and_signed_costs(self):
        rng=random.Random(741)
        for n in (2,3,4,5):
            for maximize in (False,True):
                c=[[F(rng.randrange(-50,51),8) for _ in range(n)] for _ in range(n)]
                m=model(c,maximize)
                expected=(max if maximize else min)(sum(c[i][v] for i,v in enumerate(p))
                            for p in itertools.permutations(range(n)))
                data,obj=self.solve(m)
                self.assertEqual((data['search_route'],data['bound_origin']),('assignment','structural_bound'))
                self.assertEqual(F(str(data['best_bound'])),expected+m['offset'])
                self.assertEqual(obj,expected+m['offset'])

    def test_decimal_potentials(self):
        m=model([[0,3],[2,0]])
        m['integer'] += [False,False]
        m['cost'] += [0,0]
        m['bounds'] += [(0,4),(0,4)]
        for row,_,_ in m['rows']: row.extend([0,0])
        m['rows'] += [([0,0,0,0,1,0],F(1,10),F(1,10)),
                      ([0,0,0,0,-1,1],F(1,5),F(1,5))]
        text=integer_mps(m).replace(' 0.1\n',' 0.100000000000000001\n')
        data,obj=self.solve(text)
        self.assertIn('exact matching bound',data['message'],data)
        self.assertEqual(obj,F(1,8))

    def test_cost_between_adjacent_doubles_cannot_certify_wrong_matching(self):
        m=model([[1,1],[1,1]])
        text=integer_mps(m).replace('x0 OBJ 1\n','x0 OBJ 1.000000000000000001\n')
        self.assertNotEqual(text,integer_mps(m))
        data,obj=self.solve(text)
        if 'exact matching bound' in data['message']:
            self.assertEqual(obj,F(17,8))
        self.assertEqual(data['status'],'optimal')

    def test_charge_on_continuous_variable_and_nonunit_cover(self):
        m=model([[3,0],[0,1]],True)
        m['cost'].append(1)
        m['bounds'].append((0,1))
        m['integer'].append(False)
        for r,_,_ in m['rows']:r.append(0)
        data,obj=self.solve(m)
        self.assertNotIn('exact matching bound',data['message'])
        self.assertEqual(obj,F(41,8))
        base=model([[5,0],[0,4]],True)
        text=integer_mps(base)
        modified=text.replace('x0 r0 1\n','x0 r0 0\n')
        self.assertNotEqual(text,modified)
        data,obj=self.solve(modified)
        self.assertNotIn('exact matching bound',data['message'])
        self.assertEqual(obj,F(1,8))

    def quotient_model(self, maximize=False, reverse=False):
        edges = list(itertools.product(range(4), repeat=2))
        costs = [[int((a < 2) != (b < 2)) * (-1 if maximize else 1)
                  for b in range(4)] for a in range(4)]
        m = model(costs, maximize)
        m['integer'] += [False, False]
        m['bounds'] += [(0, 1), (0, 1)]
        m['cost'] += [0, 0]
        for row, _, _ in m['rows']:
            row.extend([0, 0])
        flags = [int(a < 2 and b < 2) for a, b in edges] + [0, 0]
        m['rows'].append(([-v for v in flags], -1, None) if reverse else (flags, None, 1))
        m['rows'].append(([0]*16 + [-1, 1], F(1, 10), None))
        return m

    def test_interchangeable_vertices_with_infeasible_cheapest_matching(self):
        for maximize, reverse in itertools.product((False, True), repeat=2):
            m = self.quotient_model(maximize, reverse)
            if reverse:
                m['rows'].reverse()
            objectives = [m['offset'] + sum(m['cost'][4*a+b] for a,b in enumerate(p))
                          for p in itertools.permutations(range(4))
                          if sum(a < 2 and b < 2 for a,b in enumerate(p)) <= 1]
            with self.subTest(maximize=maximize, reverse=reverse):
                result, value = self.solve(m, '--time-limit', '10')
                self.assertEqual(result['status'], 'optimal')
                self.assertEqual(value, (max if maximize else min)(objectives))
                if result['scip_route'] != 'unavailable':
                    self.assertEqual(result['scip_route'], 'assignment_quotient')
                self.assertFalse(result['certifies_optimality'])

    def test_sub_binary64_side_signatures_remain_distinct(self):
        m = self.quotient_model()
        # All sixteen side coefficients round to the same binary64 value.
        # No two assignment vertices have the same original side column.
        coeff = [1 + F(k, 10**20) for k in range(16)] + [0, 0]
        m['rows'].append((coeff, None, 10))
        result, value = self.solve(m, '--time-limit', '10')
        self.assertEqual((result['status'], value), ('optimal', m['offset']+2))
        self.assertNotEqual(result['scip_route'], 'assignment_quotient')

    def test_quotient_respects_node_limit(self):
        m = self.quotient_model()
        result, value = self.solve(m, '--node-limit', '2')
        self.assertLessEqual(result['nodes'], 2)
        self.assertNotEqual(result['scip_route'], 'assignment_quotient')
        if result['status'] == 'optimal':
            self.assertEqual(value, m['offset']+2)

    def test_inconsistent_cycle_never_returns_a_matching(self):
        m=model([[0,0],[0,0]])
        m['cost'] += [0,0]
        m['integer'] += [False,False]
        m['bounds'] += [(0,3),(0,3)]
        for row,_,_ in m['rows']:row.extend([0,0])
        m['rows'] += [([0,0,0,0,1,-1],1,None),([0,0,0,0,-1,1],1,None)]
        data,_=self.solve(m,'--node-limit','8')
        self.assertNotIn('exact matching bound',data['message'])
        self.assertFalse(data['has_primal'])
        self.assertLessEqual(data['nodes'],8)


if __name__=='__main__':
    unittest.main()
