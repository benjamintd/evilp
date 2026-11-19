"""Verify standalone shortcut points and duals on shifted, scaled and ranged models."""
from fractions import Fraction as F
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

from test_lp_solver import mps

SOLVER, CHECKER = (str(Path(p).resolve()) for p in sys.argv[1:3])
del sys.argv[1:3]


class StandaloneProposal(unittest.TestCase):
    def solve(self, model, options=()):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            path, point, dual = [root/n for n in ('model.mps','point.sol','dual.txt')]
            path.write_text(mps(model))
            r = subprocess.run([SOLVER,path,'--relax','--solution',point,'--dual',dual,*options],
                               capture_output=True,text=True,timeout=12)
            self.assertEqual(r.returncode,0,(r.stdout,r.stderr))
            data = json.loads(r.stdout)
            self.assertEqual(data['status'],'optimal')
            c = subprocess.run([CHECKER,'--model',path,'--solution',point,'--dual',dual,'--relax',
                                '--lin-tol','1e-7','--obj-tol','1e-7'],capture_output=True,text=True)
            self.assertEqual(c.returncode,0,(c.stdout,c.stderr))
            check = json.loads(c.stdout)
            self.assertTrue(check['checked_feasible'])
            self.assertEqual(check['bound_status'],'valid')
            return data,F(check['objective']),F(check['certified_bound'])

    def test_signed_transforms_and_objective_direction(self):
        n=259
        cost=[2,-3,1]+[0]*(n-3)
        bounds=[(None,2),(None,None),(-1,4)]+[(0,1)]*(n-3)
        rows=[]
        def row(entries,lo,hi):
            a=[0]*n
            for j,v in entries:a[j]=v
            rows.append((a,lo,hi))
        row([(0,2),(1,2)],4,4)
        row([(1,8),(2,8)],8,24)
        for j in range(3,n):row([(j,1)],0,1)
        for maximize in (False,True):
            m=dict(cost=cost,bounds=bounds,rows=rows,maximize=maximize)
            fast,x,y=self.solve(m)
            old,ox,oy=self.solve(m,['--glpk-basis'])
            self.assertEqual((x,y),(ox,oy))
            self.assertEqual(fast['glpk_basis_status'],'source_verified_point_and_bound')
            self.assertNotEqual(old['glpk_basis_status'],fast['glpk_basis_status'])

    def test_rounded_cost_does_not_hide_material_gap(self):
        def padded(n):
            a=[1,1]+[0]*(n-2)
            return dict(cost=[1,F('1.00000000000000001')]+[0]*(n-2),
                        bounds=[(0,10**15)]*2+[(0,1)]*(n-2),
                        rows=[(a,10**15,10**15)],maximize=True)
        expected=F('1.00000000000000001')*10**15
        model=padded(1024)
        for options in ([], ['--glpk-basis']):
            with self.subTest(options=options):
                r,x,y=self.solve(model, options)
                self.assertEqual((x,y),(expected,expected))
                self.assertNotEqual(r['glpk_basis_status'],'source_verified_point_and_bound')
        # One original equation and 8,192 implicit bound equations exceed
        # the new sparse limit. The ambiguous basis must still fail closed.
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'oversized.mps'
            path.write_text(mps(padded(8192)))
            run=subprocess.run([SOLVER,path,'--relax','--glpk-basis'],capture_output=True,text=True)
            self.assertEqual(run.returncode,4,run.stdout)
            rejected=json.loads(run.stdout)
            self.assertEqual(rejected['status'],'numerical_failure')
            self.assertFalse(rejected['has_primal'])
            self.assertGreater(rejected['source_pricing_attempts'],0)


if __name__=='__main__':
    unittest.main()
