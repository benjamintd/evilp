#!/usr/bin/env python3
"""Independent exact feasibility/oracle and isolation contracts for primal diving."""
import itertools
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F
from test_mip_solver import integer_mps, mip_oracle
from test_lp_solver import inequalities
DRIVER=str(pathlib.Path(sys.argv.pop(1)).resolve())
CHECKER=str(pathlib.Path(sys.argv.pop(1)).resolve())
SOLVER=str(pathlib.Path(sys.argv.pop(1)).resolve())
COUNTS={'calls':0,'candidates':0,'independent_checks':0,'oracle_cases':0,'source_points':0}
EVIDENCE=[]
BACKTRACK={'cost':[1,0],'bounds':[(0,1),(0,F(3,5))],
           'rows':[([1,1],1,None)],'integer':[1,0]}
COUPLED={'cost':[1,0],'bounds':[(0,1),(0,10**6)],
         'rows':[([1,0],F(99999999,10**8),None),([10**6,-1],0,0)],'integer':[1,0]}
class Diving(unittest.TestCase):
    def call(self,model,mode='default',root=None):
        if root is None:
            with tempfile.TemporaryDirectory(prefix='dive contract ') as td:
                return self.call(model,mode,pathlib.Path(td))
        model_path=root/'model.mps';model_path.write_text(integer_mps(model));point_path=root/'point.sol'
        run=subprocess.run([DRIVER,str(model_path),str(point_path),mode],capture_output=True,text=True,timeout=20)
        self.assertEqual(run.returncode,0,(run.stdout,run.stderr));self.assertFalse(run.stderr,run.stderr)
        result=json.loads(run.stdout);COUNTS['calls']+=1
        self.assertEqual(result['root_status'],'optimal',result)
        self.assertTrue(result['inputs_unchanged'],result);self.assertTrue(result['main_session_unchanged'],result)
        self.assertLessEqual(result['steps'],128,result)
        self.assertLessEqual(result['backtracks'],8,result)
        if result['candidate']:
            COUNTS['candidates']+=1
            text=point_path.read_text();point={line.split()[0]:F(line.split()[1]) for line in text.splitlines() if not line.startswith('=')};x=[point[f'x{j}'] for j in range(len(model['cost']))]
            for j,flag in enumerate(model['integer']):
                if flag:self.assertEqual(x[j].denominator,1,(model,x,result))
            for a,b in inequalities(model):self.assertLessEqual(sum(F(c)*v for c,v in zip(a,x))-b,F(1,10**7),(model,x,result))
            checked=subprocess.run([CHECKER,'--model',str(model_path),'--solution',str(point_path),'--lin-tol','1e-7','--int-tol','0','--obj-tol','1e-7'],capture_output=True,text=True,timeout=10)
            self.assertEqual(checked.returncode,0,(model,result,checked.stdout,checked.stderr));check=json.loads(checked.stdout)
            self.assertTrue(check['checked_feasible']);self.assertFalse(check['certifies_optimality']);COUNTS['independent_checks']+=1
            COUNTS['source_points']+=1
            optimum=mip_oracle(model);self.assertIsNotNone(optimum,model)
            objective=F(model.get('offset',0))+sum(F(c)*v for c,v in zip(model['cost'],x))
            if model.get('maximize'):self.assertLessEqual(objective,optimum+F(1,10**7))
            else:self.assertGreaterEqual(objective,optimum-F(1,10**7))
        else:self.assertFalse(point_path.exists(),result)
        EVIDENCE.append({'mode':mode,'model':model,'result':result})
        return result
    def test_sibling_backtracking_and_cold_control(self):
        for mode in ['default','cold']:
            r=self.call(BACKTRACK,mode);self.assertTrue(r['candidate']);self.assertEqual(r['backtracks'],1);self.assertEqual(r['lp_solves'],2);self.assertEqual(r['objective'],1)
        r=self.call(BACKTRACK,'no-backtracks');self.assertFalse(r['candidate']);self.assertEqual(r['status'],'not_found')
    def test_final_fixing_moves_coupled_continuous_coordinate(self):
        for mode in ['default','cold']:
            r=self.call(COUPLED,mode);self.assertTrue(r['candidate']);self.assertEqual(r['final_resolves'],1);self.assertEqual(r['objective'],1)
    def test_abstention_and_resource_controls_keep_main_session(self):
        for mode,expected in [('tiny-budget','time_limit'),('row-cap','resource_limit'),('invalid-root','invalid_input'),('lp-abstain','lp_abstain')]:
            r=self.call(BACKTRACK,mode);self.assertFalse(r['candidate']);self.assertEqual(r['status'],expected)
        r=self.call(BACKTRACK,'one-step');self.assertFalse(r['candidate']);self.assertEqual(r['steps'],1)
    def test_exact_binary_and_mixed_enumeration(self):
        rng=random.Random(127991)
        for k in range(48):
            n=4 if k<24 else 3;mixed=k>=24;cols=n+int(mixed)
            cost=[rng.randrange(1,8) for _ in range(n)]+([F(rng.randrange(1,5),2)] if mixed else [])
            bounds=[(0,1)]*n+([(0,3)] if mixed else [])
            rows=[]
            for _ in range(3):
                a=[rng.randrange(0,5) for _ in range(n)]+([1] if mixed else [])
                if not any(a):a[0]=1
                upper=sum(a[:n])+int(mixed)*3
                rhs=F(rng.randrange(1,2*upper+1),2)
                rows.append((a,rhs,None))
            model={'cost':cost,'bounds':bounds,'rows':rows,'integer':[1]*n+[0]*int(mixed)}
            if k%6==0:
                model['cost']=[-c for c in model['cost']];model['maximize']=True
            optimum=mip_oracle(model);self.assertIsNotNone(optimum,model);COUNTS['oracle_cases']+=1
            for mode in ['default','cold']:
                with self.subTest(k=k,mode=mode):self.call(model,mode)
    def test_large_integer_value_abstains_without_changing_search_session(self):
        model={'cost':[1,F(1,2**67)],'bounds':[(0,1),(0,None)],
               'rows':[([2**67,1],2**65,None)],'integer':[1,1]}
        r=self.call(model)
        self.assertFalse(r['candidate']);self.assertEqual(r['status'],'lp_abstain')
        self.assertEqual(r['message'],'Integer dive value exceeds the exactly spaced range')
    def test_bb_hook_finds_incumbent_before_branching(self):
        with tempfile.TemporaryDirectory(prefix='dive bb contract ') as td:
            root=pathlib.Path(td);model=root/'model.mps';model.write_text(integer_mps(BACKTRACK))
            results=[]
            for mode,args in [('enabled',[]),('disabled',['--no-root-diving'])]:
                point=root/(mode+'.sol');run=subprocess.run([SOLVER,str(model),'--mip','--no-root-cuts','--no-presolve','--node-limit','1','--time-limit','10','--solution',str(point),*args],capture_output=True,text=True,timeout=15)
                self.assertIn(run.returncode,[0,4],(run.stdout,run.stderr));j=json.loads(run.stdout);results.append(j)
                if mode=='enabled':
                    self.assertTrue(j['has_primal'],j);self.assertEqual(j['root_dive_incumbents'],1);self.assertGreater(j['root_dive_lp_solves'],0)
                    check=subprocess.run([CHECKER,'--model',str(model),'--solution',str(point),'--int-tol','0'],capture_output=True,text=True,timeout=10);self.assertEqual(check.returncode,0,check.stdout)
                else:self.assertFalse(j['has_primal'],j);self.assertEqual(j['root_dive_attempts'],0)
                self.assertEqual(j['status'],'node_limit',j)
            self.assertEqual(results[0]['nodes'],results[1]['nodes'])
            self.assertEqual(results[0]['lp_solves'],results[1]['lp_solves'])
            self.assertEqual(results[0]['lp_iterations'],results[1]['lp_iterations'])
    def test_bb_failure_abstains_without_node_pruning(self):
        parity={'cost':[0],'bounds':[(0,1)],'rows':[([2],1,1)],'integer':[1]}
        with tempfile.TemporaryDirectory(prefix='dive bb failure ') as td:
            root=pathlib.Path(td);model=root/'model.mps';model.write_text(integer_mps(parity))
            runs=[]
            for args in [[],['--no-root-diving']]:
                r=subprocess.run([SOLVER,str(model),'--mip','--no-root-cuts','--no-presolve','--node-limit','1','--time-limit','10',*args],capture_output=True,text=True,timeout=15);self.assertEqual(r.returncode,4,r.stdout);j=json.loads(r.stdout);self.assertEqual(j['status'],'node_limit');self.assertFalse(j['has_primal']);self.assertEqual(j['pruned_nodes'],0);runs.append(j)
            self.assertEqual(runs[0]['lp_solves'],runs[1]['lp_solves']);self.assertEqual(runs[0]['best_bound'],runs[1]['best_bound'])
            self.assertEqual(runs[0]['root_dive_status'],'not_found');self.assertEqual(runs[0]['root_dive_backtracks'],1)
if __name__=='__main__':
    result=unittest.main(exit=False)
    print(json.dumps(COUNTS,sort_keys=True))
    if result.result.wasSuccessful():raise SystemExit(0)
    raise SystemExit(1)
