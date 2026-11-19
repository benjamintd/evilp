"""Exact dependence negative controls and guarded pricing against exact optima."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,subprocess,tempfile,copy
from test_lp_solver import mps,oracle
from test_pricing_integration import beale,triangular,weighted,random_models,cancellation_models
p=argparse.ArgumentParser()
for name in ['solver','checker','driver']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();solver,checker,driver=(str(x.resolve()) for x in [a.solver,a.checker,a.driver])
base={'cost':[1,2],'bounds':[(0,2),(0,2)],'rows':[([1,1],1,1),([2,2],2,2)]}
models=[(base,[2,-1],1),(base,[0,0],0),(base,[1,1],0)]
bad=copy.deepcopy(base);bad['rows'][1]=([2,2],F('2.000000000000000001'),F('2.000000000000000001'));models.append((bad,[2,-1],0))
bad=copy.deepcopy(base);bad['rows'][1]=([2,F('2.000000000000000001')],2,2);models.append((bad,[2,-1],0))
bad=copy.deepcopy(base);bad['rows'][1]=([2,2],None,2);models.append((bad,[2,-1],0))
checked=dependencies=0
with tempfile.TemporaryDirectory(prefix='adaptive-dependency-') as t:
 root=Path(t)
 for i,(model,y,expected) in enumerate(models):
  path=root/f'dep{i}.mps';path.write_text(mps(model))
  for mode in ['normal','unknown','expired']:
   out=subprocess.run([driver,str(path),mode,*map(str,y)],capture_output=True,text=True,timeout=10)
   assert out.returncode==0,(out.stdout,out.stderr)
   assert int(out.stdout)==(expected if mode=='normal' else 0),(i,mode,out.stdout)
   dependencies+=1
 cases=[beale(),*[triangular(n) for n in range(2,6)],weighted(8),weighted(32),weighted(12,True),*cancellation_models(),*list(random_models())]
 for i,model in enumerate(cases):
  expected=F(5**(i+1)) if 1<=i<=4 else F(-len(model['cost'])) if i in [5,6] else F(-12) if i==7 else oracle(model)
  for k,flags in enumerate([extra+fast+crash for crash in [[],['--crash-basis']] for fast in [[],['--fast-residual']] for extra in [[],['--backend','dense'],['--no-scaling'],['--bland-ratio','--no-updates']]]):
   path=root/f'lp{i}-{k}.mps';sol=root/f'lp{i}-{k}.sol';path.write_text(mps(model))
   out=subprocess.run([solver,str(path),'--relax','--adaptive-pricing','--solution',str(sol),'--time-limit','10','--iteration-limit','2000',*flags],capture_output=True,text=True,timeout=15)
   assert out.returncode==0,(i,k,out.stdout,out.stderr)
   result=json.loads(out.stdout);assert result['status']=='optimal' and result['has_primal']
   assert abs(F(str(result['objective']))-expected)<=F(1,10**7),(i,k,result,expected)
   out=subprocess.run([checker,'--model',str(path),'--solution',str(sol),'--relax','--lin-tol','1e-7','--obj-tol','1e-7'],capture_output=True,text=True,timeout=10)
   assert out.returncode==0,(i,k,out.stdout,out.stderr)
   result=json.loads(out.stdout);assert result['checked_feasible'] and abs(F(result['objective'])-expected)<=F(1,10**7)
   checked+=1
print(json.dumps({'passed':True,'source_dependency_controls':dependencies,'independently_checked_adaptive_points':checked}))
