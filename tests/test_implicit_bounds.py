#!/usr/bin/env python3
"""Source-proved row omission and warm reanchoring against exact LP oracles."""
import argparse,copy,json,subprocess,tempfile
from pathlib import Path
from fractions import Fraction as F
from test_lp_solver import mps,oracle,decimal

ap=argparse.ArgumentParser();ap.add_argument('--solver',type=Path,required=True);ap.add_argument('--checker',type=Path,required=True);ap.add_argument('--driver',type=Path,required=True);args=ap.parse_args()
solver,checker,driver=(p.resolve() for p in [args.solver,args.checker,args.driver])
base={'cost':[3,2],'bounds':[(0,5),(0,5)],'rows':[([2,1],None,7)],'maximize':True}
both=[[(0,5),(0,5)],[(2,5),(0,5)],[(2,2),(0,5)],[(0,5),(0,5)],[(0,5),(0,2)],[(0,5),(0,5)]]
cases=[('both',base,both,True),('large',base,both,True)]
packing={'cost':[3,2],'bounds':[(0,1),(0,1)],'rows':[([1,1],None,1)],'maximize':True}
packing_steps=[[(0,1),(0,1)],[(0,F(1,2)),(0,1)],[(0,1),(0,1)],[(0,1),(1,1)],[(0,1),(0,1)]]
cases.append(('packing',packing,packing_steps,False))
for name,coefficient,rhs,lower in [('almost_one',F('0.99999999999999999'),1,0),
                                    ('large_rhs',1,F('1.00000000000000001'),0),
                                    ('negative',1,1,-1)]:
    model=copy.deepcopy(packing)
    model['rows']=[([1,coefficient],None,rhs)]
    model['bounds'][1]=(lower,1)
    cases.append((name,model,[model['bounds']],False))
upper={'cost':[3,2],'bounds':[(None,5),(0,5)],'rows':[([-1,1],None,3),([1,1],None,7)],'maximize':True}
cases.append(('upper',upper,[[(None,u),(0,5)] for u in [5,3,5,1,5]],True))
lower=copy.deepcopy(base);lower['bounds'][0]=(0,None)
cases.append(('lower',lower,[[(x,None),(0,5)] for x in [0,2,3,0,1,0]],True))
ranged={'cost':[3,2],'bounds':[(-2,4),(-3,5)],'rows':[([2,1],F(-3,2),4),([1,-1],-3,4)],'maximize':True}
cases.append(('ranged',ranged,[[(-2,4),(-3,5)],[(0,4),(-3,5)],[(2,2),(-3,5)],[(-2,4),(-3,5)],[(-2,4),(-1,5)]],False))
uncertain=copy.deepcopy(base);q=F('0.100000000000000001');uncertain['bounds'][0]=(q,5)
cases.append(('uncertain',uncertain,[[(q,5),(0,5)],[(F(1,2),5),(0,5)],[(1,1),(0,5)],[(q,5),(0,5)]],False))

def render(model,large):
    text=mps(model)
    if large:
        text=text.replace('COLUMNS\n',''.join(f' L pad{i}\n' for i in range(199))+'COLUMNS\n')
        text=text.replace('RHS\n',''.join(f' x{i} OBJ 0\n' for i in range(2,300))+'RHS\n')
        text=text.replace('ENDATA\n',''.join(f' FX b x{i} 0\n' for i in range(2,300))+'ENDATA\n')
    return text

def bound_token(x,upper):return ('inf' if upper else '-inf') if x is None else decimal(x)

def verify(root,name,model,steps,keep,dense,noscale,forget=False,strict=False):
    root.mkdir();source=root/'model.mps';source.write_text(render(model,name=='large'));path=root/'steps.txt'
    path.write_text(''.join(' '.join(bound_token(v,u) for pair in step for u,v in enumerate(pair))+'\n' for step in steps))
    out=root/'points';out.mkdir()
    run=subprocess.run([str(driver),str(source),str(path),str(out),str(keep),str(dense),str(noscale),str(int(forget))],capture_output=True,text=True,timeout=60)
    assert run.returncode==0,(name,keep,dense,noscale,run.stdout,run.stderr)
    records=[json.loads(x) for x in run.stdout.splitlines()];assert len(records)==2*len(steps)
    for record in records:
        step=record['step'];view=copy.deepcopy(model);view['bounds']=steps[step];expected=oracle(view);assert expected is not None
        current=root/f'view{step}-{record["cold"]}.mps';current.write_text(render(view,name=='large'))
        point=out/f'step{step}-cold{record["cold"]}.sol'
        checked=subprocess.run([str(checker),'--model',str(current),'--solution',str(point),'--relax','--lin-tol','1e-7','--obj-tol','1e-7'],capture_output=True,text=True,timeout=30)
        assert checked.returncode==0,(name,record,checked.stdout,checked.stderr)
        result=json.loads(checked.stdout);assert result['checked_feasible']
        assert abs(F(result['objective'])-expected)<=F(1,10**7),(name,record,result,expected)
        assert abs(F(str(record['objective']))-expected)<=F(1,10**7)
        if strict and step and not record['cold']:assert record['warm']==1,(name,record)
        if step==0:
            wanted=0 if keep else 1 if name=='uncertain' or forget else 2
            assert record['implicit']==wanted,(name,record,wanted)
            assert record['packing']==(2 if name=='packing' and not keep else 0),(name,record)
    return len(records)

count=0
with tempfile.TemporaryDirectory(prefix='implicit-bound-contract-') as temporary:
    root=Path(temporary)
    for name,model,steps,strict in cases:
        for keep in [0,1]:
            for dense in [0,1]:
                count+=verify(root/f'{name}-k{keep}-d{dense}',name,model,steps,keep,dense,0,strict=strict)
        if name in ['both','large']:
            for dense in [0,1]:count+=verify(root/f'{name}-unscaled-d{dense}',name,model,steps,0,dense,1,strict=strict)
    for dense in [0,1]:count+=verify(root/f'unknown-d{dense}','unknown',base,[both[0]],0,dense,0,forget=True)
    for label, delta, rows in [('same',0,1),('distinct',F('0.000000000000000001'),2)]:
        model={'cost':[0,1], 'bounds':[(0,1),(0,1)],
               'rows':[([1,1],1,1),([1,1+delta],1,1)]}
        source=root/f'duplicate-{label}.mps'
        source.write_text(mps(model))
        sol=root/f'duplicate-{label}.sol'
        command=[str(solver),str(source),'--relax','--packing-bounds','--solution',str(sol)]
        solved=subprocess.run(command,capture_output=True,text=True,timeout=10)
        data=json.loads(solved.stdout)
        assert data['standard_rows']==rows,(label,data)
        assert solved.returncode==0 and data['status']=='optimal' and data['objective']==0,(label,data)
        checked=subprocess.run([str(checker),'--model',str(source),'--solution',str(sol),'--relax'],
                               capture_output=True,text=True,timeout=10)
        assert checked.returncode==0 and json.loads(checked.stdout)['checked_feasible']
    # The compact and explicit forms represent the same original LP, including
    # decimal bounds whose binary representatives are insufficient evidence.
    for name,model,steps,_ in cases:
        source=root/f'cold-{name}.mps';source.write_text(render(model,name=='large'))
        for keep in [False,True]:
            run=subprocess.run([str(solver),str(source),'--relax','--time-limit','10',*(['--keep-bound-rows'] if keep else [])],capture_output=True,text=True,timeout=30)
            assert run.returncode==0,(name,run.stdout,run.stderr)
            data=json.loads(run.stdout);assert data['status']=='optimal'
            assert data['implicit_bound_rows']==(0 if keep else 1 if name=='uncertain' else 2)
print(json.dumps({'passed':True,'independently_checked_session_points':count,'direct_form_controls':12,'original_source_and_exact_Fraction_oracle':True}))
