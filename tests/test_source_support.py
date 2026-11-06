#!/usr/bin/env python3
"""Outward support bounds checked against exact correlated interval corners."""
from pathlib import Path
from fractions import Fraction as Q
import json,itertools,subprocess,random,sys
binary=Path(sys.argv[1]).resolve();rng=random.Random(20261106);rows=[]
unit=subprocess.run([str(binary)],capture_output=True,text=True,timeout=30)
assert unit.returncode==0,(unit.returncode,unit.stdout,unit.stderr)
assert json.loads(unit.stdout)=={'passed':True,'checks':27}
for i in range(240):
 x=rng.choice([0.,.25,.5,1.]);y=rng.choice([-2.,-1.,1.,2.]);a=1.;da=rng.choice([0.,2.**-50,2.**-45]);dc=2.**-48;db=rng.choice([0.,2.**-50]);bound=1. if i%3!=1 else 1e15
 lo,hi=-bound,bound;curlo,curhi=lo,hi
 if i%3==2:curlo,curhi=.125,.875;x=.5
 # Model representatives match their source enclosures; bound overlays differ.
 q=[a-da,a+da,a*y-dc,a*y+dc,x,y,x-db,x+db,x-db,x+db,lo,lo,hi,hi,curlo,curhi,lo,hi,x,x]
 rows.append(q)
for i in range(24):
 x=1e15;y=(-1.)**i;a=1.;da=2.**-50;b=x;db=.125
 rows.append([a-da,a+da,a*y,a*y,x,y,b-db,b+db,b-db,b+db,0.,0.,2e15,2e15,0.,2e15,0.,2e15,b,b])

def gap(q):
 a0,a1,c0,c1,x,y,bl0,bl1,bu0,bu1,l0,l1,u0,u1,cl,cu,sl,su,rl,ru=map(Q,q)
 lows=(l0,l1) if cl==sl else (cl,cl);ups=(u0,u1) if cu==su else (cu,cu)
 worst=Q(0)
 for a,c,bl,bu,lo,hi in itertools.product((a0,a1),(c0,c1),(bl0,bl1),(bu0,bu1),lows,ups):
  g=abs(y*(a*x-(bl if y>0 else bu)))
  r=c-a*y
  if r:g+=abs(r*(x-(lo if r>0 else hi)))
  worst=max(worst,g)
 return worst
text=''.join(' '.join(repr(float(x)) for x in q)+'\n' for q in rows);p=subprocess.run([str(binary),'--stream'],input=text,capture_output=True,text=True,timeout=30);assert p.returncode==0,(p.returncode,p.stderr)
actual=[tuple(map(int,l.split())) for l in p.stdout.splitlines()];assert len(actual)==len(rows)
records=[];accepted=rejected=abstained=0
for i,(q,(av,bad)) in enumerate(zip(rows,actual)):
 worst=gap(q);safe=av==1 and bad==0
 if safe:assert worst<=Q(1,10**7),(i,q,worst,av,bad);accepted+=1
 elif av==1:rejected+=1
 else:abstained+=1
 # Small finite boxes are deliberately safely within the absolute threshold.
 if i<240 and i%3!=1:assert safe,(i,q,worst,av,bad)
 records.append({'id':i,'input':q,'maximum_exact_corner_gap':str(worst),'available':av,'bad':bad,'accepted':safe})
assert accepted>=160 and rejected>=80 and not abstained,(accepted,rejected,abstained)
j={'passed':True,'cases':len(rows),'exact_corners_per_case':64,'accepted':accepted,'rejected':rejected,'abstained':abstained,'oracle':'Exact Fraction exhaustive endpoints, retaining shared coefficient correlation between row activity and column reduced cost. Absolute1e-7 acceptance; no epsilon clipping.','records':records}
print(json.dumps({k:v for k,v in j.items() if k!='records'}))
