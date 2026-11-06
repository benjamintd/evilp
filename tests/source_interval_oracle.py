from fractions import Fraction as F
import math,random,subprocess,json,sys
rng=random.Random(916083)
def decimal(q):
 d=q.denominator;n=0
 while d%2==0:n+=1;d//=2
 assert d==1
 text=str(abs(q.numerator)*5**n).zfill(n+1)
 return ('-' if q<0 else '')+(text if not n else text[:-n]+'.'+text[-n:])
known=['0.5','1.0','1e0','1D+3','9007199254740992.0','18446744073709551616','-0.125','1.25e2','-7.8125e-3']
texts=known+['0e999999','0.0e-999999','1e-9999','-1e9999','0.100000000000000001','2.000000000000000001']
for i in range(2500):
 q=F(rng.randrange(-(1<<60),1<<60),1<<rng.randrange(0,1051));texts.append(decimal(q))
for i in range(2500):texts.append(f'{rng.randrange(-10**18,10**18)}e{rng.randrange(-340,310)}')
run=subprocess.run([sys.argv[1]],input='\n'.join(texts)+'\n',text=True,capture_output=True,timeout=60);assert run.returncode==0 and not run.stderr,run.stderr
lines=run.stdout.splitlines();assert len(lines)==len(texts);points=0
for i,(text,line) in enumerate(zip(texts,lines)):
 q=F(text.replace('D','e').replace('d','e'));rounded,lo,hi=map(float.fromhex,line.split())
 assert not math.isnan(lo) and not math.isnan(hi) and lo<=hi,(text,line)
 if math.isfinite(lo):assert F.from_float(lo)<=q,(text,line)
 if math.isfinite(hi):assert F.from_float(hi)>=q,(text,line)
 if lo==hi:
  assert math.isfinite(lo) and F.from_float(lo)==q,(text,line)
  points+=1
 if i<len(known):assert lo==hi==rounded,(text,line)
print(json.dumps({'passed':True,'source_number_cases':len(texts),'proved_points':points,'known_decimal_scientific_points':len(known),'false_point_claims':0}))
