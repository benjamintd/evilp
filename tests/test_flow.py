"""Independent gate enumeration with integral augmenting-path feasibility."""
from collections import deque
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import random
import subprocess
import sys
import tempfile
import unittest

from test_mip_solver import integer_mps

DRIVER, SOLVER, CHECKER = (str(Path(a).resolve()) for a in sys.argv[1:4])
del sys.argv[1:4]


def make_model(demand, edges, maximize=False):
    n = len(edges)
    rows = []
    for v, d in enumerate(demand):
        a = [int(head == v) - int(tail == v) for tail, head, _, _ in edges] + [0]*n
        rows.append((a, d, d))
    for i, (_, _, cap, _) in enumerate(edges):
        a = [0]*(2*n)
        a[i], a[n+i] = 1, -cap
        rows.append((a, None, 0))
    return dict(rows=rows, bounds=[(0,e[2]) for e in edges]+[(0,1)]*n,
                cost=[0]*n+[(-1 if maximize else 1)*e[3] for e in edges],
                integer=[False]*n+[True]*n, maximize=maximize, offset=F(19,8))


def feasible(demand, edges, mask):
    n = len(demand)
    source, sink = n, n+1
    cap = [[0]*(n+2) for _ in range(n+2)]
    for i,(u,v,bound,_) in enumerate(edges):
        if mask & (1<<i):
            cap[u][v] += bound
    if sum(demand): return False
    required = 0
    for i,b in enumerate(demand):
        if b < 0: cap[source][i] = -b
        if b > 0: cap[i][sink] = b; required += b
    while required:
        parent = {source:source}
        queue = deque([source])
        while queue and sink not in parent:
            u = queue.popleft()
            for v,c in enumerate(cap[u]):
                if c>0 and v not in parent:
                    parent[v]=u
                    queue.append(v)
        if sink not in parent: return False
        f=required
        v=sink
        while v!=source:
            u=parent[v];f=min(f,cap[u][v]);v=u
        v=sink
        while v!=source:
            u=parent[v];cap[u][v]-=f;cap[v][u]+=f;v=u
        required-=f
    return True


def oracle(demand,edges):
    result=None
    for mask in range(1<<len(edges)):
        price=sum(e[3] for i,e in enumerate(edges) if mask & (1<<i))
        if (result is None or price<result) and feasible(demand,edges,mask): result=price
    return result


class Flow(unittest.TestCase):
    def invoke(self,data,limit=100000,normal=False):
        with tempfile.TemporaryDirectory() as temp:
            m,p=Path(temp)/'network.mps',Path(temp)/'point.sol'
            m.write_text(integer_mps(data))
            args=[SOLVER,m,'--mip','--solution',p] if normal else [DRIVER,m,p,str(limit)]
            proc=subprocess.run(args,capture_output=True,text=True,timeout=10)
            self.assertEqual(proc.returncode,0,(proc.stdout,proc.stderr))
            r=json.loads(proc.stdout);value=None
            if p.exists():
                check=subprocess.run([CHECKER,'--model',m,'--solution',p,'--exact'],
                                     capture_output=True,text=True,timeout=10)
                self.assertEqual(check.returncode,0,(check.stdout,check.stderr))
                value=F(json.loads(check.stdout)['objective'])
            return r,value

    def test_series_parallel_networks_against_cut_independent_oracle(self):
        rng=random.Random(8991)
        for case in range(64):
            # Start at a triangle and add vertices on previous edges to keep
            # treewidth <= 2. Subgraphs, reversed arcs and parallel arcs vary.
            n=rng.randint(3,6)
            pairs=[(0,1),(0,2),(1,2)]
            for v in range(3,n):
                a,b=rng.choice(pairs);pairs.extend([(a,v),(b,v)])
            if case%2:pairs.append(rng.choice(pairs))
            edges=[];d=[0]*n
            for u,v in pairs:
                if rng.randrange(2):u,v=v,u
                cap=rng.randint(1,7); flow=rng.randrange(cap+1)
                edges.append((u,v,cap,rng.randrange(10)))
                d[u]-=flow;d[v]+=flow
            if case%5==0:d[0]+=10;d[-1]-=10
            best=oracle(d,edges)
            with self.subTest(case=case):
                r,x=self.invoke(make_model(d,edges,case%2))
                self.assertEqual(bool(r['solved']),best is not None)
                self.assertEqual(x,None if best is None else F(19,8)+(-best if case%2 else best))

    def test_multicomponent_wide_and_blocked_networks(self):
        edges=[(0,1,1000000000,5),(0,1,1000000000,3),(2,3,1000000000,2)]
        model=make_model([-999999999,999999999,-1,1],edges)
        r,v=self.invoke(model,normal=True)
        self.assertEqual((r['status'],v),('optimal',F(19,8)+5))
        self.assertGreater(r['flow_eliminations'],0)
        self.assertFalse(r['certifies_optimality'])
        r,v=self.invoke(model,limit=1)
        self.assertEqual((r['solved'],v),(0,None))
        self.assertLessEqual(r['eliminations'],1)
        # K4 has no allowable elimination.
        edges=[(a,b,3,1) for a,b in itertools.combinations(range(4),2)]
        self.assertEqual(self.invoke(make_model([-1,-1,1,1],edges))[0]['solved'],0)

    def test_few_terminal_directed_networks(self):
        rng=random.Random(669)
        for i in range(30):
            demand=[-4,1,2,1]
            if i%2: demand=[-d for d in demand]
            edges=[]
            for a,b in itertools.combinations(range(4),2):
                if i%2:a,b=b,a
                edges.append((a,b,4,rng.randrange(8)))
                if rng.randrange(2):edges.append((b,a,4,rng.randrange(8)))
            best=oracle(demand,edges)
            r,v=self.invoke(make_model(demand,edges))
            self.assertIsNotNone(best)
            self.assertEqual(v,F(19,8)+best)
            self.assertEqual(r['states'],28)
        # An arc too small for the aggregate supply invalidates the tree reduction.
        edges[0]=(edges[0][0],edges[0][1],1,0)
        self.assertEqual(self.invoke(make_model(demand,edges))[0]['states'],0)
        r,v=self.invoke(make_model(demand,[(a,b,4,c) for a,b,_,c in edges]),limit=10)
        self.assertEqual((r['solved'],v),(0,None))

    def test_decimal_bounds_coefficients_costs_and_missing_gates(self):
        delta=F(1,10**21)
        for kind in ['demand','coefficient','gate','upper','cost','negative','flow_cost']:
            m=make_model([-1,1],[(0,1,3,2)])
            if kind=='demand':m['rows'][0]=(m['rows'][0][0],-1+delta,-1+delta)
            elif kind=='coefficient':m['rows'][0][0][0]=-1+delta
            elif kind=='gate':m['rows'][-1][0][1]=-3+delta
            elif kind=='upper':m['bounds'][0]=(0,3-delta)
            elif kind=='cost':m['cost'][1]+=delta
            elif kind=='negative':m['cost'][1]=-1
            else:m['cost'][0]=1
            with self.subTest(kind=kind):
                self.assertEqual(self.invoke(m)[0]['solved'],0)

if __name__=='__main__':
    unittest.main()
