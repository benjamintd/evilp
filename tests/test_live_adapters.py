#!/usr/bin/env python3
"""Optional native integration test: all three engines, known MIN/MAX optima."""
import argparse,csv,fractions,hashlib,json,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='results/adapter-regressions');args=ap.parse_args()
    fixtures=ROOT/'work/live-fixtures';fixtures.mkdir(parents=True,exist_ok=True)
    original=(ROOT/'tests/fixtures/tiny.mps').read_text()
    suite=fixtures/'suite.tsv'
    with suite.open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['id','path','sha256','reference_status','reference_objective'])
        for sense,obj in [('MAX','8.03'),('MIN','5.03')]:
            m=fixtures/(sense.lower()+'.mps');m.write_text(original.replace(' MAX',' '+sense))
            w.writerow([sense.lower(),str(m),hashlib.sha256(m.read_bytes()).hexdigest(),'=opt=',obj])
    subprocess.run([str(ROOT/'build/milpbench'),'--suite',str(suite),'--out',args.out,'--limit','5','--seed','1'],cwd=ROOT,check=True)
    rows=[json.loads(x) for x in (ROOT/args.out/'runs.jsonl').read_text().splitlines()]
    assert len(rows)==6
    for r in rows:
        assert r['claimed_status']=='optimal',r
        assert r['checker']['checked_feasible'],r
        assert fractions.Fraction(r['checker']['objective'])==fractions.Fraction(r['reference_objective']),r
    print('All six native MIN/MAX adapter checks passed.')
if __name__=='__main__':main()
