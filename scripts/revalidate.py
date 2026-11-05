#!/usr/bin/env python3
"""Rerun independent validation without changing historical solver runs."""
import argparse,hashlib,json,pathlib,subprocess,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('--checker',default='build/milpcheck');args=ap.parse_args()
    directory=pathlib.Path(args.directory).resolve();checker=pathlib.Path(args.checker).resolve();digest=sha(checker)
    rows=[json.loads(x) for x in (directory/'runs.jsonl').read_text().splitlines() if x.strip()]
    counts={}
    with open(directory/'revalidated.jsonl','x') as f:
        for i,r in enumerate(rows):
            prefix=directory/r['instance']/r['solver']/('seed-'+str(r['seed']))
            sol=prefix/'solution.sol';model=pathlib.Path(r['model_path'])
            if not model.exists():model=ROOT/'data/miplib2017/instances'/(r['instance']+'.mps.gz')
            if sha(model)!=r['model_sha256']:raise SystemExit('Model changed: '+r['instance'])
            r['original_validation']={k:r[k] for k in ['checker_sha256','checker','checker_exit_code','checker_wall_s']}
            elapsed=0;code=None
            if sol.exists():
                start=time.monotonic()
                try:
                    p=subprocess.run([str(checker),'--model',str(model),'--solution',str(sol),'--format',r['solver']],capture_output=True,text=True,timeout=120)
                    code=p.returncode;c=json.loads(p.stdout)
                except (subprocess.TimeoutExpired,json.JSONDecodeError):
                    c={'status':'unavailable','checked_feasible':False}
                elapsed=time.monotonic()-start
            else:c={'status':'unavailable','checked_feasible':False}
            r.update(checker_sha256=digest,checker=c,checker_exit_code=code,checker_wall_s=elapsed)
            r['revalidation']={'checker_sha256':digest,'solver_times_unchanged':True}
            logfile=prefix/'solver.log'
            if r['solver']=='highs' and logfile.exists() and 'Parser reached timeout' in logfile.read_text(errors='replace'):
                r['original_claimed_status']=r['claimed_status'];r['claimed_status']='input_time_limit'
            (prefix/'checker-revalidated.json').write_text(json.dumps({'checker_sha256':digest,'exit_code':code,'checker':c},indent=2)+'\n')
            f.write(json.dumps(r)+'\n');f.flush();counts[c['status']]=counts.get(c['status'],0)+1
            if (i+1)%60==0:print(i+1,'revalidated',flush=True)
    if sha(checker)!=digest:raise SystemExit('Checker changed during validation')
    (directory/'revalidation.json').write_text(json.dumps({'checker_sha256':digest,'runs':len(rows),'outcomes':counts,'solver_times_unchanged':True},indent=2)+'\n')
    print(json.dumps(counts,indent=2))
if __name__=='__main__':main()
