#!/usr/bin/env python3
import argparse,concurrent.futures,csv,hashlib,json,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--checker',default='build/milpcheck');ap.add_argument('--out',default='results/import-audit');args=ap.parse_args()
    checker=(ROOT/args.checker).resolve();rows=list(csv.DictReader(open(ROOT/'data/miplib2017/suite.tsv'),delimiter='\t'))
    out=ROOT/args.out;out.mkdir(parents=True,exist_ok=True);results=[]
    def one(r):
        p=subprocess.run([str(checker),'--model',str(ROOT/r['path']),'--inspect'],capture_output=True,text=True,timeout=120)
        return {'instance':r['id'],'exit_code':p.returncode,**json.loads(p.stdout)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool,open(out/'models.jsonl','w') as f:
        for i,r in enumerate(pool.map(one,rows)):
            results.append(r);f.write(json.dumps(r)+'\n');f.flush()
            if r['exit_code'] or (i+1)%40==0:print(i+1,r['instance'],r['status'],flush=True)
    s={'instances':len(results),'supported':sum(r['exit_code']==0 for r in results),'checker_sha256':hashlib.sha256(checker.read_bytes()).hexdigest(),'errors':[r for r in results if r['exit_code']]}
    (out/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2));return bool(s['errors'])
if __name__=='__main__':raise SystemExit(main())
