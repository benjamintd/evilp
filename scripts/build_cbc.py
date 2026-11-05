#!/usr/bin/env python3
"""Build the pinned Cbc dependency tree in the workspace, then select it."""
import csv,hashlib,json,os,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
    lock=json.loads((ROOT/'config/cbc-source.lock.json').read_text())
    work=ROOT/'work';base=work/'cbc-build';base.mkdir(parents=True,exist_ok=True)
    coinbrew=work/'coinbrew-pinned'
    url='https://raw.githubusercontent.com/coin-or/coinbrew/'+lock['coinbrew_commit']+'/coinbrew'
    subprocess.run(['curl','-fLsS',url,'-o',str(coinbrew)],check=True)
    def run(args,log,env=None):
        with open(base/log,'w') as f:
            p=subprocess.run(args,cwd=base,stdout=f,stderr=subprocess.STDOUT,env=env)
        if p.returncode:raise SystemExit((base/log).read_text()[-6000:])
    run(['bash',str(coinbrew),'fetch','Cbc@'+lock['cbc_tag'],'--no-prompt','--no-third-party'],'fetch.log')
    for name,s in lock['repositories'].items():
        actual=subprocess.check_output(['git','-C',str(base/name),'rev-parse','HEAD'],text=True).strip()
        if actual!=s['commit']:raise SystemExit('Pinned dependency mismatch: '+name)
    env=os.environ.copy();flags='-O3'
    if os.uname().sysname=='Darwin':
        sdk=subprocess.check_output(['xcrun','--show-sdk-path'],text=True).strip()
        headers=pathlib.Path(sdk)/'usr/include/c++/v1'
        if headers.exists():flags+=' -isystem '+str(headers)
    env['CXXFLAGS']=flags
    prefix=ROOT/'tools/cbc-native'
    run(['bash',str(coinbrew),'build','Cbc','--off-line','--no-prompt','--parallel-jobs','4','--tests','main','--no-third-party','--verbosity','3','--prefix',str(prefix),'--enable-cbc-parallel','--enable-static','--disable-shared'],'build.log',env)
    binary=prefix/'bin/cbc';digest=hashlib.sha256(binary.read_bytes()).hexdigest()
    rows=list(csv.DictReader(open(ROOT/'config/solvers.tsv'),delimiter='\t'))
    for r in rows:
        if r['solver']=='cbc':r.update(binary=str(binary.relative_to(ROOT)),sha256=digest,version='2.10.13 (local source build)')
    with open(ROOT/'config/solvers.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
    provenance={**lock,'binary_sha256':digest,'CXXFLAGS':flags,'coinbrew_sha256':hashlib.sha256(coinbrew.read_bytes()).hexdigest(),'flags':'static COIN libraries; cbc-parallel enabled; tests main'}
    (ROOT/'tools/cbc-source-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    path=ROOT/'tools/provenance.json'
    selected=json.loads(path.read_text()) if path.exists() else {}
    selected['cbc']={'version':'2.10.13 (local source build)','binary':str(binary.relative_to(ROOT)),'binary_sha256':digest,'source_lock':'config/cbc-source.lock.json','source_provenance':'tools/cbc-source-provenance.json'}
    path.write_text(json.dumps(selected,indent=2)+'\n')
    print('Selected local Cbc:',binary)
if __name__=='__main__':main()
