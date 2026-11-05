#!/usr/bin/env python3
"""Download pinned native baselines and MIPLIB; no global installation."""
import argparse, concurrent.futures, hashlib, json, os, pathlib, shutil, subprocess, tarfile, zipfile
ROOT = pathlib.Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def fetch(url,path,expected=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        part=path.with_suffix(path.suffix+'.part')
        subprocess.run(['curl','--fail','--location','--retry','3','--silent','--show-error',url,'-o',str(part)],check=True)
        part.replace(path)
    digest=sha(path)
    if expected and digest!=expected: raise RuntimeError(f'Checksum mismatch: {path}')
    return digest
def solvers():
    spec=json.loads((ROOT/'config/baseline-assets.json').read_text())
    if os.uname().sysname!='Darwin' or os.uname().machine!='arm64':
        raise SystemExit('Pinned binary bundle is Darwin arm64; supply native binaries in config/solvers.tsv on other hosts.')
    rows=['solver\tbinary\tsha256\tversion']
    provenance={}
    def one(item):
        name,s=item
        native=ROOT/'tools/cbc-native/bin/cbc'
        native_prov=ROOT/'tools/cbc-source-provenance.json'
        if name=='cbc' and native.exists() and native_prov.exists():
            p=json.loads(native_prov.read_text())
            if sha(native)==p['binary_sha256']:
                return name,native.resolve(),{**s,'version':'2.10.13 (local source build)','source_provenance':'tools/cbc-source-provenance.json'}
        archive=ROOT/'work/downloads'/s['url'].split('/')[-1]
        fetch(s['url'],archive,s['sha256'])
        dest=ROOT/'tools'/name
        if not (dest/'.complete').exists():
            dest.mkdir(parents=True,exist_ok=True)
            with tarfile.open(archive) as tf: tf.extractall(dest,filter='data')
            (dest/'.complete').touch()
        candidates=[p for p in dest.rglob(name) if p.is_file() and p.parent.name=='bin']
        if not candidates: raise RuntimeError(f'No {name} executable in {dest}')
        binary=candidates[0].resolve()
        return name,binary,s
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for name,binary,s in sorted(pool.map(one,spec['solvers'].items())):
            rows.append(f'{name}\t{binary.relative_to(ROOT)}\t{sha(binary)}\t{s["version"]}')
            provenance[name]={**s,'binary':str(binary.relative_to(ROOT)),'binary_sha256':sha(binary)}
            print('Installed',name,s['version'],flush=True)
    (ROOT/'config/solvers.tsv').write_text('\n'.join(rows)+'\n')
    (ROOT/'tools/provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    if not (ROOT/'tools/cbc-native/bin/cbc').exists():
        import sys
        subprocess.run([sys.executable,str(ROOT/'scripts/build_cbc.py')],check=True)
def miplib():
    root=ROOT/'data/miplib2017'; root.mkdir(parents=True,exist_ok=True)
    urls={'benchmark.zip':'https://miplib.zib.de/downloads/benchmark.zip',
          'benchmark-v2.test':'https://miplib.zib.de/downloads/benchmark-v2.test',
          'miplib2017-v37.solu':'https://miplib.zib.de/downloads/miplib2017-v37.solu'}
    lockpath=ROOT/'config/miplib.lock.json'
    previous=json.loads(lockpath.read_text()) if lockpath.exists() else {}
    sources={}
    for name,url in urls.items():
        path=root/name
        digest=fetch(url,path,previous.get('sources',{}).get(name,{}).get('sha256'))
        sources[name]={'url':url,'sha256':digest,'bytes':path.stat().st_size}
        print('Fetched',name,flush=True)
    names=[line.strip() for line in (root/'benchmark-v2.test').read_text().splitlines() if line.strip()]
    if len(names)!=240 or len(set(names))!=240: raise RuntimeError('Expected 240 unique MIPLIB benchmark-v2 files')
    instances=root/'instances'; instances.mkdir(exist_ok=True)
    with zipfile.ZipFile(root/'benchmark.zip') as z:
        members={pathlib.PurePosixPath(n).name:n for n in z.namelist() if n.endswith('.mps.gz')}
        for name in names:
            target=instances/name
            if not target.exists():
                with z.open(members[name]) as src,open(target,'wb') as dst: shutil.copyfileobj(src,dst)
    refs={}
    for line in (root/'miplib2017-v37.solu').read_text().splitlines():
        fields=line.split()
        if len(fields)>=2: refs[fields[1]]=(fields[0],fields[2] if len(fields)>2 else '-')
    rows=['id\tpath\tsha256\treference_status\treference_objective']
    records=[]
    for name in names:
        ident=name.removesuffix('.mps.gz'); path=instances/name
        status,obj=refs.get(ident,('=unk=','-'))
        digest=sha(path)
        rows.append(f'{ident}\t{path.relative_to(ROOT)}\t{digest}\t{status}\t{obj}')
        records.append({'id':ident,'sha256':digest,'reference_status':status,'reference_objective':obj,
                        'source_url':'https://miplib.zib.de/WebData/instances/'+name,
                        'metadata_url':'https://miplib.zib.de/instance_details_'+ident+'.html',
                        'license':'unresolved-per-instance; fetch-only, no redistribution'})
    lock={'dataset':'MIPLIB 2017 benchmark v2','reference':'v37','count':len(records),'sources':sources,'instances':records}
    if previous and previous!=lock: raise RuntimeError('Dataset differs from lock; explicit reviewed lock update required')
    lockpath.write_text(json.dumps(lock,indent=2)+'\n')
    (root/'suite.tsv').write_text('\n'.join(rows)+'\n')
    print('Imported and hashed all 240 instances',flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('target',choices=['solvers','miplib','all']); args=ap.parse_args()
    if args.target in ('solvers','all'): solvers()
    if args.target in ('miplib','all'): miplib()
