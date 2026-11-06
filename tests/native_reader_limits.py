from pathlib import Path
import subprocess,sys,json
reader=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(exist_ok=False)
def model(rows=200,cols=300,entries=None,bounds=None):
 text=['NAME limits','ROWS',' N OBJ']+[f' L R{i}' for i in range(rows)]+['COLUMNS']
 if entries is None:text+=[' X0 OBJ -1 R0 2']+[f' X{j} OBJ 0' for j in range(1,cols)]
 else:text+=entries
 text+=['RHS']+[f' RHS R{i} 3' for i in range(rows)]+['BOUNDS']
 text+=bounds if bounds is not None else [' UP B X0 16']+[f' FX B X{j} 0' for j in range(1,cols)]
 return '\n'.join(text+['ENDATA',''])
cases=[]
cases.append(('allowed',model(),True,True,False))
cases.append(('rowcap',model(rows=2049),True,False,False))
cases.append(('colcap',model(cols=4097),True,False,False))
entries=[' X0 OBJ -1']+[' X0 R0 0.1']*100001+[f' X{j} OBJ 0' for j in range(1,300)]
cases.append(('recordcap',model(entries=entries),True,False,False))
long='1.'+'0'*1023+'1';cases.append(('tokencap',model(entries=[' X0 OBJ -1 R0 '+long]+[f' X{j} OBJ 0' for j in range(1,300)]),True,False,False))
zero='0e4097';cases.append(('exponentcap',model(entries=[' X0 OBJ -1 R0 '+zero]+[f' X{j} OBJ 0' for j in range(1,300)]),True,False,False))
lexeme='0.'+'1234567890'*7
entries=[f' X{j} R{i} {lexeme}' for j in range(300) for i in range(300)]
cases.append(('aggregate_bits',model(rows=300,entries=entries),True,False,False))
lexeme='1.'+'0'*1020+'1'
entries=[f' X{j} R{i} {lexeme}' for j in range(100) for i in range(170)]
cases.append(('textcap',model(cols=100,entries=entries),True,False,False))
cases.append(('unsupported',model().replace('RHS\n','SOS\n'),False,False,True))
cases.append(('post_snapshot_overflow',model(entries=[' X0 OBJ -1 R0 1e308',' X0 R0 1e308']+[f' X{j} OBJ 0' for j in range(1,300)]),False,None,False))
records=[]
for name,text,ok,native,unsupported in cases:
 path=out/(name+'.mps');path.write_text(text);run=subprocess.run([str(reader),str(path)],capture_output=True,text=True,timeout=90);assert not run.returncode and not run.stderr,(name,run.stdout,run.stderr);j=json.loads(run.stdout);assert j['read_ok']==ok and (native is None or j['native_exact']==native) and j['unsupported']==unsupported,(name,j);assert not j['global_exact'] or not ok,(name,j);records.append(dict(case=name,**j));(out/(name+'.checked.json')).write_text(run.stdout)
summary={'passed':True,'controls':len(records),'records':records};(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
