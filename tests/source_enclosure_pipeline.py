import json,os,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
producer=Path(sys.argv[1]).resolve();checker=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]).resolve();out.mkdir(exist_ok=False)
cases=[('integer','2','3','0',None,False,None,1),('dyadic','1.0','1.5','0.5',None,False,None,1),('scientific','2D0','3e0','0.0',None,False,None,1),('range','2','3','0','2.5',False,None,1),('uncertain_coefficient','2.000000000000000001','3','0',None,False,None,0),('uncertain_rhs','2','3.000000000000000001','0',None,False,None,0),('uncertain_zero','2','3','0',None,True,None,0),('uncertain_source_bound','2','3','0.1',None,False,None,0),('binary_overlay','2','3','0.1',None,False,('.25','16'),1)]
records=[]
for name,coef,rhs,lo,rng,zero,overlay,count in cases:
 text=['NAME large','ROWS',' N OBJ',' L R0']+[f' L R{i}' for i in range(1,200)]+['COLUMNS'," M0 'MARKER' 'INTORG'",f' X0 OBJ -1 R0 {coef}'," M1 'MARKER' 'INTEND'"]+[f' X{j} OBJ 0' for j in range(1,300)]
 if zero:
  at=text.index(" M1 'MARKER' 'INTEND'");text[at:at]=[' X0 R199 0.1',' X0 R199 -0.100000000000000001']
 text+=['RHS',f' RHS R0 {rhs}']+[f' RHS R{i} 1.0' for i in range(1,200)]
 if rng:text+=['RANGES',f' RNG R0 {rng}']
 text+=['BOUNDS',f' LO B X0 {lo}',' UP B X0 16.0']+[f' FX B X{j} 0' for j in range(1,300)]+['ENDATA','']
 model=out/(name+'.mps');model.write_text('\n'.join(text))
 for dense,scale in [(0,0),(1,0),(0,1)]:
  prefix=out/f'{name}-{dense}-{scale}';cmd=[str(producer),str(model),str(prefix),str(dense),str(scale)]+(list(overlay)+['0'] if overlay else [])
  env=os.environ.copy();env['LAB_TEST_POINT_SOURCE_ONLY']='1'
  run=subprocess.run(cmd,text=True,capture_output=True,timeout=60,env=env);prefix.with_suffix('.process.json').write_text(json.dumps({'argv':cmd,'code':run.returncode,'stdout':run.stdout,'stderr':run.stderr},indent=2));assert run.returncode==0,(name,run.stdout,run.stderr)
  j=json.loads(run.stdout);assert j['cuts']==count,(name,j);assert not run.stderr,run.stderr
  for cut in j['published']:
   checked=subprocess.run([str(checker),'--model',str(model),'--gmi',cut['certificate']],text=True,capture_output=True,timeout=30);assert checked.returncode==0,(name,checked.stdout,checked.stderr);q=json.loads(checked.stdout);assert q['certifies_cut_validity'] and not q['certifies_original_lp_infeasibility']
   Path(cut['certificate']+'.checked.json').write_text(checked.stdout)
  records.append({'case':name,'dense':dense,'no_scaling':scale,'status':j['status'],'cuts':j['cuts'],'controls':j['controls']})
summary={'passed':True,'large_source_rows':200,'large_source_columns':300,'profiles':len(records),'independently_valid_cut_proofs':sum(x['cuts'] for x in records),'unavailable_uncertain_controls':sum(x['cuts']==0 for x in records),'records':records};(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
