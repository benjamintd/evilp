#!/usr/bin/env python3
"""Bounded native-only rational provenance, ownership and independent proofs."""
import argparse,json,subprocess,sys,tempfile
from pathlib import Path

def checked(argv):
 run=subprocess.run([str(x) for x in argv],text=True,capture_output=True,timeout=180)
 if run.returncode:raise AssertionError((argv,run.returncode,run.stdout,run.stderr))
 return run.stdout

if __name__=='__main__':
 parser=argparse.ArgumentParser()
 for name in ['owner','zero','reader','producer','checker']:parser.add_argument('--'+name,type=Path,required=True)
 args=parser.parse_args();here=Path(__file__).resolve().parent
 with tempfile.TemporaryDirectory(prefix='milplab-native-rationals-') as temporary:
  root=Path(temporary);large=root/'large'
  result=json.loads(checked([sys.executable,here/'native_rational_pipeline.py',args.producer.resolve(),args.checker.resolve(),large]))
  assert result['passed'] and result['independently_valid_cut_proofs']==27
  source=large/'uncertain_zero.mps';proof=root/'owner.gmi'
  owner=json.loads(checked([args.owner.resolve(),source,proof]))
  assert owner['passed'] and owner['native_zero_coefficients']==0 and owner['global_exact_primal_attempts']==0
  validation=json.loads(checked([args.checker.resolve(),'--model',source,'--gmi',proof]))
  assert validation['certifies_cut_validity'] and not validation['certifies_original_lp_infeasibility'] and not validation['certifies_original_mip_infeasibility'] and not validation['certifies_full_tree']
  checked([args.zero.resolve(),source])
  limits=json.loads(checked([sys.executable,here/'native_reader_limits.py',args.reader.resolve(),root/'limits']))
  assert limits['passed'] and limits['controls']==10
 print(json.dumps({'passed':True,'large_source_proofs':27,'owned_model_proofs':1,'source_limit_controls':10,'source_zero_support_scales':2,'global_exact_primal_attempts_on_native_only':0}))
