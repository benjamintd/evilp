#!/usr/bin/env python3
"""Record the complete evaluator/configuration revision without changing files."""
import hashlib,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
paths=sorted([ROOT/'CMakeLists.txt',*ROOT.glob('src/*'),*ROOT.glob('include/*'),*ROOT.glob('config/*'),*ROOT.glob('tests/*.py'),*ROOT.glob('tests/fixtures/*'),*ROOT.glob('scripts/*.py')])
lock={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
revision=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True).stdout.strip() or None
data={'schema':1,'git_revision':revision,'files':lock}
if len(sys.argv)!=2:raise SystemExit('usage: python3 scripts/freeze.py OUTPUT_JSON')
pathlib.Path(sys.argv[1]).write_text(json.dumps(data,indent=2)+'\n')
