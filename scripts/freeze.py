#!/usr/bin/env python3
"""Record the complete evaluator/configuration revision without changing files."""
import argparse,hashlib,json,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[1]
paths=sorted([ROOT/'CMakeLists.txt',*ROOT.glob('src/*'),*ROOT.glob('include/*'),*ROOT.glob('config/*'),*ROOT.glob('tests/*.py'),*ROOT.glob('tests/*.c'),*ROOT.glob('tests/fixtures/*'),*ROOT.glob('scripts/*.py')])
lock={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
ap=argparse.ArgumentParser();ap.add_argument('output');ap.add_argument('--revision');args=ap.parse_args()
revision=subprocess.run(['git','rev-parse',args.revision or 'HEAD'],cwd=ROOT,capture_output=True,text=True).stdout.strip() or None
if args.revision:
    if not revision:raise SystemExit('Unknown Git revision')
    names=subprocess.check_output(['git','ls-tree','-r','--name-only',revision],cwd=ROOT,text=True).splitlines()
    names=[n for n in names if n=='CMakeLists.txt' or n.startswith(('src/','include/','config/','tests/','scripts/'))]
    lock={n:hashlib.sha256(subprocess.check_output(['git','show',revision+':'+n],cwd=ROOT)).hexdigest() for n in names}
data={'schema':1,'git_revision':revision,'source':'git_revision' if args.revision else 'working_tree','files':lock}
pathlib.Path(args.output).write_text(json.dumps(data,indent=2)+'\n')
