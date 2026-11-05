#!/usr/bin/env python3
import hashlib,json,os,pathlib,subprocess,sys,tempfile,unittest
BIN=str(pathlib.Path(sys.argv.pop(1)).resolve());CHECK=str(pathlib.Path(sys.argv.pop(1)).resolve())
FIX=pathlib.Path(__file__).parent/'fixtures'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
class Harness(unittest.TestCase):
    def setup_case(self,root,body):
        fake=root/'fake highs';fake.write_text('#!/usr/bin/env python3\n'+body);fake.chmod(0o755)
        suite=root/'suite.tsv';suite.write_text(f'id\tpath\tsha256\treference_status\treference_objective\ntiny\t{FIX/"tiny.mps"}\t{sha(FIX/"tiny.mps")}\t=opt=\t8.03\n')
        config=root/'solvers.tsv';config.write_text(f'solver\tbinary\tsha256\tversion\nhighs\t{fake}\t{sha(fake)}\ttest\n')
        return [BIN,'--suite',str(suite),'--solvers',str(config),'--checker',CHECK,'--out',str(root/'result'),'--limit','0.1','--grace','0.1']
    def test_valid_solution_and_provenance(self):
        with tempfile.TemporaryDirectory(prefix='milp lab ') as td:
            root=pathlib.Path(td)
            body='import pathlib,sys\np=pathlib.Path(sys.argv[sys.argv.index("--options_file")+1])\ns=[x.split(" = ",1)[1] for x in p.read_text().splitlines() if x.startswith("solution_file")][0]\npathlib.Path(s).write_text("Model status\\nOptimal\\n# Primal solution values\\nFeasible\\nObjective 8.03\\n# Columns 3\\nx 1\\ny 0\\nz .3\\n# Rows 0\\n")\nprint("Model status : Optimal")\n'
            args=self.setup_case(root,body);args[args.index('--limit')+1]='5';args[args.index('--grace')+1]='1'
            p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
            lines=(root/'result/runs.jsonl').read_text().splitlines();self.assertEqual(len(lines),1)
            r=json.loads(lines[0]);self.assertTrue(r['checker']['checked_feasible'],r);self.assertEqual(r['claimed_status'],'optimal')
            self.assertEqual(r['model_sha256'],sha(FIX/'tiny.mps'));self.assertEqual(r['checker_sha256'],sha(CHECK));self.assertEqual(r['runner_sha256'],sha(BIN))
            self.assertNotEqual(subprocess.run(args,capture_output=True).returncode,0)
    def test_watchdog_and_checksum_failure(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td);args=self.setup_case(root,'import time\ntime.sleep(20)\n')
            p=subprocess.run(args,capture_output=True,text=True,timeout=5);self.assertEqual(p.returncode,0,p.stderr)
            r=json.loads((root/'result/runs.jsonl').read_text());self.assertTrue(r['watchdog_timeout']);self.assertFalse(r['checker']['checked_feasible'])
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td);args=self.setup_case(root,'print("hi")\n')
            (root/'fake highs').write_text('changed')
            self.assertNotEqual(subprocess.run(args,capture_output=True).returncode,0)
    def test_seed_zero_rejected(self):
        p=subprocess.run([BIN,'--seed','0'],capture_output=True);self.assertNotEqual(p.returncode,0)
unittest.main()
