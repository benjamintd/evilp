#!/usr/bin/env python3
import hashlib,json,os,pathlib,signal,subprocess,sys,tempfile,time,unittest
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
            self.assertEqual(r['reported_solver_status'],'optimal')
            self.assertFalse(r['relaxation']);self.assertIsNone(r['dual_file'])
            self.assertFalse(json.loads((root/'result/environment.json').read_text())['relaxation'])
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
    def test_highs_nonzero_exit_preserves_report_without_overriding_process_error(self):
        for report,expected in [('Model status : Optimal','optimal'),('No terminal report','unknown')]:
            with tempfile.TemporaryDirectory() as td:
                root=pathlib.Path(td)
                body='import sys\nprint('+repr(report)+')\nprint("WARNING: Row names are not present")\nsys.exit(1)\n'
                args=self.setup_case(root,body)
                args[args.index('--limit')+1]='5';args[args.index('--grace')+1]='1'
                p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
                r=json.loads((root/'result/runs.jsonl').read_text())
                self.assertEqual(r['claimed_status'],'process_error',r)
                self.assertEqual(r['reported_solver_status'],expected,r);self.assertEqual(r['exit_code'],1)
    def test_interrupt_terminates_active_run(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td);args=self.setup_case(root,'import time\ntime.sleep(20)\n')
            args[args.index('--limit')+1]='20'
            p=subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
            deadline=time.monotonic()+5
            log=root/'result/tiny/highs/seed-1/solver.log'
            while not log.exists() and time.monotonic()<deadline:time.sleep(.02)
            self.assertTrue(log.exists())
            p.send_signal(signal.SIGINT);out,err=p.communicate(timeout=5)
            self.assertEqual(p.returncode,130,err)
            r=json.loads((root/'result/runs.jsonl').read_text());self.assertEqual(r['claimed_status'],'interrupted')
    def test_cbc_max_objective_flag(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td)
            body='import pathlib,sys\nassert "-max" in sys.argv\np=pathlib.Path(sys.argv[sys.argv.index("-solution")+1])\np.write_text("Optimal - objective value 8.03\\n0 x 1 0\\n1 y 0 0\\n2 z .3 0\\n")\nprint("Result - Optimal solution found")\n'
            args=self.setup_case(root,body);args[args.index('--limit')+1]='5';args[args.index('--grace')+1]='1'
            cfg=root/'solvers.tsv';cfg.write_text(cfg.read_text().replace('\nhighs\t','\ncbc\t'))
            p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
            r=json.loads((root/'result/runs.jsonl').read_text())
            self.assertEqual(r['checker']['objective'],'803/100',r);self.assertIn('-max',r['argv'])
    def test_ambiguous_and_lp_statuses_are_not_milp_certificates(self):
        cases=[('highs','Solving report\\n  Status            Unbounded or infeasible','infeasible_or_unbounded'),
               ('scip','SCIP Status : problem is solved [infeasible or unbounded]','infeasible_or_unbounded'),
               ('cbc','Result - Linear relaxation unbounded','lp_unbounded')]
        for solver,message,status in cases:
            with tempfile.TemporaryDirectory() as td:
                root=pathlib.Path(td);args=self.setup_case(root,'print("'+message+'")\n')
                args[args.index('--limit')+1]='5';args[args.index('--grace')+1]='1'
                cfg=root/'solvers.tsv';cfg.write_text(cfg.read_text().replace('\nhighs\t','\n'+solver+'\t'))
                p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
                r=json.loads((root/'result/runs.jsonl').read_text())
                self.assertEqual(r['claimed_status'],status,r);self.assertFalse(r['checker']['checked_feasible'])
    def test_cbc_continuous_terminal_status_and_mip_root_separation(self):
        cases=[('Optimal - objective value 0.33333333','optimal'),
               ('Optimal objective 0.3333333333 - 0 iterations time 0.002, Presolve 0.00','optimal'),
               ('Optimal - objective value unavailable','unknown'),
               ('Optimal - objective value 0.33333333\nCbc0010I After 0 nodes, 1 on tree','unknown'),
               ('Continuous objective value is 0.333333\nOptimal objective 0.3333333333 - 0 iterations time 0.002','unknown'),
               ('Optimal - objective value 0.33333333\nResult - Stopped on time limit','time_limit')]
        for log,status in cases:
            with tempfile.TemporaryDirectory() as td:
                root=pathlib.Path(td);args=self.setup_case(root,'print('+repr(log)+')\n')
                args[args.index('--limit')+1]='5';args[args.index('--grace')+1]='1'
                cfg=root/'solvers.tsv';cfg.write_text(cfg.read_text().replace('\nhighs\t','\ncbc\t'))
                p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
                r=json.loads((root/'result/runs.jsonl').read_text());self.assertEqual(r['claimed_status'],status,r)
    def lab_case(self,root,body):
        args=self.setup_case(root,body)
        args[args.index('--limit')+1]='5';args[args.index('--grace')+1]='1'
        cfg=root/'solvers.tsv';cfg.write_text(cfg.read_text().replace('\nhighs\t','\nlab\t'))
        suite=root/'suite.tsv';suite.write_text(suite.read_text().replace('\t8.03\n','\t8.53\n'))
        return args
    def test_lab_relaxation_solution_and_dual_validation(self):
        with tempfile.TemporaryDirectory(prefix='milp lab LP ') as td:
            root=pathlib.Path(td)
            body='import json,pathlib,sys\nassert "--relax" in sys.argv\nassert "--time-limit" in sys.argv\nassert not any(a in ("--seed","--random-seed") for a in sys.argv[1:])\ns=pathlib.Path(sys.argv[sys.argv.index("--solution")+1])\nd=pathlib.Path(sys.argv[sys.argv.index("--dual")+1])\ns.write_text("=obj= 8.53\\nx .5\\ny 1\\nz .3\\n")\nd.write_text("CAP -1.5\\nFIXZ -.1\\n")\nprint(json.dumps({"message":"a string containing \\\"status\\\":\\\"infeasible\\\"", "nested":{"status":"unbounded"},"status":"optimal"}))\n'
            args=self.lab_case(root,body)+['--relax']
            p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
            r=json.loads((root/'result/runs.jsonl').read_text())
            self.assertEqual(r['solver'],'lab');self.assertTrue(r['relaxation'])
            self.assertEqual(r['claimed_status'],'optimal',r)
            self.assertTrue(r['checker']['checked_feasible'],r)
            self.assertTrue(r['checker']['relaxation']);self.assertEqual(r['checker']['objective'],'853/100')
            self.assertEqual(r['checker']['bound_status'],'valid',r)
            self.assertTrue(pathlib.Path(r['dual_file']).exists())
            self.assertEqual(r['argv'][2],'--relax');self.assertNotIn('--seed',r['argv'])
            self.assertTrue(json.loads((root/'result/environment.json').read_text())['relaxation'])
            self.assertIn('reference_objectives=LP',(root/'result/tiny/lab/seed-1/settings.txt').read_text())
    def test_lab_mode_requires_explicit_relaxation_and_excludes_native_adapters(self):
        for mode in ['missing_relax','native_relax','mixed_relax']:
            with tempfile.TemporaryDirectory() as td:
                root=pathlib.Path(td)
                args=self.lab_case(root,'print("should not run")\n')
                cfg=root/'solvers.tsv'
                if mode=='native_relax':
                    cfg.write_text(cfg.read_text().replace('\nlab\t','\nhighs\t'))
                elif mode=='mixed_relax':
                    cfg.write_text(cfg.read_text()+cfg.read_text().splitlines()[1].replace('lab\t','highs\t',1)+'\n')
                if mode!='missing_relax':args+=['--relax']
                p=subprocess.run(args,capture_output=True,text=True)
                self.assertEqual(p.returncode,2,p.stderr);self.assertIn('--relax',p.stderr)
                self.assertFalse((root/'result').exists())
    def test_lab_structured_nonzero_statuses_are_preserved(self):
        cases=[('numerical_failure',4),('unsupported',3),('input_error',2),('time_limit',4),
               ('iteration_limit',4),('infeasible',0),('unbounded',0)]
        for status,code in cases:
            with tempfile.TemporaryDirectory() as td:
                root=pathlib.Path(td)
                body='import json,sys\nprint(json.dumps({"status":'+repr(status)+',"message":"optimal"}))\nsys.exit('+str(code)+')\n'
                args=self.lab_case(root,body)+['--relax']
                p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
                r=json.loads((root/'result/runs.jsonl').read_text())
                self.assertEqual(r['claimed_status'],status,r);self.assertEqual(r['exit_code'],code)
                self.assertFalse(r['checker']['checked_feasible']);self.assertIsNone(r['dual_file'])
    def test_lab_signal_overrides_previously_printed_status(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td)
            args=self.lab_case(root,'import os,signal\nprint(\'{"status":"optimal"}\',flush=True)\nos.kill(os.getpid(),signal.SIGTERM)\n')+['--relax']
            p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
            r=json.loads((root/'result/runs.jsonl').read_text())
            self.assertEqual(r['claimed_status'],'process_error',r);self.assertEqual(r['signal'],signal.SIGTERM)
            self.assertEqual(r['reported_solver_status'],'optimal',r)
    def test_lab_does_not_match_status_words_in_messages_or_malformed_json(self):
        messages=['{"message":"optimal", "nested":{"status":"optimal"}}',
                  '{"status":"optimal"} trailing', '{"status":"optimal","status":"infeasible"}',
                  '{"status":"optimal\\u0000garbage"}', '{"status":"almost_optimal"}']
        for text in messages:
            with tempfile.TemporaryDirectory() as td:
                root=pathlib.Path(td);args=self.lab_case(root,'print('+repr(text)+')\n')+['--relax']
                p=subprocess.run(args,capture_output=True,text=True);self.assertEqual(p.returncode,0,p.stderr)
                r=json.loads((root/'result/runs.jsonl').read_text());self.assertEqual(r['claimed_status'],'unknown',r)
unittest.main()
