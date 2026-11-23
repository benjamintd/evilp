#!/usr/bin/env python3
"""Acceptance/scoring contracts, independent of solver speed or log syntax."""
import json,pathlib,sys,tempfile,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
from report import checker_outcome, report
from run_evaluation import assess_evidence, classify
def record(name,solver,status,obj,checked,sense='max'):
    return {'instance':name,'solver':solver,'seed':1,'claimed_status':status,
      'checker_exit_code':0,'reported_dual_bound':'999999999999999','reference_status':'=opt=','reference_objective':'8.03',
      'wall_s':2,'checker_wall_s':.1,'checker':{'status':'valid' if checked else 'invalid' if obj else 'no_solution','checked_feasible':checked,'relaxation':False,'objective':obj,'objective_sense':sense}}
class Report(unittest.TestCase):
    def test_report_gate_and_replay_evidence_agree(self):
        model = {'sha256': '1' * 64, 'reference_status': '=opt=', 'reference_objective': '1'}
        row = dict(record('scope', 'lab', 'optimal', '1', True), mode='mip', relaxation=False,
                   reference_objective='1', model_sha256=model['sha256'], checker_exit_code=0)
        row['checker'].update(bound_status='valid', certified_bound='1', certifies_optimality=True)
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            (root/'environment.json').write_text(json.dumps({'time_limit_s':2, 'seed':1}))
            for scope, code, bound, feasible, proof in [
                    (True, 0, '1', False, False), (False, 1, '1', False, False),
                    (False, 0, '2', True, False), (False, 0, '1', True, True)]:
                row['checker']['relaxation'] = scope
                row['checker']['certified_bound'] = bound
                row['checker_exit_code'] = code
                facts = assess_evidence(row, row, relaxation=False)
                case = classify(model, row)
                self.assertTrue(all(case[k] == v for k, v in facts.items()))
                (root/'runs.jsonl').write_text(json.dumps(row)+'\n')
                summary = report(root)['summary']['lab']
                self.assertEqual(summary['checked_feasible'], int(feasible))
                self.assertEqual(summary['certified_optima'], int(proof))
                self.assertEqual(summary['reached_reference'], int(feasible))
    def test_historical_checked_witness_blocks_a_current_best_known_claim(self):
        from test_benchmark_gate import fixture
        _, current, previous, _ = fixture()
        row, witness = current[0], previous[0]
        row.update(reference_status='=best=', reported_dual_bound='1', checker_wall_s=0.1)
        row['checker']['objective_sense'] = 'min'
        witness.update(claimed_status='time_limit', campaign_id='historical')
        witness['checker'].update(objective='0.9996', objective_sense='min')
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root/'environment.json').write_text(json.dumps({'time_limit_s': 5, 'seed': 1}))
            (root/'runs.jsonl').write_text(json.dumps(row))
            archive = root/'archive.jsonl'
            archive.write_text(json.dumps(witness))
            data = report(root, witness_paths=[archive])
            self.assertEqual(data['attempts'], 1)
            self.assertEqual(data['summary']['lab']['witness_contradictions'], 1)
            self.assertEqual(data['attempt_results'][0]['contradicting_witness']['campaign_id'], 'historical')
            self.assertIn('reported_bound_excludes_checked_witness', (root/'summary.html').read_text())
            self.assertFalse(data['attempt_results'][0]['independently_certified_optimality'])

    def test_unavailable_from_older_evaluators(self):
        for code, elapsed, expected in [(None, 0, 'absent'), (None, 120, 'error'),
                                        (64, .1, 'error'), (4, .1, 'unavailable')]:
            with self.subTest(code=code, elapsed=elapsed):
                row = {'checker': {'status': 'unavailable'}, 'checker_exit_code': code,
                       'checker_wall_s': elapsed}
                self.assertEqual(checker_outcome(row), expected)
    def test_reference_and_failure_accounting(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td)
            rows=[record('a','highs','optimal','803/100',True),record('b','highs','optimal','803/100',False),
                  record('c','scip','time_limit',None,False),record('d','cbc','time_limit','903/100',True)]
            (root/'runs.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
            (root/'environment.json').write_text(json.dumps({'seed':1,'time_limit_s':2}))
            data=report(root,prefix='contract')
            self.assertEqual(data['summary']['highs']['attempts'],2)
            self.assertEqual(data['summary']['highs']['reached_reference'],1)
            self.assertEqual(data['summary']['highs']['invalid_candidates'],1)
            self.assertEqual(data['summary']['scip']['attempts'],1)
            self.assertEqual(data['summary']['scip']['checked_feasible'],0)
            self.assertEqual(data['summary']['cbc']['reference_contradictions'],1)
            self.assertTrue((root/'contract.html').exists())
            self.assertFalse((root/'summary.html').exists())
    def test_lp_reference_match_is_distinct_from_certificate(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td)
            r=record('lp','lab','optimal','803/100',True)
            r['checker'].update({'relaxation':True,'bound_status':'valid','certified_bound':'8','certifies_optimality':False})
            (root/'runs.jsonl').write_text(json.dumps(r)+'\n')
            (root/'environment.json').write_text(json.dumps({'seed':1,'time_limit_s':2,'relaxation':True}))
            data=report(root)
            self.assertTrue(data['relaxation'])
            self.assertEqual(data['summary']['lab']['reached_reference'],1)
            self.assertEqual(data['summary']['lab']['certified_lp_bounds'],1)
            self.assertEqual(data['summary']['lab']['certified_optima'],0)
            self.assertIn('Continuous LP relaxations',(root/'summary.html').read_text())
    def test_solution_status_refinement_does_not_hide_reference_contradiction(self):
        with tempfile.TemporaryDirectory() as td:
            root=pathlib.Path(td)
            r=record('cbc-claim','cbc','infeasible',None,False)
            r.update({'reported_solver_status':'infeasible_or_unbounded','reported_solution_status':'infeasible'})
            (root/'runs.jsonl').write_text(json.dumps(r)+'\n')
            (root/'environment.json').write_text(json.dumps({'seed':1,'time_limit_s':2}))
            data=report(root)
            self.assertEqual(data['summary']['cbc']['reference_contradictions'],1)
            self.assertIn('solution: infeasible',(root/'summary.html').read_text())
unittest.main()
