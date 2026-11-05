#!/usr/bin/env python3
"""Acceptance/scoring contracts, independent of solver speed or log syntax."""
import json,pathlib,sys,tempfile,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
from report import report
def record(name,solver,status,obj,checked,sense='max'):
    return {'instance':name,'solver':solver,'seed':1,'claimed_status':status,
      'reported_dual_bound':'999999999999999','reference_status':'=opt=','reference_objective':'8.03',
      'wall_s':2,'checker_wall_s':.1,'checker':{'status':'valid' if checked else 'invalid' if obj else 'no_solution','checked_feasible':checked,'objective':obj,'objective_sense':sense}}
class Report(unittest.TestCase):
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
            r['checker'].update({'bound_status':'valid','certified_bound':'8','certifies_optimality':False})
            (root/'runs.jsonl').write_text(json.dumps(r)+'\n')
            (root/'environment.json').write_text(json.dumps({'seed':1,'time_limit_s':2,'relaxation':True}))
            data=report(root)
            self.assertTrue(data['relaxation'])
            self.assertEqual(data['summary']['lab']['reached_reference'],1)
            self.assertEqual(data['summary']['lab']['certified_lp_bounds'],1)
            self.assertEqual(data['summary']['lab']['certified_optima'],0)
            self.assertIn('Continuous LP relaxations',(root/'summary.html').read_text())
unittest.main()
