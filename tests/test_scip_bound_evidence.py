"""Keep the backend bound observable through solve, measurement and replay."""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SOLVE, BENCH, CHECK = map(Path, sys.argv[1:4])
del sys.argv[1:4]
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from run_evaluation import assess_evidence


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class BoundEvidence(unittest.TestCase):
    def check_case(self, root, maximize, status, missing=False):
        source = (ROOT / 'tests/fixtures/near_tied_costs_order.mps').read_text()
        if maximize:
            source = source.replace('ROWS\n', 'OBJSENSE\n MAX\nROWS\n').replace('COST 0.022', 'COST -0.022')
        model = root / 'model.mps'
        model.write_text(source)
        objective = -22000 if maximize else 22000
        expected = None if missing else objective - Fraction(1, 1024)
        env = dict(os.environ, TEST_SCIP_STATUS=status, TEST_SCIP_BOUND='missing' if missing else 'finite')
        suite = root / 'suite.tsv'
        suite.write_text('id\tpath\tsha256\treference_status\treference_objective\n'
                         f'model\t{model}\t{sha(model)}\t=opt=\t{objective}\n')
        solvers = root / 'solvers.tsv'
        solvers.write_text(f'solver\tbinary\tsha256\tversion\nlab\t{SOLVE}\t{sha(SOLVE)}\ttest\n')
        runs = root / 'runs'
        subprocess.run([BENCH, '--mip', '--suite', suite, '--solvers', solvers, '--checker', CHECK,
                        '--out', runs, '--limit', '10'], env=env, check=True, capture_output=True)
        attempt = runs / 'model/lab/seed-1'
        cli, = [json.loads(line) for line in (attempt / 'solver.log').read_text().splitlines()
                if line.startswith('{"schema"')]
        record = json.loads((attempt / 'run.json').read_text())
        self.assertEqual((cli['scip_route'], cli['status'], cli['objective']), ('general', status, objective))
        self.assertEqual(record['reported_solver_status'], status)
        self.assertTrue(record['checker']['checked_feasible'])
        self.assertEqual(Fraction(record['checker']['objective']), objective)

        subprocess.run([sys.executable, ROOT / 'scripts/revalidate.py', runs, '--checker', CHECK],
                       check=True, capture_output=True)
        subprocess.run([sys.executable, ROOT / 'scripts/report.py', runs, '--runs', 'revalidated.jsonl'],
                       check=True, capture_output=True)
        replay = json.loads((runs / 'revalidated.jsonl').read_text())
        report, = json.loads((runs / 'summary.json').read_text())['attempt_results']
        for bound in [cli['best_bound'], record['reported_dual_bound'], replay['reported_dual_bound'],
                      replay['evaluation']['reported_bound'], report['reported_bound']]:
            self.assertEqual(None if bound is None else Fraction(bound), expected)
        self.assertEqual(report['checked_objective'], str(objective))
        self.assertEqual(report['claimed_status'], status)
        self.assertFalse(report['independently_certified_optimality'])
        for observation in (cli, record, replay, replay['evaluation'], report):
            self.assertEqual(observation['search_route'], 'general_scip')
            self.assertEqual(observation['closure_reason'], 'backend_status' if status == 'optimal' else 'none')
            self.assertEqual(observation['bound_origin'], 'none' if missing else 'backend_dual')
            self.assertEqual(observation['bound_scope'], 'none' if missing else 'original_mip')
        expected_gap = str((objective - expected) / expected) if not missing and not maximize else None
        self.assertEqual(report['reported_relative_gap'], expected_gap)
        conflicts = ['reported_bound_excludes_checked_witness'] if maximize and not missing else []
        self.assertEqual(report['witness_conflicts'], conflicts)
        self.assertEqual(replay['evaluation']['witness_conflicts'], conflicts)
        if not missing and not maximize:
            self.assertGreater(Fraction(report['checked_objective']) - Fraction(report['reported_bound']), 0)
            for damage in ({'watchdog_timeout': True}, {'signal': 9}, {'exit_code': 1},
                           {'status_conflict': True}, {'bound_origin': None}):
                facts = assess_evidence({**record, **damage}, {}, relaxation=False)
                self.assertIsNone(facts['reported_relative_gap'], damage)

    def test_public_bound(self):
        for maximize, status, missing in [(False, 'optimal', False), (True, 'optimal', False),
                                          (True, 'time_limit', False), (False, 'time_limit', True)]:
            with self.subTest(maximize=maximize, status=status, missing=missing), tempfile.TemporaryDirectory() as tmp:
                self.check_case(Path(tmp), maximize, status, missing)

    def test_ill_conditioned_limit_has_no_global_bound(self):
        source = ROOT / 'tests/fixtures/tiny_row_pricing.mps'
        with tempfile.TemporaryDirectory() as tmp:
            sol = Path(tmp) / 'point.sol'
            run = subprocess.run([SOLVE, source, '--mip', '--solution', sol], capture_output=True,
                                 env={**os.environ, 'TEST_SCIP_STATUS': 'time_limit'})
            self.assertEqual(run.returncode, 4, run.stderr)
            result = json.loads(run.stdout)
            self.assertEqual((result['status'], result['search_route']), ('time_limit', 'general_scip'))
            self.assertIsNone(result['best_bound'])
            self.assertEqual((result['bound_origin'], result['bound_scope']), ('none', 'none'))
            self.assertTrue(result['has_primal'])
            check = json.loads(subprocess.check_output([CHECK, '--model', source, '--solution', sol, '--exact']))
            self.assertGreaterEqual(Fraction(check['objective']), Fraction(1, 10**16))


if __name__ == '__main__':
    unittest.main()
