#!/usr/bin/env python3
"""Exercise live checking and replay of retained runner evidence."""
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

BENCH, SOLVE, CHECK = [pathlib.Path(arg).resolve() for arg in sys.argv[1:4]]
del sys.argv[1:4]
ROOT = pathlib.Path(__file__).resolve().parents[1]
MODEL = ROOT / 'tests/fixtures/tiny.mps'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Revalidation(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='milp replay ')
        self.addCleanup(tmp.cleanup)
        self.root = pathlib.Path(tmp.name)

    def run_solver(self, relax, solver=SOLVE, checker=CHECK, model=MODEL):
        suite = self.root / 'suite.tsv'
        reference = '8.53' if relax else '8.03'
        suite.write_text('id\tpath\tsha256\treference_status\treference_objective\n'
                         f'tiny\t{model}\t{sha(model)}\t=opt=\t{reference}\n')
        config = self.root / 'solvers.tsv'
        config.write_text(f'solver\tbinary\tsha256\tversion\nlab\t{solver}\t{sha(solver)}\ttest\n')
        original = self.root / 'original'
        subprocess.run([BENCH, '--suite', suite, '--solvers', config, '--checker', checker,
                        '--relax' if relax else '--mip', '--out', original, '--limit', '5'],
                       capture_output=True, check=True)
        self.runs = original.rename(self.root / 'moved')
        self.source = self.runs / 'runs.jsonl'
        return json.loads(self.source.read_bytes())

    def executable(self, name, body, interpreter='/bin/sh'):
        path = self.root / name
        path.write_text(f'#!{interpreter}\n{body}\n')
        path.chmod(0o755)
        return path

    def replay(self, name='revalidated', checker=CHECK, expected_exit=0):
        before = self.source.read_bytes()
        p = subprocess.run([sys.executable, ROOT / 'scripts/revalidate.py', self.runs,
                            '--checker', checker, '--name', name], capture_output=True)
        self.assertEqual(p.returncode, expected_exit, p.stderr)
        self.assertEqual(before, self.source.read_bytes())
        filename = name + '.jsonl'
        subprocess.run([sys.executable, ROOT / 'scripts/report.py', self.runs, '--runs',
                        filename, '--prefix', name + '-report'], check=True, capture_output=True)
        row = json.loads((self.runs / filename).read_bytes())
        summary = json.loads((self.runs / (name + '-report.json')).read_bytes())['summary']['lab']
        return row, summary

    def check_replay(self, original):
        for legacy in (False, True):
            if legacy:
                old = dict(original)
                for key in ('search_route', 'closure_reason', 'bound_origin', 'bound_scope'):
                    old.pop(key, None)
                del old['checker_argv']
                del old['evidence_sha256']
                self.source.write_text(json.dumps(old) + '\n')
            row, summary = self.replay('legacy' if legacy else 'recorded')
            if legacy:
                self.assertIsNone(row['evaluation']['bound_origin'])
                self.assertIsNone(row['evaluation']['reported_relative_gap'])
            elif original.get('mode') == 'mip' and original.get('search_route') == 'native':
                attempt = self.runs / 'tiny/lab/seed-1'
                cli = json.loads((attempt / 'solver.log').read_text().splitlines()[-1])
                report, = json.loads((self.runs / 'recorded-report.json').read_bytes())['attempt_results']
                for record in (cli, original, row, row['evaluation'], report):
                    self.assertEqual([record[k] for k in ('search_route', 'closure_reason', 'bound_origin', 'bound_scope')],
                                     ['native', 'search_exhausted', 'search_exhaustion', 'original_mip'])
            self.assertEqual(original['checker'], row['checker'])
            self.assertEqual(original['checker'], row['original_validation']['checker'])
            self.assertEqual(original['wall_s'], row['wall_s'])
            self.assertEqual(summary['checked_feasible'], int(row['checker']['checked_feasible']))
            self.assertEqual(summary['certified_lp_bounds'], int(row['checker'].get('bound_status') == 'valid'))
            self.assertEqual(summary['checker_errors'], 0)

    def test_integer_witness(self):
        row = self.run_solver(False)
        self.assertEqual(row['checker']['objective'], '803/100')
        self.assertFalse(row['checker']['relaxation'])
        self.check_replay(row)

    def test_fractional_primal_and_dual(self):
        row = self.run_solver(True)
        self.assertEqual(row['checker']['objective'], '853/100')
        self.assertEqual(row['checker']['certified_bound'], '853/100')
        self.assertTrue(row['checker']['relaxation'])
        self.check_replay(row)

    def test_dual_only(self):
        solver = self.executable('dual-producer',
                                 'from pathlib import Path\nimport sys\n'
                                 'dual = Path(sys.argv[sys.argv.index("--dual") + 1])\n'
                                 'dual.write_text("CAP -1.5\\nFIXZ -.1\\n")', '/usr/bin/env python3')
        row = self.run_solver(True, solver)
        self.assertFalse(row['checker']['checked_feasible'])
        self.assertEqual(row['checker']['certified_bound'], '853/100')
        self.assertNotIn('--solution', row['checker_argv'])
        self.check_replay(row)

    def test_infeasibility_proof_survives_relocation(self):
        model = self.root / 'contradiction.mps'
        model.write_text('NAME C\nROWS\n N O\n G A\n L B\nCOLUMNS\n x A 1 B 1\n'
                         'RHS\n rhs A 2 B 1\nENDATA\n')
        original = self.run_solver(True, model=model)
        self.assertEqual(original['claimed_status'], 'infeasible')
        self.assertIn('infeasible.farkas', original['evidence_sha256'])
        row, summary = self.replay()
        self.assertEqual(original['checker'], row['checker'])
        self.assertTrue(row['evaluation']['independently_certified_infeasibility'])
        self.assertFalse(row['evaluation']['checked_feasible'])
        self.assertEqual(summary['certified_infeasibilities'], 1)

    def test_no_evidence(self):
        self.run_solver(False, self.executable('no-files', 'exit 0'))
        row, summary = self.replay()
        self.assertEqual(row['checker']['status'], 'unavailable')
        self.assertEqual(summary['absent_evidence'], 1)
        self.assertEqual(summary['checker_errors'], 0)

    def test_formatted_checker_output_keeps_one_json_line_per_attempt(self):
        pretty = self.executable('pretty-checker', 'import json, subprocess, sys\n'
                                 f'p = subprocess.run([{str(CHECK)!r}, *sys.argv[1:]], capture_output=True)\n'
                                 'print(json.dumps(json.loads(p.stdout), indent=2))\n'
                                 'sys.exit(p.returncode)', '/usr/bin/env python3')
        original = self.run_solver(False, checker=pretty)
        self.assertEqual(len(self.source.read_text().splitlines()), 1)
        self.check_replay(original)

    def test_deleted_dual_cannot_silently_reduce_validation_scope(self):
        self.run_solver(True)
        (self.runs / 'tiny/lab/seed-1/dual.sol').unlink()
        row, summary = self.replay(expected_exit=1)
        self.assertEqual(row['checker']['status'], 'error')
        self.assertEqual(summary['checker_errors'], 1)
        self.assertEqual(summary['absent_evidence'], 0)

    def test_changed_witness_is_reported_even_if_it_is_still_feasible(self):
        original = self.run_solver(False)
        solution = self.runs / 'tiny/lab/seed-1/solution.sol'
        self.assertEqual(sha(solution), original['evidence_sha256']['solution.sol'])
        with solution.open('a') as stream:
            stream.write('\n')
        row, summary = self.replay(expected_exit=1)
        self.assertIn('SHA-256 mismatch', row['checker']['message'])
        self.assertIsNone(row['checker_argv'])
        self.assertEqual(summary['checker_errors'], 1)

    def test_malformed_checker_output_and_exit_disagreement(self):
        bad = self.executable('bad-checker', """echo '{"status":"valid"} trailing'""")
        original = self.run_solver(False, checker=bad)
        self.assertEqual(original['checker']['status'], 'error')
        for name, body, code in [('malformed', None, 0),
                                 ('disagreement', """echo '{"status":"valid"}'; exit 1""", 1),
                                 ('usage', 'exit 64', 64)]:
            if body:
                self.executable('bad-checker', body)
            row, summary = self.replay(name, bad, expected_exit=1)
            self.assertEqual(row['checker_exit_code'], code)
            self.assertEqual(row['checker']['status'], 'error')
            self.assertEqual(summary['checker_errors'], 1)
            self.assertEqual(summary['absent_evidence'], 0)


if __name__ == '__main__':
    unittest.main()
