#!/usr/bin/env python3
"""Check exported candidates and replay local LP proofs across independent tools."""
import copy
import gzip
import json
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

SOLVER = str(pathlib.Path(sys.argv.pop(1)).resolve())
CHECKER = str(pathlib.Path(sys.argv.pop(1)).resolve())
REPLAY = pathlib.Path(__file__).resolve().parents[1] / 'scripts/check_trace_certificates.py'


class FarkasPipeline(unittest.TestCase):
    def test_large_source_enclosure_export_and_typed_trace_margin(self):
        # Original-decimal .1*x>=.3 and .1*x<=.2 contradict each other.
        # Additional fixed rows exceed retained small-model provenance and
        # dense reconstruction limits without changing that independent fact.
        lines = ['NAME ENCLOSED', 'ROWS', ' N O', ' G A', ' L B']
        lines += [f' E F{i}' for i in range(127)]
        lines += ['COLUMNS', ' x O 0 A .1', ' x B .1']
        lines += [f' z{i} O 0 F{i} 1' for i in range(127)]
        lines += ['RHS', ' b A .3 B .2', 'BOUNDS', ' UP b x 3']
        lines += [f' FX b z{i} 0' for i in range(127)]
        model = '\n'.join(lines + ['ENDATA', ''])
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path, cert, trace = root/'model.mps', root/'proof.farkas', root/'trace.jsonl'
            path.write_text(model)
            p = subprocess.run([SOLVER, str(path), '--relax', '--farkas', str(cert)],
                               text=True, capture_output=True, timeout=30)
            self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
            result = json.loads(p.stdout)
            self.assertEqual(result['status'], 'infeasible', result)
            self.assertFalse(result['exact_model_available'], result)
            self.assertEqual(result['exact_infeasibilities'], 0, result)
            self.assertEqual(result['source_infeasibilities'], 1, result)
            self.assertEqual(result['farkas_margin_kind'], 'source_enclosure')
            self.assertFalse(result['certifies_infeasibility'])
            q = subprocess.run([CHECKER, '--model', str(path), '--farkas', str(cert)],
                               text=True, capture_output=True, timeout=30)
            self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
            actual = F(json.loads(q.stdout)['contradiction_margin'])
            self.assertGreater(actual, 0)
            p = subprocess.run([SOLVER, str(path), '--mip', '--no-presolve', '--trace', str(trace)],
                               text=True, capture_output=True, timeout=30)
            self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
            self.assertEqual(json.loads(p.stdout)['source_infeasibilities'], 1)
            events = [json.loads(line) for line in trace.read_text().splitlines()]
            proofs = [e for e in events if e['event'] == 'lp_infeasibility_certificate']
            self.assertEqual(len(proofs), 1)
            event = proofs[0]
            self.assertEqual(event['margin_kind'], 'source_enclosure')
            self.assertNotIn('contradiction_margin', event)
            lower = F(event['contradiction_margin_lower_bound'])
            self.assertGreater(lower, 0)
            self.assertLess(lower, actual)
            for i, variant in enumerate(['valid', 'too_large', 'zero', 'ambiguous', 'unknown']):
                changed = copy.deepcopy(event)
                if variant == 'too_large':
                    changed['contradiction_margin_lower_bound'] = str(actual + 1)
                elif variant == 'zero':
                    changed['contradiction_margin_lower_bound'] = '0'
                elif variant == 'ambiguous':
                    changed['contradiction_margin'] = str(actual)
                elif variant == 'unknown':
                    changed['margin_kind'] = 'unknown'
                replay_trace = root/f'replay-{i}.jsonl'
                replay_trace.write_text(json.dumps(changed)+'\n')
                p = subprocess.run([sys.executable, str(REPLAY), '--model', str(path),
                                    '--trace', str(replay_trace), '--checker', CHECKER,
                                    '--out', str(root/f'replay-{i}')],
                                   text=True, capture_output=True, timeout=60)
                self.assertEqual(p.returncode, 0 if variant == 'valid' else 1,
                                 (variant, p.stdout, p.stderr))
                if variant == 'valid':
                    summary = json.loads(p.stdout)
                    self.assertTrue(summary['proofs'][0]['margin_lower_bound_valid'])
                    self.assertFalse(summary['proofs'][0]['margin_matches'])
            disabled = subprocess.run([SOLVER, str(path), '--relax', '--no-exact'],
                                      text=True, capture_output=True, timeout=30)
            result = json.loads(disabled.stdout)
            self.assertEqual(result['source_infeasibility_attempts'], 0, result)
            self.assertFalse(result['has_farkas_candidate'], result)

    @staticmethod
    def free_model(a='3', b='7', lower='3', upper='1', duplicate=False):
        columns = (' x O 0 A .1\n x A .2 B .3\n' if duplicate else
                   f' x O 0 A {a}\n x B {b}\n')
        return ('NAME FREE\nROWS\n N O\n L A\n G B\nCOLUMNS\n' + columns +
                f'RHS\n b A {upper} B {lower}\nBOUNDS\n FR b x\nENDATA\n')

    @staticmethod
    def box_model(scale='1', fixed=False):
        a = F(scale)
        # These source numbers terminate exactly in decimal in every test case.
        decimal = lambda value: str(float(value))
        return ('NAME BOX\nROWS\n N O\n G R\nCOLUMNS\n'
                f' x O 0 R {decimal(3*a)}\n' +
                (f' z O 0 R {decimal(7*a)}\n' if fixed else '') +
                f'RHS\n b R {decimal((18 if fixed else 4)*a)}\nBOUNDS\n'
                ' LO b x -2\n UP b x 1\n' + (' FX b z 2\n' if fixed else '') + 'ENDATA\n')

    def export(self, root, model, extra=(), expected=0):
        root.mkdir()
        path, cert = root / 'model.mps', root / 'certificate.farkas'
        path.write_text(model)
        p = subprocess.run([SOLVER, str(path), '--relax', '--farkas', str(cert), *extra],
                           text=True, capture_output=True, timeout=20)
        self.assertEqual(p.returncode, expected, (p.stdout, p.stderr))
        result = json.loads(p.stdout)
        self.assertFalse(result['certifies_infeasibility'])
        self.assertTrue(result['solver_status_is_claim'])
        if expected == 0 and '--no-exact' not in extra:
            self.assertEqual(result['status'], 'infeasible', result)
            self.assertTrue(result['has_farkas_candidate'], result)
            self.assertGreater(result['exact_infeasibilities'], 0, result)
            self.assertTrue(cert.read_text().startswith('MILPLAB_FARKAS_1\n'))
            q = subprocess.run([CHECKER, '--model', str(path), '--farkas', str(cert)],
                               text=True, capture_output=True, timeout=20)
            self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
            checked = json.loads(q.stdout)
            self.assertTrue(checked['certifies_infeasibility'], checked)
            self.assertGreater(F(checked['contradiction_margin']), 0)
            self.assertFalse(checked['checked_feasible'])
            self.assertFalse(checked['certifies_optimality'])
        return result, path, cert

    def test_free_rational_ray_and_duplicate_original_decimals(self):
        models = [self.free_model(), self.free_model(lower='.6', upper='.3', duplicate=True)]
        with tempfile.TemporaryDirectory() as td:
            for i, model in enumerate(models):
                self.export(pathlib.Path(td) / str(i), model)

    def test_finite_boxes_fixed_variables_and_numerical_controls(self):
        controls = [(), ('--backend', 'dense'), ('--no-scaling',), ('--no-refinement',),
                    ('--bland-ratio',), ('--backend', 'dense', '--no-scaling', '--no-refinement')]
        with tempfile.TemporaryDirectory() as td:
            for i, extra in enumerate(controls):
                model = self.box_model('1e6' if i == 0 else '1', fixed=bool(i % 2))
                result, _, _ = self.export(pathlib.Path(td) / str(i), model, extra)
                self.assertEqual(result['scaling_enabled'], '--no-scaling' not in extra)
                self.assertEqual(result['refinement_enabled'], '--no-refinement' not in extra)
                self.assertEqual(result['backend'], 'dense_lu_reference'
                                 if 'dense' in extra else 'sparse_lu')

    def test_generated_original_decimal_strict_contradictions(self):
        rng = random.Random(20261108)
        with tempfile.TemporaryDirectory() as td:
            for case in range(36):
                a, b = rng.randint(1, 20), rng.randint(1, 20)
                # x <= t and x >= t + gap; decimal tokens stay exact.
                t = F(rng.randint(-8, 8), 4)
                gap = F(rng.randint(1, 8), 4)
                model = self.free_model(str(a), str(b), str(float(b*(t+gap))), str(float(a*t)))
                self.export(pathlib.Path(td) / str(case), model,
                            ('--backend', 'dense') if case % 2 else ())

    def test_no_exact_disables_candidate_export(self):
        with tempfile.TemporaryDirectory() as td:
            for i, model in enumerate([self.free_model(), self.box_model()]):
                result, _, cert = self.export(pathlib.Path(td) / str(i), model, ['--no-exact'])
                self.assertFalse(result['exact_enabled'])
                self.assertEqual(result['exact_attempts'], 0)
                self.assertFalse(result['has_farkas_candidate'])
                self.assertFalse(cert.exists())

    def test_existing_export_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path, cert = root / 'model.mps', root / 'proof.farkas'
            path.write_text(self.free_model())
            cert.write_text('historical proof\n')
            p = subprocess.run([SOLVER, str(path), '--relax', '--farkas', str(cert)],
                               text=True, capture_output=True, timeout=20)
            result = json.loads(p.stdout)
            self.assertNotEqual(p.returncode, 0, result)
            self.assertEqual(result['status'], 'input_error')
            self.assertEqual(cert.read_text(), 'historical proof\n')

    def test_constant_row_presolve_claim_requires_independent_proof(self):
        # Constant-row presolve may return before ray reconstruction. A status
        # remains a claim; a separately supplied proof still verifies exactly.
        model = ('NAME CONSTANT\nROWS\n N O\n G R\nCOLUMNS\n x O 0 R 3\n'
                 'RHS\n b R 4\nBOUNDS\n FX b x 1\nENDATA\n')
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path, candidate, independent = root / 'model.mps', root / 'candidate', root / 'independent'
            path.write_text(model)
            p = subprocess.run([SOLVER, str(path), '--relax', '--farkas', str(candidate)],
                               text=True, capture_output=True, timeout=20)
            self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
            result = json.loads(p.stdout)
            self.assertEqual(result['status'], 'infeasible')
            self.assertFalse(result['certifies_infeasibility'])
            self.assertEqual(candidate.exists(), result['has_farkas_candidate'])
            independent.write_text('MILPLAB_FARKAS_1\nR 1\n')
            for certificate in [independent] + ([candidate] if candidate.exists() else []):
                q = subprocess.run([CHECKER, '--model', str(path), '--farkas', str(certificate)],
                                   text=True, capture_output=True, timeout=20)
                self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
                self.assertTrue(json.loads(q.stdout)['certifies_infeasibility'])

    def test_feasible_near_proportional_decimal_rows_have_no_proof(self):
        # In original decimal arithmetic, .1 + .2 == .3 and x=1 is feasible.
        model = self.free_model(lower='.3', upper='.3', duplicate=True)
        controls = [(), ('--backend', 'dense'), ('--no-scaling',), ('--bland-ratio',)]
        with tempfile.TemporaryDirectory() as td:
            for i, extra in enumerate(controls):
                root = pathlib.Path(td) / str(i)
                root.mkdir()
                path, cert, witness = root / 'model.mps', root / 'proof.farkas', root / 'exact.sol'
                path.write_text(model)
                witness.write_text('x 1\n')
                q = subprocess.run([CHECKER, '--model', str(path), '--solution', str(witness), '--exact'],
                                   text=True, capture_output=True, timeout=20)
                self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
                self.assertTrue(json.loads(q.stdout)['checked_feasible'])
                p = subprocess.run([SOLVER, str(path), '--relax', '--farkas', str(cert), *extra],
                                   text=True, capture_output=True, timeout=20)
                self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
                result = json.loads(p.stdout)
                self.assertFalse(result['has_farkas_candidate'], result)
                self.assertFalse(result['certifies_infeasibility'])
                self.assertFalse(cert.exists())

    @staticmethod
    def replay_fixture():
        # z only occurs in BOUNDS; it still belongs to the complete source box.
        model = ('NAME TRACE\nOBJSENSE\n MAX\nROWS\n N O\n L R\nCOLUMNS\n'
                 ' x O 17 R .1\n x R .2\nRHS\n b O 7 R .9\nRANGES\n g R .6\n'
                 'BOUNDS\n LO b x 0\n UP b x 3\n FR b z\nENDATA\n')
        event = {'event': 'lp_infeasibility_certificate', 'lp_solve': 2,
                 'bounds': [{'column': 'x', 'lower': '-1/8', 'upper': '1/2'},
                            {'column': 'z', 'lower': None, 'upper': None}],
                 'multipliers': {'R': '1'}, 'contradiction_margin': '3/20'}
        return model, event

    def replay(self, root, model, events, extra=(), expected=0, compressed=False):
        root.mkdir()
        model_path = root / ('source.mps.gz' if compressed else 'source.mps')
        model_path.write_bytes(gzip.compress(model.encode()) if compressed else model.encode())
        trace, output = root / 'trace.jsonl', root / 'checked'
        trace.write_text(events if isinstance(events, str) else
                         ''.join(json.dumps(event) + '\n' for event in events))
        p = subprocess.run([sys.executable, str(REPLAY), '--model', str(model_path),
                            '--trace', str(trace), '--checker', CHECKER, '--out', str(output), *extra],
                           text=True, capture_output=True, timeout=30)
        self.assertEqual(p.returncode, expected, (p.stdout, p.stderr))
        summary = json.loads((output / 'summary.json').read_text())
        self.assertFalse(summary['certifies_tree_infeasibility'])
        self.assertIn('not a branch-and-bound tree proof', summary['scope'])
        self.assertEqual(summary, json.loads(p.stdout))
        return summary, output

    def test_trace_replay_preserves_original_matrix_ranges_and_exact_bounds(self):
        model, event = self.replay_fixture()
        other = copy.deepcopy(event)
        other['lp_solve'], other['bounds'][0]['lower'], other['bounds'][0]['upper'] = 3, '1/4', '1/4'
        other['contradiction_margin'] = '9/40'
        with tempfile.TemporaryDirectory() as td:
            summary, output = self.replay(pathlib.Path(td) / 'run', model,
                                           [{'event': 'unrelated'}, event, other], compressed=True)
            self.assertTrue(summary['all_lp_certificates_valid'])
            self.assertEqual(summary['proof_count'], 2)
            self.assertEqual(summary['valid_proof_count'], 2)
            self.assertEqual(summary['source_columns'], 2)
            for record in summary['proofs']:
                self.assertTrue(record['independently_valid'])
                self.assertTrue(record['margin_matches'])
                directory = output / record['directory']
                rewritten = (directory / 'model.mps').read_text()
                self.assertEqual(rewritten.split('BOUNDS\n')[0], model.split('BOUNDS\n')[0])
                self.assertIn(' FR replay z\n', rewritten)
                for name in ['certificate.farkas', 'event.json', 'checker.json',
                             'checker.stdout.txt', 'checker.stderr.txt']:
                    self.assertTrue((directory / name).exists(), name)
            first = (output / summary['proofs'][0]['directory'] / 'model.mps').read_text()
            self.assertIn(' LO replay x -0.125\n UP replay x 0.5\n', first)

    def test_trace_box_must_match_every_source_column(self):
        model, event = self.replay_fixture()
        mutations = []
        for bounds in [event['bounds'][:-1], event['bounds'] + [event['bounds'][0]],
                       event['bounds'] + [{'column': 'extra', 'lower': None, 'upper': None}]]:
            bad = copy.deepcopy(event)
            bad['bounds'] = bounds
            mutations.append(bad)
        with tempfile.TemporaryDirectory() as td:
            for i, bad in enumerate(mutations):
                summary, _ = self.replay(pathlib.Path(td) / str(i), model, [bad], expected=1)
                self.assertFalse(summary['all_lp_certificates_valid'])
                self.assertFalse(summary['proofs'][0]['independently_valid'])
                self.assertIsNotNone(summary['proofs'][0]['error'])

    def test_trace_nonterminating_bounds_and_malformed_rationals_rejected(self):
        model, event = self.replay_fixture()
        with tempfile.TemporaryDirectory() as td:
            for i, value in enumerate(['1/3', 'nan', '1/0', '2/-3', 0.5]):
                bad = copy.deepcopy(event)
                bad['bounds'][0]['upper'] = value
                summary, _ = self.replay(pathlib.Path(td) / str(i), model, [bad], expected=1)
                self.assertFalse(summary['proofs'][0]['independently_valid'])

    def test_trace_wrong_margins_and_invalid_proofs_cannot_pass(self):
        model, event = self.replay_fixture()
        cases = []
        for key, value in [('contradiction_margin', '1'), ('multipliers', {'R': '-1'}),
                           ('multipliers', {'UNKNOWN': '1'}), ('multipliers', {'R': '0'})]:
            bad = copy.deepcopy(event)
            bad[key] = value
            cases.append(bad)
        with tempfile.TemporaryDirectory() as td:
            summary, output = self.replay(pathlib.Path(td) / 'run', model, cases, expected=1)
            self.assertEqual(summary['proof_count'], 4)
            self.assertEqual(summary['valid_proof_count'], 0)
            for proof in summary['proofs']:
                self.assertIsNotNone(proof['error'])
                self.assertTrue((output / proof['directory'] / 'checker.stdout.txt').exists())

    def test_trace_duplicate_json_keys_fail_closed(self):
        model, event = self.replay_fixture()
        text = json.dumps(event).replace('"multipliers": {"R": "1"}',
                                         '"multipliers": {"R": "1", "R": "0"}')
        with tempfile.TemporaryDirectory() as td:
            summary, _ = self.replay(pathlib.Path(td) / 'run', model, text, expected=1)
            self.assertEqual(summary['proof_count'], 0)
            self.assertTrue(any('duplicate JSON key' in error for error in summary['errors']))

    def test_malformed_nonfinite_trace_records_fail_closed(self):
        model, _ = self.replay_fixture()
        with tempfile.TemporaryDirectory() as td:
            for i, text in enumerate(['{', '[]', '{"value":NaN}', '{"value":1e9999}']):
                summary, _ = self.replay(pathlib.Path(td) / str(i), model, text, expected=1)
                self.assertEqual(summary['proof_count'], 0)
                self.assertFalse(summary['all_lp_certificates_valid'])
                self.assertTrue(summary['errors'])

    def test_trace_empty_requires_explicit_allow_empty(self):
        model, _ = self.replay_fixture()
        with tempfile.TemporaryDirectory() as td:
            summary, _ = self.replay(pathlib.Path(td) / 'required', model, [], expected=1)
            self.assertFalse(summary['all_lp_certificates_valid'])
            summary, _ = self.replay(pathlib.Path(td) / 'optional', model, [], ['--allow-empty'])
            self.assertTrue(summary['all_lp_certificates_valid'])
            self.assertEqual(summary['proof_count'], 0)

    def test_replay_output_directory_must_be_fresh(self):
        model, event = self.replay_fixture()
        with tempfile.TemporaryDirectory() as td:
            _, output = self.replay(pathlib.Path(td) / 'run', model, [event])
            snapshot = {path.relative_to(output): path.read_bytes()
                        for path in output.rglob('*') if path.is_file()}
            p = subprocess.run([sys.executable, str(REPLAY), '--model', str(output.parent / 'source.mps'),
                                '--trace', str(output.parent / 'trace.jsonl'), '--checker', CHECKER,
                                '--out', str(output)], text=True, capture_output=True, timeout=20)
            self.assertEqual(p.returncode, 2, (p.stdout, p.stderr))
            self.assertEqual(snapshot, {path.relative_to(output): path.read_bytes()
                                       for path in output.rglob('*') if path.is_file()})

    def test_actual_milp_trace_certifies_local_lp_boxes_only(self):
        model = ('NAME HALF\nROWS\n N O\n E R\nCOLUMNS\n'
                 " M0 'MARKER' 'INTORG'\n x O 1 R 2\n y O 1 R 2\n M1 'MARKER' 'INTEND'\n"
                 'RHS\n b R 1\nBOUNDS\n BV b x\n BV b y\nENDATA\n')
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path, trace = root / 'source.mps', root / 'trace.jsonl'
            path.write_text(model)
            default = subprocess.run([SOLVER, str(path), '--mip'],
                                     text=True, capture_output=True, timeout=20)
            self.assertEqual(default.returncode, 0, (default.stdout, default.stderr))
            self.assertEqual(json.loads(default.stdout)['status'], 'infeasible')
            # Root propagation now closes this integer model. Keep this same
            # model, explicitly disable it to exercise local branch proofs.
            p = subprocess.run([SOLVER, str(path), '--mip', '--no-presolve', '--trace', str(trace)],
                               text=True, capture_output=True, timeout=20)
            self.assertEqual(p.returncode, 0, (p.stdout, p.stderr))
            self.assertEqual(json.loads(p.stdout)['status'], 'infeasible')
            events = [json.loads(line) for line in trace.read_text().splitlines()]
            proofs = [event for event in events if event.get('event') == 'lp_infeasibility_certificate']
            self.assertGreaterEqual(len(proofs), 1, events)
            summary, _ = self.replay(root / 'replay', model, events)
            self.assertEqual(summary['valid_proof_count'], len(proofs))
            # The root LP has the exact witness (x,y)=(1/2,0), even though its integer
            # branches are contradictory. Local LP proofs do not certify it.
            witness = root / 'root.sol'
            witness.write_text('x 1/2\ny 0\n')
            q = subprocess.run([CHECKER, '--model', str(path), '--solution', str(witness),
                                '--format', 'rational', '--relax', '--exact'],
                               text=True, capture_output=True, timeout=20)
            self.assertEqual(q.returncode, 0, (q.stdout, q.stderr))
            self.assertTrue(json.loads(q.stdout)['checked_feasible'])


if __name__ == '__main__':
    unittest.main()
