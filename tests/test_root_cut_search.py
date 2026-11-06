#!/usr/bin/env python3
"""Root cuts in actual search, with enumeration and independent proof replay."""
import argparse
import itertools
import json
from pathlib import Path
import random
import subprocess
import tempfile
import unittest


def model_text(weights, profit, capacity, extra=False):
    text = ['NAME packing', 'OBJSENSE', ' MAX', 'ROWS', ' N obj', ' L cap']
    if extra:
        text += [' L redundant']
    text += ['COLUMNS']
    for j, (a, c) in enumerate(zip(weights, profit)):
        text += [f' x{j} obj {c} cap {a}']
        if extra:
            text += [f' x{j} redundant 3']
    text += ['RHS', f' rhs cap {capacity}']
    if extra:
        text += [f' rhs redundant {3 * len(weights)}']
    text += ['BOUNDS', *(f' BV b x{j}' for j in range(len(weights))), 'ENDATA']
    return '\n'.join(text) + '\n'


class RootCutSearch(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='milplab-root-search-')
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.serial = 0

    def invoke(self, command, expected=0):
        run = subprocess.run([str(x) for x in command], text=True, capture_output=True,
                             timeout=25)
        self.assertEqual(run.returncode, expected, (command, run.stdout, run.stderr))
        return json.loads(run.stdout)

    def solve(self, text, flags=(), expect='optimal'):
        self.serial += 1
        stem = self.path / str(self.serial)
        model, solution, trace = (stem.with_suffix(s) for s in ('.mps', '.sol', '.jsonl'))
        model.write_text(text)
        result = self.invoke([SOLVER, model, '--mip', '--time-limit', '10',
                              '--trace', trace, '--solution', solution, *flags])
        self.assertEqual(result['status'], expect, result)
        self.assertFalse(result['certifies_optimality'])
        self.assertFalse(result['bound_certified'])
        events = [json.loads(line) for line in trace.read_text().splitlines()]
        cuts = [event for event in events if event['event'] == 'root_gmi_cut']
        self.assertEqual(len(cuts), result['root_cuts'])
        self.assertLessEqual(len(cuts), 4)
        ids = set()
        for index, cut in enumerate(cuts):
            self.assertEqual(cut['scope'], 'ORIGINAL_MIP_IN_BOX')
            self.assertNotIn(cut['cut_id'], ids)
            ids.add(cut['cut_id'])
            certificate = stem.with_suffix(f'.cut{index}')
            certificate.write_text(cut['certificate'])
            checked = self.invoke([CHECKER, '--model', model, '--gmi', certificate])
            self.assertEqual(checked['status'], 'valid')
            self.assertFalse(checked['certifies_original_lp_infeasibility'])
            self.assertFalse(checked['certifies_original_mip_infeasibility'])
            self.assertFalse(checked['certifies_full_tree'])
        for event in events:
            if event['event'] == 'mip_in_box_infeasibility_certificate':
                self.assertTrue(ids)
                self.assertEqual(set(event['cut_ids']), ids)
                self.assertIn('not original LP', event['scope'])
        if result['has_primal']:
            witness = self.invoke([CHECKER, '--model', model, '--solution', solution])
            self.assertTrue(witness['checked_feasible'])
        return result, events

    def test_packing_reduces_nodes_across_backends_and_cold_search(self):
        text = model_text([2, 2], [1, 1], 3)
        control, _ = self.solve(text, ['--no-presolve', '--no-root-cuts'])
        self.assertEqual(control['root_cut_status'], 'disabled')
        for flags in [[], ['--cold'], ['--backend', 'dense'], ['--no-scaling'],
                      ['--no-updates']]:
            with self.subTest(flags=flags):
                result, _ = self.solve(text, ['--no-presolve', *flags])
                self.assertEqual(result['objective'], 1)
                self.assertGreater(result['root_cuts'], 0)
                self.assertEqual(result['root_cut_lp_solves'], 1)
                self.assertLess(result['nodes'], control['nodes'])
                if '--cold' in flags:
                    self.assertEqual(result['warm_starts'], 0)

    def test_integer_parity_cut_is_separate_from_lp_infeasibility(self):
        text = ('NAME parity\nROWS\n N obj\n E eq\nCOLUMNS\n x obj 1 eq 2\n'
                'RHS\n rhs eq 1\nBOUNDS\n BV b x\nENDATA\n')
        result, events = self.solve(text, ['--no-presolve'], expect='infeasible')
        self.assertEqual(result['nodes'], 1)
        self.assertEqual(result['root_cuts'], 1)
        certificate = next(e['certificate'] for e in events if e['event'] == 'root_gmi_cut')
        model, proof = self.path / 'lp.mps', self.path / 'parity.cut'
        model.write_text(text)
        proof.write_text(certificate)
        checked = self.invoke([CHECKER, '--model', model, '--gmi', proof])
        self.assertTrue(checked['certifies_mip_in_box_infeasibility'])
        relaxed = self.invoke([SOLVER, model, '--relax', '--time-limit', '5'])
        self.assertEqual(relaxed['status'], 'optimal')
        self.assertEqual(relaxed['objective'], .5)
        control, _ = self.solve(text, ['--no-presolve', '--no-root-cuts'], expect='infeasible')
        self.assertGreater(control['nodes'], result['nodes'])

    def test_presolve_row_masks_keep_original_proof_names(self):
        result, events = self.solve(model_text([2, 2], [1, 1], 3, extra=True))
        self.assertEqual(result['objective'], 1)
        self.assertGreater(result['presolve_rows'], 0)
        self.assertGreater(result['root_cuts'], 0)
        self.assertTrue(any(e['event'] == 'root_gmi_cut' for e in events))

    def test_generated_knapsacks_match_independent_enumeration(self):
        rng = random.Random(819827)
        total_cuts = 0
        for sample in range(20):
            n = 3 + sample % 3
            weight = [rng.randrange(2, 9) for _ in range(n)]
            profit = [rng.randrange(1, 8) for _ in range(n)]
            capacity = sum(weight) // 2
            optimum = max(sum(c * x for c, x in zip(profit, point))
                          for point in itertools.product((0, 1), repeat=n)
                          if sum(a * x for a, x in zip(weight, point)) <= capacity)
            for control in [False, True]:
                with self.subTest(sample=sample, control=control):
                    flags = ['--no-presolve'] + (['--no-root-cuts'] if control else [])
                    result, _ = self.solve(model_text(weight, profit, capacity), flags)
                    self.assertEqual(result['objective'], optimum)
                    total_cuts += result['root_cuts']
        self.assertGreater(total_cuts, 5)

    def test_no_exact_abstains_and_cli_modes_remain_explicit(self):
        result, _ = self.solve(model_text([2, 2], [1, 1], 3),
                               ['--no-presolve', '--no-exact'])
        self.assertEqual(result['objective'], 1)
        self.assertEqual(result['root_cuts'], 0)
        model = self.path / 'mode.mps'
        model.write_text(model_text([2, 2], [1, 1], 3))
        invalid = self.invoke([SOLVER, model, '--relax', '--no-root-cuts'], expected=2)
        self.assertEqual(invalid['status'], 'input_error')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', required=True)
    parser.add_argument('--checker', required=True)
    args, remaining = parser.parse_known_args()
    SOLVER, CHECKER = Path(args.solver).resolve(), Path(args.checker).resolve()
    unittest.main(argv=[__file__, *remaining])
