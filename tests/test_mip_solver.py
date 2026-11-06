#!/usr/bin/env python3
"""Exact tiny-MILP oracles and native branch-and-bound contracts.

The oracle enumerates bounded integer assignments and uses independent Fraction
vertex enumeration for each remaining continuous LP. It does not reuse the C
solver's search, presolve, model reader or numerical stopping rules.
"""
import itertools
import json
import math
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

from test_lp_solver import inequalities, mps, oracle as lp_oracle

SOLVER = str(pathlib.Path(sys.argv.pop(1)).resolve()) if __name__ == '__main__' else None
CHECKER = (str(pathlib.Path(sys.argv.pop(1)).resolve())
           if __name__ == '__main__' and len(sys.argv) > 1 and not sys.argv[1].startswith('-')
           else None)


def integer_mps(model):
    """Add integer markers to the independent LP writer without changing bounds."""
    text = mps(model).splitlines()
    marked = {f'x{i}' for i, integer in enumerate(model['integer']) if integer}
    start, end = text.index('COLUMNS') + 1, text.index('RHS')
    columns = []
    active = False
    marker = 0
    for line in text[start:end]:
        name = line.split()[0]
        integer = name in marked
        if integer != active:
            columns.append(f" M{marker} 'MARKER' '{'INTORG' if integer else 'INTEND'}'")
            marker += 1
            active = integer
        columns.append(line)
    if active:
        columns.append(f" M{marker} 'MARKER' 'INTEND'")
    return '\n'.join(text[:start] + columns + text[end:]) + '\n'


def mip_oracle(model):
    """Exact optimum, or None for infeasibility, on bounded tiny test models."""
    integer = [j for j, flag in enumerate(model['integer']) if flag]
    continuous = [j for j, flag in enumerate(model['integer']) if not flag]
    domains = []
    for j in integer:
        lo, hi = model['bounds'][j]
        if lo is None or hi is None:
            raise ValueError('the enumeration oracle requires finite integer bounds')
        domains.append(range(math.ceil(F(lo)), math.floor(F(hi)) + 1))
    best = None
    for values in itertools.product(*domains):
        fixed = dict(zip(integer, values))
        constant = F(model.get('offset', 0)) + sum(F(model['cost'][j]) * x
                                                  for j, x in fixed.items())
        if continuous:
            residual = {'cost': [model['cost'][j] for j in continuous],
                        'bounds': [model['bounds'][j] for j in continuous],
                        'offset': constant, 'maximize': bool(model.get('maximize')),
                        'rows': []}
            for coeff, lo, hi in model['rows']:
                shifted = sum(F(coeff[j]) * x for j, x in fixed.items())
                residual['rows'].append(([coeff[j] for j in continuous],
                                         None if lo is None else F(lo) - shifted,
                                         None if hi is None else F(hi) - shifted))
            value = lp_oracle(residual)
        else:
            valid = True
            for coeff, lo, hi in model['rows']:
                activity = sum(F(coeff[j]) * x for j, x in fixed.items())
                if (lo is not None and activity < F(lo)) or (hi is not None and activity > F(hi)):
                    valid = False
                    break
            value = constant if valid else None
        if value is not None and (best is None or
                                  (value > best if model.get('maximize') else value < best)):
            best = value
    return best


def generated_models(count=300):
    rng = random.Random(20261107)
    for case in range(count):
        family = case % 3
        n = 3 if family == 0 or case % 7 == 0 else 2
        integer = [True] * n if family < 2 else [True] + [False] * (n - 1)
        bounds = [(F(0), F(1))] * n if family == 0 else [(F(-2), F(2))] * n
        if family == 1 and case % 5 == 0:
            bounds[0] = (F(-3, 2), F(5, 2))
        point = [F(rng.randint(math.ceil(lo), math.floor(hi))) if integer[j]
                 else F(rng.randint(int(2 * lo), int(2 * hi)), 2)
                 for j, (lo, hi) in enumerate(bounds)]
        if case % 17 == 0:
            bounds[0] = (point[0], point[0])
        rows = []
        for _ in range(rng.randint(1, 3)):
            coeff = [F(rng.randint(-4, 4), rng.choice((1, 2))) for _ in range(n)]
            activity = sum(a * x for a, x in zip(coeff, point))
            slack = F(rng.randint(0, 4), 2)
            kind = rng.choice(('L', 'G', 'E', 'R'))
            lo = activity - slack if kind in ('G', 'R') else activity if kind == 'E' else None
            hi = activity + slack if kind in ('L', 'R') else activity if kind == 'E' else None
            rows.append((coeff, lo, hi))
        if case % 11 == 0:
            # A half-integral equality deliberately separates integer and LP
            # feasibility; never infer MILP feasibility from the relaxation.
            fractional = F(math.ceil(bounds[0][0])) + F(1, 2)
            rows.append(([F(int(j == 0)) for j in range(n)], fractional, fractional))
        yield {'cost': [F(rng.randint(-5, 5), rng.choice((1, 2))) for _ in range(n)],
               'integer': integer, 'bounds': bounds, 'rows': rows,
               'offset': F(rng.randint(-12, 12), 4), 'maximize': bool(case % 2)}


class MilpSolver(unittest.TestCase):
    def invoke(self, model, extra=(), *, mode='--mip', checker=True):
        with tempfile.TemporaryDirectory(prefix='milp oracle ') as td:
            root = pathlib.Path(td)
            path = root / 'model.mps'
            path.write_text(model if isinstance(model, str) else integer_mps(model))
            solution = root / 'solution.sol'
            args = [SOLVER, str(path), mode, '--solution', str(solution), *map(str, extra)]
            p = subprocess.run(args, text=True, capture_output=True, timeout=30)
            try:
                result = json.loads(p.stdout)
            except ValueError:
                self.fail(f'non-JSON solver output: {p.stdout!r}; stderr: {p.stderr!r}')
            self.assertIn(p.returncode, (0, 2, 3, 4), (result, p.stderr))
            for key, value in result.items():
                if isinstance(value, float):
                    self.assertTrue(math.isfinite(value), (key, result))
            self.assertFalse(result.get('certifies_optimality', False), result)
            self.assertFalse(result.get('bound_certified', False), result)
            values = {}
            if solution.exists():
                for line in solution.read_text().splitlines():
                    if not line.strip() or line.lstrip().startswith('#'):
                        continue
                    name, value = line.split()
                    if name == '=obj=':
                        continue
                    self.assertNotIn(name, values)
                    values[name] = F(value)
                if checker and CHECKER:
                    q = subprocess.run([CHECKER, '--model', str(path), '--solution', str(solution)] +
                                       (['--relax'] if mode == '--relax' else []),
                                       text=True, capture_output=True, timeout=15)
                    self.assertEqual(q.returncode, 0, (model, result, q.stdout, q.stderr))
                    self.assertTrue(json.loads(q.stdout)['checked_feasible'])
            self.assertEqual(bool(result.get('has_primal')), bool(values), result)
            return p.returncode, result, values

    def assert_point(self, model, result, values):
        self.assertEqual(set(values), {f'x{j}' for j in range(len(model['cost']))}, result)
        point = [values[f'x{j}'] for j in range(len(model['cost']))]
        for j, flag in enumerate(model['integer']):
            if flag:
                self.assertLessEqual(abs(point[j] - round(point[j])), F(1, 10**7), (model, point, result))
        for a, b in inequalities(model):
            residual = sum(F(x) * y for x, y in zip(a, point)) - F(b)
            self.assertLessEqual(residual, F(1, 10**7), (model, point, a, b, result))
        value = F(model.get('offset', 0)) + sum(F(c) * x for c, x in zip(model['cost'], point))
        self.assertAlmostEqual(result['objective'], float(value), delta=1e-7 * max(1, abs(float(value))))
        return value

    def compare(self, model, extra=()):
        expected = mip_oracle(model)
        code, result, values = self.invoke(model, extra)
        self.assertEqual(code, 0, (model, result))
        self.assertEqual(result['status'], 'infeasible' if expected is None else 'optimal', (model, result))
        self.assertFalse(result.get('relaxation', True), result)
        if expected is None:
            self.assertFalse(result['has_primal'], result)
            self.assertIsNone(result.get('objective'), result)
            return result
        value = self.assert_point(model, result, values)
        self.assertAlmostEqual(float(value), float(expected), delta=1e-7 * max(1, abs(float(expected))),
                               msg=str((model, result)))
        return result

    def test_300_generated_exact_integer_and_mixed_oracles(self):
        for case, model in enumerate(generated_models()):
            with self.subTest(case=case):
                self.compare(model)

    def test_oracle_self_checks_fractional_gap_and_integer_infeasibility(self):
        model = {'cost': [1, 1], 'bounds': [(0, 1), (0, 1)], 'integer': [True, True],
                 'rows': [([2, 2], None, 3)], 'maximize': True}
        self.assertEqual(mip_oracle(model), 1)
        self.assertEqual(lp_oracle(model), F(3, 2))
        self.assertEqual(self.compare(model)['objective'], 1)
        code, relaxation, _ = self.invoke(model, mode='--relax')
        self.assertEqual(code, 0, relaxation)
        self.assertEqual(relaxation['objective'], 1.5)
        self.assertTrue(relaxation['relaxation'])
        impossible = {**model, 'rows': [([1, 0], F(1, 2), F(1, 2))]}
        self.assertIsNone(mip_oracle(impossible))
        self.compare(impossible)

    def test_infinite_integer_sides_and_new_child_bound_topology(self):
        cases = [
            ({'cost': [1], 'integer': [True], 'bounds': [(0, None)],
              'rows': [([1], F(1, 2), None)]}, F(1)),
            ({'cost': [1], 'integer': [True], 'bounds': [(None, F(-6, 5))],
              'rows': [], 'maximize': True}, F(-2)),
        ]
        for model, expected in cases:
            for extra in ((), ('--cold',)):
                with self.subTest(model=model, extra=extra):
                    code, result, values = self.invoke(model, extra)
                    self.assertEqual(code, 0, result)
                    self.assertEqual(result['status'], 'optimal', result)
                    self.assertEqual(self.assert_point(model, result, values), expected)

    def test_fractional_continuous_repair_free_fixed_and_objective_senses(self):
        cases = [
            {'cost': [-2, 1], 'integer': [True, False], 'bounds': [(0, 3), (0, 2)],
             'rows': [([1, 3], 4, None)], 'offset': F(7, 4)},
            {'cost': [2, 1], 'integer': [True, False], 'bounds': [(-2, 2), (None, None)],
             'rows': [([1, 3], 1, 1), ([0, 1], -2, 2)], 'maximize': True},
            {'cost': [-1, F(1, 10), 2], 'integer': [True, False, True],
             'bounds': [(F(-3, 2), F(5, 2)), (F(3, 10), F(3, 10)), (-1, 2)],
             'rows': [([1, 0, 1], -1, F(5, 2))], 'offset': F(-11, 4)},
            {'cost': [1, 3], 'integer': [False, False], 'bounds': [(0, 2), (0, 2)],
             'rows': [([1, 2], F(1, 2), F(5, 2))], 'maximize': True},
        ]
        for model in cases:
            with self.subTest(model=model):
                self.compare(model)

    def test_cold_dense_sparse_and_deterministic_replay(self):
        model = {'cost': [3, 2, -1], 'integer': [True, True, False],
                 'bounds': [(0, 3), (-1, 2), (0, 2)],
                 'rows': [([2, 2, 1], None, F(9, 2)), ([1, -1, 2], -1, 3)],
                 'maximize': True, 'offset': F(3, 4)}
        warm = self.compare(model)
        cold = self.compare(model, ['--cold'])
        dense = self.compare(model, ['--backend', 'dense'])
        sparse = self.compare(model, ['--backend', 'sparse'])
        replay = self.compare(model)
        for result in (cold, dense, sparse, replay):
            self.assertAlmostEqual(result['objective'], warm['objective'], delta=1e-7)
        self.assertEqual(warm.get('nodes'), replay.get('nodes'))
        self.assertEqual(warm.get('lp_iterations'), replay.get('lp_iterations'))

    def test_node_and_time_limits_do_not_claim_completion(self):
        model = {'cost': [1, 1], 'bounds': [(0, 1), (0, 1)], 'integer': [True, True],
                 'rows': [([2, 2], None, 3)], 'maximize': True}
        for extra, status in [(['--node-limit', '1'], 'node_limit'),
                              (['--time-limit', '0.000000000001'], 'time_limit')]:
            with self.subTest(extra=extra):
                code, result, _ = self.invoke(model, extra)
                self.assertEqual(code, 4, result)
                self.assertEqual(result['status'], status, result)
                self.assertFalse(result.get('certifies_optimality', False))

    def test_root_presolve_tightens_box_and_owned_row_view(self):
        # x+y<=3 and x-y>=1 with fixed y=1 imply x=2. Both source rows
        # become redundant only inside that returned box. The original model
        # is still used by assert_point and the independent external checker.
        model = {'cost': [1, F(1, 10)], 'integer': [True, False],
                 'bounds': [(0, 8), (1, 1)],
                 'rows': [([1, 1], None, 3), ([1, -1], 1, None)]}
        self.assertEqual(mip_oracle(model), F(21, 10))
        for cold in ([], ['--cold']):
            with self.subTest(cold=cold):
                enabled = self.compare(model, cold)
                disabled = self.compare(model, [*cold, '--no-presolve'])
                self.assertGreater(enabled['presolve_passes'], 0, enabled)
                self.assertGreater(enabled['presolve_tightenings'], 0, enabled)
                self.assertEqual(enabled['presolve_fixed'], 2, enabled)
                self.assertEqual(enabled['presolve_rows'], 2, enabled)
                self.assertEqual(disabled['presolve_passes'], 0, disabled)
                self.assertEqual(disabled['presolve_tightenings'], 0, disabled)
                self.assertEqual(disabled['presolve_rows'], 0, disabled)
                self.assertEqual(enabled['objective'], disabled['objective'])

    def test_root_presolve_integer_contradiction_is_not_lp_unboundedness(self):
        # The continuous interval [.5,.75] is nonempty, while the integer
        # domain is empty. Disable propagation to retain the branching path.
        model = {'cost': [1], 'integer': [True], 'bounds': [(0, 2)],
                 'rows': [([1], F(1, 2), F(3, 4))]}
        self.assertIsNone(mip_oracle(model))
        enabled = self.compare(model)
        disabled = self.compare(model, ['--no-presolve'])
        self.assertEqual(enabled['nodes'], 0, enabled)
        self.assertEqual(enabled['lp_solves'], 0, enabled)
        self.assertGreater(enabled['presolve_tightenings'], 0, enabled)
        self.assertGreater(disabled['nodes'], 0, disabled)
        self.assertGreater(disabled['lp_solves'], 0, disabled)
        self.assertEqual(disabled['branch_observations'], 0, disabled)

    def test_most_fractional_selection_learning_and_first_fractional_ablation(self):
        # Root LP is (.25,.5). Most fractionality chooses x1; the explicit
        # legacy ablation chooses x0. Two optimal children can train history;
        # infeasible children and the final fixed-integer repair cannot.
        base = {'cost': [1, 1], 'integer': [True, True], 'bounds': [(0, 2), (0, 2)],
                'rows': [([1, 0], F(1, 4), None), ([0, 1], F(1, 2), None)]}
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            signatures = []
            for offset in (0, 10**12):
                model = {**base, 'offset': offset}
                for policy, flags, expected in [('default', [], 'x1'),
                                                 ('replay', [], 'x1'),
                                                 ('first', ['--first-fractional'], 'x0')]:
                    trace = root / f'{offset}-{policy}.jsonl'
                    # Preserve the branching-policy experiment: root cuts can
                    # close these independent rounding gaps before branching.
                    result = self.compare(model, ['--no-presolve', '--no-root-cuts',
                                                   '--trace', trace, *flags])
                    events = [json.loads(line) for line in trace.read_text().splitlines()]
                    branches = [(event['variable'], event['value']) for event in events
                                if event['event'] == 'branch']
                    self.assertEqual(branches[0][0], expected, (result, events))
                    self.assertEqual(result['branch_observations'], 0 if policy == 'first' else 2, result)
                    if policy != 'first':
                        signatures.append(branches)
            self.assertTrue(all(signature == signatures[0] for signature in signatures))

    def test_node_limit_retains_a_checked_incumbent(self):
        # Either first down branch produces objective one, while a sibling
        # with LP bound1.5 remains. Hitting a limit must keep that witness.
        model = {'cost': [1, 1], 'bounds': [(0, 1), (0, 1)], 'integer': [True, True],
                 'rows': [([2, 2], None, 3)], 'maximize': True}
        code, result, values = self.invoke(model, ['--node-limit', '2'])
        self.assertEqual(code, 4, result)
        self.assertEqual(result['status'], 'node_limit', result)
        self.assertTrue(result['has_primal'], result)
        self.assertEqual(self.assert_point(model, result, values), 1)
        if result.get('best_bound') is not None:
            self.assertGreaterEqual(result['best_bound'], result['objective'])

    def test_constant_objective_offset_cannot_close_a_worse_rounded_node(self):
        # The relaxed optimum is just below x1, but rounding x to1 forces
        # t32; integer x0,t0 is better. A large constant objective offset must
        # not turn a material repair gap into an acceptable relative error.
        # Use exact candidate checks: a tolerance scaled by the1e12 offset
        # would itself hide the original32-unit wrong-optimum regression.
        scale = 2**30
        for offset in (0, 10**12):
            model = {'cost': [0, 1], 'integer': [True, False],
                     'bounds': [(0, 1), (-2, 100)],
                     'rows': [([-1, -1], None, 0), ([-scale, 1], -(scale - 32), None)],
                     'offset': offset}
            self.assertEqual(mip_oracle(model), offset)
            for extra in ((), ('--cold',)):
                with self.subTest(offset=offset, extra=extra):
                    code, result, values = self.invoke(model, extra)
                    self.assertEqual(code, 0, result)
                    self.assertEqual(result['status'], 'optimal', result)
                    self.assertEqual(result['objective'], offset, result)
                    self.assertEqual(self.assert_point(model, result, values), offset)
                    self.assertEqual(values, {'x0': F(0), 'x1': F(0)})

        # A costly fixed variable is also a constant. Merely removing the MPS
        # objective offset from a relative tolerance would leave this bug.
        fixed_constant = {'cost': [0, 1, 10**12], 'integer': [True, False, False],
                          'bounds': [(0, 1), (-2, 100), (1, 1)],
                          'rows': [([-1, -1, 0], None, 0),
                                   ([-scale, 1, 0], -(scale - 32), None)]}
        self.assertEqual(mip_oracle(fixed_constant), 10**12)
        for extra in ((), ('--cold',)):
            with self.subTest(fixed_constant=True, extra=extra):
                code, result, values = self.invoke(fixed_constant, extra)
                self.assertEqual(code, 0, result)
                self.assertEqual(result['status'], 'optimal', result)
                self.assertEqual(result['objective'], 10**12, result)
                self.assertEqual(self.assert_point(fixed_constant, result, values), 10**12)
                self.assertEqual(values, {'x0': F(0), 'x1': F(0), 'x2': F(1)})

    def test_search_trace_is_replayable_and_never_replaced(self):
        model = {'cost': [1, 1], 'bounds': [(0, 1), (0, 1)], 'integer': [True, True],
                 'rows': [([2, 2], None, 3)], 'maximize': True}
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path = root / 'model.mps'; path.write_text(integer_mps(model))
            traces = []
            for version in (1, 2):
                trace = root / f'trace-{version}.jsonl'
                p = subprocess.run([SOLVER, str(path), '--mip', '--trace', str(trace)],
                                   text=True, capture_output=True, timeout=30)
                result = json.loads(p.stdout)
                self.assertEqual(p.returncode, 0, result)
                self.assertEqual(result['status'], 'optimal', result)
                traces.append(trace.read_bytes())
            self.assertEqual(traces[0], traces[1])
            events = [json.loads(line) for line in traces[0].decode().splitlines()]
            self.assertTrue(any(event['event'] == 'branch' for event in events))
            self.assertTrue(any(event['event'] == 'incumbent' for event in events))
            old = root / 'trace-1.jsonl'
            p = subprocess.run([SOLVER, str(path), '--mip', '--trace', str(old)],
                               text=True, capture_output=True, timeout=15)
            result = json.loads(p.stdout)
            self.assertEqual(p.returncode, 2, result)
            self.assertEqual(old.read_bytes(), traces[0])

    def test_lp_unbounded_is_not_milp_unbounded(self):
        # The relaxation has an improving y ray but integer x cannot equal1/2.
        # This first engine may correctly stop with explicit LP uncertainty.
        model = {'cost': [0, -1], 'integer': [True, False],
                 'bounds': [(0, 1), (0, None)], 'rows': [([1, 0], F(1, 2), F(1, 2))]}
        code, result, _ = self.invoke(model)
        self.assertIn(result['status'], ('lp_unbounded', 'infeasible', 'unsupported', 'numerical_failure'), result)
        self.assertEqual(code, 0 if result['status'] == 'infeasible' else
                         3 if result['status'] == 'unsupported' else 4, result)
        self.assertFalse(result.get('certifies_optimality', False))

    def test_modes_and_limits_are_explicit(self):
        model = {'cost': [1], 'integer': [True], 'bounds': [(0, 2)], 'rows': []}
        for extra in (['--relax'], ['--node-limit', '0'], ['--node-limit', '-1'],
                      ['--time-limit', 'nan'], ['--backend', 'imaginary']):
            with self.subTest(extra=extra):
                code, result, _ = self.invoke(model, extra)
                self.assertEqual(code, 2, result)
                self.assertEqual(result['status'], 'input_error', result)

    def test_existing_solution_is_never_replaced(self):
        model = {'cost': [-1], 'integer': [True], 'bounds': [(0, 2)], 'rows': []}
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td)
            path = root / 'model.mps'; path.write_text(integer_mps(model))
            solution = root / 'old.sol'; solution.write_text('preserve historical witness\n')
            p = subprocess.run([SOLVER, str(path), '--mip', '--solution', str(solution)],
                               text=True, capture_output=True, timeout=15)
            result = json.loads(p.stdout)
            self.assertEqual(p.returncode, 2, result)
            self.assertEqual(result['status'], 'input_error', result)
            self.assertEqual(solution.read_text(), 'preserve historical witness\n')

    def test_huge_integer_range_is_explicitly_unsupported(self):
        models = [
            {'cost': [-1], 'integer': [True], 'bounds': [(0, 2**53 + 4)], 'rows': []},
            # Integer range restrictions also apply to values implied by rows.
            # In double arithmetic this RHS rounds down by1; accepting its LP
            # witness would violate the original rational model by a full unit.
            {'cost': [1], 'integer': [True], 'bounds': [(0, None)],
             'rows': [([1], 2**53 + 1, None)]},
        ]
        for model in models:
            with self.subTest(model=model):
                code, result, _ = self.invoke(model)
                self.assertEqual(code, 3, result)
                self.assertEqual(result['status'], 'unsupported', result)

    def test_tiny_positive_coefficient_is_not_an_infeasibility_proof(self):
        # Dropping a small but improving Phase-I reduced cost leaves its
        # artificial variable positive; that alone must not fabricate integer
        # infeasibility. The feasible optimum is the exactly spaced1e11.
        model = {'cost': [1], 'integer': [True], 'bounds': [(0, None)],
                 'rows': [([F(1, 10**11)], 1, None)]}
        for extra in ((), ('--cold',)):
            with self.subTest(extra=extra):
                code, result, values = self.invoke(model, extra)
                self.assertIn(result['status'], ('optimal', 'numerical_failure'), result)
                if result['status'] == 'optimal':
                    self.assertEqual(code, 0, result)
                    self.assertEqual(result['objective'], 10**11)
                    self.assertEqual(self.assert_point(model, result, values), 10**11)
                    self.assertEqual(values['x0'], 10**11)
                else:
                    self.assertEqual(code, 4, result)


if __name__ == '__main__':
    unittest.main()
