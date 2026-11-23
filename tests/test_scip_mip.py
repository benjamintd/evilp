"""Contract tests for the general C17 SCIP branch-and-cut and symmetry detector."""
from fractions import Fraction as F
from itertools import permutations, product
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from test_lp_solver import decimal
from test_mip_solver import generated_models, integer_mps, mip_oracle

SOLVER, CHECKER = (Path(p).resolve() for p in sys.argv[1:3])
del sys.argv[1:3]


class ScipGeneralMip(unittest.TestCase):
    def test_binary_cover_uses_exact_decimal_threshold(self):
        low = F(3333333333, 10**10)
        just_below = F(3333333333333333332, 10**19)
        just_above = 1 - 2 * just_below
        for weights, cost, expected in [
            ([low] * 7, [1] * 7, 4),
            ([F(3333333334, 10**10)] * 7, [1] * 7, 3),
            ([F(1, 100)] * 7, [1] * 7, 10),
            ([just_below, just_below, just_above] + [just_below] * 4, [1, 1, 1, 4, 4, 4, 4], 3),
        ]:
            for upper in (False, True):
                sign = -1 if upper else 1
                a = [sign * w for w in [F(5, 4), 2, *weights]] + [0] * 25
                data = dict(cost=[sign * c for c in [10, 20, *cost]] + [0] * 25,
                            bounds=[(0, 1)] * 34, integer=[True] * 9 + [False] * 25,
                            rows=[(a, None, -1) if upper else (a, 1, None)], offset=0, maximize=upper)
                with self.subTest(weights=weights, upper=upper), tempfile.TemporaryDirectory() as tmp:
                    result, model, solution = self.run_model(Path(tmp), integer_mps(data))
                    self.assertEqual(result['status'], 'optimal', result)
                    check = subprocess.run([CHECKER, '--model', model, '--solution', solution, '--exact'],
                                           capture_output=True, text=True, check=True)
                    self.assertEqual(F(json.loads(check.stdout)['objective']), sign * expected)

    def test_binary_cover_respects_continuous_columns_and_both_row_sides(self):
        a = [2] + [F(2, 5)] * 4 + [0] * 29
        for integer, upper, expected in [(False, None, F(1, 2)), (True, F(5, 4), 3)]:
            data = dict(cost=[1] * 5 + [0] * 29, integer=[integer] + [True] * 4 + [False] * 29,
                        bounds=[(0, 1)] * 34, rows=[(a, 1, upper)], offset=0, maximize=False)
            with self.subTest(integer=integer), tempfile.TemporaryDirectory() as tmp:
                result, model, solution = self.run_model(Path(tmp), integer_mps(data))
                self.assertEqual(result['status'], 'optimal', result)
                check = subprocess.run([CHECKER, '--model', model, '--solution', solution, '--exact'],
                                       capture_output=True, text=True, check=True)
                self.assertEqual(F(json.loads(check.stdout)['objective']), expected)

    def test_uniform_binary_cover_rounds_the_source_ratio(self):
        for n, unit, required in [(7, F(3333333333, 10**10), 4), (2, F(4999999999, 10**10), None)]:
            data = dict(cost=[1] * n + [0] * 32, integer=[True] * n + [False] * 32,
                        bounds=[(0, 1)] * (n + 32), offset=0, maximize=False,
                        rows=[([unit] * n + [0] * 32, 1, None)])
            with self.subTest(n=n), tempfile.TemporaryDirectory() as tmp:
                result, model, solution = self.run_model(Path(tmp), integer_mps(data))
                if required is None:
                    self.assertEqual(result['status'], 'infeasible', result)
                else:
                    self.assertEqual(result['status'], 'optimal', result)
                    check = subprocess.run([CHECKER, '--model', model, '--solution', solution, '--exact'],
                                           capture_output=True, text=True, check=True)
                    self.assertEqual(F(json.loads(check.stdout)['objective']), required)

    def test_signed_integer_row_rounding_against_enumeration(self):
        for unit, rhs, upper, domain in [
            (F(3333333334, 10**10), 1, True, (0, 1)),
            (F(3333333334, 10**10), -1, False, (-2, 2)),
            (F(3333333333, 10**10), -1, True, (-2, 2)),
            (F(9999999999999999999, 10**19), 1, False, (0, 1)),
        ]:
            signs = [1, -1, 1, 1]
            values = [sum(s * v for s, v in zip(signs, x))
                      for x in product(range(domain[0], domain[1] + 1), repeat=len(signs))]
            expected = (max(v for v in values if unit * v <= rhs) if upper else
                        min(v for v in values if unit * v >= rhs))
            a = [s * unit for s in signs] + [0] * 32
            data = dict(cost=signs + [0] * 32, integer=[True] * 4 + [False] * 32,
                        bounds=[domain] * 4 + [(0, 1)] * 32, offset=0, maximize=upper,
                        rows=[(a, None, rhs) if upper else (a, rhs, None)])
            with self.subTest(unit=unit, upper=upper, domain=domain), tempfile.TemporaryDirectory() as tmp:
                result, model, solution = self.run_model(Path(tmp), integer_mps(data))
                self.assertEqual(result['status'], 'optimal', result)
                check = subprocess.run([CHECKER, '--model', model, '--solution', solution, '--exact'],
                                       capture_output=True, text=True, check=True)
                self.assertEqual(F(json.loads(check.stdout)['objective']), expected)

    def test_nonuniform_integer_row_lattice_against_enumeration(self):
        # A bound can land exactly on an attainable activity in binary64 while
        # lying strictly on the other side in the original decimal model.
        cases = [
            ([F(1, 5), F(2, 5), F(3, 5)], 1),
            ([F(-1, 8), F(3, 16), F(5, 32)], F(-3, 32)),
            ([F(-3, 100), F(21, 200), F(-3, 40)], F(3, 40)),
        ]
        for (a, boundary), upper, digits in product(cases, (True, False), (10, 19)):
            rhs = boundary + (-1 if upper else 1) * F(1, 10**digits)
            cost = [int(c * 800) for c in a]
            feasible = []
            for x in product(range(-2, 3), repeat=3):
                activity = sum(c * y for c, y in zip(a, x))
                if (activity <= rhs if upper else activity >= rhs):
                    feasible.append(sum(c * y for c, y in zip(cost, x)))
            expected = (max if upper else min)(feasible)
            coeff = a + [0] * 32
            model = dict(cost=cost + [0] * 32, integer=[True] * 3 + [False] * 32,
                         bounds=[(-2, 2)] * 3 + [(0, 1)] * 32, offset=0, maximize=upper,
                         rows=[(coeff, None, rhs) if upper else (coeff, rhs, None)])
            with self.subTest(a=a, upper=upper, digits=digits), tempfile.TemporaryDirectory() as tmp:
                result, original, witness = self.run_model(Path(tmp), integer_mps(model))
                self.assertEqual((result['status'], result['objective']), ('optimal', expected))
                checked = subprocess.run([CHECKER, '--exact', '--model', original, '--solution', witness],
                                         capture_output=True, text=True, check=True)
                self.assertEqual(F(json.loads(checked.stdout)['objective']), expected)

    def run_model(self, directory, text, options=(), env=None):
        model = directory / 'test.mps'
        sol = directory / 'point.sol'
        model.write_text(text)
        r = subprocess.run([SOLVER, model, '--mip', '--solution', sol, '--time-limit', '10',
                            *options], capture_output=True, text=True, timeout=20, env=env)
        self.assertIn(r.returncode, (0, 4), (r.stdout, r.stderr))
        self.assertNotIn('runtime error:', r.stderr)
        return json.loads(r.stdout), model, sol

    def test_small_reduced_cost_on_large_transport_flow(self):
        # With equal supplies and demands the LP has an optimal permutation.
        # Every permutation remains feasible with an integer first column entry.
        increments = [[0, 61, 0], [22, 84, 61], [61, 38, 61]]
        costs = [F(22000000 + v, 10**9) for row in increments for v in row]
        expected = 10000 * min(sum(costs[3*i+j] for i, j in enumerate(order))
                               for order in permutations(range(3)))
        for algorithm in ('default', 'primal'):
            with self.subTest(algorithm=algorithm), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                lines = ['NAME TRANSPORT', 'ROWS', ' N COST']
                rows = [f'{side}{i}' for side in ('r', 's') for i in range(3)]
                lines.extend(f' E {row}' for row in rows)
                lines.append('COLUMNS')
                for i, j in product(range(3), repeat=2):
                    lines.extend([f' x{i}{j} COST {decimal(costs[3*i+j])} r{i} 1',
                                  f' x{i}{j} s{j} 1'])
                lines.extend(f' z{k} COST 0' for k in range(33))
                lines.append('RHS')
                lines.extend(f' rhs {row} 10000' for row in rows)
                lines.extend(['BOUNDS', ' LI bd x00 0', 'ENDATA'])
                env = os.environ.copy()
                if algorithm == 'primal':
                    settings = directory / 'primal.set'
                    settings.write_text('lp/initalgorithm = p\n')
                    env['LAB_SCIP_SET'] = str(settings)
                result, model, sol = self.run_model(directory, '\n'.join(lines) + '\n', env=env)
                self.assertEqual((result['search_route'], result['status']), ('general_scip', 'optimal'))
                checked = json.loads(subprocess.check_output([CHECKER, '--model', model, '--solution', sol, '--exact']))
                self.assertEqual(F(checked['objective']), expected)
                self.assertAlmostEqual(result['best_bound'], float(expected), places=9)

    def test_near_tied_prices_in_equivalent_formulations(self):
        # a+b >= 10^6 implies the objective is at least 10^6 times the
        # cheaper price. The integer b alone attains that lower bound.
        fixture = (Path(__file__).parent / 'fixtures/near_tied_costs_order.mps').read_text()
        prices = [' b COST 0.022 DEMAND 1', ' a COST 0.0220000004 DEMAND 1']
        for sign, scale, offset, swap in product((1, -1), (1, 1000), (0, 17), (False, True)):
            lines = [line.replace(price, decimal(sign * scale * F(price)))
                     for line, price in zip(prices, ['0.022', '0.0220000004'])]
            text = fixture.replace('\n'.join(prices), '\n'.join(reversed(lines) if swap else lines))
            text = text.replace('ROWS\n', 'OBJSENSE\n ' + ('MAX' if sign < 0 else 'MIN') + '\nROWS\n')
            text = text.replace('RHS\n', f'RHS\n rhs COST {-offset}\n')
            for route, options in [('general', []), ('disabled', ['--no-scip'])]:
                with self.subTest(sign=sign, scale=scale, offset=offset, swap=swap, route=route), \
                        tempfile.TemporaryDirectory() as tmp:
                    result, model, sol = self.run_model(Path(tmp), text, options)
                    self.assertEqual((result['status'], result['scip_route']), ('optimal', route), result)
                    check = subprocess.run([CHECKER, '--model', model, '--solution', sol, '--exact'],
                                           capture_output=True, text=True, check=True)
                    for value in [json.loads(check.stdout)['objective'], result['objective'], result['best_bound']]:
                        self.assertEqual((F(str(value)) - offset) / (sign * scale), 22000)

    def test_distinct_prices_below_the_strict_threshold_fall_back(self):
        source = (Path(__file__).parent / 'fixtures/near_tied_costs_order.mps').read_text()
        text = source.replace('0.0220000004', '0.0220000000000005')
        with tempfile.TemporaryDirectory() as tmp:
            result, model, sol = self.run_model(Path(tmp), text)
            self.assertEqual(result['scip_route'], 'general_fallback')
            self.assertEqual(result['scip_fallback_reason'], 'lossy_model_import')
            self.assertEqual((result['search_route'], result['bound_origin']), ('native', 'search_exhaustion'))
            self.assertEqual((result['status'], result['best_bound']), ('optimal', 22000))
            check = subprocess.run([CHECKER, '--model', model, '--solution', sol, '--exact'],
                                   capture_output=True, text=True, check=True)
            self.assertEqual(F(json.loads(check.stdout)['objective']), 22000)

    def test_relative_gap_stop_preserves_the_open_bound(self):
        # A large offset makes the all-zero feasible independent set reach the
        # requested relative gap while its objective is still improvable.
        n = 36
        rng = random.Random(240)
        rows = []
        for i in range(n):
            for j in range(i):
                if rng.random() < 0.2:
                    row = [0] * n
                    row[i] = row[j] = 1
                    rows.append((row, None, 1))
        data = dict(cost=[1 + rng.randrange(10) for _ in range(n)], bounds=[(0, 1)] * n,
                    integer=[True] * n, rows=rows, offset=10**9, maximize=True)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            model, settings, sol = root/'gap.mps', root/'gap.set', root/'gap.sol'
            model.write_text(integer_mps(data))
            settings.write_text('limits/gap = 1e-6\nlimits/absgap = 0\n')
            run = subprocess.run([SOLVER, model, '--mip', '--solution', sol],
                                 env={**os.environ, 'LAB_SCIP_SET': str(settings)},
                                 capture_output=True, text=True, timeout=20)
            self.assertEqual(run.returncode, 4, run.stderr)
            result = json.loads(run.stdout)
            self.assertEqual((result['status'], result['scip_route']), ('gap_limit', 'general'))
            self.assertFalse(result['certifies_optimality'])
            check = subprocess.run([CHECKER, '--model', model, '--solution', sol, '--exact'],
                                   capture_output=True, text=True, check=True)
            primal = F(json.loads(check.stdout)['objective'])
            dual = F(str(result['best_bound']))
            self.assertGreater(dual, primal)
            self.assertLessEqual((dual - primal) / primal, F(1, 10**6))

    def test_small_coefficient_and_nearly_dependent_equalities(self):
        fixtures = Path(__file__).parent / 'fixtures'
        small = (fixtures / 'small_coefficient_scip.mps').read_text()
        parallel = (fixtures / 'near_parallel_scip.mps').read_text()
        scaled_small = (small.replace(' G R\n', ' L R\n')
                        .replace('O 1 R 1', 'O 1 R -2').replace('y R 1e-10', 'y R -2e-10')
                        .replace('b R 1', 'b R -2'))
        reordered = (parallel.replace(' E A\n E B\n', ' E B\n E A\n')
                     .replace(' x B 1\n', ' x B -2\n')
                     .replace('B 1.000000000014551915228366851806640625',
                              'B -2.00000000002910383045673370361328125')
                     .replace('B .0000019073486328125', 'B -.000003814697265625'))
        for name, text in [('small', small), ('scaled', scaled_small), ('parallel', parallel),
                           ('reordered', reordered)]:
            with self.subTest(model=name), tempfile.TemporaryDirectory() as tmp:
                result, path, point = self.run_model(Path(tmp), text)
                self.assertEqual(result['status'], 'optimal', result)
                checker = subprocess.run([CHECKER, '--model', path, '--solution', point, '--exact'],
                                         capture_output=True, text=True, check=True)
                self.assertEqual(F(json.loads(checker.stdout)['objective']), 0)
                self.assertFalse(result['certifies_optimality'])

    def test_small_objective_over_large_domain(self):
        original = (Path(__file__).parent / 'fixtures/small_coefficient_scip.mps').read_text()
        for exponent in (10, 16):
            text = (original.replace('x O 1 R 1', 'x O 0 R 1')
                    .replace('y R 1e-10', f'y O -1e-{exponent}')
                    .replace('10000000000', '1000000000000'))
            with self.subTest(exponent=exponent), tempfile.TemporaryDirectory() as tmp:
                result, path, point = self.run_model(Path(tmp), text)
                self.assertEqual(result['status'], 'optimal', result)
                checked = subprocess.run([CHECKER, '--model', path, '--solution', point, '--exact'],
                                         capture_output=True, check=True, text=True)
                self.assertEqual(F(json.loads(checked.stdout)['objective']), -F(10**12, 10**exponent))
                if exponent == 16:
                    self.assertEqual(result['scip_fallback_reason'], 'lossy_model_import')

    def test_source_equalities_beyond_backend_precision(self):
        original = (Path(__file__).parent / 'fixtures/near_parallel_scip.mps').read_text()
        source_coefficient = '1.000000000014551915228366851806640625'
        for exponent in (50, 54):
            coefficient = '1.' + str(5 ** exponent).zfill(exponent)
            text = original.replace(source_coefficient, coefficient)
            for flags in ([], ['--no-scip']):
                with self.subTest(exponent=exponent, flags=flags), tempfile.TemporaryDirectory() as tmp:
                    result, path, point = self.run_model(Path(tmp), text, flags)
                    self.assertEqual((result['status'], result['objective']), ('optimal', 0))
                    check = subprocess.run([CHECKER, '--exact', '--model', path, '--solution', point],
                                           capture_output=True, check=True, text=True)
                    self.assertTrue(json.loads(check.stdout)['checked_feasible'])
            if exponent == 54:
                with tempfile.TemporaryDirectory() as tmp:
                    result, _, _ = self.run_model(Path(tmp), text, ['--no-exact'])
                    self.assertEqual(result['status'], 'numerical_failure', result)
            # A redundant inequality makes the matrix rectangular. If the exact
            # recovery cannot find a point, the floating cutoff must stay open.
            rectangular = text.replace(' E B\n', ' E B\n L C\n').replace(
                ' x B 1\n', ' x B 1 C 1\n').replace(' b A 0', ' b C 0\n b A 0')
            with tempfile.TemporaryDirectory() as tmp:
                result, path, _ = self.run_model(Path(tmp), rectangular)
                witness = Path(tmp) / 'known.sol'
                y = 2 ** (exponent - 19)
                witness.write_text(f'x {-y}\ny {y}\n')
                subprocess.run([CHECKER, '--exact', '--model', path, '--solution', witness],
                               capture_output=True, check=True)
                self.assertIn(result['status'], ('optimal', 'numerical_failure'))

    def test_padded_generated_models_against_exact_oracle(self):
        for case, model in enumerate(generated_models(16)):
            optimum = mip_oracle(model)
            padding = 36 - len(model['cost'])
            padded = {**model, 'cost': model['cost'] + [0] * padding,
                      'bounds': model['bounds'] + [(0, 0)] * padding,
                      'integer': model['integer'] + [False] * padding,
                      'rows': [(a + [0] * padding, lo, hi) for a, lo, hi in model['rows']]}
            with tempfile.TemporaryDirectory() as tmp:
                result, path, point = self.run_model(Path(tmp), integer_mps(padded))
                # Every eleventh generated case has a half-integral singleton equality.
                self.assertEqual(result['search_route'], 'presolve' if case % 11 == 0 else 'general_scip')
                if optimum is None:
                    self.assertEqual(result['status'], 'infeasible', result)
                else:
                    self.assertEqual(result['status'], 'optimal', result)
                    check = subprocess.run([CHECKER, '--model', path, '--solution', point],
                                           capture_output=True, text=True, check=True)
                    self.assertLessEqual(abs(F(json.loads(check.stdout)['objective']) - optimum), F(1, 10**7))
                    self.assertLessEqual(abs(F(str(result['best_bound'])) - optimum), F(1, 10**7))

    def test_observing_or_disabling_native_heuristics_preserves_backend(self):
        text = (Path(__file__).parent / 'fixtures/small_coefficient_scip.mps').read_text()
        for options in ([], ['--no-root-diving'], ['--no-root-cuts'], ['--trace']):
            with tempfile.TemporaryDirectory() as tmp:
                trace = Path(tmp) / 'events.jsonl'
                args = [*options, trace] if options == ['--trace'] else options
                result, _, _ = self.run_model(Path(tmp), text, args)
                self.assertEqual(result['scip_route'], 'general', result)
                self.assertEqual((result['status'], result['objective']), ('optimal', 0))
                if trace.exists():
                    events = [json.loads(line)['event'] for line in trace.read_text().splitlines()]
                    self.assertIn('scip_mip_start', events)
                    self.assertIn('scip_mip_optimum', events)
        with tempfile.TemporaryDirectory() as tmp:
            result, _, _ = self.run_model(Path(tmp), text, ['--no-scip'])
            self.assertEqual(result['scip_route'], 'disabled')
            self.assertEqual((result['status'], result['objective']), ('optimal', 0))

    def test_mixed_integer_min_max_and_offset(self):
        # 36 columns: 18 binary + 18 continuous, with linking big-M rows and offset.
        rng = random.Random(1024)
        for maximize in (False, True):
            n_bin = 18
            n = 36
            costs = [rng.randint(-4, 6) for _ in range(n_bin)] + [F(rng.randint(-3, 5), 2) for _ in range(n_bin)]
            bounds = [(0, 1)] * n_bin + [(0, 4)] * n_bin
            integer = [True] * n_bin + [False] * n_bin
            rows = []
            # Each continuous y_i <= 4 * x_i
            for i in range(n_bin):
                a = [0] * n
                a[i] = -4
                a[n_bin + i] = 1
                rows.append((a, None, 0))
            # Budget on binary sum <= 7 and >= 3
            rows.append(([1] * n_bin + [0] * n_bin, 3, 7))
            # Coupled continuous demand >= 10
            rows.append(([0] * n_bin + [1] * n_bin, 10, 20))
            m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows,
                     offset=F(11, 4), maximize=maximize)
            # For fixed binaries this is a box with a single sum constraint and
            # integer endpoints, so an integral continuous optimum exists.
            best = max if maximize else min
            dp = {(0, 0): F(0)}
            for i in range(n_bin):
                after = dict(dp)
                for (opened, units), value in dp.items():
                    if opened == 7:
                        continue
                    for y in range(5):
                        state = opened + 1, units + y
                        score = value + costs[i] + y * costs[n_bin + i]
                        after[state] = best(after.get(state, score), score)
                dp = after
            expected = m['offset'] + best(v for (k, u), v in dp.items()
                                          if 3 <= k <= 7 and 10 <= u <= 20)
            with self.subTest(maximize=maximize), tempfile.TemporaryDirectory() as tmp:
                text = integer_mps(m)
                v, model, sol = self.run_model(Path(tmp), text)
                self.assertEqual(v['status'], 'optimal', v)
                self.assertIn('SCIP branch-and-cut', v['message'])
                self.assertFalse(v['certifies_optimality'])
                chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                     capture_output=True, text=True, timeout=15)
                self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
                c = json.loads(chk.stdout)
                self.assertTrue(c['checked_feasible'], c)
                self.assertEqual(F(c['objective']), expected)

    def test_symmetric_packing_exercises_c17_symmetry(self):
        # 6 identical bins x 8 items = 48 binary vars + 6 bin-open vars = 54 vars.
        # Highly symmetric (S_6 on bins + item twins) -> exercises scip_symmetry.c!
        bins, items = 6, 8
        weights = [4, 4, 4, 4, 3, 3, 3, 3]
        cap = 7
        n = bins * items + bins
        costs = [0] * (bins * items) + [1] * bins
        bounds = [(0, 1)] * n
        integer = [True] * n
        rows = []
        # Each item assigned to exactly 1 bin
        for it in range(items):
            a = [0] * n
            for b in range(bins):
                a[b * items + it] = 1
            rows.append((a, 1, 1))
        # Each bin capacity <= cap * open_b
        for b in range(bins):
            a = [0] * n
            for it in range(items):
                a[b * items + it] = weights[it]
            a[bins * items + b] = -cap
            rows.append((a, None, 0))
        # Add 1 continuous slack variable so binary_sat abstains and SCIP+symmetry runs.
        costs.append(0)
        bounds.append((0, 0))
        integer.append(False)
        rows = [(a + [0], lo, hi) for a, lo, hi in rows]
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, model, sol = self.run_model(Path(tmp), integer_mps(m))
            self.assertEqual(v['status'], 'optimal', v)
            self.assertIn('SCIP branch-and-cut', v['message'])
            chk = subprocess.run([CHECKER, '--exact', '--model', model, '--solution', sol],
                                 capture_output=True, text=True, timeout=15)
            self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
            c = json.loads(chk.stdout)
            self.assertTrue(c['checked_feasible'], c)
            # 4 items of weight 4 + 4 items of weight 3, cap 7 => min 4 bins (each 4+3=7).
            self.assertEqual(F(c['objective']), F(4))

    def test_infeasible_mixed_integer(self):
        n = 36
        costs = [1] * n
        bounds = [(0, 1)] * 35 + [(0, 10)]
        integer = [True] * 35 + [False]
        rows = [
            ([1] * 35 + [0], 10, None),
            ([1] * 35 + [1], None, 9),
        ]
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, _, _ = self.run_model(Path(tmp), integer_mps(m))
            self.assertEqual(v['status'], 'infeasible', v)

    def test_feasjump_pure_feasibility_and_mixed_polish(self):
        # Mixed-integer feasibility + optimization with linked general integers and continuous vars
        # that exercises scip_feasjump.c (both zero-objective feasibility and coordinate polish).
        rng = random.Random(2026)
        n_int = 24
        n_cont = 16
        n = n_int + n_cont
        bounds = [(0, 3)] * n_int + [(0, 5)] * n_cont
        integer = [True] * n_int + [False] * n_cont
        # Construct a planted feasible point so feasibility is guaranteed.
        x_star = [rng.randint(0, 3) for _ in range(n_int)] + [rng.randint(0, 5) for _ in range(n_cont)]
        rows = []
        for _ in range(28):
            a = [rng.randint(-2, 3) if rng.random() < 0.35 else 0 for _ in range(n)]
            if not any(a):
                a[0] = 1
            rhs = sum(ai * xi for ai, xi in zip(a, x_star))
            rows.append((a, rhs - 2, rhs + 2))
        for zero_obj in (True, False):
            costs = [0] * n if zero_obj else [rng.randint(-3, 5) for _ in range(n)]
            m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
            with self.subTest(zero_obj=zero_obj), tempfile.TemporaryDirectory() as tmp:
                v, model, sol = self.run_model(Path(tmp), integer_mps(m))
                self.assertEqual(v['status'], 'optimal', v)
                self.assertIn('SCIP branch-and-cut', v['message'])
                chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                     capture_output=True, text=True, timeout=15)
                self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
                c = json.loads(chk.stdout)
                self.assertTrue(c['checked_feasible'], c)
                if zero_obj:
                    self.assertEqual(F(c['objective']), F(0))
                else:
                    # Must be at least as good as the planted feasible solution.
                    planted_obj = sum(ci * xi for ci, xi in zip(costs, x_star))
                    self.assertLessEqual(F(c['objective']), F(planted_obj))

    def test_parallel_row_presolver_merge_and_infeasibility(self):
        # Mixed-integer model with scalar-multiple parallel rows (positive and negative ratios):
        # Exercises scip_presol_parallel.c both for feasible merging and infeasible cutoff.
        n = 36
        costs = [1 + (j % 5) for j in range(n)]
        bounds = [(0, 1)] * 24 + [(0, 4)] * 12
        integer = [True] * 24 + [False] * 12
        base_row = [2 if j % 3 == 0 else (-1 if j % 3 == 1 else 3) for j in range(n)]
        # Feasible parallel pair: base_row >= 10, 3*base_row <= 36 (i.e. base_row <= 12), -2*base_row <= -22 (i.e. base_row >= 11)
        rows_feas = [
            (base_row, 10, None),
            ([3 * a for a in base_row], None, 36),
            ([-2 * a for a in base_row], None, -22),
            ([1] * 24 + [0] * 12, 4, 16),
        ]
        m_feas = dict(cost=costs, bounds=bounds, integer=integer, rows=rows_feas, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, model, sol = self.run_model(Path(tmp), integer_mps(m_feas))
            self.assertEqual(v['status'], 'optimal', v)
            chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                 capture_output=True, text=True, timeout=15)
            self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
            self.assertTrue(json.loads(chk.stdout)['checked_feasible'])
        # Infeasible parallel pair: base_row >= 15 and -2*base_row >= -24 (i.e. base_row <= 12)
        rows_inf = [
            (base_row, 15, None),
            ([-2 * a for a in base_row], -24, None),
        ]
        m_inf = dict(cost=costs, bounds=bounds, integer=integer, rows=rows_inf, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            v, _, _ = self.run_model(Path(tmp), integer_mps(m_inf))
            self.assertEqual(v['status'], 'infeasible', v)

    def test_big_m_coefftight_presolver(self):
        # Mixed-integer model with deliberately loose big-M coefficients on both
        # <= and >= rows (positive and negative binary coefficients) so
        # scip_presol_coefftight.c fires in SCIP_STAGE_PRESOLVING without errors
        # and preserves exact optimality.
        n_bin = 18
        n_cont = 18
        n = n_bin + n_cont
        costs = [5 + (i % 4) for i in range(n_bin)] + [-2 - (i % 3) for i in range(n_cont)]
        bounds = [(0, 1)] * n_bin + [(0, 3)] * n_cont
        integer = [True] * n_bin + [False] * n_cont
        rows = []
        # 1) sum_{k} y_k - 1000 * x_i <= 0 (loose negative big-M on <= row)
        for i in range(0, n_bin, 2):
            a = [0] * n
            a[i] = -1000
            a[n_bin + i] = 1
            a[n_bin + i + 1] = 1
            rows.append((a, None, 0))
        # 2) 1000 * x_i + y_i + y_{i+1} <= 1003 (loose positive big-M on <= row)
        for i in range(1, n_bin, 2):
            a = [0] * n
            a[i] = 1000
            a[n_bin + i] = 1
            a[n_bin + (i - 1)] = 1
            rows.append((a, None, 1003))
        # 3) 1000 * x_i + y_i >= 2 (loose positive big-M on >= row)
        for i in range(0, n_bin, 3):
            a = [0] * n
            a[i] = 1000
            a[n_bin + i] = 1
            rows.append((a, 2, None))
        # 4) -1000 * x_i + y_i >= -998 (loose negative big-M on >= row)
        for i in range(1, n_bin, 3):
            a = [0] * n
            a[i] = -1000
            a[n_bin + i] = 1
            rows.append((a, -998, None))
        # Budget row
        rows.append(([1] * n_bin + [0] * n_cont, 3, 8))
        m = dict(cost=costs, bounds=bounds, integer=integer, rows=rows, offset=0, maximize=False)
        with tempfile.TemporaryDirectory() as tmp:
            text = integer_mps(m)
            v, model, sol = self.run_model(Path(tmp), text)
            self.assertEqual(v['status'], 'optimal', v)
            self.assertIn('SCIP branch-and-cut', v['message'])
            chk = subprocess.run([CHECKER, '--model', model, '--solution', sol],
                                 capture_output=True, text=True, timeout=15)
            self.assertEqual(chk.returncode, 0, (chk.stdout, chk.stderr))
            c = json.loads(chk.stdout)
            self.assertTrue(c['checked_feasible'], c)
            # Fixing each binary pair leaves two box-bounded continuous
            # variables and an integer upper bound on their sum. Enumerate
            # their integral vertices and combine independent pairs by budget.
            dp = {0: 0}
            for i in range(0, n_bin, 2):
                after = {}
                for x0, x1, y0, y1 in product(range(2), range(2), range(4), range(4)):
                    if y0 + y1 > min(1000 * x0, 1003 - 1000 * x1):
                        continue
                    if any(j % 3 == 0 and 1000 * x + y < 2 or
                           j % 3 == 1 and y - 1000 * x < -998
                           for j, x, y in [(i, x0, y0), (i + 1, x1, y1)]):
                        continue
                    cost = sum(a * b for a, b in zip((x0, x1, y0, y1),
                               (costs[i], costs[i+1], costs[n_bin+i], costs[n_bin+i+1])))
                    for used, value in dp.items():
                        count = used + x0 + x1
                        after[count] = min(after.get(count, value + cost), value + cost)
                dp = after
            self.assertEqual(F(c['objective']), min(v for k, v in dp.items() if 3 <= k <= 8))


if __name__ == '__main__':
    unittest.main()
