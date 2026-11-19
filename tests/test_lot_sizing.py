"""Compare line production and shared openings to enumeration on the original MILP."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import random
import subprocess
import sys
import tempfile
import unittest

from test_mip_solver import integer_mps

SOLVER, CHECKER = (str(Path(s).resolve()) for s in sys.argv[1:3])
del sys.argv[1:3]


def instance(rng, sizes, maximize=False):
    sites = sum(sizes)
    cols = 4 * sites - 2 * len(sizes)
    groups = max(sizes)
    rows = [[0] * cols for _ in range(2 * sites + groups)]
    bounds, sides, cost, integer = [], [], [], []
    gate_groups = []
    demands = [rng.randint(1, 6) for _ in range(sites)]
    rates = [rng.randrange(8) for _ in range(sites)]
    setups = [rng.randrange(16) for _ in range(sites)]
    sides.extend((d, d) for d in demands)
    sides.extend([(None, 0)] * sites + [(None, 1)] * groups)
    first = 0
    arc_costs = []
    for n in sizes:
        cap = sum(demands[first:first+n])
        labels = list(range(groups))
        rng.shuffle(labels)
        right, left = [], []
        for i in range(n):
            v = first + i
            group = labels[i]
            gate_groups.append(group)
            j = len(bounds)
            rows[v][j] = rows[sites+v][j] = 1
            rows[sites+v][j+1] = -cap
            rows[2*sites+group][j+1] = 1
            bounds.extend([(0, cap), (0, 1)])
            cost.extend([rates[v], setups[v]])
            integer.extend([False, True])
            if i+1 < n:
                j += 2
                rows[v][j], rows[v+1][j] = -1, 1
                rows[v][j+1], rows[v+1][j+1] = 1, -1
                right.append(rng.randrange(5))
                left.append(rng.randrange(5))
                cost.extend([right[-1], left[-1]])
                bounds.extend([(0, cap)] * 2)
                integer.extend([False] * 2)
        arc_costs.append((right, left))
        first += n
    # Enumerate every legal opening pattern; independently route each demand
    # unit to its cheapest available supplier for this choice.
    best = None
    for bits in product((0, 1), repeat=sites):
        groups_open = [g for g, bit in zip(gate_groups, bits) if bit]
        if len(set(groups_open)) != len(groups_open):
            continue
        val = sum(f * b for f, b in zip(setups, bits))
        first = 0
        for n, (right, left) in zip(sizes, arc_costs):
            suppliers = [j for j in range(n) if bits[first+j]]
            if not suppliers:
                break
            for i in range(n):
                val += demands[first+i] * min(
                    rates[first+j] + (sum(right[j:i]) if j<i else sum(left[i:j]))
                    for j in suppliers)
            first += n
        else:
            best = val if best is None else min(best, val)
    offset = F(-17, 8)
    assert best is not None
    m = dict(cost=[-c if maximize else c for c in cost], bounds=bounds, integer=integer,
             rows=[(row, lo, hi) for row, (lo, hi) in zip(rows, sides)],
             offset=offset, maximize=maximize)
    return m, offset + (-best if maximize else best)


class LotSizing(unittest.TestCase):
    def run_case(self, m, options=()):
        with tempfile.TemporaryDirectory() as tmp:
            model, point = (Path(tmp)/p for p in ('original.mps', 'witness.sol'))
            model.write_text(integer_mps(m))
            r = subprocess.run([SOLVER, model, '--mip', '--solution', point, *options],
                               capture_output=True, text=True, timeout=15)
            self.assertIn(r.returncode, (0, 4), (r.stdout, r.stderr))
            self.assertNotIn('runtime error:', r.stderr)
            status = json.loads(r.stdout)
            if not point.exists():
                return status, None
            # The integer path reconstructs an exact point. A limited fallback
            # search may export a continuous LP incumbent at the normal tolerance.
            flags = ['--exact'] if 'Line production' in status['message'] else []
            c = subprocess.run([CHECKER, *flags, '--model', model, '--solution', point],
                               capture_output=True, text=True, timeout=10)
            self.assertEqual(c.returncode, 0, c.stdout)
            return status, F(json.loads(c.stdout)['objective'])

    def test_exhaustive_openings(self):
        rng = random.Random(12391)
        for i in range(40):
            sizes = ([1, 3], [2, 3], [4, 3, 3], [6])[i % 4]
            m, want = instance(rng, sizes, maximize=bool(i % 2))
            with self.subTest(i=i):
                r, value = self.run_case(m)
                self.assertIn('Line production', r['message'])
                self.assertEqual((r['status'], value), ('optimal', want))

    def test_capacity_and_decimal_guards(self):
        for mode in ['supply', 'transport', 'setup', 'coefficient']:
            m, _ = instance(random.Random(35), [3, 3])
            if mode == 'supply':
                m['bounds'][0] = (0, 1)
            elif mode == 'transport':
                m['bounds'][2] = (0, 1)
            elif mode == 'setup':
                m['cost'][1] = F(10000000000000001, 10**16)
            else:
                # Keep the irrelevant tiny source discrepancy in the
                # at-most-one coefficient, where rounding would change it.
                for row, _, hi in m['rows']:
                    if hi == 1 and row[1]:
                        row[1] = F(10000000000000001, 10**16)
            r, _ = self.run_case(m, ['--node-limit', '12'])
            self.assertNotIn('Line production', r['message'])

    def test_exhaustion_keeps_status_open(self):
        m, expected = instance(random.Random(42), [4, 4, 4])
        r, value = self.run_case(m, ['--node-limit', '1'])
        self.assertNotIn('Line production', r['message'])
        if r['status'] == 'optimal':
            self.assertEqual(value, expected)


if __name__ == '__main__':
    unittest.main()
