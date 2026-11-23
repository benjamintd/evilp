"""Check all separated cliques against exact feasible arrival intervals."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from test_mip_solver import integer_mps

DRIVER = Path(sys.argv.pop(1)).resolve()


class PathCliques(unittest.TestCase):
    def separate(self, model):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / 'paths.mps'
            src.write_text(integer_mps(model))
            run = subprocess.run([DRIVER, src], check=True, capture_output=True, text=True)
            return [tuple(int(name[1:]) for name in line.split())
                    for line in run.stdout.splitlines() if line]

    def test_cliques_preserve_all_feasible_binary_choices(self):
        eps = F(1, 10**19)
        for gap, reverse, packing, peer in product((-eps, 0, eps, 2), (False, True),
                                                   (False, True), (False, True)):
            lags = [4+gap, 5, 3, 4]
            rows = []
            if packing:
                rows.extend([([1, 1, 0, 0] + [0]*5, None, 1),
                             ([0, 0, 1, 1] + [0]*5, None, 1)])
            for j, (a, b) in enumerate([(0, 2), (1, 2), (2, 3), (2, 4)]):
                co = [0]*9
                co[j], co[4+a], co[4+b] = -30, -1, 1
                if peer and j == 0:
                    co[1] = 9
                side = lags[j]-30
                rows.append(([-w for w in co], None, -side) if reverse else (co, side, None))
            model = dict(cost=[0]*9, integer=[True]*4 + [False]*5, offset=0, maximize=False,
                         bounds=[(0, 1)]*4 + [(0, 0), (0, 0), (0, 10), (7, 7), (7, 7)], rows=rows)
            with self.subTest(gap=gap, reverse=reverse, packing=packing, peer=peer):
                cuts = self.separate(model)
                for x in product((0, 1), repeat=4):
                    if packing and (sum(x[:2]) > 1 or sum(x[2:]) > 1):
                        continue
                    low = max(0, lags[0]-30*(1-x[0])-(9*x[1] if peer else 0),
                              lags[1]-30*(1-x[1]))
                    high = min(10, 7-lags[2]+30*(1-x[2]), 7-lags[3]+30*(1-x[3]))
                    if low <= high:
                        for cut in cuts:
                            self.assertLessEqual(sum(x[j] for j in cut), 1)
                if packing:
                    self.assertEqual(any(0 in c and 2 in c for c in cuts), gap > 0)

    def test_source_coefficient_just_above_one(self):
        # Both choices fit with the common time at exactly four. The first
        # constraint's positive coefficient rounds to one in binary64.
        eps = F(1, 10**19)
        model = dict(cost=[0]*5, integer=[True, True, False, False, False],
                     bounds=[(0, 1)]*2 + [(0, 0), (0, 10), (7, 7)],
                     offset=0, maximize=False, rows=[
                         ([-30, 0, -1, 1+eps, 0], -26+eps, None),
                         ([0, -30, 0, -1, 1], -27, None)])
        self.assertEqual(self.separate(model), [])


if __name__ == '__main__':
    unittest.main()
