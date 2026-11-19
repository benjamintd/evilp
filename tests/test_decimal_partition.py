"""Check exact decimal objective conversion against subset enumeration."""
import random
import unittest
from fractions import Fraction

import test_partition as base


class DecimalPartition(unittest.TestCase):
    run_case = base.Partition.run_case

    def test_signed_costs_and_repaired_prices(self):
        rng = random.Random(99236)
        for i in range(40):
            n = rng.randint(2, 7)
            masks = [1 << j for j in range(n)] + [rng.randrange(1, 1 << n) for _ in range(14)]
            costs = [Fraction(rng.randrange(-50, 60), rng.choice([4, 10, 100, 1000000]))
                     for _ in masks]
            maximize = bool(i % 2)
            optimum = base.oracle(n, masks, costs, maximize)
            with self.subTest(i=i):
                r = self.run_case(base.model(n, masks, costs, maximize), mode=i % 4)
                self.assertTrue(r['optimal'])
                self.assertEqual(sum(costs[j] for j in r['columns']), optimum)

    def test_close_costs_and_scale_limits(self):
        r = self.run_case(base.model(2, [1, 2, 3],
                                     [Fraction(1, 8), Fraction(7, 8), Fraction(999999, 1000000)]))
        self.assertEqual(r['columns'], [2])
        # The enlarged accumulator can now rank the formerly oversized cost.
        recovered = self.run_case(base.model(2, [1, 2, 3], [Fraction(1, 1000000), 2000, 1]))
        self.assertTrue(recovered['optimal'])
        self.assertEqual(recovered['columns'], [2])
        # Mixed exact denominators and an oversized scaled objective must abstain.
        for costs in ([Fraction(1, 10000000000), 0, 1],
                      [Fraction(1, 1000000), 300000, 1],
                      [Fraction('1.0000000000000000001'), Fraction(1, 2), 1]):
            with self.subTest(costs=costs):
                r = self.run_case(base.model(2, [1, 2, 3], costs))
                self.assertEqual((r['nodes'], r['point']), (0, False))


if __name__ == '__main__':
    unittest.main()
