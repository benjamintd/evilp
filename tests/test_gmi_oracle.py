#!/usr/bin/env python3
"""Independent Fraction validity and safe-rounding oracles."""
from fractions import Fraction as Q
import argparse
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import unittest

PARSER = argparse.ArgumentParser()
PARSER.add_argument("--driver", required=True, type=Path)
ARGS, REST = PARSER.parse_known_args()
COUNTS = dict(driver_calls=0, mixed_points=0, split_vertices=0,
              split_branch_certificates=0, rounding_vertices=0)


def checked_zip(*vectors):
    """Retain strict length checks on the supported Python 3.9 runtime."""
    if len({len(vector) for vector in vectors}) != 1:
        raise ValueError("oracle vector lengths differ")
    return zip(*vectors)


def call(mode, coefficients, rhs, metadata):
    lines = [f"{mode} {len(coefficients)}", str(rhs)]
    for a, extra in checked_zip(coefficients, metadata):
        if mode == "gmi":
            lines.append(f"{a} {int(extra)}")
        else:
            lower, upper = extra
            lines.append(f"{a} {'*' if lower is None else lower} "
                         f"{'*' if upper is None else upper}")
    result = subprocess.run([str(ARGS.driver.resolve())], input="\n".join(lines) + "\n",
                            text=True, capture_output=True, timeout=20, check=True)
    if result.stderr:
        raise AssertionError(result.stderr)
    COUNTS["driver_calls"] += 1
    return json.loads(result.stdout)


def build(a, beta, integer):
    result = call("gmi", a, beta, integer)
    if result["status"] != "ok":
        raise AssertionError(result)
    return tuple(map(Q, result["coefficient"])), Q(result["rhs"])


def dot(a, x):
    return sum((p * q for p, q in checked_zip(a, x)), Q(0))


def fractional(q):
    return q - q.numerator // q.denominator


def box_corners(n, upper=Q(4)):
    return itertools.product((Q(0), upper), repeat=n)


def section_vertices(a, rhs, upper=Q(4)):
    """Exact equality/box vertices without using the GMI formula."""
    n = len(a)
    result = set(x for x in box_corners(n, upper) if dot(a, x) == rhs)
    for free in range(n):
        if not a[free]:
            continue
        for fixed in box_corners(n - 1, upper):
            x = list(fixed)
            x.insert(free, Q(0))
            x[free] = (rhs - dot(a, x)) / a[free]
            if 0 <= x[free] <= upper:
                result.add(tuple(x))
    return result


class ExactGmiOracles(unittest.TestCase):
    def test_native_integrality_is_a_real_premise(self):
        a, beta, z = (Q(3, 4),), Q(1, 2), (Q(2, 3),)
        self.assertEqual(dot(a, z), beta)  # x_B=0 is integer; z is continuous.
        valid, rhs = build(a, beta, (0,))
        self.assertGreaterEqual(dot(valid, z), rhs)
        falsely_integer, rhs = build(a, beta, (1,))
        self.assertLess(dot(falsely_integer, z), rhs)

    def test_known_negative_and_mixed(self):
        a = tuple(map(Q, ("-5/4", "-2/3", "-7/5", "2/7")))
        gamma, f = build(a, Q(7, 3), (1, 1, 0, 0))
        self.assertEqual(f, Q(1, 3))
        self.assertEqual(gamma, tuple(map(Q, ("1/8", "1/3", "7/10", "2/7"))))
        self.assertEqual(build((), Q(-7, 3), ()), ((), Q(2, 3)))
        for beta in (Q(-100), Q(0), Q(27)):
            self.assertEqual(call("gmi", a, beta, (1, 1, 0, 0))["status"], "integral_beta")

    def test_mixed_integer_tableau_enumeration(self):
        rng = random.Random(90412)
        for trial in range(96):
            ni, nc = trial % 3, (trial // 3) % 3
            if not ni + nc:
                nc = 1
            integer = (1,) * ni + (0,) * nc
            a = tuple(Q(rng.randint(-13, 13), rng.randint(1, 9)) for _ in integer)
            beta = Q(rng.randint(-13, 13), rng.randint(2, 11))
            if beta.denominator == 1:
                beta += Q(1, 13)
            gamma, f = build(a, beta, integer)
            for zi in itertools.product(map(Q, range(5)), repeat=ni):
                for xb in range(-12, 13):
                    continuous = section_vertices(a[ni:], beta - xb - dot(a[:ni], zi))
                    for zc in continuous:
                        z = zi + zc
                        self.assertEqual(Q(xb) + dot(a, z), beta)
                        self.assertGreaterEqual(dot(gamma, z), f, (a, beta, integer, xb, z))
                        COUNTS["mixed_points"] += 1
        self.assertGreater(COUNTS["mixed_points"], 500)

    def test_exact_unbounded_disjunctive_certificates_and_vertices(self):
        rng = random.Random(104931)
        cases = []
        for _ in range(160):
            n = rng.randint(1, 4)
            a = tuple(Q(rng.randint(-21, 21), rng.randint(1, 13)) for _ in range(n))
            flags = tuple(rng.randrange(2) for _ in a)
            beta = Q(rng.randint(-20, 20), rng.randint(2, 17))
            if beta.denominator == 1:
                beta += Q(1, 19)
            cases.append((a, beta, flags))
        # Near-integer fractions are exact, without a numerical tolerance.
        for beta in (Q(1, 10**100), 1 - Q(1, 10**100), -1 - Q(1, 10**100)):
            cases.append(((Q(-17, 13), Q(4, 7), Q(-10**25, 3)), beta, (1, 1, 0)))
        for a, beta, flags in cases:
            gamma, f = build(a, beta, flags)
            # Absorb integer coefficient parts, then optionally one integer z_j,
            # into x_B. Its exact split is r*z>=f OR r*z<=f-1.
            r = tuple((fractional(aj) - int(fractional(aj) > f)) if integer else aj
                      for aj, integer in checked_zip(a, flags))
            for gj, rj in checked_zip(gamma, r):
                self.assertGreaterEqual(gj, 0)
                self.assertGreaterEqual(gj, rj)
                self.assertGreaterEqual(gj, -f * rj / (1 - f))
            # These two coefficientwise implications prove validity over the
            # unbounded nonnegative split branches, not just sampled points.
            COUNTS["split_branch_certificates"] += 2
            vertices = set(box_corners(len(a)))
            vertices |= section_vertices(r, f)
            vertices |= section_vertices(r, f - 1)
            for z in vertices:
                rz = dot(r, z)
                if rz >= f or rz <= f - 1:
                    self.assertGreaterEqual(dot(gamma, z), f)
                    COUNTS["split_vertices"] += 1

    def test_integer_translation_and_continuous_scaling(self):
        a = (Q(-17, 13), Q(4, 7), Q(-5, 9), Q(2, 3))
        flags = (1, 1, 0, 0)
        beta = Q(-11, 7)
        original, f = build(a, beta, flags)
        for shift in (-100000, -1, 0, 2, 100001):
            translated = (a[0] + shift, a[1] - shift, a[2], a[3])
            self.assertEqual(build(translated, beta + shift, flags), (original, f))
        for scale in (Q(1, 10**50), Q(1, 7), Q(11), Q(10**50)):
            transformed = a[:2] + tuple(scale * q for q in a[2:])
            gamma, rhs = build(transformed, beta, flags)
            self.assertEqual(rhs, f)
            self.assertEqual(gamma, original[:2] + tuple(scale * q for q in original[2:]))


class BinaryRoundingOracles(unittest.TestCase):
    def verify(self, a, b, box, result):
        self.assertEqual(result["status"], "ok")
        d = tuple(map(Q, result["coefficient_binary"]))
        rhs = Q(result["rhs_binary"])
        lower = Q(0)
        for exact, binary, (lo, hi) in checked_zip(a, d, box):
            self.assertTrue(math.isfinite(float(binary)))
            self.assertEqual(binary == 0, exact == 0)
            if exact:
                self.assertEqual(binary > 0, exact > 0)
            delta = binary - exact
            bound = lo if delta > 0 else hi if delta < 0 else Q(0)
            self.assertIsNotNone(bound)
            lower += delta * bound
        self.assertEqual(Q(result["error_lower"]), lower)
        self.assertEqual(Q(result["adjusted_rhs"]), b + lower)
        self.assertLessEqual(rhs, b + lower)
        successor = math.nextafter(float(rhs), math.inf)
        if math.isfinite(successor):
            self.assertGreater(Q(successor), b + lower)
        # All finite vertices satisfy the independently computed perturbation
        # bound; every exact feasible vertex also survives the converted cut.
        if all(lo is not None and hi is not None for lo, hi in box):
            for x in itertools.product(*[(lo, hi) for lo, hi in box]):
                self.assertGreaterEqual(dot(d, x) - dot(a, x), lower)
                if dot(a, x) >= b:
                    self.assertGreaterEqual(dot(d, x), rhs)
                COUNTS["rounding_vertices"] += 1
        return d, rhs

    def test_random_signed_boxes_and_error_compensation(self):
        rng = random.Random(59192)
        for trial in range(144):
            n = rng.randint(1, 4)
            a = tuple(Q(rng.randint(-100, 100), rng.randint(1, 37)) for _ in range(n))
            box = []
            for _ in a:
                lo = Q(rng.randint(-100, 40), rng.randint(1, 11))
                hi = lo + Q(rng.randint(0, 100), rng.randint(1, 13))
                kind = trial % 3
                box.append((lo, hi) if kind == 0 else (lo, None) if kind == 1 else (None, hi))
            b = Q(rng.randint(-40, 40), rng.randint(1, 17))
            result = call("round", a, b, box)
            self.verify(a, b, box, result)

    def test_adversarial_fixed_negative_box_and_naive_rounding(self):
        a, x, b = (Q(1, 10),), (Q(-10**20),), Q(-10**19)
        box = [(x[0], x[0])]
        result = call("round", a, b, box)
        d, rhs = self.verify(a, b, box, result)
        self.assertEqual(dot(a, x), b)
        self.assertLess(dot(d, x), b)  # A naively unchanged rhs is invalid.
        self.assertGreaterEqual(dot(d, x), rhs)
        self.assertLess(Q(result["error_lower"]), -500)
        # Large, cancelling affine offsets and a tiny nonzero coefficient.
        a = (Q(1, 3), Q(-1, 7), Q(1, 10**100))
        box = [(Q(-10**60), Q(10**60)), (Q(10**50), Q(10**50)), (Q(0), Q(10**90))]
        self.verify(a, Q(1, 10**80), box, call("round", a, Q(1, 10**80), box))

    def test_unbounded_sides_exact_zeros_overflow_and_underflow(self):
        for a in (Q(1, 3), Q(-1, 3)):
            self.assertEqual(call("round", [a], Q(1, 3), [(None, None)])["status"],
                             "unbounded_rounding_error")
            for box in ([(Q(0), None)], [(None, Q(0))]):
                self.verify([a], Q(1, 3), box, call("round", [a], Q(1, 3), box))
        self.verify([Q(1, 2), Q(0)], Q(1, 3), [(None, None)] * 2,
                    call("round", [Q(1, 2), Q(0)], Q(1, 3), [(None, None)] * 2))
        tiny, huge = Q(1, 2**2000), Q(2**2000)
        box = [(Q(0), None)]
        result = call("round", [tiny], Q(1, 3), box)
        d, _ = self.verify([tiny], Q(1, 3), box, result)
        self.assertEqual(d[0], Q(1, 2**1074))
        self.assertEqual(call("round", [tiny], Q(1, 3), [(None, Q(0))])["status"],
                         "unbounded_rounding_error")
        self.assertEqual(call("round", [huge], Q(1, 3), box)["status"],
                         "unbounded_rounding_error")
        self.verify([huge], Q(1, 3), [(Q(0), Q(0))],
                    call("round", [huge], Q(1, 3), [(Q(0), Q(0))]))
        self.assertEqual(call("round", [], -huge, [])["status"], "nonfinite_rounding")
        self.verify([], huge, [], call("round", [], huge, []))


if __name__ == "__main__":
    suite = unittest.main(argv=[sys.argv[0]] + REST, exit=False)
    print(json.dumps({"oracle_counts": COUNTS}, sort_keys=True))
    sys.exit(not suite.result.wasSuccessful())
