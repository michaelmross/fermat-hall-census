#!/usr/bin/env python3
"""Enumerate the F2 family of Hall near-misses (Proposition 2 of the paper).

Members: x = t^2 +- 4w with t | 6w^2, j = 6w^2/t, and identity residual
  k_id = 8w^3 + j^2   (plus branch)   or   |8w^3 - j^2|   (minus branch).
Band membership uses the deterministic identity criterion k_id <= x^theta,
tested in exact integer arithmetic (no floating logarithms).
Dependency-free (stdlib only).

Usage:
  python family_enumerate.py --x-lo 1e9 --x-hi 1e10 --thetas 0.8,0.9
Reproduces the preregistered family counts (dec 9: 465/1568; dec 10: 832/3283).
"""
import argparse
from fractions import Fraction


def divisors(n):
    f = {}
    d, m = 2, n
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1; m //= d
        d += 1 if d == 2 else 2
    if m > 1: f[m] = f.get(m, 0) + 1
    divs = [1]
    for p, e in f.items():
        divs = [q * p ** k for q in divs for k in range(e + 1)]
    return divs


def family_best(x_lo, x_hi, theta=1.0, w_max=None):
    """Return {x: min identity residual} for all F2 members in [x_lo, x_hi).

    The result is COMPLETE for every threshold exponent <= theta.  Callers that
    count several columns must pass the largest exponent they intend to use.

    Certified cutoff (Lemma): every qualifying member has
      plus  branch:  k = 8w^3 + j^2 >= 8w^3            =>  w <= (1/8)^(1/3) x^(theta/3)
      minus branch:  k = 4w^3|2x-w|/(x+4w) >= (4/5)w^3 for w <= x
                     (w > x forces k > (4/5)x^2 > x^theta for theta <= 1, x >= 2,
                      unless w = 2x, which forces x = (t/3)^2, a square with k = 0)
                                                       =>  w <= (5/4)^(1/3) x^(theta/3)
    The minus-branch constant governs.  Corollary: no qualifying member is a
    perfect square, since x = t^2 +- 4w square requires 4w >= 2t + 1, i.e.
    w >~ sqrt(x)/2, which the cutoff excludes for all x above a small bound.
    """
    if w_max is None:
        w_max = int((1.25 * x_hi ** theta) ** (1 / 3)) + 2
    best = {}
    for w in range(1, w_max):
        w6 = 6 * w * w
        for t in divisors(w6):
            j = w6 // t
            for sgn in (+1, -1):
                x = t * t + sgn * 4 * w
                if x_lo <= x < x_hi:
                    k = 8 * w ** 3 + j * j if sgn > 0 else abs(8 * w ** 3 - j * j)
                    if k > 0 and (x not in best or k < best[x]):
                        best[x] = k
    return best


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--x-lo", type=float, required=True)
    ap.add_argument("--x-hi", type=float, required=True)
    ap.add_argument("--thetas", default="0.8,0.9")
    ap.add_argument("--w-max", type=int, default=None)
    a = ap.parse_args()

    thetas = [Fraction(tok.strip()) for tok in a.thetas.split(",")]
    theta_enum = float(max(thetas))          # enumerate deep enough for every column

    best = family_best(int(a.x_lo), int(a.x_hi), theta=theta_enum, w_max=a.w_max)
    print(f"family members in [{a.x_lo:.3g}, {a.x_hi:.3g}): {len(best)}"
          f"   (enumerated to theta = {theta_enum})")
    for th in thetas:
        p, q = th.numerator, th.denominator      # k <= x^(p/q)  <=>  k**q <= x**p
        n = sum(1 for x, k in best.items() if k ** q <= x ** p)
        print(f"  members with k_id <= x^{float(th)}: {n}")
