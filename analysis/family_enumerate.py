#!/usr/bin/env python3
"""Enumerate the F2 family of Hall near-misses (Proposition 2 of the paper).

Members: x = t^2 +- 4w with t | 6w^2, j = 6w^2/t, and identity residual
  k_id = 8w^3 + j^2   (plus branch)   or   |8w^3 - j^2|   (minus branch).
Band membership uses the deterministic identity criterion k_id <= x^theta.
Dependency-free (stdlib only).

Usage:
  python3 family_enumerate.py --x-lo 1e9 --x-hi 1e10 --thetas 0.8,0.9
Reproduces the preregistered family counts (dec 9: 465/1568; dec 10: 832/3283).
"""
import argparse

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

def family_best(x_lo, x_hi, w_max=None):
    """Dict {x: minimal identity residual k_id} over the family, x in [x_lo, x_hi)."""
    if w_max is None:
        w_max = 10 * int((x_hi ** 0.9 / 8) ** (1 / 3)) + 50
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
    best = family_best(int(a.x_lo), int(a.x_hi), a.w_max)
    print(f"family members in [{a.x_lo:.3g}, {a.x_hi:.3g}): {len(best)}")
    for th in (float(t) for t in a.thetas.split(",")):
        n = sum(1 for x, k in best.items() if k <= x ** th)
        print(f"  members with k_id <= x^{th}: {n}")
