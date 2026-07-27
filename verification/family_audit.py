#!/usr/bin/env python
"""
Independent audit of the rung-2 Hall near-miss family enumeration.

Shares no code with family_enumerate.py.  Enumerates over (w, t | 6w^2),
forms x = t^2 +/- 4w, and computes k(x) DIRECTLY from the nearest square
(exact integer arithmetic, no identity assumed, no floating logarithms).

Certified w-bound (proved):
    plus  branch  k = 8w^3 + j^2 >= 8w^3          =>  w <= (1/8)^(1/3) x^(t/3) = 0.500 x^(t/3)
    minus branch  k = 4w^3|2x-w|/(x+4w)
                  w <= x  =>  k >= (4/5)w^3       =>  w <= (5/4)^(1/3) x^(t/3) = 1.077 x^(t/3)
                  w >  x  =>  k >  (4/5)x^2 > x^t unless w = 2x, which forces x = (t/3)^2 (a square, k=0)
W_SCALE below is deliberately generous so the bound itself is tested, not assumed.
"""
import sys
from math import isqrt

W_SCALE = 3.0          # margin over the certified 1.077
THETAS = {"0.8": (4, 5), "0.9": (9, 10)}   # k <= x^(num/den) <=> k**den <= x**num


def spf_sieve(n):
    s = list(range(n + 1))
    for i in range(2, isqrt(n) + 1):
        if s[i] == i:
            for j in range(i * i, n + 1, i):
                if s[j] == j:
                    s[j] = i
    return s


def factor(n, spf):
    f = {}
    while n > 1:
        p = spf[n]
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    return f


def divisors(fac):
    ds = [1]
    for p, e in fac.items():
        ds = [d * p**i for d in ds for i in range(e + 1)]
    return ds


def k_exact(x):
    """k(x) = min over y of |y^2 - x^3|, computed from the nearest square. 0 iff x is a square."""
    c = x**3
    y = isqrt(c)
    return min(c - y * y, (y + 1) ** 2 - c)


def le_pow(k, x, num, den):
    """exact test  k <= x^(num/den)   <=>   k**den <= x**num"""
    return k**den <= x**num


def enumerate_families(lo, hi):
    w_max = int(W_SCALE * hi ** (0.9 / 3)) + 2
    spf = spf_sieve(w_max + 1)
    t_max = isqrt(hi) + 2

    def fac6w2(w):
        f = {p: 2 * e for p, e in factor(w, spf).items()}
        f[2] = f.get(2, 0) + 1
        f[3] = f.get(3, 0) + 1
        return f

    hits = {th: {} for th in THETAS}     # theta -> {x: (w, t, j, branch, k)}
    maxw = {th: 0 for th in THETAS}
    pairs = 0

    for w in range(1, w_max + 1):
        for t in divisors(fac6w2(w)):
            if t > t_max:
                continue
            j = 6 * w * w // t
            for sgn, tag in ((1, "+"), (-1, "-")):
                x = t * t + sgn * 4 * w
                if x < lo or x >= hi:
                    continue
                pairs += 1
                k = k_exact(x)
                if k == 0:
                    continue                      # x is a perfect square: excluded
                for th, (num, den) in THETAS.items():
                    if le_pow(k, x, num, den):
                        hits[th][x] = (w, t, j, tag, k)
                        if w > maxw[th]:
                            maxw[th] = w
    return hits, maxw, w_max, pairs


if __name__ == "__main__":
    lo = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**10
    hi = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**11

    hits, maxw, w_max, pairs = enumerate_families(lo, hi)

    print(f"range [{lo:.6g}, {hi:.6g})   w scanned to {w_max}   (w,t) pairs landing in range: {pairs}")
    for th in sorted(THETAS):
        d = hits[th]
        nplus = sum(1 for v in d.values() if v[3] == "+")
        nminus = len(d) - nplus
        cert = 1.077 * hi ** (float(th) / 3)
        print(f"  theta={th}: {len(d):6d} distinct x   (+branch {nplus}, -branch {nminus})"
              f"   max qualifying w = {maxw[th]}   certified bound = {cert:.0f}"
              f"   plus-branch-only bound = {0.5*hi**(float(th)/3):.0f}")
        above = sum(1 for v in d.values() if v[0] > 0.5 * hi ** (float(th) / 3))
        print(f"           members with w above the plus-branch bound: {above}")
