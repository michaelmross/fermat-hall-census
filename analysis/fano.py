#!/usr/bin/env python3
"""Dispersion diagnostics for the Hall near-miss shrinking-target counts.

Measures the index of dispersion (Fano factor) of the family-subtracted counts
of {x : k(x) <= x^theta} over equal-model-mass blocks, together with the block
autocorrelation, a Bartlett long-run-variance inflation factor, and a
moving-block bootstrap for the aggregate count.

The Poisson surrogate used in the census assumes Var(N) = E(N), i.e. F = 1.
This script tests that assumption at values of s = x^(theta-3/2) where events
are dense enough to measure it, so that the extrapolation to the census cells
rests on a measured F(s) rather than an assumption.

Run from the repository root:

  python analysis/fano.py 7 0.9,1.0
  python analysis/fano.py 6 0.9,1.0

Reproduces the table of Section 6.  A sub-range run at theta = 0.9 over
[1e8, 4e8) supplies the third point; see the --x-lo/--x-hi options.
"""
import argparse, os, random, statistics as st, sys
from fractions import Fraction
from math import isqrt, sqrt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
from family_enumerate import family_best


def scan(lo, hi, thetas):
    """Exact k(x) for every x in [lo,hi); return {theta: [x with k(x) <= x^theta]}."""
    hits = {th: [] for th in thetas}
    pq = [(th, th.numerator, th.denominator) for th in thetas]
    for x in range(lo, hi):
        c = x * x * x
        y = isqrt(c)
        k = c - y * y
        k2 = (y + 1) * (y + 1) - c
        if k2 < k:
            k = k2
        if 0 < k <= x:                      # exact theta=1 gate; k=0 iff x square
            for th, p, q in pq:
                if k ** q <= x ** p:
                    hits[th].append(x)
    return hits


def mass(x, th, lo):
    """cumulative uniform-model mass integral_lo^x u^(theta-3/2) du"""
    a = float(th) - 0.5
    return (x ** a - lo ** a) / a


def blocks_equal_mass(xs, lo, hi, th, B):
    """bin x-values into B blocks of equal uniform-model mass"""
    tot = mass(hi, th, lo)
    counts = [0] * B
    for x in xs:
        b = int(B * mass(x, th, lo) / tot)
        counts[min(b, B - 1)] += 1
    return counts


def autocorr(v, lag):
    n = len(v); m = sum(v) / n
    num = sum((v[i] - m) * (v[i + lag] - m) for i in range(n - lag))
    den = sum((vi - m) ** 2 for vi in v)
    return num / den if den else 0.0


def bartlett_lrv(v, L):
    """long-run variance with Bartlett kernel; returns LRV and inflation vs gamma_0"""
    n = len(v); m = sum(v) / n
    g = lambda l: sum((v[i] - m) * (v[i + l] - m) for i in range(n - l)) / n
    g0 = g(0)
    lrv = g0 + 2 * sum((1 - l / (L + 1)) * g(l) for l in range(1, L + 1))
    return lrv, (lrv / g0 if g0 else float('nan'))


def moving_block_bootstrap(v, block_len, reps=4000, seed=0):
    rng = random.Random(seed)
    n = len(v); nb = n // block_len
    starts = n - block_len + 1
    tots = []
    for _ in range(reps):
        t = 0
        for _ in range(nb):
            s = rng.randrange(starts)
            t += sum(v[s:s + block_len])
        tots.append(t * n / (nb * block_len))
    return st.mean(tots), st.pstdev(tots)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("decade", type=int,
                    help="decade d; scans [10^d, 10^(d+1)) unless --x-lo/--x-hi given")
    ap.add_argument("thetas", default="0.9,1.0", nargs="?")
    ap.add_argument("--blocks", type=int, default=0)
    ap.add_argument("--x-lo", type=float, default=None)
    ap.add_argument("--x-hi", type=float, default=None)
    a = ap.parse_args()

    lo = int(a.x_lo) if a.x_lo else 10 ** a.decade
    hi = int(a.x_hi) if a.x_hi else 10 ** (a.decade + 1)
    thetas = [Fraction(t.strip()) for t in a.thetas.split(",")]

    hits = scan(lo, hi, thetas)
    fam = family_best(lo, hi, theta=float(max(thetas)))

    print(f"decade [{lo:.3g}, {hi:.3g})")
    for th in thetas:
        p, q = th.numerator, th.denominator
        famx = {x for x, k in fam.items() if k ** q <= x ** p}
        t = 2                                   # rung 1: x = t^2 +/- 1, t even
        while t * t <= hi + 1:
            for x in (t * t - 1, t * t + 1):
                if lo <= x < hi:
                    c = x * x * x; y = isqrt(c)
                    kk = min(c - y * y, (y + 1) ** 2 - c)
                    if 0 < kk and kk ** q <= x ** p:
                        famx.add(x)
            t += 2
        xs = [x for x in hits[th] if x not in famx]
        E = mass(hi, th, lo)
        s_mid = (10 ** (a.decade + 0.5)) ** (float(th) - 1.5)

        B = a.blocks or max(20, min(400, len(xs) // 25))
        cnt = blocks_equal_mass(xs, lo, hi, th, B)
        m = sum(cnt) / B
        var = st.pvariance(cnt)
        F = var / m
        seF = sqrt(2.0 / (B - 1))
        r = [autocorr(cnt, l) for l in (1, 2, 3)]
        L = max(1, int(B ** (1 / 3)))
        lrv, infl = bartlett_lrv(cnt, L)
        bl = max(2, B // 20)
        bmean, bsd = moving_block_bootstrap(cnt, bl)

        print(f"\n  theta={float(th)}  log10 s={__import__('math').log10(s_mid):.2f}   "
              f"non-family hits={len(xs)}  (family removed: {len(hits[th])-len(xs)})")
        print(f"    blocks B={B}   mean/block={m:.2f}   var/block={var:.2f}")
        print(f"    Fano F = {F:.3f} +/- {seF:.3f}        (Poisson surrogate assumes F = 1)")
        print(f"    autocorr r1={r[0]:+.3f} r2={r[1]:+.3f} r3={r[2]:+.3f}   "
              f"Bartlett(L={L}) LRV inflation = {infl:.3f}")
        print(f"    aggregate sd: Poisson={sqrt(sum(cnt)):.1f}   "
              f"block-bootstrap={bsd:.1f}   ratio={bsd/sqrt(sum(cnt)):.3f}")
