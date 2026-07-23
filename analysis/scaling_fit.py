#!/usr/bin/env python3
"""Regenerate Empirical Law 1: the depletion scaling D(s) = c * s^gamma.

Reads the Hall census state.json, subtracts the exact F2 family per decade,
computes the non-family deficit against the uniform model, and fits
log10 D vs log10 s over the seven cells with model expectation >= ~100
(theta=0.9, decades 6-9; theta=0.8, decades 7-9).

Usage: python3 scaling_fit.py path/to/state.json
Expected output: gamma = 0.190 +/- 0.013, c = 0.86, R^2 = 0.98.
"""
import json, math, sys
from family_enumerate import family_best

THETAS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
CELLS = [(0.9, d) for d in (6, 7, 8, 9)] + [(0.8, d) for d in (7, 8, 9)]

def model(dec, th):
    p = th - 1.5
    return ((10 ** (dec + 1)) ** (p + 1) - (10 ** dec) ** (p + 1)) / (p + 1)

def main(path):
    state = {int(k): v for k, v in json.load(open(path)).items()}
    fam_cache = {d: family_best(10 ** d, 10 ** (d + 1)) for d in {d for _, d in CELLS}}
    pts = []
    print(f"{'theta':>6} {'decade':>7} {'family':>7} {'obs':>7} {'model':>9} {'D':>7}")
    for th, d in CELLS:
        obs = state[d][THETAS.index(th)]
        fam = sum(1 for x, k in fam_cache[d].items() if k <= x ** th)
        m = model(d, th)
        D = (m - (obs - fam)) / m
        s = (d + 0.5) * (th - 1.5)
        pts.append((s, D))
        print(f"{th:>6} {d:>7} {fam:>7} {obs:>7} {m:>9.1f} {D:>+7.1%}")
    xs = [p[0] for p in pts]; ys = [math.log10(p[1]) for p in pts]
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((a - mx) ** 2 for a in xs)
    slope = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sxx
    inter = my - slope * mx
    resid = [b - (slope * a + inter) for a, b in zip(xs, ys)]
    se = math.sqrt(sum(r * r for r in resid) / (n - 2) / sxx)
    r2 = 1 - sum(r * r for r in resid) / sum((b - my) ** 2 for b in ys)
    print(f"\nD(s) = c * s^gamma:  gamma = {slope:.3f} +/- {se:.3f},"
          f"  c = {10 ** inter:.2f},  R^2 = {r2:.3f}")

if __name__ == "__main__":
    main(sys.argv[1])
