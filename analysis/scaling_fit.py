#!/usr/bin/env python3
"""Regenerate Empirical Law 1: the depletion scaling D(s) = c * s^gamma.

Reads the Hall census state.json, subtracts the exact F2 family per decade,
computes the non-family deficit with per-cell Poisson uncertainties, and
reports both the unweighted fit (canonical: gamma = 0.191 +/- 0.013, c = 0.86,
R^2 = 0.98) and the Poisson-weighted fit (gamma = 0.21 +/- 0.05,
chi2/dof = 0.04; same-decade cells across columns share events, so the
effective chi2 is optimistic).

Usage: python3 scaling_fit.py path/to/state.json
"""
import json, math, sys
from pathlib import Path
# locate analysis/family_enumerate.py relative to this file, wherever it was copied
_here = Path(__file__).resolve().parent
for _p in [_here] + list(_here.parents):
    if (_p / "analysis" / "family_enumerate.py").exists():
        sys.path.insert(0, str(_p / "analysis")); break
    if (_p / "family_enumerate.py").exists():
        sys.path.insert(0, str(_p)); break
from family_enumerate import family_best

THETAS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
CELLS = [(0.9, d) for d in (6, 7, 8, 9)] + [(0.8, d) for d in (7, 8, 9)]

def model(dec, th):
    p = th - 1.5
    return ((10 ** (dec + 1)) ** (p + 1) - (10 ** dec) ** (p + 1)) / (p + 1)

def fit(pts, weighted):
    xs = [p[0] for p in pts]; ys = [math.log10(p[1]) for p in pts]
    ws = ([1 / ((p[2] / (p[1] * math.log(10))) ** 2) for p in pts]
          if weighted else [1.0] * len(pts))
    W = sum(ws)
    mx = sum(w * a for w, a in zip(ws, xs)) / W
    my = sum(w * b for w, b in zip(ws, ys)) / W
    sxx = sum(w * (a - mx) ** 2 for w, a in zip(ws, xs))
    slope = sum(w * (a - mx) * (b - my) for w, a, b in zip(ws, xs, ys)) / sxx
    inter = my - slope * mx
    if weighted:
        se = math.sqrt(1 / sxx)
        stat = sum(w * (b - (slope * a + inter)) ** 2
                   for w, a, b in zip(ws, xs, ys)) / (len(xs) - 2)
        tag = "chi2/dof"
    else:
        resid = [b - (slope * a + inter) for a, b in zip(xs, ys)]
        se = math.sqrt(sum(r * r for r in resid) / (len(xs) - 2) / sxx)
        stat = 1 - sum(r * r for r in resid) / sum((b - my) ** 2 for b in ys)
        tag = "R^2"
    return slope, se, 10 ** inter, stat, tag

def main(path):
    state = {int(k): v for k, v in json.load(open(path)).items()}
    fam_cache = {d: family_best(10 ** d, 10 ** (d + 1)) for d in {d for _, d in CELLS}}
    pts = []
    print(f"{'theta':>6} {'decade':>7} {'family':>7} {'obs':>7} {'model':>9}"
          f" {'D':>7} {'+/-':>6}")
    for th, d in CELLS:
        obs = state[d][THETAS.index(th)]
        fam = sum(1 for x, k in fam_cache[d].items() if k <= x ** th)
        m = model(d, th)
        D = (m - (obs - fam)) / m
        sD = math.sqrt(obs) / m           # Poisson on the observed count
        pts.append(((d + 0.5) * (th - 1.5), D, sD))
        print(f"{th:>6} {d:>7} {fam:>7} {obs:>7} {m:>9.1f} {D:>+7.1%} {sD:>6.1%}")
    for weighted in (False, True):
        g, se, c, stat, tag = fit(pts, weighted)
        name = "weighted" if weighted else "unweighted"
        print(f"{name:>10}: gamma = {g:.3f} +/- {se:.3f},"
              f"  c = {c:.2f},  {tag} = {stat:.2f}")

if __name__ == "__main__":
    main(sys.argv[1])
