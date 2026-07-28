#!/usr/bin/env python3
"""Regenerate the scaling analysis of the Hall near-miss census.

Run from the repository root:

    python analysis/scaling_fit.py

Emits, in order, every number reported in Sections 5 and 8 of the note:

  1. the cell table  (observed, uniform expectation, family count, D, sigma, s)
  2. the matched-pair collapse test              (Section 5.1)
  3. the seven-, nine- and eleven-cell fits      (Section 5.2)
  4. the theta = 1 column and its shortfall      (Section 5.3)
  5. within-column local exponents               (Section 8)
  6. decade-10 extrapolations and standardized residuals

Family counts are recomputed from analysis/family_enumerate.py rather than
hard-coded, so a change to the enumerator propagates here.  Observed counts are
read from the census state file if given, else taken from the published table.
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_enumerate import family_best

THETAS = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

PUBLISHED = {
    0.8: [2, 5, 10, 18, 49, 95, 222, 468, 929, 2050, 4180],
    0.9: [3, 6, 15, 54, 133, 360, 957, 2458, 6188, 15640, 39691],
    1.0: [3, 14, 51, 187, 604, 1981, 6381, 20526, 64908, 206233, 653324],
}
DISPERSION = 0.866          # measured aggregate sd / Poisson sd, Section 6
RUNG1 = {9: 68377, 10: 216227}


def E(th, d):
    a = th - 0.5
    return math.log(10) if abs(a) < 1e-12 else (10 ** ((d + 1) * a) - 10 ** (d * a)) / a


def log10s(th, d):
    return (th - 1.5) * (d + 0.5)


def load_obs(path):
    if not path or not os.path.exists(path):
        print(f"[state file {path!r} not found; using the published table]")
        return PUBLISHED
    try:
        st = json.load(open(path))
    except Exception as ex:
        print(f"[could not read {path}: {ex}; using the published table]")
        return PUBLISHED
    obs = {th: [] for th in (0.8, 0.9, 1.0)}
    for d in range(11):
        row = st.get(str(d))
        for th in obs:
            obs[th].append(row[THETAS.index(th)] if row else None)
    return obs


def fit(cells, weighted=False):
    xs = [c["ls"] for c in cells]
    ys = [math.log10(c["D"]) for c in cells]
    w = [(math.log(10) * c["D"] / c["sig"]) ** 2 if weighted else 1.0 for c in cells]
    S = sum(w); Sx = sum(a * b for a, b in zip(w, xs)); Sy = sum(a * b for a, b in zip(w, ys))
    Sxx = sum(a * b * b for a, b in zip(w, xs)); Sxy = sum(a * b * c for a, b, c in zip(w, xs, ys))
    den = S * Sxx - Sx * Sx
    m = (S * Sxy - Sx * Sy) / den
    b = (Sxx * Sy - Sx * Sxy) / den
    res = [y - (m * x + b) for x, y in zip(xs, ys)]
    n = len(cells)
    s2 = sum(a * r * r for a, r in zip(w, res)) / (n - 2) if n > 2 else 0.0
    sm = math.sqrt(s2 * S / den) if n > 2 else 0.0
    ybar = sum(a * y for a, y in zip(w, ys)) / S
    sst = sum(a * (y - ybar) ** 2 for a, y in zip(w, ys))
    ssr = sum(a * r * r for a, r in zip(w, res))
    return m, sm, 10 ** b, 1 - ssr / sst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state", nargs="?", default="data/hall/state_full.json")
    a = ap.parse_args()
    obs = load_obs(a.state)

    fam = {}
    for d in range(4, 11):
        best = family_best(10 ** d, 10 ** (d + 1), theta=1.0)
        fam[d] = {}
        for th, (p, q) in ((0.8, (4, 5)), (0.9, (9, 10)), (1.0, (1, 1))):
            fam[d][th] = sum(1 for x, k in best.items() if k ** q <= x ** p)

    print("1. Cells (family-subtracted depletion)")
    print(f"   {'cell':>14} {'obs':>8} {'E':>10} {'fam':>7} {'D %':>7} {'sig %':>6} "
          f"{'log10 s':>8}  eligible")
    cells = {}
    for th in (0.8, 0.9):
        for d in range(4, 11):
            o = obs[th][d]
            if o is None:
                continue
            e = E(th, d); f = fam[d][th]
            D = 1 - (o - f) / e
            sig = math.sqrt(o) / e
            ok = e >= 100 and d <= 9
            cells[(th, d)] = dict(th=th, d=d, ls=log10s(th, d), D=D, sig=sig,
                                  obs=o, E=e, fam=f, elig=ok)
            print(f"   (t={th}, d={d:2d}) {o:8d} {e:10.1f} {f:7d} {100*D:7.2f} "
                  f"{100*sig:6.2f} {log10s(th, d):8.2f}  {'yes' if ok else '-'}")

    print("\n2. Matched-pair collapse test (Section 5.1)")
    for d8 in range(5, 10):
        c8 = cells.get((0.8, d8))
        if not c8 or not c8["elig"]:
            continue
        near = [cells[(0.9, d)] for d in range(4, 10)
                if (0.9, d) in cells and cells[(0.9, d)]["elig"]
                and abs(cells[(0.9, d)]["ls"] - c8["ls"]) <= 0.30]
        for c9 in near:
            diff = 100 * (c8["D"] - c9["D"])
            sd = 100 * math.hypot(c8["sig"], c9["sig"])
            print(f"   log s {c8['ls']:.2f} / {c9['ls']:.2f}:  "
                  f"D(0.8)={100*c8['D']:6.2f} +/- {100*c8['sig']:.2f}   "
                  f"D(0.9)={100*c9['D']:6.2f} +/- {100*c9['sig']:.2f}   "
                  f"ratio {c8['D']/c9['D']:.2f}   {abs(diff)/sd:.2f} sigma")

    SEVEN = [(0.9, d) for d in (6, 7, 8, 9)] + [(0.8, d) for d in (7, 8, 9)]
    NINE = SEVEN + [(0.8, 5), (0.8, 6)]
    ELEVEN = NINE + [(0.9, 4), (0.9, 5)]

    print("\n3. Fits and selection sensitivity (Section 5.2)")
    gammas = []
    for name, keys in (("seven-cell (superseded)", SEVEN),
                       ("nine-cell", NINE),
                       ("eleven-cell (frozen rule)", ELEVEN)):
        sel = [cells[k] for k in keys]
        g, sg, c, r2 = fit(sel)
        gw, sgw, _, _ = fit(sel, weighted=True)
        gammas.append(g)
        print(f"   {name:26s} n={len(sel):2d}  gamma={g:.3f} +/- {sg:.3f}  "
              f"c={c:.3f}  R2={r2:.3f}   weighted gamma={gw:.3f} +/- {sgw:.3f}")
    print(f"   spread across cell sets: {max(gammas)-min(gammas):.3f}"
          f"  (against a within-fit error of 0.013 on the seven-cell value)")

    print("\n4. The theta = 1 column (Section 5.3)")
    g11, _, c11, _ = fit([cells[k] for k in ELEVEN])
    for d in (9, 10):
        o = obs[1.0][d]
        if o is None:
            continue
        e = E(1.0, d); r1 = RUNG1.get(d, 0); r2 = fam[d][1.0]
        ls = log10s(1.0, d)
        Dfit = c11 * 10 ** (g11 * ls)
        need = o - e * (1 - Dfit)
        print(f"   decade {d}: observed {o}  E={e:.0f}  rung1={r1}  rung2={r2}")
        print(f"     D from enumerated families: {100*(1-(o-r1-r2)/e):.2f}%   "
              f"D from the eleven-cell fit at log s={ls:.2f}: {100*Dfit:.2f}%")
        print(f"     family count the fit would require: {need:.0f}"
              f"   -> unenumerated shortfall {need-r1-r2:.0f}")

    print("\n5. Within-column local exponents, theta = 0.9 (Section 8)")
    prev = None
    for d in range(4, 10):
        c1, c2 = cells[(0.9, d)], cells[(0.9, d + 1)]
        dy = math.log10(c2["D"]) - math.log10(c1["D"])
        dx = c2["ls"] - c1["ls"]
        sy = math.hypot(c1["sig"] / (c1["D"] * math.log(10)),
                        c2["sig"] / (c2["D"] * math.log(10)))
        g, sg = dy / dx, sy / abs(dx)
        print(f"   d{d} -> d{d+1}:  gamma_loc = {g:.3f} +/- {sg:.3f}")
        if d == 9 and prev:
            diff = g - prev[0]; sd = math.hypot(sg, prev[1])
            print(f"   steepening over the previous step: {diff:+.3f} +/- {sd:.3f}"
                  f"  ->  {abs(diff)/sd:.2f} sigma"
                  f"  ({abs(diff)/(sd*DISPERSION):.2f} sigma dispersion-corrected)")
        prev = (g, sg)

    col = [cells[(0.9, d)] for d in (6, 7, 8, 9)]
    gc, sgc, cc, _ = fit(col)
    print(f"   theta=0.9 column alone, decades 6-9: gamma = {gc:.3f} +/- {sgc:.3f}")

    print("\n6. Decade-10 extrapolations")
    print(f"   {'fit':26s} {'D pred':>8} {'cell':>9} {'z (Poisson)':>12} {'z (calib)':>10}")
    print(f"   {'registered (as committed)':26s} {5.00:7.2f}% {39190:9d} "
          f"{(obs[0.9][10]-39190)/math.sqrt(39190):+12.2f} "
          f"{(obs[0.9][10]-39190)/math.sqrt(39190)/DISPERSION:+10.2f}")
    for name, g, c in (("seven-cell", *fit([cells[k] for k in SEVEN])[::2][:2]),
                       ("nine-cell", *fit([cells[k] for k in NINE])[::2][:2]),
                       ("eleven-cell", *fit([cells[k] for k in ELEVEN])[::2][:2]),
                       ("theta=0.9 column", gc, cc)):
        ls = log10s(0.9, 10)
        Dp = c * 10 ** (g * ls)
        cell = E(0.9, 10) * (1 - Dp) + fam[10][0.9]
        z = (obs[0.9][10] - cell) / math.sqrt(cell)
        print(f"   {name:26s} {100*Dp:7.2f}% {cell:9.0f} {z:+12.2f} {z/DISPERSION:+10.2f}")


if __name__ == "__main__":
    main()
