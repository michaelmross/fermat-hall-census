#!/usr/bin/env python3
"""
hall_analysis.py -- corrected theta-band table for the Hall census (reads state.json).

Model: E[count of |k| <= x^theta, x in decade] = sum_x min(1, x^(theta-3/2)),
computed exactly (integral with capping; theta=1/2 is the logarithmic column).

Corrections applied / structure flagged:
  F1 (theta=1.0): the first-order family x = t^2 +- 1, t even, k = (3/4)t^2 + O(1),
      contributing 2 * #(even t with t^2 in decade) events. Subtracted exactly.
  F2 (theta=0.8): second-order layer x = t^2 + u, t even, u = 0 mod 4,
      u ~ sqrt(8jt/3), k ~ c * x^(3/4) with quantized c. Detected empirically
      (k/x^0.75 clustering; x = 0 mod 4 density doubling); not yet in closed form,
      so it is FLAGGED (excess reported) rather than subtracted.

Usage: python3 hall_analysis.py state.json [x_max]
"""

import json, math, sys

THETAS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def expected_cell(lo, hi, th):
    """sum_{x=lo}^{hi} min(1, x^(th-3/2)), continuous approximation with capping."""
    p = th - 1.5
    x_cap = 1.0 if p >= 0 else 1.0  # x^p <= 1 iff x >= 1 (p<0): cap only affects x=1
    # p < 0 throughout; capping matters only for x where x^p >= 1, i.e. x = 1
    e = 0.0
    a, b = max(lo, 2), hi
    if a <= b:
        if abs(p + 1) < 1e-12:
            e += math.log((b + 1) / a)
        else:
            e += ((b + 1) ** (p + 1) - a ** (p + 1)) / (p + 1)
    if lo <= 1 <= hi:
        e += 1.0
    return e


def family1_count(lo, hi):
    """#events from x = t^2 +- 1, t even, within [lo, hi] (k = 3x/4 <= x needs x >= 4)."""
    n = 0
    t = 2
    while t * t - 1 <= hi:
        for x in (t * t - 1, t * t + 1):
            if lo <= x <= hi and x >= 44:   # below x=44, k=(3t^2+4)/4 > x fails theta=1
                n += 1
        t += 2
    return n


def main():
    state = {int(k): v for k, v in json.load(open(sys.argv[1])).items()}
    x_max = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** (max(state) )
    print(f"{'decade':>12} " + " ".join(f"{'th='+str(t):>16}" for t in THETAS))
    tot_o = [0.0] * len(THETAS)
    tot_e = [0.0] * len(THETAS)
    for dec in sorted(state):
        lo, hi = 10 ** dec, min(10 ** (dec + 1) - 1, x_max)
        if lo > x_max or sum(state[dec]) == 0:
            continue
        cells = []
        for j, th in enumerate(THETAS):
            o = state[dec][j]
            e = expected_cell(lo, hi, th)
            note = ""
            if th == 1.0:
                f1 = family1_count(lo, hi)
                o_adj = o - f1
                z = (o_adj - e) / math.sqrt(e)
                cells.append(f"{o}-{f1}={o_adj}/{e:.0f}")
                tot_o[j] += o_adj
            else:
                z = (o - e) / math.sqrt(e) if e > 0 else float("nan")
                flag = "*" if th == 0.8 and z > 2 else " "
                cells.append(f"{o:>6}/{e:>7.1f}{flag}({z:+.1f})")
                tot_o[j] += o
            tot_e[j] += e
        print(f"[1e{dec},1e{dec+1}) " + " ".join(f"{c:>16}" for c in cells))
    print(f"{'TOTALS':>12} " + " ".join(
        f"{o:>7.0f}/{e:>7.1f}" for o, e in zip(tot_o, tot_e)))
    print(f"{'z':>12} " + " ".join(
        f"{(o - e) / math.sqrt(e):>+15.2f}" for o, e in zip(tot_o, tot_e)))
    print("\ntheta=1.0 column shown as observed - F1_family = adjusted/expected")
    print("* = theta=0.8 cells with z > 2: the F2 (second-order, k ~ c x^(3/4)) layer, "
          "detected but not yet subtracted (no closed form)")
    print("theta=0.5 is the Hall column: log density, deficit direction is the "
          "conjecture-flavored signal; totals row gives the current verdict")


if __name__ == "__main__":
    main()
