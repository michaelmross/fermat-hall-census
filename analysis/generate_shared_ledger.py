#!/usr/bin/env python3
r"""Authoritative shared numerical ledger for the Fermat--Hall companion papers.

Computes the raw density-model expectation over the 478 canonical anchors
a^m <= 1e16 (m >= 7), per orientation and per scope, as the EXACT DISCRETE
SUM the model defines -- one Bernoulli trial per integer x with success
probability (1/2) t^{-1/2}, t = |x^3 +- a^m| -- not as a continuum integral.

Rigor: heads (first XHEAD terms past every endpoint, and all of orientation 3)
are summed exactly with integer t (no float overflow: all |t| < 2^63); tails
use the decreasing-function bracket
    int_{A+1}^{B+1} f  <=  sum_{x=A+1}^{B} f(x)  <=  int_{A}^{B} f
whose midpoint is taken and half-width accumulated as the reported bound.
Quadrature is mpmath at 30 dps on a smooth integrand away from all branch
points; its error is far below the bracket width.

Known prior defect this ledger supersedes: earlier values (33.09; 33.36 with
a finer integrator) were quadratures of the continuum integral, which
acquires ~0.38 of spurious mass at the branch points x^3 = a^m of
orientations 2 and 3 (pure integral: 33.42). Refinement therefore moved the
value away from the model value. The discrete sum is 33.0399 +- 6e-6.

Emits: fermat_ledger.json (machine), fermat_ledger.tex (\newcommand macros +
a booktabs table body for \input), and a printed summary.
"""
import json
import numpy as np
import mpmath as mp
from fractions import Fraction

mp.mp.dps = 30
CEIL, M_MIN, XSCAN, XHEAD = 10**16, 7, 10**9, 10**5


def icbrt(n):
    x = int(round(n ** (1 / 3)))
    while x ** 3 > n:
        x -= 1
    while (x + 1) ** 3 <= n:
        x += 1
    return x


def canonical_anchors(ceil=CEIL, m_min=M_MIN):
    seen = {}
    a = 2
    while a ** m_min <= ceil:
        v, m = a ** m_min, m_min
        while v <= ceil:
            if v not in seen or m > seen[v][1]:
                seen[v] = (a, m)
            v *= a
            m += 1
        a += 1
    return seen


def ell(a):
    r, n, p, fac = Fraction(1), a, 2, set()
    while p * p <= n:
        while n % p == 0:
            fac.add(p)
            n //= p
        p += 1
    if n > 1:
        fac.add(n)
    for p in fac:
        r *= Fraction(1, 2) if p == 2 else Fraction(p - 1, p)
    return r


def head_sum_exact(s, x0, x1, sign):
    if x1 < x0:
        return 0.0
    tot, x, step = 0.0, x0, 4_000_000
    while x <= x1:
        xe = min(x1, x + step - 1)
        xs = np.arange(x, xe + 1, dtype=np.int64)
        t = xs * xs * xs + np.int64(sign) * np.int64(s)
        assert t.min() >= 1
        tot += float(np.sum(0.5 / np.sqrt(t.astype(np.float64))))
        x = xe + 1
    return tot


def sum3_exact(s, xmax):
    if xmax < 1:
        return 0.0
    xs = np.arange(1, xmax + 1, dtype=np.int64)
    t = np.int64(s) - xs * xs * xs
    assert t.min() >= 1
    return float(np.sum(0.5 / np.sqrt(t.astype(np.float64))))


def gint(s, sign, A, B):
    return float(mp.quad(lambda x: 0.5 / mp.sqrt(x ** 3 + sign * s), [A, B]))


def orient_mass(s, sign, xmin, xcap):
    head_end = min(xcap, xmin + XHEAD - 1) if xcap else xmin + XHEAD - 1
    v, err = head_sum_exact(s, xmin, head_end, sign), 0.0
    if xcap is None or xcap > head_end:
        lo = gint(s, sign, head_end + 1, mp.inf if xcap is None else xcap + 1)
        hi = gint(s, sign, head_end, mp.inf if xcap is None else xcap)
        v += 0.5 * (lo + hi)
        err = 0.5 * (hi - lo)
    return v, err


def main():
    anchors = canonical_anchors()
    acc = dict(E1a=0.0, E1s=0.0, E2a=0.0, E2s=0.0, E3=0.0,
               erra=0.0, errs=0.0, La=0.0, Ls=0.0)
    for s in sorted(anchors):
        a, _m = anchors[s]
        la = float(ell(a))
        v1a, e1a = orient_mass(s, +1, 1, None)
        v1s, e1s = orient_mass(s, +1, 1, XSCAN)
        xmin = icbrt(s) + 1
        v2a, e2a = orient_mass(s, -1, xmin, None)
        v2s, e2s = orient_mass(s, -1, xmin, XSCAN)
        v3 = sum3_exact(s, icbrt(s - 1))
        acc["E1a"] += v1a; acc["E1s"] += v1s
        acc["E2a"] += v2a; acc["E2s"] += v2s
        acc["E3"] += v3
        acc["erra"] += e1a + e2a; acc["errs"] += e1s + e2s
        acc["La"] += la * (v1a + v2a + v3)
        acc["Ls"] += la * (v1s + v2s + v3)

    E_all = acc["E1a"] + acc["E2a"] + acc["E3"]
    E_scan = acc["E1s"] + acc["E2s"] + acc["E3"]
    led = {
        "anchors": len(anchors),
        "o1_all": acc["E1a"], "o1_scanned": acc["E1s"],
        "o2_all": acc["E2a"], "o2_scanned": acc["E2s"],
        "o3_forced": acc["E3"],
        "Eraw_all": E_all, "Eraw_scanned": E_scan,
        "Eraw_tail": E_all - E_scan,
        "err_all": acc["erra"], "err_scanned": acc["errs"],
        "ell_all": acc["La"], "ell_scanned": acc["Ls"],
        "avg_ell_scanned": acc["Ls"] / E_scan,
        "coverage_pct": 100.0 * E_scan / E_all,
        "conventions": "sum over integer x; y>=1; canonical anchors; "
                       "scanned = x<=1e9 on orientations 1,2; "
                       "orientation 3 forced complete",
    }
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "fermat_ledger.json"), "w") as fh:
        json.dump(led, fh, indent=2)

    def f4(x):
        return f"{x:.4f}"

    with open(os.path.join(here, "fermat_ledger.tex"), "w") as fh:
        fh.write("% generated by generate_shared_ledger.py -- do not edit\n")
        for name, val in (
            ("LedgerAnchors", str(led["anchors"])),
            ("LedgerErawAll", f4(E_all)),
            ("LedgerErawScanned", f4(E_scan)),
            ("LedgerErawTail", f4(led["Eraw_tail"])),
            ("LedgerEllAll", f4(led["ell_all"])),
            ("LedgerEllScanned", f4(led["ell_scanned"])),
            ("LedgerAvgEll", f"{led['avg_ell_scanned']:.4f}"),
            ("LedgerCoverage", f"{led['coverage_pct']:.2f}"),
            ("LedgerErrAll", f"{led['err_all']:.1e}"),
        ):
            fh.write(f"\\newcommand{{\\{name}}}{{{val}}}\n")
        fh.write("% table body: orientation & all-x & scanned (x<=1e9)\n")
        fh.write("\\newcommand{\\LedgerTableBody}{%\n")
        fh.write(f"$x^3 + a^m = y^2$ & {f4(acc['E1a'])} & {f4(acc['E1s'])} \\\\\n")
        fh.write(f"$x^3 - a^m = y^2$ & {f4(acc['E2a'])} & {f4(acc['E2s'])} \\\\\n")
        fh.write(f"$x^3 + y^2 = a^m$ & {f4(acc['E3'])} & {f4(acc['E3'])} \\\\\n")
        fh.write("\\midrule\n")
        fh.write(f"total & {f4(E_all)} & {f4(E_scan)} \\\\\n")
        fh.write(f"$\\ell(a)$-weighted & {f4(acc['La'])} & {f4(acc['Ls'])} \\\\\n")
        fh.write("}\n")

    for k, v in led.items():
        print(f"{k}: {v}")
    print(f"wrote {os.path.join(here, 'fermat_ledger.json')}")
    print(f"wrote {os.path.join(here, 'fermat_ledger.tex')}")


if __name__ == "__main__":
    main()
