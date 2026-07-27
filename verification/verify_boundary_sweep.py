#!/usr/bin/env python3
"""Regression test for the cube-root boundary defect (paper Sec. 2.3).

The pre-repair scanner computed the negative orientation's lower bound as
int(round(s**(1/3))) with a strict inequality, so for every anchor whose cube
root has fractional part exceeding 1/2 the first admissible candidate
x = ceil(s^(1/3)) was never tested. This check does three things:

  1. exercises the failure mode directly -- for all 478 anchors it asserts
     that the scanner's own bound function returns exactly ceil(s^(1/3)),
     computed against an exact integer reference. The pre-repair bound fails
     this on precisely the 196 affected anchors.
  2. re-enumerates those 196 skipped candidates and re-tests each exactly,
     confirming none yields a square (so no solution was missed).
  3. diffs the enumeration against the committed artifact
     data/fc23m/boundary_sweep_196.json.

Usage (from repo root):
    python verification/verify_boundary_sweep.py
    python verification/verify_boundary_sweep.py --write   # regenerate artifact
"""
import argparse
import json
import os
import sys
from math import isqrt

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, "scanners")

ARTIFACT = "data/fc23m/boundary_sweep_196.json"


def icbrt(n):
    """Exact integer cube root (floor), by float seed plus integer correction."""
    x = round(n ** (1 / 3))
    while x ** 3 > n:
        x -= 1
    while (x + 1) ** 3 <= n:
        x += 1
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s-max", default="1e16")
    ap.add_argument("--m-min", type=int, default=7)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    from fc23m_scan import gen_anchors, min_cube_base   # the repaired bound

    s_max = int(float(args.s_max))
    anchors = gen_anchors(s_max, args.m_min)
    print(f"anchors: {len(anchors)}")

    # (1) the regression proper: the scanner's bound must equal ceil(s^(1/3))
    bad = []
    for s, a, m in anchors:
        if min_cube_base(s) != icbrt(s) + 1:
            bad.append((a, m, s, min_cube_base(s), icbrt(s) + 1))
    print(f"bound function correct on {len(anchors) - len(bad)}/{len(anchors)} anchors")
    for a, m, s, got, want in bad[:10]:
        print(f"  {a}^{m}: bound {got}, expected {want}")

    # (2) re-enumerate and re-test the formerly skipped candidates
    recs, squares = [], 0
    for s, a, m in anchors:
        c = icbrt(s)
        if c ** 3 == s:
            continue
        if 8 * s > (2 * c + 1) ** 3:          # frac(s^(1/3)) > 1/2, exactly
            x = c + 1
            t = x ** 3 - s
            y = isqrt(t)
            sq = y * y == t
            squares += sq
            recs.append({"a": a, "m": m, "anchor": str(s), "x": x,
                         "residual": str(t), "is_square": sq})
    print(f"boundary candidates: {len(recs)}   yielding a square: {squares}")

    doc = {
        "artifact": "boundary_sweep_196",
        "purpose": "Candidates skipped by the pre-repair lower bound in the "
                   "negative orientation x^3 - a^m = y^2 of fc23m_scan.py. The "
                   "bound was int(round(s**(1/3))) with a strict inequality, so "
                   "for every anchor whose cube root has fractional part "
                   "exceeding 1/2 the first admissible candidate "
                   "x = ceil(s^(1/3)) was excluded.",
        "selection_rule": "anchors a^m <= 1e16, m >= 7, not perfect cubes, with "
                          "8*a^m > (2*floor((a^m)^(1/3)) + 1)^3",
        "anchor_ceiling": "1e16",
        "anchors_total": len(anchors),
        "candidates": len(recs),
        "squares_found": squares,
        "regenerate": "python3 verification/verify_boundary_sweep.py --write",
        "records": recs,
    }
    blob = json.dumps(doc, indent=1, sort_keys=False) + "\n"

    if args.write:
        with open(ARTIFACT, "w", newline="\n") as f:
            f.write(blob)
        print(f"wrote {ARTIFACT}")
        return 0

    # (3) diff against the committed artifact
    with open(ARTIFACT) as f:
        committed = f.read()
    same = committed == blob
    print(f"artifact matches recomputation: {same}")

    if bad or squares or not same:
        print("\n== BOUNDARY REGRESSION FAILED ==")
        return 1
    print(f"\n== VERIFIED: bound correct on all {len(anchors)} anchors; "
          f"{len(recs)} formerly skipped candidates re-tested, none a square ==")
    return 0


if __name__ == "__main__":
    sys.exit(main())
