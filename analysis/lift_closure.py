#!/usr/bin/env python3
"""Closure report for the {2,3,m} census (paper Sec. 2).

Every improper solution of this family is a sixth-power lift of a primitive
seed (paper Lemma 1: p | gcd(x,y,a) forces p^2|x, p^3|y, p^6|a^m, because
v_p(a^m) >= m >= 7). The lift set of each seed is therefore exactly
enumerable inside the coverage region, so the improper population is a
closure certificate: a predicted count with no free parameters. A shortfall
is a coverage defect, not a statistical fluctuation.

Reports:
  1. census totals by orientation (determined arithmetically, not from labels)
  2. the lift lemma check -- every improper solution must descend
  3. the primitive-seed decomposition and largest towers
  4. the Catalan tower closure test: 2^3 + 1 = 3^2 lifts to
     (2d^2)^3 + d^6 = (3d^3)^2, admissible iff d^6 <= s_max and d^6 has
     exponent >= m_min, which forces d to be a perfect power

Usage (from repo root):
    python3 analysis/lift_closure.py data/fc23m/hits.jsonl
    python3 analysis/lift_closure.py --bands 0,1e4,1e8,1e12,1e16 data/fc23m/*.jsonl
"""
import argparse
import glob
import json
import sys
from collections import defaultdict
from math import gcd

ORIENTATIONS = ("x^3 + a^m = y^2", "x^3 - a^m = y^2", "x^3 + y^2 = a^m")


def vp(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def prime_factors(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out.append(n)
    return out


def orientation(x, y, s):
    if x ** 3 + s == y * y:
        return ORIENTATIONS[0]
    if x ** 3 - s == y * y:
        return ORIENTATIONS[1]
    if x ** 3 + y * y == s:
        return ORIENTATIONS[2]
    return None


def lift_depth(x, y, s):
    """Largest d with d^2|x, d^3|y, d^6|s."""
    d = 1
    for p in prime_factors(gcd(x, y)):
        d *= p ** min(vp(x, p) // 2, vp(y, p) // 3, vp(s, p) // 6)
    return d


def perfect_powers_upto(lim):
    """Integers 2..lim that are b^e with e >= 2."""
    out = set()
    b = 2
    while b * b <= lim:
        v = b * b
        while v <= lim:
            out.add(v)
            v *= b
        b += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hits", nargs="+")
    ap.add_argument("--bands", default="0,1e4,1e8,1e12,1e16")
    ap.add_argument("--s-max", default="1e16")
    ap.add_argument("--m-min", type=int, default=7)
    args = ap.parse_args()

    s_max = int(float(args.s_max))
    edges = [int(float(t)) for t in args.bands.split(",")]
    bands = list(zip(edges, edges[1:]))

    files = []
    for pat in args.hits:
        files.extend(sorted(glob.glob(pat)) or [pat])

    recs, seen = [], set()
    for path in files:
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            key = (r["phase"], int(r["a"]), int(r["m"]), int(r["x"]), int(r["y"]))
            if key not in seen:
                seen.add(key)
                recs.append(r)

    # --- 1. orientation totals ---
    tot = defaultdict(lambda: [0, 0, 0])
    unmatched = 0
    for r in recs:
        a, m, x, y = int(r["a"]), int(r["m"]), int(r["x"]), int(r["y"])
        o = orientation(x, y, a ** m)
        if o is None:
            unmatched += 1
            print(f"UNMATCHED IDENTITY: {r}", file=sys.stderr)
            continue
        tot[o][0] += 1
        tot[o][1 if r["proper"] else 2] += 1
    print(f"{'orientation':22s} {'all':>6s} {'proper':>7s} {'improper':>9s}")
    T = [0, 0, 0]
    for o in ORIENTATIONS:
        n, p, i = tot[o]
        T = [T[0] + n, T[1] + p, T[2] + i]
        print(f"{o:22s} {n:6d} {p:7d} {i:9d}")
    print(f"{'total':22s} {T[0]:6d} {T[1]:7d} {T[2]:9d}")

    # --- 2/3. lift lemma and seed decomposition ---
    cores = defaultdict(list)
    nonlift = []
    band_counts = [0] * len(bands)
    for r in recs:
        if r["proper"]:
            continue
        a, m, x, y = int(r["a"]), int(r["m"]), int(r["x"]), int(r["y"])
        s = a ** m
        d = lift_depth(x, y, s)
        if d == 1:
            nonlift.append(r)
            continue
        cores[(r["phase"], x // d ** 2, y // d ** 3, s // d ** 6)].append((x, y, a, m))
        for bi, (lo, hi) in enumerate(bands):
            if lo <= s < hi or (bi == len(bands) - 1 and s == hi):
                band_counts[bi] += 1
                break

    print(f"\nlift lemma: {len(nonlift)} of {T[2]} improper solutions fail to descend "
          f"(must be 0)")
    for r in nonlift[:5]:
        print("  NON-LIFT:", r.get("equation", r))
    print(f"primitive seeds: {len(cores)}")
    top = sorted(cores.items(), key=lambda kv: -len(kv[1]))[:5]
    for (ph, u, v, w), members in top:
        print(f"  seed ({u}, {v}, {w}) [{ph}]: {len(members)} lifts")
    print("improper by anchor band: " +
          "  ".join(f"[{lo:.0e},{hi:.0e}): {n}" for (lo, hi), n in zip(bands, band_counts)))

    # --- 4. Catalan tower closure test ---
    d_lim = int(round(s_max ** (1 / 6)))
    while (d_lim + 1) ** 6 <= s_max:
        d_lim += 1
    while d_lim ** 6 > s_max:
        d_lim -= 1
    predicted = {d for d in perfect_powers_upto(d_lim)}
    observed = {int(x) ** 0 for x in []}
    obs = set()
    for (ph, u, v, w), members in cores.items():
        if (u, v, w) == (2, 3, 1):
            for (x, y, a, m) in members:
                obs.add(round((x // 2) ** 0.5))
    print(f"\nCatalan tower (2^3 + 1 = 3^2), lifts (2d^2)^3 + d^6 = (3d^3)^2:")
    print(f"  d <= {d_lim} with d a perfect power: {len(predicted)} predicted")
    print(f"  observed in census: {len(obs)}")
    missing = sorted(predicted - obs)
    if missing:
        print(f"  MISSING d values: {missing[:20]}")
    extra = sorted(obs - predicted)
    if extra:
        print(f"  UNEXPECTED d values: {extra[:20]}")

    fail = bool(nonlift or unmatched or missing or extra)
    print("\n== CLOSURE " + ("FAILED ==" if fail else "VERIFIED =="))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
