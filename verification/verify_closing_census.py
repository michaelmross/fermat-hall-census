#!/usr/bin/env python3
"""Independent re-census of the anchor-largest orientation x^3 + y^2 = a^m,
diffed against the committed hits file.

This shares no code with scanners/fc23m_scan.py: anchors are re-enumerated,
every admissible x is tested by exact integer square root with no congruence
sieve, and the result is compared record-for-record. It is deliberately the
dumbest possible implementation, because its purpose is to disagree with the
clever one if the clever one is wrong. It found a 65-solution coverage gap
that every internal check had passed over.

Runtime: a few seconds at the 1e16 ceiling.

Usage (from repo root):
    python verification/verify_closing_census.py
    python verification/verify_closing_census.py --hits data/fc23m/hits.jsonl --s-max 1e16
"""
import argparse
import json
import math
import sys
from math import gcd, isqrt

KNOWN_CLOSING = {
    (7, 13, 2, 9),
    (1414, 2213459, 65, 7),
    (9262, 15312283, 113, 7),
}


def gen_anchors(s_max, m_min=7):
    best = {}
    m_cap = int(math.log2(s_max)) + 1
    for m in range(m_min, m_cap + 1):
        a = 2
        while (v := a ** m) <= s_max:
            if v not in best or best[v][1] < m:
                best[v] = (a, m)
            a += 1
    return sorted((v, a, m) for v, (a, m) in best.items())


def icbrt(n):
    x = round(n ** (1 / 3))
    while x ** 3 > n:
        x -= 1
    while (x + 1) ** 3 <= n:
        x += 1
    return x


def enumerate_closing(s_max, m_min=7):
    out = []
    for s, a, m in gen_anchors(s_max, m_min):
        for x in range(1, icbrt(s - 1) + 1):
            r = s - x ** 3
            y = isqrt(r)
            if y * y == r and y > 0:
                proper = gcd(x, y) == 1 and gcd(x, a) == 1 and gcd(y, a) == 1
                out.append((x, y, a, m, proper))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default="data/fc23m/hits.jsonl")
    ap.add_argument("--s-max", default="1e16")
    ap.add_argument("--m-min", type=int, default=7)
    args = ap.parse_args()

    s_max = int(float(args.s_max))
    truth = enumerate_closing(s_max, args.m_min)
    tset = {t[:4] for t in truth}
    n_proper = sum(1 for t in truth if t[4])
    print(f"independent enumeration: {len(truth)} solutions "
          f"({n_proper} proper, {len(truth) - n_proper} improper)")

    found = set()
    for line in open(args.hits):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        a, m, x, y = int(r["a"]), int(r["m"]), int(r["x"]), int(r["y"])
        if x ** 3 + y * y == a ** m:
            found.add((x, y, a, m))
    print(f"ledger ({args.hits}): {len(found)} closing-orientation records")

    missing = sorted(tset - found)
    extra = sorted(found - tset)
    for k in missing:
        tag = " <-- KNOWN COPRIME SOLUTION" if k in KNOWN_CLOSING else ""
        print(f"  MISSING: {k[0]}^3 + {k[1]}^2 = {k[2]}^{k[3]}{tag}")
    for k in extra:
        print(f"  EXTRA:   {k[0]}^3 + {k[1]}^2 = {k[2]}^{k[3]}")

    if missing or extra:
        print(f"\n== DISAGREEMENT: {len(missing)} missing, {len(extra)} extra ==")
        return 1
    known_ok = KNOWN_CLOSING <= tset
    print(f"\n== AGREEMENT: {len(truth)} solutions, record-for-record ==")
    print(f"   all three known closing-orientation coprime solutions present: {known_ok}")
    return 0 if known_ok else 1


if __name__ == "__main__":
    sys.exit(main())
