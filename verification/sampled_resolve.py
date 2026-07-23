#!/usr/bin/env python3
"""Tier-2 verification for the {3,3,m} census: re-solve sampled anchors
independently and confirm every solution appears in the committed hits file.

Usage: python3 sampled_resolve.py --hits ../data/beal33m/hits_1e30.jsonl \
           --s-lo 1e24 --s-hi 1e30 --n 60 [--seed 30]
"""
import argparse, json, random, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scanners"))
from beal33m_scan import solve_anchor

def is_perfect_power(a):
    for m in range(2, a.bit_length() + 1):
        r = round(a ** (1.0 / m))
        for c in (r - 1, r, r + 1):
            if c >= 2 and c ** m == a:
                return True
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--s-lo", type=float, required=True)
    ap.add_argument("--s-hi", type=float, required=True)
    ap.add_argument("--n", type=int, default=60, help="samples per exponent")
    ap.add_argument("--seed", type=int, default=30)
    a = ap.parse_args()
    his = {json.loads(l)["equation"] for l in open(a.hits)}
    random.seed(a.seed)
    sample = []
    m = 4
    while 2 ** m <= a.s_hi:
        lo = max(2, int(a.s_lo ** (1 / m)) + 1)
        hi = int(a.s_hi ** (1 / m))
        if hi > lo:
            for _ in range(a.n):
                b = random.randint(lo, hi)
                if not is_perfect_power(b) and a.s_lo < b ** m <= a.s_hi:
                    sample.append((b ** m, b, m))
        m += 1
    # smooth (family-rich) bases at m=4
    for base in (720720, 1081080, 2162160, 4324320, 8648640, 12252240, 24504480):
        for mult in (1, 2, 3, 5, 7):
            b = base * mult
            if a.s_lo < b ** 4 <= a.s_hi and not is_perfect_power(b):
                sample.append((b ** 4, b, 4))
    sample = list({s: (s, b, m) for s, b, m in sample}.values())
    found = []
    solve_emit = lambda ph, aa, mm, x, w, eq, pr: found.append(eq)
    model = [0.0]
    for s, b, m in sample:
        solve_anchor(s, b, m, solve_emit, model)
    missing = [eq for eq in found if eq not in his]
    print(f"anchors sampled: {len(sample)} | solutions re-derived: {len(found)}"
          f" | missing from hits file: {len(missing)}")
    print("sampled completeness:", "PASS" if not missing else f"FAIL {missing[:5]}")
    sys.exit(0 if not missing else 1)

if __name__ == "__main__":
    main()
