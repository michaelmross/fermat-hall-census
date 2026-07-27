#!/usr/bin/env python3
"""Verify the deep-exponent regularity of the {3,3,m} census (paper Sec. 3).

Claim: every solution with m >= 49 belongs to one of two one-parameter
skeleton families, and every admissible value of each is realized exactly
once:

    2x^3 = 2^(3j+1)          x = 2^j,   m = 3j+1 = 49, 52, ..., 97
    x^3 + (2x)^3 = 3^(3k+2)  x = 3^k,   m = 3k+2 = 50, 53, 56, 59, 62

The 3-family exits at m = 62 because 3^65 > 10^30: a window artifact of base
growth, not an arithmetic asymmetry. Identities are recomputed from scratch;
the record's own labels are not trusted.

Usage (from repo root, after unzipping the census):
    unzip -o data/beal33m/hits_1e30.zip -d data/beal33m/
    python verification/verify_deep_skeletons.py
"""
import argparse
import json
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default="data/beal33m/hits_1e30.jsonl")
    ap.add_argument("--m-min", type=int, default=49)
    ap.add_argument("--s-max", default="1e30")
    args = ap.parse_args()
    s_max = int(float(args.s_max))

    deep = []
    with open(args.hits) as f:
        for line in f:
            line = line.strip()
            if line and int(json.loads(line)["m"]) >= args.m_min:
                deep.append(json.loads(line))

    exp2 = set(range(args.m_min, 98, 3))
    exp2 = {m for m in exp2 if m % 3 == 1 and 2 ** m <= s_max}
    exp3 = {m for m in range(args.m_min, 98) if m % 3 == 2 and 3 ** m <= s_max}
    expected = len(exp2) + len(exp3)
    print(f"deep records (m >= {args.m_min}): {len(deep)}   expected {expected} "
          f"({len(exp2)} base-2 + {len(exp3)} base-3)")

    seen = {2: set(), 3: set()}
    errors = []
    for r in deep:
        a, m = int(r["a"]), int(r["m"])
        u, v = sorted((int(r["x"]), int(r["y"])))
        if u ** 3 + v ** 3 != a ** m:
            errors.append(f"identity fails: {r.get('equation', r)}")
            continue
        if r["proper"]:
            errors.append(f"proper flag set on deep record: {r}")
        if a == 2:
            if m % 3 != 1 or u != v or u != 2 ** ((m - 1) // 3):
                errors.append(f"not 2-skeleton: {r}")
                continue
        elif a == 3:
            if m % 3 != 2 or u != 3 ** ((m - 2) // 3) or v != 2 * u:
                errors.append(f"not 3-skeleton: {r}")
                continue
        else:
            errors.append(f"deep solution with unexpected base {a}: {r}")
            continue
        if m in seen[a]:
            errors.append(f"duplicate m={m} for base {a}")
        seen[a].add(m)

    if seen[2] != exp2:
        errors.append(f"base-2 exponents {sorted(seen[2])} != expected {sorted(exp2)}")
    if seen[3] != exp3:
        errors.append(f"base-3 exponents {sorted(seen[3])} != expected {sorted(exp3)}")

    for e in errors:
        print("  " + e)
    if errors:
        print("== DEEP-SKELETON CHECK FAILED ==")
        return 1
    print(f"== VERIFIED: {len(deep)} deep solutions = "
          f"{len(exp2)} x (2^j, 2^j, 2^(3j+1)) + {len(exp3)} x (3^k, 2*3^k, 3^(3k+2)) ==")
    return 0


if __name__ == "__main__":
    sys.exit(main())
