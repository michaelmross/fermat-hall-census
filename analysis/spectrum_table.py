#!/usr/bin/env python3
"""Regenerate the {3,3,m} exponent-spectrum and provenance numbers (paper Sec. 3).

Usage: python3 spectrum_table.py path/to/hits.jsonl
"""
import json, sys
from collections import Counter

def main(path):
    hits = [json.loads(l) for l in open(path)]
    md = Counter(h['m'] for h in hits)
    assert sum(h['proper'] for h in hits) == 0, "coprime hit present -- investigate"
    assert not any(m % 3 == 0 for m in md), "3|m hit present -- contradicts Euler"
    bins = [(4, 6), (7, 9), (10, 12), (13, 24), (25, 48), (49, 97)]
    print(f"{'m range':>12} {'m=1(3)':>9} {'m=2(3)':>9} {'m=0(3)':>9}")
    tot = {0: 0, 1: 0, 2: 0}
    for lo, hi in bins:
        row = {0: 0, 1: 0, 2: 0}
        for m, v in md.items():
            if lo <= m <= hi:
                row[m % 3] += v; tot[m % 3] += v
        print(f"{lo:>5} -- {hi:<4} {row[1]:>9} {row[2]:>9} {row[0]:>9}")
    print(f"{'total':>12} {tot[1]:>9} {tot[2]:>9} {tot[0]:>9}")
    big = max(hits, key=lambda h: int(h['y']) ** 3)
    print(f"\nsolutions: {len(hits)}  |  max height: 2^{(int(big['y'])**3).bit_length()}"
          f"  |  largest: {big['equation']}")

if __name__ == "__main__":
    main(sys.argv[1])
