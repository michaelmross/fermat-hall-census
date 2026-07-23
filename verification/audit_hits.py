#!/usr/bin/env python3
"""Post-run audit for fc23m_scan hits files (paper Sec. 2).

Usage: python3 audit_hits.py hits.jsonl [more_hits.jsonl ...]
"""
import json, sys
from collections import Counter

paths = sys.argv[1:] or ["fc23m_out/hits.jsonl"]
hits, seen = [], set()
for p in paths:
    for l in open(p):
        h = json.loads(l)
        if h["equation"] not in seen:
            seen.add(h["equation"]); hits.append(h)
print("unique records:", len(hits))
proper = [h for h in hits if h["proper"]]
novel = [h for h in proper if not h["known"]]
print("proper:", len(proper), "| known:", sum(h["known"] for h in hits),
      "| proper-and-NOT-known:", len(novel))
for h in proper:
    print("  ", h["equation"], "(known)" if h["known"] else "  <<< INVESTIGATE")
for h in hits:
    a, m, x, y, s = h["a"], h["m"], h["x"], int(h["y"]), h["a"] ** h["m"]
    ok = {"A": x**3 + y*y == s, "B+": x**3 + s == y*y, "B-": x**3 - s == y*y}[h["phase"]]
    assert ok, f"RE-VERIFICATION FAILED: {h}"
print("exact re-verification: all", len(hits), "records check")
print("m (all):", dict(sorted(Counter(h["m"] for h in hits).items())))
sys.exit(0 if not novel else 1)
