#!/usr/bin/env python3
"""Merge Hall census state.json files from disjoint x-range workers.

Usage: python3 merge_states.py merged_state.json worker1/state.json worker2/state.json ...
Ranges must be disjoint (overlaps would double-count); ledgers/hits files
concatenate trivially (cat) since they are append-only per-range records.
"""
import json, sys

out, paths = sys.argv[1], sys.argv[2:]
merged = {}
for p in paths:
    for dec, row in json.load(open(p)).items():
        if dec in merged:
            merged[dec] = [a + b for a, b in zip(merged[dec], row)]
        else:
            merged[dec] = list(row)
json.dump({k: merged[k] for k in sorted(merged, key=int)}, open(out, "w"))
print(f"merged {len(paths)} states -> {out}")
for d in sorted(merged, key=int):
    print(f"  decade {d}: {merged[d]}")
