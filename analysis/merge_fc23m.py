#!/usr/bin/env python3
"""Merge {2,3,m} run directories into the canonical census hits file.

The census is the union of two runs -- the main scan and the gap-fill that
closed the anchor-band hole in it (paper Sec. 2.3) -- so no single scanner
invocation reproduces it. This script performs the union, deduplicating on the
full record key and re-verifying every identity from scratch, and refuses to
write an output that does not contain all seven known coprime solutions.

Usage (from repo root):
    python analysis/merge_fc23m.py data/fc23m/run1/hits.jsonl \\
        data/fc23m/gapfill/hits.jsonl --out data/fc23m/hits.jsonl

    python analysis/merge_fc23m.py ... --dry-run     # report, write nothing
"""
import argparse
import json
import sys
from collections import Counter

KNOWN = {
    frozenset({2 ** 7, 17 ** 3, 71 ** 2}),
    frozenset({17 ** 7, 76271 ** 3, 21063928 ** 2}),
    frozenset({43 ** 8, 96222 ** 3, 30042907 ** 2}),
    frozenset({33 ** 8, 1549034 ** 2, 15613 ** 3}),
    frozenset({9262 ** 3, 15312283 ** 2, 113 ** 7}),
    frozenset({1414 ** 3, 2213459 ** 2, 65 ** 7}),
    frozenset({7 ** 3, 13 ** 2, 2 ** 9}),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hits", nargs="+", help="per-run hits.jsonl files")
    ap.add_argument("--out", default="data/fc23m/hits.jsonl")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    seen, recs, dups = set(), [], 0
    for path in args.hits:
        n = 0
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            key = (r["phase"], int(r["a"]), int(r["m"]), int(r["x"]), int(r["y"]))
            if key in seen:
                dups += 1
                continue
            seen.add(key)
            recs.append((line, r))
            n += 1
        print(f"{path}: {n} records")

    # re-verify every identity from the record's own fields
    bad = []
    present = set()
    for _, r in recs:
        a, m, x, y = int(r["a"]), int(r["m"]), int(r["x"]), int(r["y"])
        s = a ** m
        ok = {"A": x ** 3 + y * y == s,
              "B+": x ** 3 + s == y * y,
              "B-": x ** 3 - s == y * y}.get(r["phase"])
        if not ok:
            bad.append(r)
        trip = frozenset({x ** 3, y * y, s})
        if trip in KNOWN:
            present.add(trip)

    phases = Counter(r["phase"] for _, r in recs)
    proper = sum(1 for _, r in recs if r["proper"])
    print(f"\nmerged: {len(recs)}   duplicates skipped: {dups}")
    print(f"by phase: {dict(phases)}")
    print(f"proper: {proper}   known present: {len(present)}/7")

    if bad:
        print(f"\nIDENTITY FAILURES: {len(bad)}")
        for r in bad[:5]:
            print("  ", r.get("equation", r))
        return 1
    if len(present) != 7:
        missing = KNOWN - present
        print(f"\nREFUSING TO WRITE: only {len(present)}/7 known coprime solutions "
              f"present ({len(missing)} missing). A merge that loses a known "
              f"solution indicates a coverage gap in one of the inputs; run "
              f"verification/coverage_manifest_fc23m.py on the run ledgers.")
        return 1

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return 0
    with open(args.out, "w") as f:
        f.write("\n".join(line for line, _ in recs) + "\n")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
