#!/usr/bin/env python3
"""Artifact audit for the committed {2,3,m} census (paper Sec. 2).

This audits the *artifact*, not the instrument. A scanner self-test certifies
that the program rediscovers known solutions when it is run; it says nothing
about whether the committed hits file contains them. The coverage gap of
Sec. 2.3 passed every instrument-level check while 41% of the census was
absent from the ledger, so the completeness assertions below are not
redundant with `fc23m_scan.py --selftest`.

Known-solution membership and every identity are recomputed from the record's
own numeric fields; the `known` and `proper` flags written by the scanner are
checked against that recomputation rather than trusted.

Usage (from repo root):
    python verification/audit_hits.py data/fc23m/hits.jsonl
    python verification/audit_hits.py --any run1/hits.jsonl   # partial run

Without --any the file must be the full committed census: 855 records,
orientation totals 412 / 285 / 158, all seven known coprime solutions present.
Exit status 0 iff every check passes.
"""
import argparse
import json
import sys
from collections import Counter
from math import gcd

# The seven known coprime {2,3,m} Fermat-Catalan solutions, as unordered
# triples of power values. The other three of the ten known solutions have
# signatures outside this family.
KNOWN_TRIPLES = {
    frozenset({2 ** 7, 17 ** 3, 71 ** 2}),
    frozenset({17 ** 7, 76271 ** 3, 21063928 ** 2}),
    frozenset({43 ** 8, 96222 ** 3, 30042907 ** 2}),
    frozenset({33 ** 8, 1549034 ** 2, 15613 ** 3}),
    frozenset({9262 ** 3, 15312283 ** 2, 113 ** 7}),
    frozenset({1414 ** 3, 2213459 ** 2, 65 ** 7}),
    frozenset({7 ** 3, 13 ** 2, 2 ** 9}),
}

CENSUS_RECORDS = 855
CENSUS_ORIENTATIONS = {"B+": 412, "B-": 285, "A": 158}

failures = []


def check(cond, msg):
    """Record a failure without aborting, so one run reports every problem."""
    if not cond:
        failures.append(msg)
    return cond


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hits", nargs="*", default=["data/fc23m/hits.jsonl"],
                    help="hits.jsonl file(s); default is the committed census")
    ap.add_argument("--any", action="store_true",
                    help="audit an arbitrary run: skip the census-completeness "
                         "assertions (record count, orientation totals, all "
                         "seven known solutions) and check only internal "
                         "consistency")
    args = ap.parse_args()

    hits, seen = [], set()
    dups = 0
    for path in args.hits:
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            h = json.loads(line)
            key = (h["phase"], int(h["a"]), int(h["m"]), int(h["x"]), int(h["y"]))
            if key in seen:
                dups += 1
                continue
            seen.add(key)
            hits.append(h)
    print(f"files: {', '.join(args.hits)}")
    print(f"unique records: {len(hits)}" + (f"  ({dups} duplicates skipped)" if dups else ""))

    # --- recompute identity, coprimality and known-ness from the record ------
    present = set()
    proper = []
    for h in hits:
        a, m, x, y = int(h["a"]), int(h["m"]), int(h["x"]), int(h["y"])
        s = a ** m
        ok = {"A": x ** 3 + y * y == s,
              "B+": x ** 3 + s == y * y,
              "B-": x ** 3 - s == y * y}.get(h["phase"])
        check(ok, f"RE-VERIFICATION FAILED: {h.get('equation', h)}")

        is_proper = gcd(gcd(x, y), a) == 1
        check(is_proper == bool(h["proper"]),
              f"proper flag disagrees with gcd: {h.get('equation', h)}")

        trip = frozenset({x ** 3, y * y, s})
        is_known = trip in KNOWN_TRIPLES
        check(is_known == bool(h["known"]),
              f"known flag disagrees with the known list: {h.get('equation', h)}")
        if is_known:
            present.add(trip)
        if is_proper:
            proper.append((h, is_known))

    novel = [h for h, k in proper if not k]
    print(f"proper: {len(proper)} | known present: {len(present)}/7 | "
          f"proper-and-NOT-known: {len(novel)}")
    for h, is_known in proper:
        print("   ", h["equation"], "(known)" if is_known else "  <<< INVESTIGATE")
    print(f"exact re-verification: all {len(hits)} records recomputed")
    print("m (all):", dict(sorted(Counter(int(h["m"]) for h in hits).items())))

    orient = Counter(h["phase"] for h in hits)
    print("by orientation:", dict(orient))

    # --- census completeness (skipped with --any) ---------------------------
    if args.any:
        print("\n--any: census-completeness assertions skipped")
    else:
        check(len(present) == 7,
              f"ONLY {len(present)}/7 KNOWN SOLUTIONS PRESENT IN HITS FILE -- "
              f"the artifact is incomplete even if the scanner passes its "
              f"self-test; run verification/coverage_manifest_fc23m.py")
        check(len(hits) == CENSUS_RECORDS,
              f"record count {len(hits)} != {CENSUS_RECORDS}")
        got = {k: orient.get(k, 0) for k in CENSUS_ORIENTATIONS}
        check(got == CENSUS_ORIENTATIONS,
              f"orientation totals changed: {got} != {CENSUS_ORIENTATIONS}")

    check(not novel, f"{len(novel)} coprime solution(s) not on the known list "
                     f"-- verify independently before celebrating")

    if failures:
        print(f"\n== AUDIT FAILED: {len(failures)} problem(s) ==")
        for f in failures[:20]:
            print("  " + f)
        if len(failures) > 20:
            print(f"  ... and {len(failures) - 20} more")
        return 1
    print("\n== AUDIT PASSED ==")
    return 0


if __name__ == "__main__":
    sys.exit(main())
