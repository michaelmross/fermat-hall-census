#!/usr/bin/env python3
r"""Record-for-record diff: independent {2,3,m} scan vs. the census ledger.

Compares the canonical record set emitted by scan_23m_independent.py
against the repository's solution ledger, at the level of individual
(orientation, x, y, anchor) tuples rather than totals. Totals agreeing is
necessary but not sufficient -- two implementations can agree on 855 while
disagreeing on which 855.

Inputs
  --scan  scan23m_1e+09.json  (or the .csv emitted alongside it)
  --ledger PATH               the repository's hits/solutions file

Ledger format: JSON Lines (one object per solution, fields a, m, x, y --
values may be strings) or CSV/TSV with a header. In both cases the anchor
is recomputed as a**m (or read from an anchor column) and the ORIENTATION
IS DETERMINED BY EXACT ARITHMETIC rather than trusted from a label: a
record is o1 if x^3 + s = y^2, o2 if x^3 - s = y^2, o3 if x^3 + y^2 = s.
Any record satisfying none of the three is reported as malformed. Records
on o1/o2 with x beyond the scan's xcap are outside the comparison's scope
and are excluded, with a count reported.

Exit status 0 iff the two sets are identical.

Usage:
  diff_23m_ledger.py --scan scan23m_1e+09.json --ledger data/hits_23m.csv
"""
import argparse, csv, json, sys
from collections import Counter


def load_scan(path):
    if path.endswith(".json"):
        d = json.load(open(path))
        return {(o, int(x), int(y), int(s)) for o, x, y, s in d["records"]}
    recs = set()
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        o, x, y, s = line.split(",")
        recs.add((o.strip(), int(x), int(y), int(s)))
    return recs


ORIENT_CSV = {"o1": "o1", "1": "o1", "o2": "o2", "2": "o2", "o3": "o3",
              "3": "o3", "plus": "o1", "minus": "o2"}


def classify(x, y, s):
    """Orientation by exact test -- independent of how the ledger labels or
    spells it. Returns None if the record satisfies no orientation."""
    if x ** 3 + s == y * y:
        return "o1"
    if x ** 3 - s == y * y:
        return "o2"
    if x ** 3 + y * y == s:
        return "o3"
    return None


def parse_ledger_jsonl(path, xcap=None):
    """JSON Lines: one object per solution, fields a, m, x, y (any of which
    may be strings). The anchor is recomputed as a**m and the orientation is
    determined by exact arithmetic, so no reliance on 'equation', 'phase',
    or field ordering."""
    recs, skipped, malformed, outscope = set(), 0, [], 0
    for ln, line in enumerate(open(path), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            malformed.append((ln, "unparseable JSON"))
            continue
        if not isinstance(row, dict) or not {"a", "m", "x", "y"} <= set(row):
            skipped += 1          # metadata / provenance / summary lines
            continue
        try:
            a, m = int(row["a"]), int(row["m"])
            x, y = int(row["x"]), int(row["y"])
        except (TypeError, ValueError):
            malformed.append((ln, "non-integer a/m/x/y"))
            continue
        s = a ** m
        o = classify(x, y, s)
        if o is None:
            malformed.append((ln, f"satisfies no orientation: "
                                  f"x={x} y={y} s={a}^{m}"))
            continue
        if xcap and o in ("o1", "o2") and x > xcap:
            outscope += 1         # beyond the scan's sieve ceiling
            continue
        recs.add((o, x, y, s))
    if skipped:
        print(f"  (ledger: {skipped} non-solution lines skipped)")
    if outscope:
        print(f"  (ledger: {outscope} records with x > xcap on o1/o2, "
              f"outside the scan's scope, excluded)")
    if malformed:
        print(f"  LEDGER PROBLEMS: {len(malformed)}")
        for ln, why in malformed[:10]:
            print(f"    line {ln}: {why}")
    return recs


def parse_ledger_csv(path, xcap=None):
    recs = set()
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if "\t" in sample.splitlines()[0] else ","
        rdr = csv.DictReader((l for l in fh if not l.startswith("#")),
                             delimiter=delim)
        cols = {c.lower().strip(): c for c in rdr.fieldnames or []}

        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None
        cx, cy = col("x"), col("y")
        cs = col("anchor", "s", "value", "n")
        ca, cm = col("a", "base"), col("m", "exponent")
        co = col("orientation", "orient", "o")
        if cx is None or cy is None or (cs is None and (ca is None or cm is None)):
            raise SystemExit(f"ledger: cannot find x/y and anchor columns; "
                             f"found {list(cols)} -- edit parse_ledger_csv()")
        for row in rdr:
            x, y = int(row[cx]), int(row[cy])
            s = int(row[cs]) if cs else int(row[ca]) ** int(row[cm])
            o = classify(x, y, s)
            if o is None and co:
                o = ORIENT_CSV.get(str(row[co]).strip().lower())
            if o is None:
                raise SystemExit(f"ledger: record satisfies no orientation: "
                                 f"x={x} y={y} s={s}")
            if xcap and o in ("o1", "o2") and x > xcap:
                continue
            recs.add((o, x, y, s))
    return recs


def parse_ledger(path, xcap=None):
    with open(path) as fh:
        head = fh.readline().lstrip()
    if head.startswith("{") or path.endswith((".jsonl", ".ndjson")):
        return parse_ledger_jsonl(path, xcap)
    return parse_ledger_csv(path, xcap)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--max-show", type=int, default=25)
    ns = ap.parse_args()

    mine = load_scan(ns.scan)
    xcap = None
    if ns.scan.endswith(".json"):
        xcap = json.load(open(ns.scan)).get("xcap")
    theirs = parse_ledger(ns.ledger, xcap)
    only_s = sorted(mine - theirs, key=lambda r: (r[3], r[0], r[1]))
    only_l = sorted(theirs - mine, key=lambda r: (r[3], r[0], r[1]))

    print(f"independent scan: {len(mine)} records")
    print(f"census ledger:    {len(theirs)} records")
    print(f"scan-only:   {len(only_s)}")
    print(f"ledger-only: {len(only_l)}")

    for tag, rows in (("SCAN ONLY  ", only_s), ("LEDGER ONLY", only_l)):
        for r in rows[:ns.max_show]:
            print(f"  {tag}: {r}")
        if len(rows) > ns.max_show:
            print(f"  ... and {len(rows) - ns.max_show} more")

    if only_s or only_l:
        # localize: which anchors and orientations carry the disagreement
        anchors = Counter(r[3] for r in only_s + only_l)
        orients = Counter(r[0] for r in only_s + only_l)
        print(f"disagreement spans {len(anchors)} anchors; "
              f"by orientation: {dict(orients)}")
        print("worst anchors:", anchors.most_common(10))
        sys.exit(1)

    print("IDENTICAL: every record agrees, record for record.")
    sys.exit(0)


if __name__ == "__main__":
    main()
