#!/usr/bin/env python3
r"""Record-for-record diff: independent {3,3,m} verifier vs. the census ledger.

The verifier's JSON keeps only summaries; the complete record set lives in
its checkpoint file (one JSON line per chunk, records in the last field).
This reads the checkpoint -- so no re-run is needed -- and compares tuples
(orientation, x, y_or_z, anchor) against the repository's ledger.

As in the {2,3,m} diff, THE ORIENTATION IS DETERMINED BY EXACT ARITHMETIC
rather than trusted from either side's labels:
    sum   iff  x^3 + y^3 = s        (recorded with x <= y)
    diff  iff  z^3 - x^3 = s        (recorded as (x, z))
Any ledger record satisfying neither is reported as malformed, with its
line number, instead of being silently dropped.

Ledger input: JSON Lines (fields a, m, x and one of y/z/yz -- values may be
strings) or CSV/TSV with a header. Non-solution lines are skipped and
counted.

Usage:
  diff_33m_ledger.py --ckpt verify_33m_1e+30.ckpt --ledger data/fc33m/hits.jsonl
  diff_33m_ledger.py --ckpt ... --ledger ... --max-show 25
"""
import argparse, csv, json, sys
from collections import Counter


def load_ckpt(path):
    recs = set()
    nchunks = 0
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        chunk = json.loads(line)
        nchunks += 1
        for r in chunk[-1]:
            recs.add((r[0], int(r[1]), int(r[2]), int(r[3])))
    print(f"verifier checkpoint: {nchunks} chunks, {len(recs)} distinct records")
    return recs


def classify(x, y, s):
    if x ** 3 + y ** 3 == s:
        return "sum"
    if y ** 3 - x ** 3 == s:
        return "diff"
    if x ** 3 - y ** 3 == s:          # ledger may store (z, x) order
        return "diff_swapped"
    return None


def normalize(x, y, s):
    """return the verifier's canonical tuple, or None"""
    c = classify(x, y, s)
    if c == "sum":
        return ("sum", min(x, y), max(x, y), s)
    if c == "diff":
        return ("diff", x, y, s)
    if c == "diff_swapped":
        return ("diff", y, x, s)
    return None


def parse_ledger(path):
    with open(path) as fh:
        head = fh.readline().lstrip()
    if head.startswith("{") or path.endswith((".jsonl", ".ndjson")):
        return _jsonl(path)
    return _csv(path)


def _jsonl(path):
    recs, skipped, malformed = set(), 0, []
    for ln, line in enumerate(open(path), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            malformed.append((ln, "unparseable JSON"))
            continue
        if not isinstance(row, dict):
            skipped += 1
            continue
        ykey = next((k for k in ("y", "z", "yz") if k in row), None)
        if ykey is None or not {"x"} <= set(row):
            skipped += 1
            continue
        try:
            x, y = int(row["x"]), int(row[ykey])
            if "a" in row and "m" in row:
                s = int(row["a"]) ** int(row["m"])
            else:
                s = int(row[next(k for k in ("anchor", "s", "n") if k in row)])
        except (TypeError, ValueError, StopIteration):
            malformed.append((ln, "missing or non-integer x / y / anchor"))
            continue
        t = normalize(x, y, s)
        if t is None:
            malformed.append((ln, f"satisfies neither orientation: "
                                  f"x={x} y={y} s={s}"))
            continue
        recs.add(t)
    if skipped:
        print(f"  (ledger: {skipped} non-solution lines skipped)")
    if malformed:
        print(f"  LEDGER PROBLEMS: {len(malformed)}")
        for ln, why in malformed[:10]:
            print(f"    line {ln}: {why}")
    return recs


def _csv(path):
    recs, malformed = set(), []
    with open(path, newline="") as fh:
        sample = fh.read(4096); fh.seek(0)
        delim = "\t" if "\t" in sample.splitlines()[0] else ","
        rdr = csv.DictReader((l for l in fh if not l.startswith("#")),
                             delimiter=delim)
        cols = {c.lower().strip(): c for c in rdr.fieldnames or []}
        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None
        cx = col("x")
        cy = col("y", "z", "yz")
        cs = col("anchor", "s", "value", "n")
        ca, cm = col("a", "base"), col("m", "exponent")
        if cx is None or cy is None or (cs is None and (ca is None or cm is None)):
            raise SystemExit(f"ledger: need x, y/z, and anchor (or a and m); "
                             f"found {list(cols)}")
        for i, row in enumerate(rdr, 2):
            x, y = int(row[cx]), int(row[cy])
            s = int(row[cs]) if cs else int(row[ca]) ** int(row[cm])
            t = normalize(x, y, s)
            if t is None:
                malformed.append((i, f"satisfies neither orientation: "
                                     f"x={x} y={y} s={s}"))
                continue
            recs.add(t)
    if malformed:
        print(f"  LEDGER PROBLEMS: {len(malformed)}")
        for ln, why in malformed[:10]:
            print(f"    line {ln}: {why}")
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--max-show", type=int, default=25)
    ns = ap.parse_args()

    mine = load_ckpt(ns.ckpt)
    theirs = parse_ledger(ns.ledger)
    only_v = sorted(mine - theirs, key=lambda r: (r[3], r[0], r[1]))
    only_l = sorted(theirs - mine, key=lambda r: (r[3], r[0], r[1]))

    print(f"census ledger:      {len(theirs)} records")
    print(f"verifier-only: {len(only_v)}")
    print(f"ledger-only:   {len(only_l)}")
    for tag, rows in (("VERIFIER ONLY", only_v), ("LEDGER ONLY  ", only_l)):
        for r in rows[:ns.max_show]:
            print(f"  {tag}: {r}")
        if len(rows) > ns.max_show:
            print(f"  ... and {len(rows) - ns.max_show} more")

    if only_v or only_l:
        anchors = Counter(r[3] for r in only_v + only_l)
        orients = Counter(r[0] for r in only_v + only_l)
        print(f"disagreement spans {len(anchors)} anchors; "
              f"by orientation: {dict(orients)}")
        print("worst anchors:", anchors.most_common(10))
        sys.exit(1)
    print("IDENTICAL: every record agrees, record for record.")
    sys.exit(0)


if __name__ == "__main__":
    main()
