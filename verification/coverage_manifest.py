#!/usr/bin/env python3
"""Derive and verify the Hall census coverage certificate from worker checkpoints.

The authoritative range of a worker is the ``config`` block inside its
``checkpoint.json``, not its directory name. This script reads every worker
checkpoint under a root directory, verifies that the ranges

  * are individually complete   (next_x == x_max + 1),
  * are pairwise disjoint,
  * tile their span contiguously with no gaps,

and writes a COVERAGE.md table derived from the checkpoints themselves. If a
directory name encodes a range (e.g. ``x_4.0e10-5.5e10``), the encoded range is
checked against the checkpoint and a mismatch is a failure: names become
checked claims rather than assertions.

Exit status is 0 only if every check passes, so this is safe to run in CI.

Usage:
    python3 coverage_manifest.py data/hall/dec10                 # verify + print
    python3 coverage_manifest.py data/hall/dec10 --write         # also write COVERAGE.md
    python3 coverage_manifest.py data/hall/dec10 --expect 1e10 1e11
"""
import argparse
import json
import re
import sys
from pathlib import Path

NAME_RANGE = re.compile(r'([0-9.]+e[0-9]+)\s*[-_]\s*([0-9.]+e[0-9]+)', re.I)


def find_workers(root):
    """Every directory under root containing a checkpoint.json, sorted by x_min."""
    workers = []
    for ck_path in sorted(Path(root).rglob('checkpoint.json')):
        ck = json.loads(ck_path.read_text())
        cfg = ck.get('config')
        if not cfg or 'x_min' not in cfg or 'x_max' not in cfg:
            print(f"  ! {ck_path}: no config block (pre-v1.2 checkpoint?) -- skipped")
            continue
        workers.append({
            'dir': ck_path.parent,
            'name': ck_path.parent.name,
            'x_min': int(cfg['x_min']),
            'x_max': int(cfg['x_max']),
            'next_x': int(ck['next_x']),
        })
    workers.sort(key=lambda w: w['x_min'])
    return workers


def check(workers, expect=None):
    """Run every coverage check. Returns (ok, list_of_failure_strings)."""
    fail = []
    if not workers:
        return False, ["no worker checkpoints found"]

    for w in workers:
        # completeness: the checkpoint must sit exactly one past the ceiling
        if w['next_x'] != w['x_max'] + 1:
            done = w['next_x'] - w['x_min']
            span = w['x_max'] - w['x_min'] + 1
            fail.append(f"{w['name']}: incomplete, next_x={w['next_x']} "
                        f"({done / span:.1%} of [{w['x_min']}, {w['x_max']}])")
        # directory name, if it encodes a range, must agree with the checkpoint
        m = NAME_RANGE.search(w['name'])
        if m:
            lo, hi = float(m.group(1)), float(m.group(2))
            if abs(lo - w['x_min']) > max(1.0, 1e-9 * lo) or \
               abs(hi - w['x_max']) > max(1.0, 1e-9 * hi):
                fail.append(f"{w['name']}: directory name claims [{lo:.6g}, {hi:.6g}] "
                            f"but checkpoint says [{w['x_min']}, {w['x_max']}]")

    # pairwise disjoint and contiguous
    for a, b in zip(workers, workers[1:]):
        if b['x_min'] <= a['x_max']:
            fail.append(f"{a['name']} and {b['name']} overlap: "
                        f"[..., {a['x_max']}] meets [{b['x_min']}, ...]")
        elif b['x_min'] != a['x_max'] + 1:
            fail.append(f"gap between {a['name']} and {b['name']}: "
                        f"[{a['x_max'] + 1}, {b['x_min'] - 1}] uncovered")

    if expect:
        lo, hi = int(float(expect[0])), int(float(expect[1]))
        if workers[0]['x_min'] != lo:
            fail.append(f"span starts at {workers[0]['x_min']}, expected {lo}")
        if workers[-1]['x_max'] != hi:
            fail.append(f"span ends at {workers[-1]['x_max']}, expected {hi}")

    return not fail, fail


def manifest(workers):
    lo, hi = workers[0]['x_min'], workers[-1]['x_max']
    total = hi - lo + 1
    lines = [
        "# Coverage certificate",
        "",
        "Derived from the worker checkpoints by `verification/coverage_manifest.py`;",
        "regenerate with that script rather than editing by hand. The range of a",
        "worker is authoritative in its `checkpoint.json`, not in its directory name.",
        "",
        f"**Span covered: [{lo:,}, {hi:,}]** &mdash; {total:,} values of $x$,",
        f"in {len(workers)} disjoint contiguous ranges, each checkpoint terminating",
        "exactly one past its ceiling.",
        "",
        "| worker | x_min | x_max | values | complete |",
        "|---|---:|---:|---:|:-:|",
    ]
    for w in workers:
        n = w['x_max'] - w['x_min'] + 1
        done = "yes" if w['next_x'] == w['x_max'] + 1 else f"no ({w['next_x']})"
        lines.append(f"| `{w['name']}` | {w['x_min']:,} | {w['x_max']:,} | {n:,} | {done} |")
    lines += ["",
              "Checks performed: per-worker completeness, pairwise disjointness,",
              "contiguity of the tiling, and (where a directory name encodes a range)",
              "agreement between the name and the checkpoint.",
              ""]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="directory containing the worker subdirectories")
    ap.add_argument("--write", action="store_true",
                    help="write COVERAGE.md into the root directory")
    ap.add_argument("--expect", nargs=2, metavar=("X_MIN", "X_MAX"),
                    help="also assert the total span, e.g. --expect 1e10 1e11")
    a = ap.parse_args()

    workers = find_workers(a.root)
    ok, fail = check(workers, a.expect)

    for w in workers:
        flag = "OK " if w['next_x'] == w['x_max'] + 1 else "INC"
        print(f"  [{flag}] {w['name']:<24} [{w['x_min']:>15,}, {w['x_max']:>15,}]")
    if workers:
        print(f"  span: [{workers[0]['x_min']:,}, {workers[-1]['x_max']:,}] "
              f"in {len(workers)} ranges")

    if ok:
        print("COVERAGE CERTIFICATE: PASS")
        if a.write:
            out = Path(a.root) / "COVERAGE.md"
            out.write_text(manifest(workers))
            print(f"wrote {out}")
    else:
        print("COVERAGE CERTIFICATE: FAIL")
        for f in fail:
            print(f"  - {f}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
