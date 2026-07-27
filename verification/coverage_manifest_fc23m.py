#!/usr/bin/env python3
"""Verify that the {2,3,m} census covers its full claimed rectangle.

Coverage here is two-dimensional -- (anchor set) x (cube-base range) -- and the
failure mode this audit exists to catch is a run whose anchor set changed
partway through, leaving anchors scanned for only part of the x-range. Tiling
the x-range alone does not detect it.

For each run ledger the anchor subset of a Phase-B block is reconstructed from
the recorded job count (jobs = 2 x |anchors|, one job per sign) together with
the run's declared lower anchor bound. Phase-A coverage is read from the
Phase-A records. The union over all runs must equal the full rectangle.

Usage (from repo root):
    python3 verification/coverage_manifest_fc23m.py \
        --run data/fc23m/ledger.jsonl:0 \
        --run data/fc23m/gapfill/ledger.jsonl:1e14 \
        --s-max 1e16 --x-max 1e9

Each --run is FILE:S_MIN. Exit status 0 iff the rectangle is fully covered.
"""
import argparse
import json
import math
import sys


def gen_anchors(s_max, m_min=7):
    """Distinct a^m <= s_max, m >= m_min, deduplicated to maximal exponent."""
    best = {}
    m_cap = int(math.log2(s_max)) + 1
    for m in range(m_min, m_cap + 1):
        a = 2
        while (v := a ** m) <= s_max:
            if v not in best or best[v] < m:
                best[v] = m
            a += 1
    return sorted(best)


def merge(intervals):
    """Merge a list of inclusive integer intervals."""
    if not intervals:
        return []
    out = []
    for lo, hi in sorted(intervals):
        if out and lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True,
                    help="LEDGER:S_MIN, repeatable")
    ap.add_argument("--s-max", default="1e16")
    ap.add_argument("--x-max", default="1e9")
    ap.add_argument("--m-min", type=int, default=7)
    args = ap.parse_args()

    s_max = int(float(args.s_max))
    x_max = int(float(args.x_max))
    anchors = gen_anchors(s_max, args.m_min)
    idx = {v: i for i, v in enumerate(anchors)}
    print(f"anchor set: {len(anchors)} values a^m <= {s_max:.3g}, m >= {args.m_min}")

    phase_a = set()                       # anchor indices with Phase A done
    xcov = {i: [] for i in range(len(anchors))}
    problems = []                         # coverage failures (fatal)
    warnings = []                         # provenance anomalies (non-fatal)

    for spec in args.run:
        path, _, smin_s = spec.rpartition(":")
        s_min = int(float(smin_s))
        band = [v for v in anchors if v > s_min]
        nblocks = 0
        jobsets = set()
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("phase") == "A":
                n = int(r["anchors"])
                if n > len(band):
                    problems.append(f"{path}: Phase A claims {n} anchors, "
                                    f"band above {s_min:.3g} holds {len(band)}")
                for v in band[:n]:
                    phase_a.add(idx[v])
            elif r.get("phase") == "B":
                n = int(r["jobs"]) // 2
                jobsets.add(n)
                lo, hi = int(r["x_lo"]), int(r["x_hi"])
                for v in band[:n]:
                    xcov[idx[v]].append((lo, hi))
                nblocks += 1
        print(f"{path}: s_min={s_min:.3g}  band={len(band)} anchors  "
              f"B-blocks={nblocks}  distinct anchor counts in B={sorted(jobsets)}")
        if len(jobsets) > 1:
            warnings.append(f"{path}: anchor set CHANGED mid-run {sorted(jobsets)} "
                            f"-- blocks before the change cover fewer anchors")

    # --- Phase A completeness ---
    missing_a = [anchors[i] for i in range(len(anchors)) if i not in phase_a]
    print(f"\nPhase A: {len(phase_a)}/{len(anchors)} anchors covered")
    if missing_a:
        problems.append(f"Phase A missing {len(missing_a)} anchors, "
                        f"smallest {missing_a[0]:.4g}, largest {missing_a[-1]:.4g}")

    # --- Phase B completeness ---
    bad_b = []
    for i, v in enumerate(anchors):
        m = merge(xcov[i])
        if m != [[1, x_max]]:
            bad_b.append((v, m))
    print(f"Phase B: {len(anchors) - len(bad_b)}/{len(anchors)} anchors covered "
          f"for all x in [1, {x_max:.3g}]")
    for v, m in bad_b[:10]:
        got = ", ".join(f"[{a},{b}]" for a, b in m) or "nothing"
        problems.append(f"anchor {v:.6g}: x covered {got}, expected [1,{x_max:.3g}]")
    if len(bad_b) > 10:
        problems.append(f"... and {len(bad_b) - 10} further anchors with x gaps")

    if warnings:
        print("\n== provenance warnings (not coverage failures) ==")
        for w in warnings:
            print("  " + w)
        print("  -> the rectangle may still be complete via a later gap-fill run;")
        print("     coverage verdict below is what governs.")
    if problems:
        print("\n== COVERAGE INCOMPLETE ==")
        for p in problems:
            print("  " + p)
        return 1
    print("\n== COVERAGE COMPLETE: every anchor x every x in [1, %g] ==" % x_max)
    return 0


if __name__ == "__main__":
    sys.exit(main())