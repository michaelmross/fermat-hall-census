#!/usr/bin/env python3
r"""Independent scanner for the {2,3,m} sieve orientations.

Second implementation sharing no code or sieve design with the census
scanner: quadratic-residue lookup tables (one gather per modulus) instead
of a CRT wheel, with moduli disjoint from the original's (43200 = 2^6.3^3.5^2
then 49, 121, 169, 289, 361, then primes 23..61; the original used 5040 then
primes to 41). Survivors are confirmed by exact integer isqrt.

Orientations:
  o1: x^3 + s = y^2, x = 1..XCAP
  o2: x^3 - s = y^2, x = icbrt(s)+1..XCAP
  o3: x^3 + y^2 = s, x = 1..icbrt(s-1)   (forced complete; exact, no sieve)

Output: JSON with all records, per-orientation totals, coprime/proper
records, the Catalan-tower count (seeds 2^3+1=3^2 lifts with d a perfect
power, d^6 <= ceiling), and a weighted-descent audit (every improper record
must descend per Lemma lift).

Usage: scan_23m_independent.py [--xcap 1e9] [--procs N] [--out file.json]
"""
import argparse, json, math, sys
from math import isqrt, gcd
import numpy as np
from multiprocessing import Pool

M1 = 43200  # 2^6 * 3^3 * 5^2
M2 = [49, 121, 169, 289, 361, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
CHUNK = 5_000_000


def icbrt(n):
    x = int(round(n ** (1 / 3)))
    while x ** 3 > n:
        x -= 1
    while (x + 1) ** 3 <= n:
        x += 1
    return x


def canonical_anchors(ceil=10**16, m_min=7):
    seen = {}
    a = 2
    while a ** m_min <= ceil:
        v, m = a ** m_min, m_min
        while v <= ceil:
            if v not in seen or m > seen[v][1]:
                seen[v] = (a, m)
            v *= a
            m += 1
        a += 1
    return seen


def sq_residues(m):
    return np.unique(np.array([(y * y) % m for y in range(m)], dtype=np.int64))


_SQSETS = {m: None for m in [M1] + M2}


def sqmask(m):
    if _SQSETS[m] is None:
        mask = np.zeros(m, dtype=bool)
        mask[sq_residues(m)] = True
        _SQSETS[m] = mask
    return _SQSETS[m]


def anchor_table(s, sign, m):
    """table[r] = is (r^3 + sign*s) a square residue mod m"""
    r = np.arange(m, dtype=np.int64)
    t = (r * r % m * r + np.int64(sign) * np.int64(s % m)) % m
    return sqmask(m)[t]


def scan_orientation(s, sign, x0, x1):
    """exact solutions of x^3 + sign*s = y^2 for x in [x0, x1]"""
    if x1 < x0:
        return []
    tabs = [(m, anchor_table(s, sign, m)) for m in [M1] + M2]
    hits = []
    x = x0
    while x <= x1:
        xe = min(x1, x + CHUNK - 1)
        xs = np.arange(x, xe + 1, dtype=np.int64)
        surv = xs
        for m, tab in tabs:
            if surv.size == 0:
                break
            surv = surv[tab[surv % m]]
        for xv in surv.tolist():
            t = xv * xv * xv + sign * s
            r = isqrt(t)
            if r * r == t:
                hits.append((int(xv), int(r)))
        x = xe + 1
    return hits


def scan_o3_exact(s):
    """x^3 + y^2 = s: exact, forced complete"""
    hits = []
    for xv in range(1, icbrt(s - 1) + 1):
        t = s - xv ** 3
        r = isqrt(t)
        if r * r == t and r >= 1:
            hits.append((xv, r))
    return hits


def do_anchor(args):
    s, a, m, xcap = args
    recs = []
    for xv, yv in scan_orientation(s, +1, 1, xcap):
        recs.append(("o1", xv, yv, s))
    for xv, yv in scan_orientation(s, -1, icbrt(s) + 1, xcap):
        recs.append(("o2", xv, yv, s))
    for xv, yv in scan_o3_exact(s):
        recs.append(("o3", xv, yv, s))
    return s, a, m, recs


def weighted_descent_ok(x, y, s, a):
    """Lemma lift: improper solution must admit d>1 with d^2|x, d^3|y, d^6|s"""
    g = gcd(gcd(x, y), a)
    if g == 1:
        return True  # proper; nothing to check
    # find maximal d per the lemma
    d = 1
    for p in set(_factor(g)):
        al = _vp(x, p); be = _vp(y, p); ga = _vp(s, p)
        e = min(al // 2, be // 3, ga // 6)
        d *= p ** e
    return d > 1


def _factor(n):
    f, d = [], 2
    while d * d <= n:
        while n % d == 0:
            f.append(d); n //= d
        d += 1
    if n > 1:
        f.append(n)
    return f


def _vp(n, p):
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xcap", default="1e9")
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--out", default=None)
    ap.add_argument("--anchors", default=None,
                    help="comma-separated anchor values (default: all 478)")
    ns = ap.parse_args()
    xcap = int(float(ns.xcap))
    meta = canonical_anchors()
    keys = sorted(meta) if not ns.anchors else [int(v) for v in ns.anchors.split(",")]
    jobs = [(s, meta[s][0], meta[s][1], xcap) for s in keys]

    # checkpoint/resume: one JSON line per completed anchor; rerunning the
    # same command skips finished anchors, so interruption costs at most
    # one anchor per worker.
    ckpt = (ns.out or f"scan23m_{xcap:.0e}.json") + ".ckpt"
    done = {}
    try:
        with open(ckpt) as fh:
            for line in fh:
                rec = json.loads(line)
                done[rec[0]] = rec
        if done:
            print(f"resuming: {len(done)}/{len(jobs)} anchors in {ckpt}",
                  file=sys.stderr)
    except FileNotFoundError:
        pass
    todo = [j for j in jobs if j[0] not in done]
    def _save(res):
        s_, a_, m_, recs_ = res
        with open(ckpt, "a") as fh:
            fh.write(json.dumps([s_, a_, m_, recs_]) + "\n")
    if ns.procs > 1:
        with Pool(ns.procs) as pool:
            for i, res in enumerate(pool.imap_unordered(do_anchor, todo)):
                _save(res)
                if (i + 1) % 10 == 0:
                    print(f"  ...{len(done)+i+1}/{len(jobs)} anchors",
                          file=sys.stderr, flush=True)
    else:
        for i, j in enumerate(todo):
            _save(do_anchor(j))
            if (i + 1) % 10 == 0:
                print(f"  ...{len(done)+i+1}/{len(jobs)} anchors",
                      file=sys.stderr, flush=True)
    results = []
    with open(ckpt) as fh:
        for line in fh:
            s_, a_, m_, recs_ = json.loads(line)
            results.append((s_, a_, m_, [tuple(r) for r in recs_]))
    results.sort(key=lambda r: r[0])

    tot = {"o1": 0, "o2": 0, "o3": 0}
    proper = []
    allrecs = []
    descent_fail = []
    catalan = 0
    for s, a, m, recs in results:
        for o, xv, yv, _ in recs:
            tot[o] += 1
            allrecs.append((o, xv, yv, s))
            g = gcd(gcd(xv, yv), a)
            if g == 1:
                proper.append((o, xv, yv, f"{a}^{m}"))
            elif not weighted_descent_ok(xv, yv, s, a):
                descent_fail.append((o, xv, yv, s))
            # Catalan lift: o1, x = 2d^2, y = 3d^3, s = d^6
            if o == "o1":
                d6 = s
                d = round(d6 ** (1 / 6))
                for dd in (d - 1, d, d + 1):
                    if dd >= 2 and dd ** 6 == s and xv == 2 * dd * dd and yv == 3 * dd ** 3:
                        catalan += 1
    n = sum(tot.values())
    print(f"xcap {xcap:.0e}: anchors {len(keys)}, solutions {n} "
          f"(o1 {tot['o1']}, o2 {tot['o2']}, o3 {tot['o3']})")
    print(f"proper (coprime) solutions: {len(proper)}")
    for p in proper:
        print("  PROPER:", p)
    print(f"Catalan-tower lifts found: {catalan}")
    print(f"weighted-descent failures: {len(descent_fail)}", descent_fail[:5])
    out = ns.out or f"scan23m_{xcap:.0e}.json"
    with open(out, "w") as fh:
        json.dump({"xcap": xcap, "totals": tot, "n": n,
                   "proper": proper, "catalan_lifts": catalan,
                   "descent_failures": descent_fail,
                   "records": allrecs}, fh)
    print("written", out)


if __name__ == "__main__":
    main()
