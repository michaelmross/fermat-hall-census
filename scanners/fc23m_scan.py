#!/usr/bin/env python3
"""
fc23m_scan.py -- Fermat-Catalan anomaly scanner for the signature family {2,3,m}, m >= 7.

Strategy: the high power a^m ("anchor") is enumerated exactly -- there are only a few
hundred values a^m <= 1e14 with m >= 7 -- and the square/cube pair is recovered by a
two-stage congruence sieve plus exact integer square test.

Orientations covered (all solutions of the {2,3,m} generalized Fermat equation are one
of these, classified by where the anchor sits in the size ordering):

  Phase A  (anchor is the LARGEST term):   x^3 + y^2 = a^m      -> x <= (a^m)^(1/3), tiny range, complete
  Phase B+ (anchor smaller than both):     x^3 + a^m = y^2      -> sieve x up to X_MAX
  Phase B- (anchor smaller than both):     x^3 - a^m = y^2      -> sieve x up to X_MAX

Auditable coverage claim per completed run (recorded in the ledger):
  every coprime solution x^3 (+/-) a^m = y^2 or x^3 + y^2 = a^m with
  a^m <= S_MAX, m >= M_MIN, and (for Phase B) x <= X_MAX.

Status of the family: (2,3,7) resolved by Poonen-Schaefer-Stoll; (2,3,8),(2,3,9) by
Bruin; (2,3,10) by Brown; (2,3,11) under GRH by Freitas-Naskrecki-Stoll. Hits with
m <= 10 must reproduce known solutions (built-in validation). A verified pairwise-
coprime hit with m >= 11 not on the known list would be an eleventh Fermat-Catalan
solution. Do not expect one; expect an auditable empty ledger.

Usage:
  python3 fc23m_scan.py --selftest
  python3 fc23m_scan.py --s-max 1e14 --x-max 1e8 --out ./fc23m_out
  (re-run the same command to resume from checkpoint; --reset to start over)

Dependencies: numpy. gmpy2 used automatically if present (~2-3x on exact tests).
"""

import argparse, json, math, os, sys, time
import numpy as np

try:
    from gmpy2 import mpz, isqrt as _isqrt, gcd as _gcd
    def isqrt(n): return int(_isqrt(mpz(n)))
    def gcd(a, b): return int(_gcd(mpz(a), mpz(b)))
    BACKEND = "gmpy2"
except ImportError:
    isqrt = math.isqrt
    gcd = math.gcd
    BACKEND = "stdlib"

M1 = 5040                                  # stage-1 combined modulus (2^4 * 3^2 * 5 * 7)
P2 = [11, 13, 17, 19, 23, 29, 31, 37, 41]  # stage-2 prime moduli

# The known coprime {2,3,m} Fermat-Catalan solutions, keyed by their three power values.
KNOWN = {
    frozenset({2**7, 17**3, 71**2}),
    frozenset({17**7, 76271**3, 21063928**2}),
    frozenset({43**8, 96222**3, 30042907**2}),
    frozenset({33**8, 1549034**2, 15613**3}),
    frozenset({9262**3, 15312283**2, 113**7}),
    frozenset({1414**3, 2213459**2, 65**7}),
    frozenset({7**3, 13**2, 2**9}),
}
N_KNOWN = len(KNOWN)  # 7 in-family; the other three of the ten have signatures outside {2,3,m}


def square_table(mod):
    t = np.zeros(mod, dtype=bool)
    y = np.arange(mod, dtype=np.int64)
    t[(y * y) % mod] = True
    return t


SQ_M1 = square_table(M1)
SQ_P2 = {p: square_table(p) for p in P2}


def gen_anchors(s_max, m_min=7, s_min=0):
    """All values s_min < a^m <= s_max with a >= 2, m >= m_min, deduped to maximal exponent.

    NOTE: dedup must see ALL representations, so enumerate from a=2 regardless of s_min
    and band-filter at the end (e.g. 2^14 in a band must not resurface as 4^7)."""
    best = {}
    m_cap = int(math.log2(s_max)) + 1
    for m in range(m_min, m_cap + 1):
        a = 2
        while (v := a ** m) <= s_max:
            if v not in best or best[v][1] < m:
                best[v] = (a, m)
            a += 1
    return sorted((v, a, m) for v, (a, m) in best.items() if v > s_min)


def job_tables(s, sign):
    """Allowed x-residue tables for t = x^3 + sign*s being a square, stage 1 and stage 2."""
    x0 = np.arange(M1, dtype=np.int64)
    t1 = SQ_M1[(x0 * x0 * x0 + sign * (s % M1)) % M1]
    t2 = {}
    for p in P2:
        xp = np.arange(p, dtype=np.int64)
        t2[p] = SQ_P2[p][(xp * xp * xp + sign * (s % p)) % p]
    return t1, t2


def emit_hit(out, phase, a, m, x, y, counters):
    s = a ** m
    terms = {"A": (x**3, y*y, s), "B+": (x**3, s, y*y), "B-": (s, y*y, x**3)}[phase]
    assert terms[0] + terms[1] == terms[2], "hit failed re-verification"
    g = gcd(gcd(x, y), a)
    rec = {
        "phase": phase, "a": a, "m": m, "x": x, "y": str(y),
        "equation": {"A": f"{x}^3 + {y}^2 = {a}^{m}",
                     "B+": f"{x}^3 + {a}^{m} = {y}^2",
                     "B-": f"{x}^3 - {a}^{m} = {y}^2"}[phase],
        "proper": g == 1,
        "known": frozenset(terms) in KNOWN,
        "signature": [2, 3, m],
    }
    counters["hits"] += 1
    if rec["proper"] or rec["known"]:
        flag = "KNOWN" if rec["known"] else "NEW"
        print(f"  [{flag}] {rec['equation']}")
    if not rec["known"] and rec["proper"]:
        print("  *** COPRIME HIT NOT ON KNOWN LIST -- verify independently, then breathe. ***")
    with open(os.path.join(out, "hits.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")


def phase_a(anchors, out, counters):
    """Complete search of x^3 + y^2 = a^m for every enumerated anchor (x <= s^(1/3))."""
    t0 = time.time()
    for s, a, m in anchors:
        cb = int(round(s ** (1 / 3))) + 2
        x = np.arange(1, cb + 1, dtype=np.int64)
        if s < 2**62:                       # int64-safe fast path
            t = s - x * x * x
            t = t[t > 0]
            x = x[: len(t)]
            cand = x[SQ_M1[t % M1]]
        else:                               # rare: fall back to exact loop
            cand = x
        for xv in map(int, cand):
            t = s - xv ** 3
            if t <= 0:
                continue
            counters["exact"] += 1
            y = isqrt(t)
            if y * y == t and y > 0:
                emit_hit(out, "A", a, m, xv, y, counters)
    return time.time() - t0


def phase_b_block(jobs, x_lo, x_hi, out, counters):
    """Sieve one block of x for every (anchor, sign) job."""
    x = np.arange(x_lo, x_hi, dtype=np.int64)
    xm1 = x % M1
    xmp = {p: x % p for p in P2}
    for s, a, m, sign, t1, t2 in jobs:
        mask = t1[xm1]
        if sign < 0:                        # need x^3 > s
            xmin = isqrt(isqrt(s)) if False else int(round(s ** (1 / 3)))
            mask &= x > xmin
        idx = np.flatnonzero(mask)
        counters["s1"] += len(idx)
        for p in P2:
            if len(idx) == 0:
                break
            idx = idx[t2[p][xmp[p][idx]]]
        for xv in map(int, x[idx]):
            t = xv ** 3 + sign * s
            if t <= 0:
                continue
            counters["exact"] += 1
            y = isqrt(t)
            if y * y == t and y > 0:
                emit_hit(out, "B+" if sign > 0 else "B-", a, m, xv, y, counters)


def ledger(out, **kw):
    kw["ts"] = round(time.time(), 1)
    with open(os.path.join(out, "ledger.jsonl"), "a") as f:
        f.write(json.dumps(kw) + "\n")


def run(s_max, x_max, block, m_min, out, reset=False, quiet=False, s_min=0):
    os.makedirs(out, exist_ok=True)
    ck_path = os.path.join(out, "checkpoint.json")
    config = {"s_max": s_max, "s_min": s_min, "m_min": m_min}
    ck = {"phase_a_done": False, "next_x": 1, "config": config}
    if os.path.exists(ck_path) and not reset:
        ck = json.load(open(ck_path))
        if ck.get("config") != config:
            sys.exit(f"[fc23m] REFUSING TO RESUME: checkpoint in {out} was written with config "
                     f"{ck.get('config')} but this run requests {config}.\n"
                     f"  Resuming across configs silently loses coverage (Phase A and low-x Phase B "
                     f"are skipped for anchors the old run never saw).\n"
                     f"  Either use --reset, a fresh --out, or run the missing band as a gap-fill, e.g.:\n"
                     f"    --s-min <old s_max> --s-max <new s_max> --x-max <old x_max> --out <new dir>")
    anchors = gen_anchors(s_max, m_min, s_min)
    print(f"[fc23m] backend={BACKEND}  anchors={len(anchors)} ({s_min:.3g} < a^m <= {s_max:.3g}, m >= {m_min})"
          f"  x_max={x_max:.3g}  resume_x={ck['next_x']}")
    counters = {"hits": 0, "exact": 0, "s1": 0}

    if not ck["phase_a_done"]:
        secs = phase_a(anchors, out, counters)
        ck["phase_a_done"] = True
        json.dump(ck, open(ck_path, "w"))
        ledger(out, phase="A", s_min=s_min, s_max=s_max, anchors=len(anchors), exact=counters["exact"],
               hits=counters["hits"], secs=round(secs, 2))
        print(f"[fc23m] phase A complete ({secs:.1f}s): anchor-largest orientation exhausted.")

    jobs = []
    for s, a, m in anchors:
        for sign in (+1, -1):
            t1, t2 = job_tables(s, sign)
            jobs.append((s, a, m, sign, t1, t2))

    x = ck["next_x"]
    while x <= x_max:
        hi = min(x + block, x_max + 1)
        t0 = time.time()
        c0 = dict(counters)
        phase_b_block(jobs, x, hi, out, counters)
        secs = time.time() - t0
        ledger(out, phase="B", x_lo=x, x_hi=hi - 1, jobs=len(jobs),
               s1_survivors=counters["s1"] - c0["s1"],
               exact_tests=counters["exact"] - c0["exact"],
               hits=counters["hits"] - c0["hits"], secs=round(secs, 2))
        if not quiet:
            rate = (hi - x) / max(secs, 1e-9)
            print(f"[fc23m] B: x in [{x}, {hi-1}]  {secs:.1f}s  ({rate:.2e} x/s across {len(jobs)} jobs)")
        x = hi
        ck["next_x"] = x
        json.dump(ck, open(ck_path, "w"))

    print(f"[fc23m] DONE. Coverage: all {{2,3,m}} solutions with {s_min:.3g} < a^m <= {s_max:.3g}, m >= {m_min},"
          f" and (phase B) x <= {x_max:.3g}. hits={counters['hits']}"
          f"  (ledger + hits in {out}/)")
    return counters


def selftest():
    import tempfile
    out = tempfile.mkdtemp(prefix="fc23m_selftest_")
    print("=== selftest: S_MAX=2.5e14, X_MAX=1e5 -- must rediscover all 7 known {2,3,m} solutions ===")
    run(s_max=2.5e14, x_max=100_000, block=100_000, m_min=7, out=out, reset=True)
    hits = [json.loads(l) for l in open(os.path.join(out, "hits.jsonl"))]
    found_known = sum(1 for h in hits if h["known"])
    novel_proper = [h for h in hits if not h["known"] and h["proper"]]
    print(f"=== selftest: {found_known}/{N_KNOWN} known solutions rediscovered; "
          f"{len(novel_proper)} unexpected coprime hits ===")
    ok = found_known == N_KNOWN and not novel_proper
    print("=== selftest:", "PASS" if ok else "FAIL", "===")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--s-max", type=float, default=1e14)
    ap.add_argument("--s-min", type=float, default=0, help="anchor band lower bound (gap-fill runs)")
    ap.add_argument("--x-max", type=float, default=1e7)
    ap.add_argument("--block", type=int, default=1_000_000)
    ap.add_argument("--m-min", type=int, default=7)
    ap.add_argument("--out", default="./fc23m_out")
    ap.add_argument("--reset", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    run(a.s_max, int(a.x_max), a.block, a.m_min, a.out, a.reset, a.quiet, a.s_min)
