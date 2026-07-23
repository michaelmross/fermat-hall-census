#!/usr/bin/env python3
"""
beal33m_scan.py -- EXACT solver for the Beal-proper signature family {3,3,m}, m >= 4.

Unlike the {2,3,m} scanner, nothing here is searched. Both orientations factor:

  D33:  x^3 + a^m = z^3   =>  d := z - x  divides a^m,  and  3x^2 + 3dx + d^2 = a^m / d
  S33:  x^3 + y^3 = a^m   =>  e := x + y  divides a^m,  and  3x^2 - 3ex + e^2 = a^m / e

so for each anchor s = a^m the COMPLETE solution set is: for each divisor d|s, test
disc = 12(s/d) - 3d^2 for being a positive perfect square t^2, then x = (t -+ 3d)/6
when integral. Divisor sets of prime-power-ish anchors are tiny; the whole family up
to s <= 1e16 resolves in seconds, with cube heights up to ~s^(3/2) and no x-cutoff.

Coverage claim per run: EVERY integer solution of x^3 +- y^3 = a^m with a^m <= S_MAX,
m >= M_MIN, x,y >= 1 -- unconditionally complete, since d <= (4s)^(1/3) is forced.

Theoretical status (why observed proper count should be 0):
  3 | m reduces to Fermat's cubic (Euler): no coprime solutions.
  (3,3,4), (3,3,5): resolved by Bruin, no nontrivial coprime solutions.
  (3,3,p): resolved for large prime ranges (Kraus; Chen--Siksek and successors).
A verified pairwise-coprime hit here would contradict published theorems for most m
-- treat any such output as a bug until proven otherwise, then as a very loud event.

Usage:
  python3 beal33m_scan.py --selftest
  python3 beal33m_scan.py --s-max 1e16 --out ./beal33m_out
"""

import argparse, json, math, os, sys, time
from math import isqrt, gcd


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def divisors_of_power(a, m):
    """Divisors of a^m from the factorization of a."""
    divs = [1]
    for p, e in factor(a).items():
        divs = [d * p ** k for d in divs for k in range(m * e + 1)]
    return divs


def gen_anchors(s_max, m_min=4):
    best = {}
    for m in range(m_min, int(math.log2(s_max)) + 2):
        a = 2
        while (v := a ** m) <= s_max:
            if v not in best or best[v][1] < m:
                best[v] = (a, m)
            a += 1
    return sorted((v, a, m) for v, (a, m) in best.items())


def solve_anchor(s, a, m, emit, model):
    """Complete solution set of x^3 + y^3 = s and z^3 - x^3 = s. Also accumulates the
    random-model expectation: P(disc is a square) ~ 1/(2 sqrt(disc)), times 1/6 for
    the integrality of x."""
    for d in divisors_of_power(a, m):
        disc = 12 * (s // d) - 3 * d * d
        if disc < 0:
            continue
        if disc > 0:
            model[0] += 2 / (12 * math.sqrt(disc))  # both orientations share the disc; disc=0 is deterministic
        t = isqrt(disc)
        if t * t != disc:
            continue
        # D33: x = (t - 3d)/6, z = x + d
        if (t - 3 * d) % 6 == 0:
            x = (t - 3 * d) // 6
            if x >= 1:
                z = x + d
                assert z ** 3 - x ** 3 == s
                emit("D33", a, m, x, z, f"{x}^3 + {a}^{m} = {z}^3", gcd(x, z) == 1)
        # S33: x = (3e - t)/6, y = e - x  (minus root only; x <= y)
        if (3 * d - t) % 6 == 0:
            x = (3 * d - t) // 6
            y = d - x
            if 1 <= x <= y:
                assert x ** 3 + y ** 3 == s
                emit("S33", a, m, x, y, f"{x}^3 + {y}^3 = {a}^{m}", gcd(x, y) == 1)


def run(s_max, m_min, out, quiet=False):
    os.makedirs(out, exist_ok=True)
    anchors = gen_anchors(int(s_max), m_min)
    print(f"[beal33m] anchors: {len(anchors)} (a^m <= {s_max:.3g}, m >= {m_min}); "
          f"exact divisor enumeration, no x-cutoff")
    hits, proper_hits, model = [], [], [0.0]
    t0 = time.time()

    def emit(phase, a, m, x, w, eq, proper):
        rec = {"phase": phase, "a": a, "m": m, "x": x, "y": str(w), "equation": eq,
               "proper": proper, "known": False, "signature": [3, 3, m]}
        hits.append(rec)
        if proper:
            proper_hits.append(rec)
            print(f"  [PROPER] {eq}")
            print("  *** COPRIME (3,3,m) SOLUTION: contradicts published theorems for most m."
                  " Assume bug first. ***")
        elif not quiet:
            pass  # improper hits logged to file only

    for s, a, m in anchors:
        solve_anchor(s, a, m, emit, model)
    secs = time.time() - t0

    with open(os.path.join(out, "hits.jsonl"), "w") as f:
        for h in hits:
            f.write(json.dumps(h) + "\n")
    summary = {"s_max": s_max, "m_min": m_min, "anchors": len(anchors),
               "hits": len(hits), "proper_hits": len(proper_hits),
               "model_expected": round(model[0], 4), "secs": round(secs, 2),
               "coverage": f"complete for all (3,3,m) solutions with a^m <= {s_max:.3g}, m >= {m_min}"}
    with open(os.path.join(out, "ledger.jsonl"), "a") as f:
        f.write(json.dumps(summary) + "\n")
    print(f"[beal33m] DONE in {secs:.1f}s: hits={len(hits)} (improper families), "
          f"proper={len(proper_hits)}, raw model expectation={model[0]:.3f}")
    print(f"[beal33m] coverage: {summary['coverage']}")
    return summary, hits


def selftest():
    import tempfile
    out = tempfile.mkdtemp(prefix="beal33m_selftest_")
    print("=== selftest 1: brute-force cross-check, all x,y,z <= 300 vs anchors <= 1e8 ===")
    S = 10 ** 8
    anchor_vals = {s for s, _, _ in gen_anchors(S, 4)}
    brute = set()
    for x in range(1, 301):
        for w in range(x, 301):
            if x ** 3 + w ** 3 in anchor_vals:
                brute.add(("S", x, w, x ** 3 + w ** 3))
        for z in range(x + 1, 301):
            if z ** 3 - x ** 3 in anchor_vals:
                brute.add(("D", x, z, z ** 3 - x ** 3))
    _, hits = run(S, 4, out, quiet=True)
    mine = {("S" if h["phase"] == "S33" else "D", h["x"], int(h["y"]), h["a"] ** h["m"])
            for h in hits}
    mine_in_range = {h for h in mine if max(h[1], h[2]) <= 300}
    missed, extra = brute - mine_in_range, mine_in_range - brute
    print(f"brute-force solutions: {len(brute)} | solver (restricted to range): {len(mine_in_range)}"
          f" | missed: {len(missed)} | spurious: {len(extra)}")
    ok1 = not missed and not extra
    print("=== selftest 2: canary improper family x=y=2^j (2*2^(3j) = 2^(3j+1)) present ===")
    canaries = {(2, 2), (4, 4), (8, 8)}   # anchors 2^4, 2^7, 2^10
    found = {(h["x"], int(h["y"])) for h in hits if h["phase"] == "S33"}
    ok2 = canaries <= found
    print("canaries found:", sorted(canaries & found))
    ok3 = not any(h["proper"] for h in hits)
    print("proper hits (must be 0):", sum(h['proper'] for h in hits))
    print("=== selftest:", "PASS" if (ok1 and ok2 and ok3) else "FAIL", "===")
    return 0 if (ok1 and ok2 and ok3) else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--s-max", type=float, default=1e16)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--out", default="./beal33m_out")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    run(a.s_max, a.m_min, a.out)
