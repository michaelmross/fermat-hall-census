#!/usr/bin/env python3
r"""Standalone verifier for the {3,3,m} exact census.

Independent reimplementation from Lemma 2 (the divisor identity) alone;
shares no code with the census scanner. From the ceiling it re-derives the
canonical anchor set, factors each base (smallest-prime-factor sieve),
enumerates admissible divisors (d^3 <= 4s), and solves the divisor
quadratics in exact integer arithmetic:

  sum orientation   x^3 + y^3 = s :  e = x+y | s,  3x^2 - 3ex + e^2 = s/e,
      x = (3e -+ sqrt(12 s/e - 3 e^2))/6, recorded once under 1 <= x <= y
      (discriminant zero gives the x = y case);
  diff orientation  z^3 - x^3 = s :  d = z-x | s,  3x^2 + 3dx + d^2 = s/d,
      x = (sqrt(12 s/d - 3 d^2) - 3d)/6,  x >= 1.

Reports totals by orientation, coprime count (any coprime solution is
printed loudly), exponent spectrum by m mod 3, deep-tail (m >= 49) records,
the exact maximum-height record, and per-block SHA-256 hashes over the
canonical serialization of records in sorted anchor order, for localizing
any disagreement with the census ledger to an anchor block.

Ledger diff: pass --ledger FILE where FILE has one record per line,
"orient,x,yz,s" with orient in {sum,diff} (adapt load_ledger() to the
repository's actual format); the verifier reports records on either side
only.

Usage:
  verify_33m_census.py CEIL [--procs N] [--blocks N] [--ledger FILE]
  e.g.  verify_33m_census.py 1e30 --procs 8
"""
__version__ = "4"
__features__ = "base-range chunking + checkpoint/resume + worker heartbeats"

import sys, math, hashlib, json, argparse
from math import isqrt
from multiprocessing import Pool


def canonical_anchors(ceil, m_min=4):
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


def spf_sieve(n):
    """smallest prime factor for 0..n (array('i'): ~4 bytes/entry, so the
    worst chunk at ceiling 1e30, a_max ~ 3.2e7, costs ~126 MB per worker)"""
    from array import array
    spf = array('i', range(n + 1))
    i = 2
    while i * i <= n:
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf


def factor_spf(n, spf):
    f = {}
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        f[p] = e
    return f


def divisors_upto(fac_pow, bound):
    divs = [1]
    for p, e in fac_pow.items():
        new = []
        for d in divs:
            v = d
            for _ in range(e + 1):
                if v > bound:
                    break
                new.append(v)
                v *= p
        divs = new
    return sorted(divs)


def solve_anchor(s, a, m, spf=None):
    fa = factor_spf(a, spf) if spf else _factor_trial(a)
    fac_pow = {p: e * m for p, e in fa.items()}
    dbound = int(round((4 * s) ** (1 / 3))) + 2
    while dbound ** 3 > 4 * s:
        dbound -= 1
    while (dbound + 1) ** 3 <= 4 * s:
        dbound += 1
    recs = []
    for d in divisors_upto(fac_pow, dbound):
        q = s // d
        disc = 12 * q - 3 * d * d
        if disc < 0:
            continue
        r = isqrt(disc)
        if r * r != disc:
            continue
        for num in (3 * d - r, 3 * d + r):        # sum orientation
            if num > 0 and num % 6 == 0:
                x = num // 6
                y = d - x
                if 1 <= x <= y:
                    recs.append(("sum", x, y, s))
        num = r - 3 * d                            # diff orientation
        if num > 0 and num % 6 == 0:
            recs.append(("diff", num // 6, num // 6 + d, s))
    return sorted(set(recs))


def _factor_trial(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


_SPF = {"arr": None, "n": 0}

def _is_perfect_power(a):
    if a < 4:
        return False
    b = a.bit_length()
    for k in range(2, b + 1):
        r = round(a ** (1.0 / k))
        for rr in (r - 1, r, r + 1):
            if rr >= 2 and rr ** k == a:
                return True
    return False

def _work(args):
    a_lo, a_hi, ceil = args
    # anchors of this chunk: non-perfect-power bases in [a_lo, a_hi),
    # all exponents m >= 4 with a^m <= ceil (canonical by construction)
    amax_global = int(ceil ** 0.25) + 2
    while amax_global ** 4 > ceil:
        amax_global -= 1
    if _SPF["arr"] is None or _SPF["n"] < amax_global:
        print(f"[worker] building SPF sieve to {amax_global} "
              f"(one-time, ~1-2 min)...", file=sys.stderr, flush=True)
        _SPF["arr"] = spf_sieve(amax_global)
        _SPF["n"] = amax_global
        print("[worker] sieve ready", file=sys.stderr, flush=True)
    print(f"[worker] starting bases {a_lo}..{a_hi-1}",
          file=sys.stderr, flush=True)
    spf = _SPF["arr"]
    keys = []
    meta = {}
    for a in range(a_lo, min(a_hi, amax_global + 1)):
        if _is_perfect_power(a):
            continue
        v, m = a ** 4, 4
        while v <= ceil:
            keys.append(v)
            meta[v] = (a, m)
            v *= a
            m += 1
    keys.sort()
    if not keys:
        return a_lo, a_lo, {"sum": 0, "diff": 0}, {0: 0, 1: 0, 2: 0}, [], \
               (0.0, None), [], hashlib.sha256().hexdigest(), 0, []
    out, h = [], hashlib.sha256()
    tot = {"sum": 0, "diff": 0}
    spec = {0: 0, 1: 0, 2: 0}
    coprime, best, deep = [], (0, None), []
    for s in keys:
        a, m = meta[s]
        for rec in solve_anchor(s, a, m, spf):
            orient, x, yz, _ = rec
            tot[orient] += 1
            spec[m % 3] += 1
            hh = max(x ** 3, yz ** 3, s)
            if hh > best[0]:
                best = (hh, (orient, x, yz, f"{a}^{m}"))
            if m >= 49:
                deep.append((orient, x, yz, f"{a}^{m}"))
            if math.gcd(x, yz) == 1 and math.gcd(x, a) == 1:
                coprime.append((rec, f"{a}^{m}"))
            h.update(repr((orient, x, yz, s, m)).encode())
            out.append(rec)
    bl = (math.log2(best[0]) if best[0] else 0.0, best[1])
    return a_lo, a_hi, tot, spec, coprime, bl, deep, h.hexdigest(), \
        len(keys), out


def load_ledger(path):
    recs = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            o, x, yz, s = line.split(",")
            recs.add((o.strip(), int(x), int(yz), int(s)))
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ceil")
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--blocks", type=int, default=64)
    ap.add_argument("--ledger")
    ns = ap.parse_args()
    ceil = int(float(ns.ceil))
    print(f"verify_33m_census v{__version__} ({__features__})", flush=True)

    print("partitioning base ranges...", file=sys.stderr, flush=True)
    amax = int(ceil ** 0.25) + 2
    while amax ** 4 > ceil:
        amax -= 1
    nb = max(ns.blocks, ns.procs * 8)
    size = (amax - 1) // nb + 1
    chunks = [(lo, min(lo + size, amax + 1), ceil)
              for lo in range(2, amax + 1, size)]
    print(f"bases 2..{amax}, {len(chunks)} chunks", file=sys.stderr, flush=True)

    # checkpoint/resume: each chunk's result is appended to CKPT as one
    # JSON line keyed by its first anchor; an interrupted run (Ctrl-C, kill,
    # reboot) restarts with the same command and skips completed chunks.
    ckpt = f"verify_33m_{ceil:.0e}.ckpt"
    done = {}
    try:
        with open(ckpt) as fh:
            for line in fh:
                rec = json.loads(line)
                done[rec[0]] = rec
        if done:
            print(f"resuming: {len(done)}/{len(chunks)} chunks in {ckpt}",
                  file=sys.stderr)
    except FileNotFoundError:
        pass
    todo = [c for c in chunks if c[0] not in done]
    def _save(res):
        lo, hi, t, sp, cp, b, dp, hx, na, out = res
        with open(ckpt, "a") as fh:
            fh.write(json.dumps([lo, hi, t, sp, cp,
                                 [b[0], b[1]], dp, hx, na, out]) + "\n")
    if ns.procs > 1:
        with Pool(ns.procs) as pool:
            for i, res in enumerate(pool.imap_unordered(_work, todo)):
                _save(res)
                print(f"  chunk {len(done)+i+1}/{len(chunks)} done",
                      file=sys.stderr, flush=True)
    else:
        for i, c in enumerate(todo):
            _save(_work(c))
            print(f"  chunk {len(done)+i+1}/{len(chunks)} done",
                  file=sys.stderr, flush=True)
    results = []
    with open(ckpt) as fh:
        for line in fh:
            lo, hi, t, sp, cp, b, dp, hx, na, out = json.loads(line)
            t = {k: t[k] for k in ("sum", "diff")}
            sp = {int(k): v for k, v in sp.items()}
            cp = [(tuple(r), lbl) for r, lbl in cp]
            b = (b[0], tuple(b[1]) if b[1] else None)
            dp = [tuple(r) for r in dp]
            out = [tuple(r) for r in out]
            results.append((lo, hi, t, sp, cp, b, dp, hx, na, out))
    results.sort(key=lambda r: r[0])

    tot = {"sum": 0, "diff": 0}
    spec = {0: 0, 1: 0, 2: 0}
    best = (0.0, None)
    deep, coprime, hashes, allrecs = [], [], [], []
    nanchors = 0
    for lo, hi, t, sp, cp, b, dp, hx, na, out in results:
        nanchors += na
        for k in tot:
            tot[k] += t[k]
        for k in spec:
            spec[k] += sp[k]
        coprime += cp
        deep += dp
        if b[0] > best[0]:
            best = b
        hashes.append({"base_lo": lo, "base_hi": hi, "sha256": hx})
        allrecs += out

    n = tot["sum"] + tot["diff"]
    print(f"ceiling {ceil:.0e}: anchors {nanchors}, solutions {n} "
          f"(sum {tot['sum']}, diff {tot['diff']}), coprime {len(coprime)}")
    for cp in coprime:
        print("COPRIME SOLUTION:", cp)
    if best[0]:
        print(f"max height 2^{best[0]:.4f}  record {best[1]}")
    print(f"spectrum m mod 3: 1:{spec[1]} 2:{spec[2]} 0:{spec[0]}")
    print(f"deep tail (m>=49): {len(deep)} records")
    with open(f"verify_33m_{ceil:.0e}.json", "w") as fh:
        json.dump({"ceiling": ceil, "anchors": nanchors, "solutions": n,
                   "by_orientation": tot, "coprime": len(coprime),
                   "max_height_log2": best[0] or None,
                   "max_height_record": best[1], "spectrum_mod3": spec,
                   "deep_tail": deep, "block_hashes": hashes}, fh, indent=2)

    print(f"complete; checkpoint {ckpt} can be deleted")

    if ns.ledger:
        mine = set(allrecs)
        theirs = load_ledger(ns.ledger)
        only_v = sorted(mine - theirs)
        only_l = sorted(theirs - mine)
        print(f"ledger diff: verifier-only {len(only_v)}, "
              f"ledger-only {len(only_l)}")
        for r in only_v[:20]:
            print("  VERIFIER ONLY:", r)
        for r in only_l[:20]:
            print("  LEDGER ONLY:  ", r)
        sys.exit(1 if (only_v or only_l) else 0)


if __name__ == "__main__":
    main()
