#!/usr/bin/env python3
"""
hall_census.py -- census of Hall near-misses k(x) = min |y^2 - x^3| over x <= X_MAX.

For each x, the nearest square to x^3 is checked (y = isqrt(x^3) and y+1), giving the
minimal residual k(x). Three outputs:

1. THETA-BAND CENSUS (the Hall statistic): per decade of x, observed counts of
   |k(x)| <= x^theta for theta in {0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, against the random
   model E = sum_x min(1, x^(theta - 3/2)). Hall's conjecture asserts the theta=1/2+eps
   column stays finite; the model predicts ~log-density there. This table is the
   empirical content.

2. QUALITY RECORDS: every x with Hall ratio r = sqrt(x)/|k| above --r-floor is logged
   with full data (cross-check against published record tables offline).

3. POWER CROSS-CHECK: any k(x) that is a perfect m-th power (m >= 7) with gcd checks
   is a {2,3,m} Fermat-Catalan solution and MUST already appear in the fc23m census
   (k <= x^(3/2) is far below the anchor ceiling for any reachable x). Agreement
   validates both scanners; disagreement is a bug alarm.

Canaries (selftest): x=2 gives k=1 (Catalan: 3^2-2^3=1, r=1.414...);
x=5234 gives k=17 (378661^2 - 5234^3 = 17, r ~ 4.26, the classical Hall example).

Usage:
  python3 hall_census.py --selftest
  python3 hall_census.py --x-max 1e9 --out ./hall_out      (resumable; ~tens of min)
"""

import argparse, json, math, os, sys, time

try:
    from gmpy2 import mpz, isqrt as _gisqrt
    def isqrt(n): return int(_gisqrt(mpz(n)))
    BACKEND = "gmpy2"
except ImportError:
    isqrt = math.isqrt
    BACKEND = "stdlib"

THETAS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


_PP = {"bound": 0, "table": {}}

def power_table(bound):
    """Dict {a^m: (a, m)} for all m >= 7, a >= 2, a^m <= bound (m maximal). Tiny: O(bound^(1/7))."""
    if bound > _PP["bound"]:
        t = {}
        for m in range(7, bound.bit_length() + 1):
            a = 2
            while (v := a ** m) <= bound:
                t[v] = (a, m)       # increasing m overwrites -> maximal exponent kept
                a += 1
        _PP["bound"], _PP["table"] = bound, t
    return _PP["table"]


def scan_block(x_lo, x_hi, r_floor, records, census, power_hits):
    ptab = power_table(isqrt(x_hi ** 3) + 1)
    for x in range(x_lo, x_hi):
        c = x * x * x
        y = isqrt(c)
        k1 = c - y * y
        k2 = (y + 1) * (y + 1) - c
        if k1 <= k2:
            k, ysign = k1, y
        else:
            k, ysign = k2, y + 1
        if k == 0:          # x a perfect square -> x^3 is a square; not a near-miss
            continue
        lx = math.log(x) if x > 1 else 0.0
        lk = math.log(k)
        dec = int(math.log10(x)) if x > 1 else 0
        for j, th in enumerate(THETAS):
            if lk <= th * lx or x == 1:
                census[dec][j] += 1
        r = math.sqrt(x) / k
        if r >= r_floor:
            records.append({"x": x, "y": ysign, "k": k if ysign * ysign > c else -k,
                            "r": round(r, 4)})
        pp = ptab.get(k)
        if pp:
            a, m = pp
            power_hits.append({"x": x, "y": ysign, "a": a, "m": m,
                               "equation": (f"{x}^3 + {a}^{m} = {ysign}^2" if ysign * ysign > c
                                            else f"{x}^3 - {a}^{m} = {ysign}^2"),
                               "proper": math.gcd(x, ysign) == 1})


def expected_census(x_lo, x_hi):
    """E[count of |k| <= x^theta] per decade, random model, exact sums (fast approx)."""
    exp = {}
    for dec in range(int(math.log10(max(x_lo, 1))), int(math.log10(x_hi)) + 1):
        lo, hi = max(10 ** dec, x_lo), min(10 ** (dec + 1) - 1, x_hi - 1)
        if lo > hi:
            continue
        exp[dec] = []
        for th in THETAS:
            p = th - 1.5
            if abs(p + 1) < 1e-12:          # theta = 1/2: logarithmic density
                e = math.log(hi / lo)
            else:
                e = (hi ** (p + 1) - lo ** (p + 1)) / (p + 1)
            exp[dec].append(e)
    return exp


def run(x_max, out, r_floor=0.2, block=1_000_000):
    os.makedirs(out, exist_ok=True)
    ck_path = os.path.join(out, "checkpoint.json")
    st_path = os.path.join(out, "state.json")
    ck = {"next_x": 1, "x_max": x_max}
    census = {}
    if os.path.exists(ck_path):
        ck_old = json.load(open(ck_path))
        if ck_old.get("x_max") == x_max:
            ck = ck_old
            census = {int(k): v for k, v in json.load(open(st_path)).items()}
        else:
            print("[hall] config changed; starting fresh")
    print(f"[hall] backend={BACKEND}  x_max={x_max:.3g}  resume_x={ck['next_x']}")
    x = ck["next_x"]
    while x <= x_max:
        hi = min(x + block, x_max + 1)
        for dec in range(int(math.log10(max(x, 1))), int(math.log10(hi)) + 1):
            census.setdefault(dec, [0] * len(THETAS))
        recs, pows = [], []
        t0 = time.time()
        scan_block(x, hi, r_floor, recs, census, pows)
        secs = time.time() - t0
        with open(os.path.join(out, "hall_hits.jsonl"), "a") as f:
            for r in recs:
                f.write(json.dumps(r) + "\n")
        with open(os.path.join(out, "power_hits.jsonl"), "a") as f:
            for p in pows:
                f.write(json.dumps(p) + "\n")
        with open(os.path.join(out, "ledger.jsonl"), "a") as f:
            f.write(json.dumps({"x_lo": x, "x_hi": hi - 1, "records": len(recs),
                                "power_hits": len(pows), "secs": round(secs, 2)}) + "\n")
        x = hi
        ck["next_x"] = x
        json.dump(ck, open(ck_path, "w"))
        json.dump(census, open(st_path, "w"))
        print(f"[hall] x <= {hi-1:.3g}  ({(hi - (x - block)) / max(secs,1e-9):.2e} x/s)")

    exp = expected_census(1, x_max)
    print(f"\n{'decade':>12} " + " ".join(f"{'th='+str(t):>14}" for t in THETAS))
    for dec in sorted(census):
        row = census[dec]
        erow = exp.get(dec, [float('nan')] * len(THETAS))
        cells = " ".join(f"{o:>6}/{e:>7.1f}" for o, e in zip(row, erow))
        print(f"[1e{dec},1e{dec+1}) {cells}")
    print("(cells are observed/expected counts of |k| <= x^theta)")
    return census


def selftest():
    import tempfile
    out = tempfile.mkdtemp(prefix="hall_selftest_")
    recs, pows, census = [], {0: [0]*len(THETAS)}, []
    census = {d: [0]*len(THETAS) for d in range(0, 5)}
    r_list, p_list = [], []
    scan_block(1, 6001, 0.2, r_list, census, p_list)
    by_x = {r["x"]: r for r in r_list}
    ok1 = 2 in by_x and abs(by_x[2]["k"]) == 1
    ok2 = 5234 in by_x and abs(by_x[5234]["k"]) == 17 and abs(by_x[5234]["r"] - 4.256) < 0.01
    print("canary x=2, k=1:", "PASS" if ok1 else "FAIL")
    print("canary x=5234, k=17 (classical Hall point):", "PASS" if ok2 else "FAIL")
    ok3 = all(r["y"] ** 2 - r["x"] ** 3 == r["k"] for r in r_list)
    print("exact residual verification on", len(r_list), "records:", "PASS" if ok3 else "FAIL")
    print("=== selftest:", "PASS" if (ok1 and ok2 and ok3) else "FAIL", "===")
    return 0 if (ok1 and ok2 and ok3) else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--x-max", type=float, default=1e7)
    ap.add_argument("--out", default="./hall_out")
    ap.add_argument("--r-floor", type=float, default=0.2)
    ap.add_argument("--block", type=int, default=1_000_000)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    run(int(a.x_max), a.out, a.r_floor, a.block)
