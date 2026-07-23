#!/usr/bin/env python3
"""
deficit_analysis.py -- observed vs predicted incidental-square counts for fc23m_scan runs.

Model: for a job (anchor s, sign), the value t = x^3 + sign*s is treated as a random
integer of its size, so P(t is a perfect square) ~ 1/(2*sqrt(t)). Expected counts are
exact integrals of that density over the covered (anchor, x) region:

  Phase B:  E = sum_{s,sign} int dx / (2*sqrt(x^3 + sign*s))     [substitution u=sqrt(t)
            makes the integrand du / (3*(u^2 - sign*s)^(2/3)), smooth everywhere]
  Phase A:  E = sum_s C * s^(-1/6) / 2,  C = int_0^1 (1-u^3)^(-1/2) du ~ 1.40218

The model counts ALL squares -- it does not know about the structured scalar-multiple
families (improper hits), which are certain rather than random events. Interpretation:
observed > expected at small x = the structured families; observed ~ expected within
Poisson at large x = calibration that the random model is honest where structure thins;
the proper residual (coprime, m >= 11) is the quantity whose emptiness is being measured.

Usage:
  python3 deficit_analysis.py --ledger fc23m_out/ledger.jsonl fc23m_gap/ledger.jsonl \
      --hits fc23m_out/hits.jsonl fc23m_gap/hits.jsonl \
      --segment 0:1e14 --segment 1e14:1e16 --x-max 1e9
"""

import argparse, json, math
import numpy as np


def anchors(s_min, s_max, m_min=7):
    best = {}
    for m in range(m_min, int(math.log2(s_max)) + 2):
        a = 2
        while (v := a ** m) <= s_max:
            if v not in best or best[v][1] < m:
                best[v] = (a, m)
            a += 1
    return sorted((v, a, m) for v, (a, m) in best.items() if v > s_min)


def C_phase_a(n=200001):
    v = np.linspace(0, 1, n)                       # u = 1 - v^2 substitution
    f = 2.0 / np.sqrt(3.0 - 3.0 * v * v + v ** 4)
    return float(np.trapezoid(f, v))


def e_b(s, sign, x1, x2, n=4001):
    """Expected squares of x^3 + sign*s for integer x in [x1, x2] (discrete model)."""
    if sign < 0:
        r = round(s ** (1 / 3))
        while r ** 3 <= s: r += 1
        while (r - 1) ** 3 > s: r -= 1
        x1 = max(x1, r)
    if x1 > x2:
        return 0.0
    count = x2 - x1 + 1
    if count <= 200_000:                      # exact discrete sum, python ints (no cancellation)
        return sum(0.5 / math.sqrt(x * x * x + sign * s) for x in range(x1, x2 + 1))
    # boundary strip exactly, smooth remainder by integral in x-space
    K = 10_000
    head = sum(0.5 / math.sqrt(x * x * x + sign * s) for x in range(x1, x1 + K))
    xs = np.geomspace(float(x1 + K), float(x2), n)
    t = xs ** 3 + sign * float(s)             # safe: |x^3| >> cancellation zone after the strip
    return head + float(np.trapezoid(0.5 / np.sqrt(t), xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", nargs="+", required=True)
    ap.add_argument("--hits", nargs="+", default=[])
    ap.add_argument("--segment", action="append", required=True,
                    help="anchor band 's_min:s_max', repeatable; bands must be disjoint")
    ap.add_argument("--x-max", type=float, required=True)
    ap.add_argument("--m-min", type=int, default=7)
    args = ap.parse_args()

    segs = []
    for seg in args.segment:
        lo, hi = (float(v) for v in seg.split(":"))
        segs.append((lo, hi, anchors(lo, hi, args.m_min)))
    all_anchors = [a for _, _, band in segs for a in band]
    x_max = int(args.x_max)
    C = C_phase_a()

    # ---- observed: per-decade counts MUST come from exact hit x-values.
    # Ledger blocks are wider than decades; pro-rating fabricates sub-block structure.
    hits = []
    for path in args.hits:
        hits += [json.loads(l) for l in open(path)]
    seen = set()
    hits = [h for h in hits if not (h["equation"] in seen or seen.add(h["equation"]))]
    decades = [(10 ** k, min(10 ** (k + 1), x_max)) for k in range(0, math.ceil(math.log10(x_max)))]
    obs = np.zeros(len(decades))
    obs_a = 0
    for h in hits:
        if h["phase"] == "A":
            obs_a += 1
        else:
            for i, (lo, hi) in enumerate(decades):
                if lo <= h["x"] < hi or (h["x"] == x_max and hi == x_max):
                    obs[i] += 1
    led_a = led_b = 0
    for path in args.ledger:
        for rec in map(json.loads, open(path)):
            if rec["phase"] == "A": led_a += rec["hits"]
            else: led_b += rec["hits"]
    print(f"cross-check -- ledger hits A/B: {led_a}/{led_b}; hits-file records A/B: "
          f"{obs_a}/{int(obs.sum())} (mismatch => hits/ledger files from different run sets)")

    # ---- expected ----
    exp_b = np.zeros(len(decades))
    for i, (lo, hi) in enumerate(decades):
        exp_b[i] = sum(e_b(s, sg, lo, hi) for s, _, _ in all_anchors for sg in (+1, -1))
    exp_a = sum(C * s ** (-1 / 6) / 2 for s, _, _ in all_anchors)

    print(f"anchors: {len(all_anchors)} across {len(segs)} band(s); x <= {x_max:.3g}; C_A = {C:.5f}")
    print(f"\n{'x band':>22} {'observed':>9} {'expected':>9} {'obs/exp':>8} {'Poisson z':>10}")
    for i, (lo, hi) in enumerate(decades):
        z = (obs[i] - exp_b[i]) / math.sqrt(max(exp_b[i], 1e-12))
        print(f"[{lo:>9.0e}, {hi:>9.0e}] {obs[i]:>9.0f} {exp_b[i]:>9.2f} "
              f"{obs[i] / max(exp_b[i], 1e-12):>8.2f} {z:>+10.1f}")
    z_a = (obs_a - exp_a) / math.sqrt(exp_a)
    print(f"{'Phase A (all)':>22} {obs_a:>9.0f} {exp_a:>9.2f} {obs_a / exp_a:>8.2f} {z_a:>+10.1f}")
    tot_o, tot_e = obs.sum() + obs_a, exp_b.sum() + exp_a
    print(f"{'TOTAL':>22} {tot_o:>9.0f} {tot_e:>9.2f} {tot_o / tot_e:>8.2f}")

    if hits:
        proper = [h for h in hits if h["proper"]]
        print(f"\nhits composition ({len(hits)} unique): proper={len(proper)} "
              f"known={sum(h['known'] for h in hits)} "
              f"proper-not-known={sum(1 for h in proper if not h['known'])}")
        from collections import Counter
        cnt = Counter(h["m"] for h in hits if not h["proper"])
        print("improper by m:", dict(sorted(cnt.items())))

    # ---- the deficit statement: residual expectation beyond coverage ----
    tail = sum(2.0 / math.sqrt(x_max) for _ in all_anchors)   # int_{x_max}^inf x^{-3/2} dx per sign
    print(f"\nresidual expected squares beyond x = {x_max:.0e} (all anchors, both signs): {tail:.4f}")
    print("of which coprime with m >= 11 (the live target) is a small fraction of that.")


if __name__ == "__main__":
    main()
