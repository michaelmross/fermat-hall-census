#!/usr/bin/env python3
"""Regenerate every table and numerical value in the coprimality note.

  M. M. Ross, "Coprimality Density and the Proper-Solution Deficit in the
  Generalized Fermat Family {2,3,m}"

Stages
------
  skeleton   enumerate anchors, derive the sixth-power-free curve classes,
             write skeleton.json and skeleton_classes.txt (input to the
             PARI/GP rank census, analysis/skeleton_ranks.gp)
  local      the coprimality density l(a), the exact 2- and 3-adic densities
             with their stabilization certificate, the good-prime decay
             table, and the truncation study E(300)/E(1000)/E(3000)
  ranks      the rank decomposition of the diagnostic mass (requires
             skeleton_ranks.txt from the GP stage)
  all        every stage, in order, plus coprimality_results.json for the
             plotting script

Usage
-----
    python coprimality_pipeline.py all
    python coprimality_pipeline.py skeleton          # then run the .gp script
    python coprimality_pipeline.py ranks --ranks skeleton_ranks.txt

Dependencies: numpy, sympy. Runtime for `all`: a few minutes.
"""
import argparse
import json
import math
import sys
from collections import defaultdict

import numpy as np
from sympy import factorint, primerange

# --- census parameters (must match the companion census) ---------------------
S_MAX = 10 ** 16          # anchor ceiling a^m
M_MIN = 7                 # exponent floor for the {2,3,m} family
X_MAX = 10 ** 9           # cube-base ceiling of the census coverage
K2, K3 = 10, 7            # working levels for the 2- and 3-adic densities
CUTOFFS = (300, 1000, 3000)

# The seven known coprime solutions, as (anchor base, exponent, orientation
# sign): sign +1 for x^3 + a^m = y^2 and the closing orientation, -1 for
# x^3 - a^m = y^2. Used only for the calibration statistics.
KNOWN = [(2, 7, +1), (2, 9, +1), (17, 7, +1), (65, 7, +1),
         (113, 7, +1), (43, 8, +1), (33, 8, -1)]


# --- anchors -----------------------------------------------------------------
def anchors(s_max=S_MAX, m_min=M_MIN):
    """Distinct values a^m <= s_max, deduplicated to maximal exponent."""
    best = {}
    for m in range(m_min, int(math.log2(s_max)) + 2):
        a = 2
        while (v := a ** m) <= s_max:
            if v not in best or best[v][1] < m:
                best[v] = (a, m)
            a += 1
    return sorted((v, a, m) for v, (a, m) in best.items())


def sixth_power_free(a, m, sign):
    """The class k0 = sign * prod p^((m v_p(a)) mod 6) of y^2 = x^3 + sign a^m."""
    core = 1
    for p, e in factorint(a).items():
        core *= p ** ((m * e) % 6)
    return sign * core


# --- raw density model -------------------------------------------------------
def e_raw_branch(s, sign, x_max=X_MAX, n=3001):
    """Sum of 1/(2 sqrt(x^3 + sign s)) over the admissible x <= x_max."""
    x_lo = 1
    if sign < 0:                       # need x^3 > s: exact integer bound
        r = round(s ** (1 / 3))
        while r ** 3 <= s:
            r += 1
        while (r - 1) ** 3 > s:
            r -= 1
        x_lo = r
    if x_lo > x_max:
        return 0.0
    head_hi = min(x_lo + 10000, x_max + 1)
    total = sum(0.5 / math.sqrt(x * x * x + sign * s) for x in range(x_lo, head_hi))
    if head_hi <= x_max:
        xs = np.geomspace(float(head_hi), float(x_max), n)
        total += float(np.trapezoid(0.5 / np.sqrt(xs ** 3 + float(sign * s)), xs))
    return total


def e_raw_closing(s):
    """Third orientation x^3 + y^2 = a^m, whose cube base is forced.

    Discrete sum of 1/(2 sqrt(s - x^3)) over 1 <= x <= floor((s-1)^(1/3)),
    replacing the continuum closed form C_A s^(-1/6)/2, which collects mass
    on the branch-point interval containing no integer (the third defect of
    the census paper's Appendix B)."""
    r = round(s ** (1 / 3))
    while r ** 3 >= s:
        r -= 1
    while (r + 1) ** 3 < s:
        r += 1
    xs = np.arange(1, r + 1, dtype=np.int64)
    diff = np.int64(s) - xs ** 3
    return float(np.sum(0.5 / np.sqrt(diff.astype(np.float64))))


# --- local factors -----------------------------------------------------------
def ell(a):
    """Coprimality density: prod over p | a of (1 - 1/p), with 1/2 at p = 2."""
    f = 1.0
    for p in factorint(a):
        f *= 0.5 if p == 2 else (1 - 1 / p)
    return f


def _sqrt_counts(pk):
    """tab[r] = #{y mod pk : y^2 = r}."""
    tab = [0] * pk
    for y in range(pk):
        tab[(y * y) % pk] += 1
    return tab


_TAB = {}


def delta(s, sign, p, K):
    """Normalized local density of y^2 = x^3 + sign s modulo p^K."""
    pk = p ** K
    if (p, K) not in _TAB:
        _TAB[(p, K)] = _sqrt_counts(pk)
    tab = _TAB[(p, K)]
    sp = (sign * s) % pk
    return sum(tab[(x * x * x + sp) % pk] for x in range(pk)) / pk


def point_count_tables(z):
    """NP[p][r] = affine point count of y^2 = x^3 + r over F_p, all r at once.

    Computed by FFT cross-correlation of the square-root count with the cube
    count; verified against direct character sums on a random sample.
    """
    NP = {}
    for p in primerange(5, z + 1):
        w = np.bincount((np.arange(p) ** 2) % p, minlength=p).astype(float)
        cc = np.bincount((np.arange(p, dtype=object) ** 3 % p).astype(int),
                         minlength=p).astype(float)
        NP[p] = np.rint(np.fft.irfft(np.fft.rfft(w) * np.conj(np.fft.rfft(cc)),
                                     n=p)).astype(int)
    return NP


def verify_tables(NP, trials=40, seed=1):
    """Direct character-sum check of the FFT tables on random (p, r)."""
    import random
    random.seed(seed)
    ps = list(NP)
    for _ in range(trials):
        p = random.choice(ps)
        r = random.randrange(p)
        direct = 0
        for x in range(p):
            t = (x * x * x + r) % p
            direct += 1 if t == 0 else (2 if pow(t, (p - 1) // 2, p) == 1 else 0)
        if direct != NP[p][r]:
            raise AssertionError(f"FFT table wrong at p={p}, r={r}: "
                                 f"{NP[p][r]} != {direct}")
    return trials


def good_prime_product(s, sign, NP, z):
    """prod over good p <= z, p not dividing a, of rho_p = N_p / p."""
    total = 0.0
    for p, tab in NP.items():
        if p > z:
            continue
        if s % p == 0:
            continue
        total += math.log(tab[(sign * s) % p] / p)
    return math.exp(total)


# --- stages ------------------------------------------------------------------
def stage_skeleton(out_json="skeleton.json", out_txt="skeleton_classes.txt"):
    print("== skeleton: sixth-power-free curve classes ==")
    anc = anchors()
    skel = defaultdict(list)
    for s, a, m in anc:
        for sign in (+1, -1):
            skel[sixth_power_free(a, m, sign)].append([str(s), a, m, sign])
    ks = sorted(skel, key=abs)
    json.dump({str(k): v for k, v in skel.items()}, open(out_json, "w"))
    with open(out_txt, "w") as f:
        f.write("\n".join(str(k) for k in ks) + "\n")
    sizes = sorted((len(v) for v in skel.values()), reverse=True)
    print(f"  anchors: {len(anc)}   curve-orientations: {2 * len(anc)}")
    print(f"  distinct classes: {len(ks)}   max |k0|: {max(abs(k) for k in ks)}")
    print(f"  largest classes: {sizes[:4]}   singletons: {sizes.count(1)}")
    print(f"  wrote {out_json}, {out_txt}")
    return anc, skel


def stage_local(anc):
    print("\n== local: coprimality, p-adic densities, truncation study ==")

    # decay table for the rigid class, quoted in the note
    NP_big = point_count_tables(30000)
    print(f"  FFT point-count tables verified on "
          f"{verify_tables(NP_big)} random (p, r) samples")
    decay = [good_prime_product(1, +1, NP_big, z)
             for z in (300, 1000, 3000, 30000)]
    print("  rigid class y^2 = X^3 + 1, prod_{5<=p<=z} rho_p at "
          "z = 300, 1000, 3000, 30000:")
    print("    " + ", ".join(f"{v:.3f}" for v in decay))
    del NP_big

    NP = point_count_tables(max(CUTOFFS))
    raw = cop = 0.0
    stab_checked = stab_fail = 0
    per_anchor = {}
    totals = {z: 0.0 for z in CUTOFFS}

    for s, a, m in anc:
        L = ell(a)
        e_plus = e_raw_branch(s, +1) + e_raw_closing(s)
        e_minus = e_raw_branch(s, -1)
        raw += e_plus + e_minus
        cop += L * (e_plus + e_minus)

        d = {}
        for sign in (+1, -1):
            d2 = delta(s, sign, 2, K2) if a % 2 else 1.0
            d3 = delta(s, sign, 3, K3) if a % 3 else 1.0
            if a % 2:                       # stabilization certificate
                stab_checked += 1
                if abs(delta(s, sign, 2, K2 + 1) - d2) > 1e-12:
                    stab_fail += 1
            if a % 3:
                stab_checked += 1
                if abs(delta(s, sign, 3, K3 + 1) - d3) > 1e-12:
                    stab_fail += 1
            d[sign] = d2 * d3

        mass = {}
        for z in CUTOFFS:
            contrib = {}
            for sign, e in ((+1, e_plus), (-1, e_minus)):
                contrib[sign] = (L * e * d[sign]
                                 * good_prime_product(s, sign, NP, z))
            totals[z] += contrib[+1] + contrib[-1]
            mass[z] = contrib
        per_anchor[s] = {"a": a, "m": m,
                         "mass_1000": {1: mass[1000][+1], -1: mass[1000][-1]}}

    print(f"  stabilization: N_(K+1) = p N_K at all {stab_checked} "
          f"checked densities, {stab_fail} failures")
    print(f"  raw model total:                 {raw:.2f}")
    print(f"  x coprimality l(a):              {cop:.2f}   "
          f"(mean factor {cop / raw:.3f})")
    for z in CUTOFFS:
        print(f"  E({z}):{'':>{max(0, 5 - len(str(z)))}}"
              f"                        {totals[z]:.2f}")

    # the 6 | m layer
    lam6 = sum(v["mass_1000"][1] + v["mass_1000"][-1]
               for s, v in per_anchor.items() if v["m"] % 6 == 0)
    n6 = sum(1 for v in per_anchor.values() if v["m"] % 6 == 0)
    print(f"  lambda_6 ({n6} anchors with 6 | m): {lam6:.2f} of E(1000)")

    # discrimination: top-decile mass share vs solution-bearing anchors
    ranked = sorted(per_anchor.items(),
                    key=lambda kv: -(kv[1]["mass_1000"][1] + kv[1]["mass_1000"][-1]))
    order = {s: i + 1 for i, (s, _) in enumerate(ranked)}
    dec = max(1, round(len(ranked) / 10))
    top_mass = sum(v["mass_1000"][1] + v["mass_1000"][-1] for _, v in ranked[:dec])
    print(f"  top decile ({dec} anchors) carries "
          f"{top_mass / totals[1000]:.1%} of E(1000) mass")
    hits = 0
    for a, m, sign in KNOWN:
        s = a ** m
        r = order.get(s)
        if r is not None and r <= dec:
            hits += 1
        print(f"    {a}^{m}: rank {r} of {len(ranked)}"
              f"{'  (top decile)' if r and r <= dec else ''}")
    print(f"  solution-bearing anchors in top decile: {hits}/{len(KNOWN)} "
          f"= {hits / len(KNOWN):.1%}")

    return {"raw": raw, "coprimality": cop,
            "E": {str(z): totals[z] for z in CUTOFFS},
            "lambda6": lam6, "decay_rigid": decay,
            "top_decile_mass_share": top_mass / totals[1000],
            "top_decile_hits": hits,
            "per_anchor": {str(s): v for s, v in per_anchor.items()}}


def stage_ranks(skel_path="skeleton.json", ranks_path="skeleton_ranks.txt",
                local=None):
    print("\n== ranks: decomposition of the diagnostic mass ==")
    try:
        ranks = {}
        for line in open(ranks_path):
            parts = line.strip().split("|")
            if len(parts) == 4 and parts[1] != "ERR":
                ranks[int(parts[0])] = int(parts[1])
                if parts[1] != parts[2]:
                    print(f"  ! class {parts[0]}: rank bounds "
                          f"[{parts[1]}, {parts[2]}] do not coincide")
    except FileNotFoundError:
        print(f"  {ranks_path} not found -- run analysis/skeleton_ranks.gp first")
        return None
    skel = json.load(open(skel_path))
    k_of = {}
    for k, members in skel.items():
        for s, a, m, sign in members:
            k_of[(int(s), sign)] = int(k)

    by_rank_mass = defaultdict(float)
    by_rank_orient = defaultdict(int)
    for s_str, v in local["per_anchor"].items():
        s = int(s_str)
        for sign in (1, -1):
            r = ranks[k_of[(s, sign)]]
            by_rank_mass[r] += v["mass_1000"][sign] if isinstance(v["mass_1000"], dict) \
                else v["mass_1000"][str(sign)]
            by_rank_orient[r] += 1
    total = sum(by_rank_mass.values())
    n_class = defaultdict(int)
    for k, r in ranks.items():
        n_class[r] += 1

    print(f"  {'rank':>5} {'classes':>8} {'orientations':>13} "
          f"{'E(1000) mass':>13} {'share':>7}")
    for r in sorted(by_rank_mass):
        print(f"  {r:>5} {n_class[r]:>8} {by_rank_orient[r]:>13} "
              f"{by_rank_mass[r]:>13.2f} {by_rank_mass[r] / total:>6.1%}")
    dead = by_rank_mass[0]
    obs = 7
    print(f"  rank-0 mass (provably dead): {dead:.2f} = {dead / total:.1%} "
          f"of diagnostic mass")
    print(f"  gap E(1000) - observed = {total:.2f} - {obs} = {total - obs:.1f}; "
          f"rank-0 accounts for {dead / (total - obs):.1%} of it")
    print(f"  positive-rank mass: {total - dead:.2f} against {obs} observed")
    return {"by_rank_mass": dict(by_rank_mass),
            "by_rank_orientations": dict(by_rank_orient),
            "by_rank_classes": dict(n_class),
            "ranks": {str(k): v for k, v in ranks.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["skeleton", "local", "ranks", "all"])
    ap.add_argument("--ranks", default="skeleton_ranks.txt")
    ap.add_argument("--out", default="coprimality_results.json")
    a = ap.parse_args()

    results = {}
    if a.stage in ("skeleton", "all"):
        stage_skeleton()
    if a.stage in ("local", "all"):
        results["local"] = stage_local(anchors())
    if a.stage in ("ranks", "all"):
        if "local" not in results:
            print("stage 'ranks' needs the local masses; run 'all' instead")
            sys.exit(2)
        r = stage_ranks(ranks_path=a.ranks, local=results["local"])
        if r:
            results["ranks"] = r
    if results:
        json.dump(results, open(a.out, "w"))
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
