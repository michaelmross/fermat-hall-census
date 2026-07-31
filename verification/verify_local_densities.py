#!/usr/bin/env python3
"""Regression test of the closed-form local densities at p = 2, 3.

Proposition (stabilization lemma): for gcd(c, 6) = 1 and
N_K = #{(x, y) mod p^K : y^2 = x^3 + c}, delta_p = N_K / p^K, every K >= 3:

    delta_2 = 5/2 if c = 1 (mod 8), else 1/2
    delta_3 = 5/3 if c = +-1 (mod 9), else 2/3

and consequently N_{K+1} = p * N_K for all K >= 3.

Three phases:
  (1) brute-force the closed forms over small c at a spread of levels;
  (2) enumerate the canonical anchors a^m <= CEIL independently (must be 478);
  (3) re-run the paper's own computation -- every anchor, both signs, at the
      working levels 2^10 / 3^7 and one further level 2^11 / 3^8 -- and
      assert each count equals the closed form exactly.

Phase (3) supersedes the former stabilization check: the same counts that
were previously certified only against one further lifting level are now
certified against the proposition. Note the case tally: delta_2 applies to
the 2*236 = 472 odd-anchor curve-orientations and delta_3 to the 2*319 = 638
anchor-orientations with 3 not dividing a; the paper's "956 cases" counts
curve-orientations, each checked at its applicable primes.

Exit code 0 iff all checks pass. Pure stdlib; exact rational arithmetic.
"""

from fractions import Fraction
import sys
import time

CEIL = 10**16          # anchor ceiling of the {2,3,m} census
M_MIN = 7              # canonical exponent floor
BRUTE_C_RANGE = 200    # brute-force |c| range for phase (1)
BRUTE_K = {2: (3, 4, 5, 6, 7, 8), 3: (3, 4, 5, 6)}
WORKING_K = {2: (10, 11), 3: (7, 8)}   # paper's levels plus one


def predicted(p, c):
    if p == 2:
        return Fraction(5, 2) if c % 8 == 1 else Fraction(1, 2)
    if p == 3:
        return Fraction(5, 3) if c % 9 in (1, 8) else Fraction(2, 3)
    raise ValueError(p)


def root_table(M):
    roots = {}
    for y in range(M):
        t = y * y % M
        roots[t] = roots.get(t, 0) + 1
    return roots


def NK(p, K, c, roots=None):
    M = p ** K
    if roots is None:
        roots = root_table(M)
    return sum(roots.get((x * x * x + c) % M, 0) for x in range(M))


def canonical_anchors(ceil, m_min):
    """Distinct values a^m <= ceil with canonical exponent >= m_min, a >= 2."""
    seen = set()
    a = 2
    while a ** m_min <= ceil:
        v = a ** m_min
        while v <= ceil:
            seen.add(v)
            v *= a
        a += 1
    return sorted(seen)


def main():
    t0 = time.time()
    fails = 0

    # (1) closed forms, small c, spread of levels
    n1 = 0
    for c in range(-BRUTE_C_RANGE, BRUTE_C_RANGE + 1):
        if c == 0 or c % 2 == 0 or c % 3 == 0:
            continue
        for p, Ks in BRUTE_K.items():
            pred = predicted(p, c)
            for K in Ks:
                if Fraction(NK(p, K, c), p ** K) != pred:
                    print(f"FAIL brute p={p} c={c} K={K}")
                    fails += 1
                n1 += 1
    print(f"phase 1: {n1} brute-force cases OK" if not fails else
          f"phase 1: {fails} failures")

    # (2) independent anchor enumeration
    anchors = canonical_anchors(CEIL, M_MIN)
    print(f"phase 2: canonical anchors <= {CEIL:.0e}: {len(anchors)} "
          f"(expected 478)")
    if len(anchors) != 478:
        fails += 1

    # (3) the paper's computation, re-certified against the proposition
    tally = {2: 0, 3: 0}
    for p, Ks in WORKING_K.items():
        for K in Ks:
            roots = root_table(p ** K)
            for s in anchors:
                if s % p == 0:
                    continue
                for c in (s, -s):
                    if Fraction(NK(p, K, c, roots), p ** K) != predicted(p, c):
                        print(f"FAIL anchor p={p} K={K} c={'-' if c<0 else '+'}{s}")
                        fails += 1
                    if K == Ks[0]:
                        tally[p] += 1
    print(f"phase 3: working-level cases OK -- "
          f"delta_2: {tally[2]}, delta_3: {tally[3]}, "
          f"each also at one further level")

    print(f"FAILS: {fails}  ({time.time()-t0:.1f}s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
