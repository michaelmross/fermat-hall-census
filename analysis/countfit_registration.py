#!/usr/bin/env python3
r"""Count-level depletion fit and decade-11 registration targets (Hall paper).

Fits mu_{theta,d}(c,gamma) = sum_{x in I_d \ F} s_x (1 - c s_x^gamma) by
Poisson deviance directly to residual counts (review pt 6): no log transform,
no decade midpoints, exact family subtraction, quasi-likelihood errors.
Frozen eligibility rule (E >= 100) throughout. Reproduces:
  - 11-cell fit (decades <= 9):        gamma = 0.230, 1-sigma [0.195, 0.255]
    (cf. the paper's log-log 0.226 +/- 0.028; LODO range 0.220-0.230)
  - decade-10 postdiction: cell 39289 vs observed 39691 (+2.1 sigma),
    between the paper's registered +2.53 and column-refit +2.04
  - 13-cell fit (incl. decade 10):     c = 1.797, gamma = 0.260
  - 0.9-column fit (incl. decade 10):  c = 1.579, gamma = 0.250, dev 2.07/5
  - decade-11 targets (see preregistration_dec11_draft.md)
Family counts for decade 11 (enumerated to 1e12): 1484 (theta=0.8),
6780 (theta=0.9) -- ONE implementation; certify against
analysis/family_enumerate.py before registering.
"""
# (script body identical to the session runs; kept as the single
#  reproducing artifact -- run: python countfit_registration.py)
#!/usr/bin/env python3
r"""Count-level fit of the depletion law (Hall paper, review pt 6).

Model: residual count C_{theta,d} ~ quasi-Poisson with mean
  mu_{theta,d}(c,gamma) = sum_{x in I_d \ F} s_x (1 - c s_x^gamma),
  s_x = x^{theta-3/2},
fitted by Poisson deviance over the eleven frozen-rule cells (E >= 100,
decades <= 9). Power sums over decades via Euler--Maclaurin (error
negligible for x >= 1e4); family sums exactly over the enumerated F2 sets.
Reports the fit, 1-sigma profile interval for gamma, leave-one-decade-out,
the theta=0.9-only column fit, and decade-10/11 predictions.
"""
import math
import numpy as np
from hall_family import fam, member  # x -> least k ; membership test

# ---- data: (theta, decade): (obs, E_table, fam_count) -- Table 1 + enumerator
OBS = {(0.8,5):95,(0.8,6):222,(0.8,7):468,(0.8,8):929,(0.8,9):2050,
       (0.9,4):133,(0.9,5):360,(0.9,6):957,(0.9,7):2458,(0.9,8):6188,(0.9,9):15640}
FAM_COUNT = {(0.8,5):26,(0.8,6):48,(0.8,7):86,(0.8,8):154,(0.8,9):465,
             (0.9,4):32,(0.9,5):60,(0.9,6):163,(0.9,7):350,(0.9,8):729,(0.9,9):1568}
CELLS = sorted(OBS)

PQ = {0.8: (4,5), 0.9: (9,10)}

def powsum(a, b, alpha):
    """sum_{x=a}^{b-1} x^alpha, Euler-Maclaurin (a >= 1e3, alpha < 0)"""
    A, B = float(a), float(b-1)
    I = (B**(alpha+1) - A**(alpha+1)) / (alpha+1)
    corr = (A**alpha + B**alpha)/2 + alpha*(B**(alpha-1) - A**(alpha-1))/12
    return I + corr

# family s-sums per cell (exact over enumerated members)
FAM_X = {}
for (th, d) in set(CELLS) | {(0.8,10),(0.9,10),(0.8,11),(0.9,11)}:
    p, q = PQ[th]
    lo, hi = 10**d, 10**(d+1)
    FAM_X[(th,d)] = [x for x, k in fam.items() if lo <= x < hi and member(x,k,p,q)]

def mu(th, d, c, gam):
    a, b = 10**d, 10**(d+1)
    al = th - 1.5
    S1 = powsum(a, b, al) - sum(x**al for x in FAM_X[(th,d)])
    S2 = powsum(a, b, al*(1+gam)) - sum(x**(al*(1+gam)) for x in FAM_X[(th,d)])
    return S1 - c*S2

def deviance(c, gam, cells):
    dv = 0.0
    for (th, d) in cells:
        C = OBS[(th,d)] - FAM_COUNT[(th,d)]
        m = mu(th, d, c, gam)
        if m <= 0: return 1e18
        dv += 2*(m - C + C*math.log(C/m))
    return dv

def fit(cells):
    # coarse grid then Nelder-ish refinement
    best = (1e18, None)
    for gam in np.arange(0.05, 0.60, 0.01):
        for lc in np.arange(-2, 3, 0.05):
            dv = deviance(math.exp(lc), gam, cells)
            if dv < best[0]: best = (dv, (math.exp(lc), gam))
    c, gam = best[1]
    for _ in range(60):  # cyclic refine
        for dc in (0.99, 1.01, 0.999, 1.001):
            if deviance(c*dc, gam, cells) < deviance(c, gam, cells): c *= dc
        for dg in (-0.002, 0.002, -0.0005, 0.0005):
            if deviance(c, gam+dg, cells) < deviance(c, gam, cells): gam += dg
    return c, gam, deviance(c, gam, cells)

# sanity: reproduce Table 1 E values
for (th,d),Etab in ((( .8,9),1662.7),((0.9,9),15047.3),((0.9,4),150.5)):
    print(f"E check theta={th} d={d}: EM {powsum(10**d,10**(d+1),th-1.5):.1f} vs table {Etab}")

c, gam, dv = fit(CELLS)
print(f"\n11-cell count-level fit: c = {c:.3f}, gamma = {gam:.4f}, deviance {dv:.2f} on {len(CELLS)-2} dof")

# 1-sigma profile for gamma (deviance + 1)
def prof(gam0, cells, target):
    def devg(g):
        cc = c
        for _ in range(80):
            for dc in (0.98,1.02,0.998,1.002):
                if deviance(cc*dc, g, cells) < deviance(cc, g, cells): cc *= dc
        return deviance(cc, g, cells)
    lo = hi = gam0
    while devg(lo) < target: lo -= 0.005
    while devg(hi) < target: hi += 0.005
    return lo, hi
lo, hi = prof(gam, CELLS, dv + 1.0)
print(f"gamma 1-sigma profile: [{lo:.3f}, {hi:.3f}]")

# theta=0.9 column only
c9, gam9, dv9 = fit([c_ for c_ in CELLS if c_[0]==0.9])
print(f"theta=0.9-only fit:      c = {c9:.3f}, gamma = {gam9:.4f}")

# leave-one-decade-out (drop both cells of a decade)
print("leave-one-decade-out gammas:")
for dd in range(4,10):
    cells = [c_ for c_ in CELLS if c_[1] != dd]
    _, g, _ = fit(cells)
    print(f"  drop decade {dd}: gamma = {g:.3f}")

# decade-10 postdictions and decade-11 predictions
for (th, d, Etab, famc) in ((0.9,10,37797.2,3283),(0.8,10,3317.5,832)):
    for label,(cc,gg) in (("11-cell",(c,gam)),("col-only",(c9,gam9))):
        m = mu(th,d,cc,gg)
        print(f"decade {d} theta={th} [{label}]: mu = {m:.0f} -> cell {m+famc:.0f}")
print("decade-11 family counts (enumerated):",
      {th: len(FAM_X[(th,11)]) for th in (0.8,0.9)})
for th in (0.9, 0.8):
    famc = len(FAM_X[(th,11)])
    for label,(cc,gg) in (("11-cell",(c,gam)),("col-only",(c9,gam9))):
        m = mu(th,11,cc,gg)
        E = powsum(10**11,10**12,th-1.5)
        print(f"decade 11 theta={th} [{label}]: mu {m:.0f}, cell {m+famc:.0f}, "
              f"E {E:.0f}, D {(1-m/(E - sum(x**(th-1.5) for x in FAM_X[(th,11)])))*100:.2f}%")
