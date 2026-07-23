# Preregistration: Hall census decade [10^10, 10^11]

**Status: registered before any scan of this decade. Commit timestamp is the notary.**

## Exact family counts (Proposition 2; enumerable in advance)

Command: `python3 analysis/family_enumerate.py --x-lo 1e10 --x-hi 1e11`

| cell | F2 members |
|---|---|
| theta = 0.8 | **832** |
| theta = 0.9 | **3283** |

## Competing cell predictions

Uniform-model backgrounds: theta=0.8: 3318; theta=0.9: 37797.
Scaling-law depletion (Empirical Law 1, D(s) = 0.86 s^0.191):
D(0.8, dec 10) = 3.2%; D(0.9, dec 10) = 5.0%.

| cell | H_uniform | H_scaling |
|---|---|---|
| theta = 0.8 | 3318 + 832 = **4150** | 3212 + 832 = **4045** |
| theta = 0.9 | 37797 + 3283 = **41080** | 35907 + 3283 = **39190** |

Separation at theta = 0.9: ~1890 events vs Poisson width ~200 (~9 sigma).

## Protocol

1. This file is committed BEFORE `hall_census.py --x-max 1e11` is launched.
2. After the scan, the resulting `state.json` is committed and the decade-10
   cells are compared against the table above without adjustment.
3. Family membership is individually checkable via the identity
   (t^3 + 6wt + j)^2 - (t^2 + 4w)^3 = 8w^3 + j^2.
