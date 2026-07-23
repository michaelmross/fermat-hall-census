# fermat-hall-census

Code and data for *A Census of Generalized Fermat Solutions and Hall
Near-Misses, with a Scaling Law for the Distribution of x^(3/2) Modulo One*
(M. M. Ross, 2026). The paper's claims are exhaustive ("every solution in
region R is listed"), so this repository is part of the argument, not
supplementary material: every number in the paper maps to a command below.

## Layout

    scanners/       fc23m_scan.py, beal33m_scan.py, hall_census.py
    analysis/       deficit_analysis, hall_analysis, family_enumerate,
                    scaling_fit, spectrum_table
    verification/   tier-1 selftests, tier-2 sampled re-solving, audits,
                    preregistration records
    data/           committed evidence files + SHA-256 checksums

Dependencies: Python 3.11+, numpy. Optional: gmpy2 (2-3x on exact tests).

## Verification tiers

| tier | command | time | what it certifies |
|---|---|---|---|
| 1 | `bash verification/run_selftests.sh` | ~1 min | every instrument re-derives its known solutions (runs in CI on every push) |
| 2 | `python3 verification/sampled_resolve.py --hits data/beal33m/hits_1e30.jsonl --s-lo 1e24 --s-hi 1e30` | minutes | sampled anchors re-solved independently match the committed census |
| 3 | scanner commands below | hours | full reproduction (wall times in committed ledgers) |

## Claim-to-artifact table

| paper claim | evidence | regenerate with |
|---|---|---|
| Sec. 2: 7 known solutions, none new (a^m <= 1e16, x <= 1e9) | `data/fc23m/` | `python3 scanners/fc23m_scan.py --s-max 1e16 --x-max 1e9` (+ gap-fill band); audit: `python3 verification/audit_hits.py data/fc23m/*.jsonl` |
| Sec. 2: deficit accounting (33.36 raw mass, 99.91% exhausted) | -- | `python3 analysis/deficit_analysis.py --ledger data/fc23m/ledger*.jsonl --hits data/fc23m/hits*.jsonl --segment 0:1e14 --segment 1e14:1e16 --x-max 1e9` |
| Sec. 3, Tables 1-2: 193,776 solutions, 0 coprime, spectral gap | `data/beal33m/hits_1e30.jsonl` | `python3 scanners/beal33m_scan.py --s-max 1e30`; tables: `python3 analysis/spectrum_table.py data/beal33m/hits_1e30.jsonl` |
| Sec. 4, Table 3: theta-band census to 1e10 | `data/hall/state.json` | `python3 scanners/hall_census.py --x-max 1e10`; table: `python3 analysis/hall_analysis.py data/hall/state.json 1e10` |
| Sec. 4.2: record points, power channel | `data/hall/hall_hits.jsonl`, `power_hits.jsonl` | cross-certification: every power-channel equation must appear in `data/fc23m/hits*.jsonl` |
| Sec. 4.3: dec-9 preregistered family counts (465 / 1568) | -- | `python3 analysis/family_enumerate.py --x-lo 1e9 --x-hi 1e10` |
| Empirical Law 1: gamma = 0.190 +/- 0.013, c = 0.86, R^2 = 0.98 | -- | `python3 analysis/scaling_fit.py data/hall/state.json` |
| Sec. 4.4: dec-10 preregistration (832 / 3283; 4045/4150; 39190/41080) | `verification/preregistration_dec10.md` | committed BEFORE the dec-10 scan; git timestamp is the notary |

## Preregistration protocol

Predictions for untouched regions are committed before scanning them
(`verification/preregistration_dec10.md`). The commit-then-scan ordering is
verifiable from git history. The decade-9 test of this kind rejected both
prior hypotheses at ~7 sigma each and selected the scaling law; decade 10 is
registered and pending.

## Data integrity

`sha256sum -c data/CHECKSUMS` from the repo root. Releases are archived via
the Zenodo-GitHub integration; the paper cites the release DOI.

## Licenses

Code: MIT. Data: CC-BY 4.0 (`data/LICENSE-DATA`).
