# fermat-hall-census

Code, data, and verification infrastructure for the Fermat-Hall Census Program: A three-paper investigation by M. M. Ross (2026):

| Title | DOI |
|---|---|
| Audited Censuses of Two Generalized Fermat Families: Coverage Ledgers, Closure Certificates, and Datasets | [10.5281/zenodo.21584509](https://doi.org/10.5281/zenodo.21584509) |
| Coprimality Density and the Proper-Solution Deficit in the Generalized Fermat Family {2,3,m} | [10.5281/zenodo.21584627](https://doi.org/10.5281/zenodo.21584627) |
| The Hall Near-Miss Census to 10^11: Exact Families, a Depleted Residual Population, and a Preregistered Scaling Test | [10.5281/zenodo.21584727](https://doi.org/10.5281/zenodo.21584727) |

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21517861.svg)](https://doi.org/10.5281/zenodo.21517861)[![tier-1 selftests](https://github.com/michaelmross/fermat-hall-census/actions/workflows/selftest.yml/badge.svg)](https://github.com/michaelmross/fermat-hall-census/actions/workflows/selftest.yml)

The papers' claims are exhaustive ("every solution in region R is listed") or
preregistered ("these predictions were fixed before the scan"), so this
repository is part of the argument, not supplementary material: every number
maps to a command below, and the commit history is the notary for the
preregistration. The papers themselves are archived on Zenodo at the DOIs
above and are not duplicated here.

## Layout

    scanners/       fc23m_scan.py, beal33m_scan.py, hall_census.py
    analysis/       deficit_analysis, lift_closure, merge_fc23m, hall_analysis,
                    family_enumerate, scaling_fit, spectrum_table, merge_states,
                    coprimality_pipeline, skeleton_ranks.gp, plot_skeleton
    verification/   run_selftests.sh, audit_hits.py, sampled_resolve.py,
                    coverage_manifest.py (Hall), coverage_manifest_fc23m.py,
                    verify_closing_census.py, verify_deep_skeletons.py,
                    preregistration records
    data/           committed evidence files + SHA-256 checksums
                    fc23m/ holds the merged census hits.jsonl together with
                    run1/ and gapfill/, the two runs whose union it is

Dependencies: Python 3.11+, numpy. Optional: sympy and matplotlib (coprimality
pipeline and figure), PARI/GP 2.14+ (rank census), gmpy2 (2-3x on exact tests).

**All commands are run from the repository root.**

## Verification tiers

| tier | command | time | what it certifies |
|---|---|---|---|
| 1 | `bash verification/run_selftests.sh` | ~1 min | the whole tier-1 suite below, in order (runs in CI on every push) |
| 1 | `sha256sum -c data/CHECKSUMS` | seconds | committed evidence matches its manifest |
| 1 | `python3 verification/audit_hits.py data/fc23m/hits.jsonl` | seconds | every record re-verifies exactly; all seven known coprime solutions are present *in the artifact*; orientation totals unchanged |
| 1 | `python3 verification/coverage_manifest_fc23m.py --run data/fc23m/run1/ledger.jsonl:0 --run data/fc23m/gapfill/ledger.jsonl:1e14 --s-max 1e16 --x-max 1e9` | seconds | the {2,3,m} coverage *rectangle* -- all 478 anchors x all x in [1,1e9], both phases -- is covered by the union of the runs |
| 1 | `python3 verification/verify_closing_census.py` | seconds | the anchor-largest orientation agrees record-for-record with an independent from-scratch enumeration sharing no code with the scanner |
| 1 | `python3 analysis/lift_closure.py data/fc23m/hits.jsonl` | seconds | every improper solution descends under the sixth-power lift; the Catalan tower contains its predicted 28 members |
| 1 | `python3 verification/verify_deep_skeletons.py` | seconds | every {3,3,m} solution with m >= 49 is one of the two skeleton families, each admissible value realized exactly once |
| 1 | `python3 verification/coverage_manifest.py data/hall --expect 1e10 1e11` | seconds | the decade-10 worker ranges are complete, disjoint, and contiguous |
| 2 | `python3 verification/sampled_resolve.py --hits data/beal33m/hits_1e30.jsonl --s-lo 1e24 --s-hi 1e30` | minutes | sampled anchors re-solved independently match the committed census |
| 3 | scanner commands below | hours to days | full reproduction (wall times in the committed ledgers) |

The `{3,3,m}` census is committed compressed; unzip before tier-2 or the
deep-skeleton check: `unzip -o data/beal33m/hits_1e30.zip -d data/beal33m/`
(`run_selftests.sh` does this itself).

Tier 1 separates two things that are easy to conflate: whether the
*instrument* rediscovers known solutions when run, and whether the *artifact*
contains them. A scanner self-test certifies only the first. The coverage gap
described in Sec. 2.3 of the censuses paper passed every instrument-level
check while 41% of the census was absent from the committed ledger, so the
artifact-level audits above are not redundant with the scanner self-tests.

## Claim-to-artifact table

### Censuses paper ([21584509](https://doi.org/10.5281/zenodo.21584509))

| claim | evidence | regenerate with |
|---|---|---|
| 855 solutions: exactly the 7 known coprime solutions, 848 improper (a^m <= 1e16, x <= 1e9) | `data/fc23m/hits.jsonl` | two runs then merge -- see "Reproducing the {2,3,m} census" below |
| census totals by orientation: 412 / 285 / 158 | `data/fc23m/hits.jsonl` | `python3 analysis/lift_closure.py data/fc23m/hits.jsonl` |
| coverage: all 478 anchors x all x in [1,1e9], both phases | `data/fc23m/run1/ledger.jsonl`, `data/fc23m/gapfill/ledger.jsonl` | `python3 verification/coverage_manifest_fc23m.py --run data/fc23m/run1/ledger.jsonl:0 --run data/fc23m/gapfill/ledger.jsonl:1e14 --s-max 1e16 --x-max 1e9` |
| anchor-largest orientation certified against an independent enumeration | `data/fc23m/hits.jsonl` | `python3 verification/verify_closing_census.py` |
| every improper solution is a sixth-power lift (Lemma 1); 608 primitive seeds; 28-lift Catalan tower | `data/fc23m/hits.jsonl` | `python3 analysis/lift_closure.py data/fc23m/hits.jsonl` |
| boundary repair: 196 anchors swept, no solution missed | `scanners/fc23m_scan.py` (exact integer cube root + regression test) | `python3 scanners/fc23m_scan.py --selftest` |
| coverage repair: 347 solutions recovered, one of them a known coprime solution | `data/fc23m/run1/` (pre-repair) against `data/fc23m/hits.jsonl` | `python3 analysis/merge_fc23m.py data/fc23m/run1/hits.jsonl --dry-run` reproduces the pre-repair state and refuses it |
| deficit accounting (33.09 raw mass, 99.91% exhausted) | -- | `python3 analysis/deficit_analysis.py --ledger data/fc23m/*/ledger.jsonl --hits data/fc23m/hits.jsonl --segment 0:1e14 --segment 1e14:1e16 --x-max 1e9` |
| 193,776 solutions, 0 coprime, Euler spectral gap, two-skeleton deep tail (m >= 49) | `data/beal33m/hits_1e30.zip` | `python3 scanners/beal33m_scan.py --s-max 1e30`; tables: `python3 analysis/spectrum_table.py data/beal33m/hits_1e30.jsonl`; regularity: `python3 verification/verify_deep_skeletons.py` |

### Reproducing the {2,3,m} census

The census is the union of two runs and is not reproduced by any single
committed scanner invocation. The first run scanned anchors to 1e16 but,
having been started at a 1e14 ceiling and resumed after the ceiling was
raised, skipped the anchor-largest phase for the new anchors entirely and
covered them in the sieve phases only for x > 1e8. The gap-fill closed exactly
that rectangle. Both ledgers are committed; the coverage audit verifies their
union.

    python3 scanners/fc23m_scan.py --s-max 1e16 --x-max 1e9 --out run1
    python3 scanners/fc23m_scan.py --s-min 1e14 --s-max 1e16 --x-max 1e8 --out gapfill
    python3 analysis/merge_fc23m.py run1/hits.jsonl gapfill/hits.jsonl --out hits.jsonl

The merge re-verifies every identity and refuses to write an output missing
any of the seven known coprime solutions. A fresh single run at
`--s-max 1e16 --x-max 1e9` covers the same rectangle in one pass and should
reproduce the same 855 records; the two-run form is what the committed ledgers
document. The scanner now refuses to resume across a changed configuration,
which is what would have prevented the gap.

### Coprimality note ([21584627](https://doi.org/10.5281/zenodo.21584627))

| claim | evidence | regenerate with |
|---|---|---|
| l(a) reduces 33.1 to 19.3; E(z) = 24.6 / 26.1 / 27.6; rigid-class decay; stabilization certificate; lambda_6 = 1.30; top-decile mass 73.8% against 6/7 | -- | `python3 analysis/coprimality_pipeline.py all` |
| 670 sixth-power-free classes, ranks determined unconditionally | `data/coprimality/skeleton_ranks.txt` | `python3 analysis/coprimality_pipeline.py skeleton` then `gp -q analysis/skeleton_ranks.gp` |
| rank decomposition: 3.17 / 13.86 / 8.55 / 0.48 (12.2% of mass, 16.6% of gap) | -- | `python3 analysis/coprimality_pipeline.py all` (ranks stage) |
| Figure 1 (the Mordell skeleton) | -- | `python3 analysis/plot_skeleton.py --out skeleton_ranks_note.pdf` |

### Hall census paper ([21584727](https://doi.org/10.5281/zenodo.21584727))

| claim | evidence | regenerate with |
|---|---|---|
| theta-band census to 1e11 | `data/hall/state_full.json`, `data/hall/w1..w6/` | `python3 scanners/hall_census.py --x-min A --x-max B --out wN` per worker; merge: `python3 analysis/merge_states.py state_full.json data/hall/state.json data/hall/w*/state.json`; table: `python3 analysis/hall_analysis.py data/hall/state_full.json 1e11` |
| decade-10 coverage: 6 disjoint complete ranges over [1e10, 1e11] | `data/hall/w1..w6/checkpoint.json` | `python3 verification/coverage_manifest.py data/hall --expect 1e10 1e11 --write` |
| exact family counts (dec 9: 465 / 1568; dec 10: 832 / 3283) | -- | `python3 analysis/family_enumerate.py --x-lo 1e9 --x-hi 1e10` (and `--x-lo 1e10 --x-hi 1e11`) |
| Empirical Law: gamma = 0.191 +/- 0.013 unweighted, 0.212 +/- 0.049 weighted | -- | `python3 analysis/scaling_fit.py data/hall/state.json` |
| decade-10 verdict: no-depletion excluded (-6.85 at theta = 0.9); committed law exceeded (+2.53) | `verification/preregistration_dec10.md` (committed before the scan) + `data/hall/state_full.json` | compare the preregistration table against `python3 analysis/hall_analysis.py data/hall/state_full.json 1e11` |
| record points and power channel | `data/hall/hall_hits.jsonl`, `power_hits.jsonl` | cross-certification: every power-channel equation must appear in `data/fc23m/hits.jsonl` |

## Preregistration protocol

Predictions for untouched regions are committed before those regions are
scanned; the commit-then-scan ordering is verifiable from git history, and the
Zenodo release freezes it. Two such tests have been executed:

- **Decade 9** rejected both registered hypotheses (uniform continuum and
  constant depletion), and the scaling law was then fitted to the accumulated
  data. Fixed before scanning in the working record, but not independently
  notarized; the Hall paper says so explicitly.
- **Decade 10** (`verification/preregistration_dec10.md`, committed before any
  scanning of that decade) registered exact family counts 832 / 3283 and cell
  totals under both hypotheses. Outcome: the undepleted-continuum hypothesis
  is excluded for the second consecutive decade (-6.85 Poisson-standardized at
  theta = 0.9), while the committed law predicted the discriminating cell to
  within 1.3% but was exceeded by +2.53 standardized units. The depletion is
  real and shallower than the fitted power; the law survives as an interior
  approximation, not a global form.

## Data integrity
	
	find data -type f ! -name CHECKSUMS ! -name hits_1e30.jsonl | sort | xargs sha256sum > data/CHECKSUMS
    sha256sum -c data/CHECKSUMS

from the repository root. Everything under `data/` is covered except the
manifest itself and the decompressed working copy `data/beal33m/hits_1e30.jsonl`,
which is gitignored and whose own SHA-256 is recorded in `data/README.md` --
itself hashed here, so the chain closes. `.gitattributes` marks `data/**` as
`-text`, so the evidence bytes are identical on every platform and the
manifest verifies on Linux and Windows alike.

## Licenses

Code: MIT. Data: CC-BY 4.0 (`data/LICENSE-DATA`).
