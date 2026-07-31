# fermat-hall-census

Code, data, and verification infrastructure for the Fermat-Hall Census Program: A three-paper investigation by M. M. Ross (2026):

| Title | DOI |
|---|---|
| Audited Censuses of Two Generalized Fermat Families: Coverage Ledgers and Large Datasets | [10.5281/zenodo.21584509](https://doi.org/10.5281/zenodo.21584509) |
| Coprimality Density and the Proper-Solution Deficit in the Generalized Fermat Family {2,3,m} | [10.5281/zenodo.21584627](https://doi.org/10.5281/zenodo.21584627) |
| The Hall Near-Miss Census to 10^11: Exact Families, a Depleted Residual Population, and a Preregistered Scaling Test | [10.5281/zenodo.21584727](https://doi.org/10.5281/zenodo.21584727) |

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21517860.svg)](https://doi.org/10.5281/zenodo.21517860)[![tier-1 selftests](https://github.com/michaelmross/fermat-hall-census/actions/workflows/selftest.yml/badge.svg)](https://github.com/michaelmross/fermat-hall-census/actions/workflows/selftest.yml)

The papers' claims are exhaustive ("every solution in region R is listed") or
preregistered ("these predictions were fixed before the scan"), so this
repository is part of the argument, not supplementary material: every number
maps to a command below, and the commit history is the notary for the
preregistration. The papers are not duplicated here; the DOIs above are
reserved and resolve once the papers are deposited. The badge points to the
repository's *concept* DOI, which always resolves to the current release.

## Layout

    scanners/       fc23m_scan.py, beal33m_scan.py, hall_census.py
    analysis/       deficit_analysis, lift_closure, merge_fc23m, hall_analysis,
                    family_enumerate, hall_family, scaling_fit,
                    countfit_registration, spectrum_table, merge_states,
                    coprimality_pipeline, skeleton_ranks.gp, sha_census.gp,
                    skeleton_figure
    verification/   run_selftests.sh, audit_hits.py, sampled_resolve.py,
                    coverage_manifest.py (Hall), coverage_manifest_fc23m.py,
                    verify_closing_census.py, verify_deep_skeletons.py,
                    verify_local_densities.py, scan_23m_independent.py,
                    verify_33m_census.py, diff_23m_ledger.py,
                    diff_33m_ledger.py, mwrank_check.py,
                    results/ (independent-implementation outputs),
                    preregistration records
    data/           committed evidence files + SHA-256 checksums
                    fc23m/ holds the merged census hits.jsonl together with
                    run1/ and gapfill/, the two runs whose union it is;
                    hall/ carries its own README explaining the three
                    per-worker channels and why two of them are empty

Dependencies: Python 3.11+, numpy. Optional: sympy and matplotlib (coprimality
pipeline and figure), PARI/GP 2.14+ (rank census, Sha census), eclib /
`mwrank` (independent rank re-certification; no native Windows build -- use
WSL), gmpy2 (2-3x on exact tests).

**All commands are run from the repository root** unless noted. `python`
below means a Python 3.11+ interpreter. On systems without the unversioned
name -- most Linux distributions -- use `python3`, or work inside an activated
virtualenv, where `python` is always correct. `verification/run_selftests.sh`
picks the working name itself.

## Independent implementations

Five programs in `verification/` re-derive results from the specification
rather than checking the scanners against themselves. None shares code with
the scanner it checks; each re-derives its own anchor set from the stated
ceiling rather than reading one.

| program | what it re-derives | outcome |
|---|---|---|
| `scan_23m_independent.py` | the whole {2,3,m} solution set, using a disjoint sieve design (QR lookup tables mod 43200, prime squares, primes 23-61; the scanner uses a 5040 wheel and primes to 41) | all 855 records, written without sight of the ledger |
| `diff_23m_ledger.py` | tuple-level comparison against `data/fc23m/hits.jsonl` | IDENTICAL; 855/855, none on either side alone |
| `verify_33m_census.py` | the {3,3,m} census from the divisor identity alone, at every provenance ceiling | anchor counts, totals, orientation splits and extremal records all reproduced |
| `diff_33m_ledger.py` | tuple-level comparison against `data/beal33m/hits_1e30.jsonl` | IDENTICAL; 193,776/193,776 |
| `mwrank_check.py` | the Mordell rank census under eclib rather than PARI | 647/670 determined, 0 disagreements; the 23 exceptions are explained in `analysis/sha_census.gp` |

In every comparison the orientation is recomputed by exact arithmetic rather
than read from either side's labels, and agreement on totals is treated as
necessary but not sufficient: two implementations can agree on a count while
disagreeing about which records compose it.

What this does **not** establish: agreement between implementations tests
implementation, not specification. Both members of each pair enumerate the
same canonical anchor set under the same ceiling, so a defect in the
specification would be reproduced identically by both.

## Verification tiers

| tier | command | time | what it certifies |
|---|---|---|---|
| 1 | `bash verification/run_selftests.sh` | ~1 min | the whole tier-1 suite below, in order (runs in CI on every push) |
| 1 | `sha256sum -c data/CHECKSUMS` | seconds | committed evidence matches its manifest |
| 1 | `python verification/audit_hits.py data/fc23m/hits.jsonl` | seconds | every record re-verifies exactly; all seven known coprime solutions are present *in the artifact*; orientation totals unchanged |
| 1 | `python verification/coverage_manifest_fc23m.py --run data/fc23m/run1/ledger.jsonl:0 --run data/fc23m/gapfill/ledger.jsonl:1e14 --s-max 1e16 --x-max 1e9` | seconds | the {2,3,m} coverage *rectangle* -- all 478 anchors x all x in [1,1e9], both phases -- is covered by the union of the runs |
| 1 | `python verification/verify_closing_census.py` | seconds | the anchor-largest orientation agrees record-for-record with an independent from-scratch enumeration sharing no code with the scanner |
| 1 | `python verification/diff_23m_ledger.py --scan verification/results/scan23m_1e+09.json --ledger data/fc23m/hits.jsonl` | seconds | **all three** orientations agree record-for-record with the independent scan; exits nonzero on any asymmetric difference |
| 1 | `python analysis/lift_closure.py data/fc23m/hits.jsonl` | seconds | every improper solution descends under the sixth-power lift; the Catalan tower contains its predicted 28 members |
| 1 | `python verification/verify_deep_skeletons.py` | seconds | every {3,3,m} solution with m >= 49 is one of the two skeleton families, each admissible value realized exactly once |
| 1 | `python verification/verify_local_densities.py` | seconds | the closed-form 2- and 3-adic densities agree with direct modular counts at every anchor and orientation |
| 1 | `python verification/coverage_manifest.py data/hall --expect 1e10 1e11` | seconds | the decade-10 worker ranges are complete, disjoint, and contiguous |
| 2 | `python verification/verify_33m_census.py 1e24 --procs 4 --blocks 512` | minutes | the {3,3,m} census recomputed from the divisor identity to 1e24 (1e16 and 1e20 take seconds) |
| 2 | `python verification/sampled_resolve.py --hits data/beal33m/hits_1e30.jsonl --s-lo 1e24 --s-hi 1e30` | minutes | sampled anchors re-solved independently match the committed census |
| 2 | `python verification/mwrank_check.py --pari-ranks data/coprimality/skeleton_ranks.txt` | ~45 s | the 670 rank classes re-certified under eclib. **Exits nonzero by design**: 22 rank-zero classes are not pinned by plain 2-descent. See below. |
| 3 | `python verification/verify_33m_census.py 1e30 --procs 8 --blocks 512` | 1-2 h | the full {3,3,m} census; keep the `.ckpt` file, then `python verification/diff_33m_ledger.py --ckpt verify_33m_1e+30.ckpt --ledger data/beal33m/hits_1e30.jsonl` |
| 3 | scanner commands below | hours to days | full reproduction (wall times in the committed ledgers) |

Both long-running verifiers checkpoint automatically and resume by re-running
the identical command; keep `--blocks 512` fixed across resumes.

The `{3,3,m}` census is committed compressed; unzip before tier-2 or the
deep-skeleton check: `unzip -o data/beal33m/hits_1e30.zip -d data/beal33m/`
(`run_selftests.sh` does this itself).

**On `mwrank_check.py`'s exit status.** 647 of 670 classes are determined and
agree with PARI exactly; the other 23 return `0 <= rank <= selmer-rank = 2`.
That is not a disagreement. PARI's `ellrank` reports a nontrivial Sha[2] on
exactly 25 classes, where the 2-Selmer rank exceeds the rank by 2; both
programs compute the same Selmer rank and differ only in whether they can
attribute the excess to Sha. The 23 are those 25 less the two with rational
2-torsion (k0 = -22^3, -46^3), where eclib's second descent applies.
`gp -q analysis/sha_census.gp` reproduces the 25. Since 22 of the 23 are
rank-zero classes, the independent certification is partial exactly where it
is load bearing -- but those classes carry 0.0145 of the E(1000) diagnostic
mass, so 99.5% of the provably dead mass is certified by both implementations.

Tier 1 separates two things that are easy to conflate: whether the
*instrument* rediscovers known solutions when run, and whether the *artifact*
contains them. A scanner self-test certifies only the first. The coverage gap
described in Appendix B of the censuses paper passed every instrument-level
check while 41% of the census was absent from the committed ledger, so the
artifact-level audits above are not redundant with the scanner self-tests.

## Claim-to-artifact table

### Censuses paper ([21584509](https://doi.org/10.5281/zenodo.21584509))

| claim | evidence | regenerate with |
|---|---|---|
| 855 solutions: exactly the 7 known coprime solutions, 848 improper (a^m <= 1e16, x <= 1e9) | `data/fc23m/hits.jsonl` | two runs then merge -- see "Reproducing the {2,3,m} census" below |
| census totals by orientation: 412 / 285 / 158 | `data/fc23m/hits.jsonl` | `python analysis/lift_closure.py data/fc23m/hits.jsonl` |
| coverage: all 478 anchors x all x in [1,1e9], both phases | `data/fc23m/run1/ledger.jsonl`, `data/fc23m/gapfill/ledger.jsonl` | `python verification/coverage_manifest_fc23m.py --run data/fc23m/run1/ledger.jsonl:0 --run data/fc23m/gapfill/ledger.jsonl:1e14 --s-max 1e16 --x-max 1e9` |
| all 855 records reproduced by an independent scan with a disjoint sieve design | `verification/results/scan23m_1e+09.json`, `scan23m_1e9_records.csv` | `python verification/scan_23m_independent.py --xcap 1e9 --procs 6`, then `diff_23m_ledger.py` |
| every improper solution is a sixth-power lift (Lemma 1); 608 primitive seeds; 28-lift Catalan tower | `data/fc23m/hits.jsonl` | `python analysis/lift_closure.py data/fc23m/hits.jsonl` |
| boundary repair: 196 anchors swept, no solution missed | `scanners/fc23m_scan.py` (exact integer cube root + regression test) | `python scanners/fc23m_scan.py --selftest` |
| coverage repair: 347 solutions recovered, one of them a known coprime solution | `data/fc23m/run1/` (pre-repair) against `data/fc23m/hits.jsonl` | `python analysis/merge_fc23m.py data/fc23m/run1/hits.jsonl --dry-run` reproduces the pre-repair state and refuses it |
| deficit accounting (33.0399 raw mass over all x, 33.0097 scanned, 99.91% exhausted) | `analysis/fermat_ledger.json` | `python analysis/generate_shared_ledger.py`; cross-check with `python analysis/deficit_analysis.py --ledger data/fc23m/*/ledger.jsonl --hits data/fc23m/hits.jsonl --segment 0:1e14 --segment 1e14:1e16 --x-max 1e9` |
| 193,776 solutions, 0 coprime, Euler spectral gap, two-skeleton deep tail | `data/beal33m/hits_1e30.zip` | `python scanners/beal33m_scan.py --s-max 1e30`; tables: `python analysis/spectrum_table.py data/beal33m/hits_1e30.jsonl` |
| the deep tail is exhaustive for m >= 43 at this ceiling (a theorem, not a range regularity) | -- | `python verification/verify_deep_skeletons.py`; the independent verifier returns the same 22 records at m >= 49 from the ceiling alone |
| {3,3,m} provenance rows: anchor counts, orientation splits, exact log2 heights and extremal records at 1e16 / 1e20 / 1e24 / 1e30 | `verification/results/verify_33m_*.json` | `python verification/verify_33m_census.py CEIL --blocks 512` |
| raw model expectation 3.55 for {3,3,m} | -- | closed form in the paper; the anchor mass sum is reproduced by `verify_33m_census.py`'s anchor enumeration |

### Reproducing the {2,3,m} census

The census is the union of two runs and is not reproduced by any single
committed scanner invocation. The first run scanned anchors to 1e16 but,
having been started at a 1e14 ceiling and resumed after the ceiling was
raised, skipped the anchor-largest phase for the new anchors entirely and
covered them in the sieve phases only for x > 1e8. The gap-fill closed exactly
that rectangle. Both ledgers are committed; the coverage audit verifies their
union.

    python scanners/fc23m_scan.py --s-max 1e16 --x-max 1e9 --out run1
    python scanners/fc23m_scan.py --s-min 1e14 --s-max 1e16 --x-max 1e8 --out gapfill
    python analysis/merge_fc23m.py run1/hits.jsonl gapfill/hits.jsonl --out hits.jsonl

The merge re-verifies every identity and refuses to write an output missing
any of the seven known coprime solutions. A fresh single run at
`--s-max 1e16 --x-max 1e9` covers the same rectangle in one pass and should
reproduce the same 855 records; the two-run form is what the committed ledgers
document. The scanner now refuses to resume across a changed configuration,
which is what would have prevented the gap.

### Coprimality note ([21584627](https://doi.org/10.5281/zenodo.21584627))

The pipeline writes `analysis/skeleton.json` and
`analysis/coprimality_results.json`, which every number and the figure below
read from. Both are committed; regenerating them takes about 30 seconds.

| claim | evidence | regenerate with |
|---|---|---|
| l(a) reduces 33.0399 to 19.2497 (mean factor 0.5826); E(z) = 24.45 / 25.91 / 27.47; rigid-class decay; lambda_6 = 1.28; top-decile mass 73.7% against 6/7 | `analysis/coprimality_results.json` | `python analysis/coprimality_pipeline.py all --ranks data/coprimality/skeleton_ranks.txt` |
| closed-form 2- and 3-adic densities (5/2 or 1/2; 5/3 or 2/3) | -- | `python verification/verify_local_densities.py` |
| 670 sixth-power-free classes, ranks determined unconditionally by 2-descent | `data/coprimality/skeleton_ranks.txt` | `cd data/coprimality && python ../../analysis/coprimality_pipeline.py skeleton` then `gp -q ../../analysis/skeleton_ranks.gp` |
| rank decomposition: 3.14 / 13.71 / 8.58 / 0.48 (12.1% of mass, 16.6% of the 18.9-unit gap) | `analysis/coprimality_results.json` | `python analysis/coprimality_pipeline.py all` (ranks stage) |
| independent rank re-certification: 647/670 under eclib, 0 disagreements | `verification/results/mwrank_results.json` | `python verification/mwrank_check.py --pari-ranks data/coprimality/skeleton_ranks.txt` |
| the 23 undetermined classes are the Sha[2] classes without rational 2-torsion (25 minus 2) | -- | `cd analysis && gp -q sha_census.gp` |
| the local machinery relocates solution-bearing anchors: 2/7 in the raw top decile against 6/7 in the diagnostic's | `analysis/coprimality_results.json` | `python analysis/coprimality_pipeline.py all` |
| Figure 1 (the Mordell skeleton) | `skeleton_ranks_note.pdf` | `cd data/coprimality && python skeleton_figure.py` (writes the PDF, a PNG, and `skeleton_ranks_figure_numbers.txt` recording every caption number) |

### Hall census paper ([21584727](https://doi.org/10.5281/zenodo.21584727))

| claim | evidence | regenerate with |
|---|---|---|
| theta-band census to 1e11 | `data/hall/state_full.json`, `data/hall/w1..w6/` | `python scanners/hall_census.py --x-min A --x-max B --out wN` per worker; merge: `python analysis/merge_states.py state_full.json data/hall/state.json data/hall/w*/state.json`; table: `python analysis/hall_analysis.py data/hall/state_full.json 1e11` |
| decade-10 coverage: 6 disjoint complete ranges over [1e10, 1e11] | `data/hall/w1..w6/checkpoint.json` | `python verification/coverage_manifest.py data/hall --expect 1e10 1e11 --write` |
| exact family counts (dec 9: 465 / 1568; dec 10: 832 / 3283) | -- | `python analysis/family_enumerate.py --x-lo 1e9 --x-hi 1e10` (and `--x-lo 1e10 --x-hi 1e11`); third independent enumeration: `python analysis/hall_family.py` |
| the family-subtracted primary table (C, E, D by decade and threshold) | -- | `python analysis/hall_family.py` |
| Empirical Law under the frozen eleven-cell rule: gamma = 0.226 +/- 0.028 unweighted, 0.237 +/- 0.026 weighted; the registered seven-cell law was gamma = 0.191 +/- 0.013 | -- | `python analysis/scaling_fit.py data/hall/state.json` |
| count-level Poisson-deviance fit: gamma = 0.230, 1-sigma profile [0.195, 0.255], deviance 5.17/9, LODO 0.220-0.230; decade-10 postdiction +2.1 sigma | -- | `python analysis/countfit_registration.py` |
| decade-10 verdict: H_IU excluded (-6.85 at theta = 0.9); committed law exceeded (+2.53) | `verification/preregistration_dec10.md` (committed before the scan) + `data/hall/state_full.json` | compare the preregistration table against `python analysis/hall_analysis.py data/hall/state_full.json 1e11` |
| record points and power channel | `data/hall/hall_hits.jsonl`, `power_hits.jsonl` | cross-certification: every power-channel equation must appear in `data/fc23m/hits.jsonl`; see `data/hall/README.md` for why the per-worker channels are empty in w5 and w6 |

## Preregistration protocol

Predictions for untouched regions are committed before those regions are
scanned; the commit-then-scan ordering is verifiable from git history, and the
Zenodo release freezes it. Two such tests have been executed:

- **Decade 9** rejected both registered hypotheses (undepleted residual
  population and constant depletion), and the scaling law was then fitted to
  the accumulated data. Fixed before scanning in the working record, but not
  independently notarized; the Hall paper says so explicitly.
- **Decade 10** (`verification/preregistration_dec10.md`, committed before any
  scanning of that decade) registered exact family counts 832 / 3283 and cell
  totals under both hypotheses. Outcome: the independent-uniform null H_IU is
  excluded for the second consecutive decade
  (-6.85 Poisson-standardized at theta = 0.9), while the committed law
  predicted the discriminating cell to within 1.3% but was exceeded by +2.53
  standardized units. The depletion is real and shallower than the fitted
  power; the law survives as an interior approximation, not a global form.

No registration for decade 11 is committed. When one is, it will appear as a
single artifact alongside the decade-10 record and before any decade-11
scanning: drafts are deliberately kept out of the repository, since a history
containing several candidate cell sets would undercut the claim that one set
was fixed in advance.

## Data integrity
	
	find data -type f ! -name CHECKSUMS ! -name hits_1e30.jsonl | sort | xargs sha256sum > data/CHECKSUMS
    sha256sum -c data/CHECKSUMS

from the repository root. Everything under `data/` is covered except the
manifest itself and the decompressed working copy `data/beal33m/hits_1e30.jsonl`,
which is gitignored and whose own SHA-256 is recorded in `data/README.md` --
itself hashed here, so the chain closes. `.gitattributes` marks `data/**` as
`-text`, so the evidence bytes are identical on every platform and the
manifest verifies on Linux and Windows alike.

Outputs of the independent implementations live in `verification/results/`
with their own README recording what each certifies, and are hashed the same
way. Checkpoint files (`*.ckpt`) are working state, are gitignored, and are
regenerable by re-running the command that produced them.

## Licenses

Code: MIT. Data: CC-BY 4.0 (`data/LICENSE-DATA`).
