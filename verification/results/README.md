# verification/results — provenance

Outputs of the verification tier. These are **not** census data: every file
here was produced by an implementation that shares no code with the
scanners in `analysis/`, and exists to be compared against the census
artifacts in `data/`. Nothing downstream should read from this directory;
if a published number depends on a file here, it is in the wrong place.

## Files

| file | produced by | what it certifies |
|---|---|---|
| `scan23m_1e+09.json` | `verification/scan_23m_independent.py --xcap 1e9` | the complete {2,3,m} solution set at x ≤ 10⁹ (855 records, 7 proper, 28 Catalan lifts, 0 descent failures) |
| `scan23m_1e9_records.csv` | canonical sort of the above | diffable record set; SHA-256 in CHECKSUMS |
| `scan23m_1e9_per_anchor.csv` | canonical sort of the above | per-anchor solution counts (346 of 478 anchors carry solutions) |
| `verify_33m_1e+16.json` | `verification/verify_33m_census.py 1e16` | provenance row 1: 12,322 anchors, 798 solutions |
| `verify_33m_1e+20.json` | `… 1e20` | provenance row 2: 113,095 anchors, 3,841 solutions |
| `verify_33m_1e+24.json` | `… 1e24` | provenance row 3: 1,076,369 anchors, 18,495 solutions |
| `verify_33m_1e+30.json` | `… 1e30 --procs N --blocks 512` | provenance row 4: 32,744,685 anchors, 193,776 solutions; 512 block hashes for localizing any future disagreement |

## Diff results (record-for-record, not counts)

Both censuses agree with their ledgers as
(orientation, x, y, anchor) tuples, orientation recomputed by exact
arithmetic rather than read from either side's labels:

```
verification/diff_23m_ledger.py --scan verification/results/scan23m_1e+09.json \
    --ledger data/fc23m/hits.jsonl
    -> 855 / 855, IDENTICAL

verification/diff_33m_ledger.py --ckpt verify_33m_1e+30.ckpt \
    --ledger data/beal33m/hits_1e30.jsonl
    -> 193776 / 193776, IDENTICAL
```

## Not archived here, deliberately

- **Checkpoint files** (`*.ckpt`). Working state, hundreds of MB at the
  10³⁰ ceiling, regenerable by re-running. Gitignored. The 10³⁰ checkpoint
  was retained only until the record-level diff had run, since the
  verifier's JSON keeps summaries rather than records.
- **A duplicate dump of the 193,776 {3,3,m} records.** They are identical
  to `data/beal33m/hits_1e30.jsonl` by the diff above; committing a second
  copy adds bytes, not assurance. The block hashes in
  `verify_33m_1e+30.json` are the verifier's fingerprint of that record
  set and are what a future disagreement would be localized against.

## Scope of what these files establish

They establish that each ledger is exactly what its stated specification
yields at its stated ceiling — the property the {2,3,m} coverage gap
violated. They do not test the specification itself: both implementations
enumerate the same canonical anchor sets and the same orientations under
the same ceilings, so a specification defect would be reproduced
identically by both. They say nothing about {2,3,m} completeness above
x = 10⁹ on the two sieve orientations.
