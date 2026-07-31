# data/hall — auxiliary channels

The Hall scan writes three per-worker outputs. Only the first is a census
artifact; the other two are sparse side channels, and **both are legitimately
empty in some workers**. Do not delete zero-byte files here: an empty file is
a positive record that the worker ran and the channel fired nowhere, and
`verification/power_channel.py` globs `w*/power_hits.jsonl` and will break
without them.

## ledger.jsonl — the census accounting channel

One record per scanned block:

```json
{"x_lo": 70000000001, "x_hi": 70001000000, "records": 0, "power_hits": 0, "secs": 3.56}
```

`records` and `power_hits` count that block's contributions to the two
channels below, written on a separate code path from the channels
themselves. Summing them per worker and comparing against the channel files'
line counts is an internal consistency check that costs nothing:

| worker | x-range | ln share | ledger `records` | `hall_hits` lines | `power_hits` |
|---|---|---|---|---|---|
| w1 | [1.0, 2.5]×10¹⁰ | 0.916 | 4 | 4 | 0 |
| w2 | [2.5, 4.0]×10¹⁰ | 0.470 | 2 | 2 | 0 |
| w3 | [4.0, 5.5]×10¹⁰ | 0.318 | 1 | 1 | 0 |
| w4 | [5.5, 7.0]×10¹⁰ | 0.241 | 4 | 4 | 0 |
| w5 | [7.0, 8.5]×10¹⁰ | 0.194 | 0 | 0 | 0 |
| w6 | [8.5, 10.0]×10¹⁰ | 0.163 | 0 | 0 | 0 |

All six agree exactly.

## hall_hits.jsonl — the Hall-ratio channel (NOT the census)

```json
{"x": 10163640792, "y": 1024646265248575, "k": -142463, "r": 0.7077}
```

`k` = y² − x³ and `r` = √x/|k|, the classical Hall ratio; records are
emitted when r exceeds a fixed threshold. This is a rarity channel: 90
records across the entire scan, against 39,691 census records in decade 10
alone at θ = 0.9. It is not a subset of the census hit set and an empty
file here cannot constitute a coverage gap.

**Why the top workers are empty.** For a fixed ratio threshold r₀ the
expected count of x ∈ [A,B] with |k| ≤ √x/r₀ is ln(B/A)/r₀ — proportional
to *log* measure, hence constant per decade. The workers are sharded by
equal *linear* width, so their log shares differ by a factor of 5.6 between
w1 and w6. Observed counts against expectation proportional to log measure
give G² = 7.94 on 5 df (p = 0.16): no significant departure. The expected
counts in w5 and w6 are 0.93 and 0.78, and P(both zero) = 0.18. The only
mild outlier is w4 (4 observed, 1.15 expected), which contributes almost
the whole deviance; the zeros contribute none.

Cross-scan consistency: the decade-10 rate is 4.78 records per unit ln x
(implying r₀ ≈ 0.21), which extrapolates to 77–99 records over the full
scan range against the 90 in the merged file — so the constant-per-log-measure
model holds across the whole census, confirming the sharding and merge lost
nothing.

## power_hits.jsonl — the exact-power channel

Empty in all six workers, which is the correct global result: a hit would
require x³ = y², i.e. x a perfect square, and perfect squares are excluded
by the nonsquare rule. The channel exists so that the exclusion is observed
rather than assumed, and `verification/power_channel.py` reads these files
to assert the null.
