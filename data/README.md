# Data

Committed evidence files, one directory per instrument:

- `fc23m/` -- hits.jsonl + ledger.jsonl for the {2,3,m} census
  (main run + gap-fill; coverage a^m <= 1e16, x <= 1e9)
- `beal33m/` -- hits_1e30.jsonl + ledger.jsonl for the exact {3,3,m} census
- `hall/` -- state.json, hall_hits.jsonl, power_hits.jsonl, ledger.jsonl (x <= 1e10)

## Integrity

Generate/verify checksums from the repo root:

    sha256sum data/*/*.json data/*/*.jsonl > data/CHECKSUMS
    sha256sum -c data/CHECKSUMS

CHECKSUMS is committed; the paper cites the release DOI, and the release
freezes these hashes.

Data license: CC-BY 4.0 (see LICENSE-DATA).
