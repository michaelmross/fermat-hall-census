# Data

Committed evidence files, one directory per instrument:

- `fc23m/` -- hits.jsonl + ledger.jsonl for the {2,3,m} census, with
  `run1/` and `gapfill/` subdirectories and boundary_sweep_196.json
  (coverage a^m <= 1e16, x <= 1e9)
- `beal33m/` -- hits_1e30.zip (compressed) + ledger.jsonl for the exact
  {3,3,m} census
- `hall/` -- state.json, state_dec10.json, state_full.json,
  merged_state.json, hall_hits.jsonl, power_hits.jsonl, ledger.jsonl,
  plus per-worker directories w1/ .. w6/ (x <= 1e11)
- `coprimality/` -- skeleton.json, coprimality_results.json,
  skeleton_classes.txt, skeleton_ranks.txt,
  skeleton_ranks_figure_numbers.txt, the figure PDF/PNG, and
  skeleton_figure.py

## Integrity

Generate/verify checksums from the repository root:

    sha256sum $(git ls-files data | grep -v '^data/CHECKSUMS$') > data/CHECKSUMS
    sha256sum -c data/CHECKSUMS

Run these in Git Bash; its sha256sum writes the `hash *path` form that the
CI manifest check expects. The manifest covers every tracked file under
data/ except CHECKSUMS itself, so regenerate it after editing any file in
this directory -- including this README.

CHECKSUMS is committed; the papers cite the release DOI, and the release
freezes these hashes.

These hashes are valid only while data/** -text stays in .gitattributes.

## Decompressed working copy

The {3,3,m} census is committed compressed. Unzip before tier-2 or any
analysis that reads the .jsonl:

    unzip -o data/beal33m/hits_1e30.zip -d data/beal33m/

The decompressed file is gitignored and therefore absent from CHECKSUMS.
Its SHA-256 is recorded here so the chain covers both the committed
artifact and the file the scripts actually read:

    SHA256(data/beal33m/hits_1e30.jsonl) = PLACEHOLDER_SHA

Data license: CC-BY 4.0 (see LICENSE-DATA).
