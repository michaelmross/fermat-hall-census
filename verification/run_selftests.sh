#!/usr/bin/env bash
# Tier-1 verification. All scripts run from the repository root. ~1 minute.
set -e
cd "$(dirname "$0")/.."

# Interpreter: Linux CI has python3; Windows ships a stub named python3 that
# exists on PATH but only advertises the Microsoft Store, so test that it
# actually executes rather than merely resolving.
if python3 -c '' >/dev/null 2>&1; then
    PY=python3
elif python -c '' >/dev/null 2>&1; then
    PY=python
else
    echo "No working Python interpreter found (tried python3, python)." >&2
    exit 1
fi
echo "interpreter: $PY ($($PY --version 2>&1))"

echo "== fc23m_scan selftest (7 known Fermat-Catalan solutions) =="
$PY scanners/fc23m_scan.py --selftest
echo "== beal33m_scan selftest (brute-force cross-check + canaries) =="
$PY scanners/beal33m_scan.py --selftest
echo "== hall_census selftest (Catalan point, classical Hall point) =="
$PY scanners/hall_census.py --selftest

echo "== fc23m artifact audit (hits file must contain all 7 known solutions) =="
$PY verification/audit_hits.py data/fc23m/hits.jsonl

echo "== fc23m coverage rectangle (478 anchors x [1,1e9], both phases) =="
$PY verification/coverage_manifest_fc23m.py \
    --run data/fc23m/run1/ledger.jsonl:0 \
    --run data/fc23m/gapfill/ledger.jsonl:1e14 \
    --s-max 1e16 --x-max 1e9

echo "== fc23m closing orientation vs independent enumeration =="
$PY verification/verify_closing_census.py --hits data/fc23m/hits.jsonl

echo "== fc23m lift closure (Lemma 1, seeds, Catalan tower) =="
$PY analysis/lift_closure.py data/fc23m/hits.jsonl

echo "== beal33m deep-skeleton regularity =="
# Python's zipfile rather than unzip(1): Git Bash on Windows has no unzip.
$PY -c "import zipfile; zipfile.ZipFile('data/beal33m/hits_1e30.zip').extractall('data/beal33m')"
$PY verification/verify_deep_skeletons.py

echo "== ALL TIER-1 SELFTESTS PASS =="