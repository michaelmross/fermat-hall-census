#!/usr/bin/env bash
# Tier-1 verification: every instrument re-derives its known solutions. ~1 minute.
set -e
cd "$(dirname "$0")/.."
echo "== fc23m_scan selftest (7 known Fermat-Catalan solutions) =="
python3 scanners/fc23m_scan.py --selftest
echo "== beal33m_scan selftest (brute-force cross-check + canaries) =="
python3 scanners/beal33m_scan.py --selftest
echo "== hall_census selftest (Catalan point, classical Hall point) =="
python3 scanners/hall_census.py --selftest
echo "== ALL TIER-1 SELFTESTS PASS =="
