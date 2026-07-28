#!/usr/bin/env python3
"""Assemble the preregistration timestamp chain for the Hall census paper.

Run from the repository root:

    python verification/timestamps.py

Collects, and checks the ordering of:
  1. the local git record for verification/preregistration_dec10.md
  2. the local git record for the decade-10 data artifacts
  3. the GitHub push timestamps for those commits (third-party attested)
  4. the Zenodo version deposits and their creation dates

Only (3) and (4) are evidence to anyone but the author: committer dates are
settable via GIT_COMMITTER_DATE and prove nothing on their own.  Dependency
free; needs network access for steps 3 and 4.
"""
import json
import subprocess
import sys
import urllib.error
import urllib.request

REPO = "michaelmross/fermat-hall-census"
CONCEPT = "21517860"
PREREG = "verification/preregistration_dec10.md"
DATA_PATHS = ["data/hall/w1", "data/hall/w2", "data/hall/w3",
              "data/hall/w4", "data/hall/w5", "data/hall/w6",
              "data/hall/state_dec10.json", "data/hall/state_full.json"]


def git(*args):
    try:
        out = subprocess.run(["git"] + list(args), capture_output=True, text=True)
        return out.stdout.strip()
    except FileNotFoundError:
        return ""


def get_json(url, accept="application/json"):
    req = urllib.request.Request(url, headers={"Accept": accept,
                                               "User-Agent": "hall-census-timestamps"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def section(t):
    print("\n" + t)
    print("-" * len(t))


print("Preregistration timestamp chain for", REPO)

section("1. Local git record: preregistration")
prereg = git("log", "--follow", "--date=iso-strict",
             "--format=%H|%ad|%cd|%s", "--", PREREG)
if not prereg:
    print(f"  no git history found for {PREREG}")
prereg_shas = []
for line in prereg.splitlines():
    h, ad, cd, subj = line.split("|", 3)
    prereg_shas.append(h)
    print(f"  {h}")
    print(f"    authored  {ad}")
    print(f"    committed {cd}")
    print(f"    subject   {subj}")

section("2. Local git record: decade-10 data first added")
data_first = git("log", "--diff-filter=A", "--date=iso-strict",
                 "--format=%H|%cd|%s", "--", *DATA_PATHS)
data_shas = []
for line in data_first.splitlines():
    h, cd, subj = line.split("|", 2)
    data_shas.append(h)
    print(f"  {h}  {cd}  {subj}")
if not data_shas:
    print("  none found -- check the paths in DATA_PATHS against your tree")

section("3. GitHub push timestamps (third-party attested)")
pushed = {}
try:
    events = get_json(f"https://api.github.com/repos/{REPO}/events?per_page=100",
                      accept="application/vnd.github+json")
    for e in events:
        if e.get("type") != "PushEvent":
            continue
        when = e["created_at"]
        for c in e["payload"].get("commits", []):
            pushed.setdefault(c["sha"], when)
        head = e["payload"].get("head")
        if head:
            pushed.setdefault(head, when)
    print(f"  {len(pushed)} commits seen in the last 100 push events")
    for label, shas in (("preregistration", prereg_shas), ("decade-10 data", data_shas)):
        for h in shas:
            w = pushed.get(h)
            print(f"    {label:16s} {h[:12]}  pushed {w if w else 'NOT IN RECENT EVENTS'}")
except urllib.error.HTTPError as ex:
    print(f"  GitHub API error {ex.code} -- events are retained ~90 days only")
except Exception as ex:
    print(f"  GitHub API unavailable: {ex}")

section("4. Zenodo deposits")
try:
    vers = get_json(f"https://zenodo.org/api/records/{CONCEPT}/versions?size=100")
    hits = vers.get("hits", {}).get("hits", [])
    print(f"  concept {CONCEPT}: {len(hits)} version(s)")
    for h in sorted(hits, key=lambda d: d.get("created", "")):
        print(f"    {str(h.get('id')):12s} created {str(h.get('created'))[:19]}"
              f"  {h.get('metadata', {}).get('version', '')}"
              f"  {h.get('metadata', {}).get('title', '')[:50]}")
except urllib.error.HTTPError as ex:
    print(f"  Zenodo returned {ex.code} for concept {CONCEPT}"
          f" -- 404 means restricted or unpublished; 406 means a bad Accept header")
except Exception as ex:
    print(f"  Zenodo unavailable: {ex}")

section("Ordering check")
print("  The claim the paper can make depends on which of these holds:")
print("   (a) prereg pushed  <  data pushed        -> publicly timestamped, prospective")
print("   (b) a Zenodo version created after the")
print("       prereg commit and before the scan    -> immutably notarized (strongest)")
print("  Committer dates alone establish neither.")
