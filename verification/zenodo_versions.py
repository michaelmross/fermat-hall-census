#!/usr/bin/env python3
"""List the Zenodo deposits for the repository archive, with diagnostics.

    python3 verification/zenodo_versions.py [recid]

Fetches the record, then follows the versions link the API returns rather than
constructing one. On an HTTP error it prints the response body, which is where
Zenodo explains what it objected to.
"""
import json
import sys
import urllib.error
import urllib.request

RECID = sys.argv[1] if len(sys.argv) > 1 else "21517860"


def get(url):
    req = urllib.request.Request(
        url, headers={"Accept": "application/json",
                      "User-Agent": "hall-census-zenodo/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r)
    except urllib.error.HTTPError as ex:
        body = ex.read().decode("utf-8", "replace")[:600]
        print(f"  HTTP {ex.code} for {url}")
        print(f"  body: {body}")
        return None
    except Exception as ex:
        print(f"  request failed for {url}: {ex}")
        return None


def show(rec, indent="  "):
    md = rec.get("metadata", {}) or {}
    print(f"{indent}id={rec.get('id')}  conceptrecid={rec.get('conceptrecid')}"
          f"  created={str(rec.get('created'))[:19]}")
    print(f"{indent}  version={md.get('version','-')}  doi={rec.get('doi', md.get('doi','-'))}")
    print(f"{indent}  title={str(md.get('title','-'))[:70]}")


print(f"Zenodo record {RECID}")
rec = get(f"https://zenodo.org/api/records/{RECID}")
if rec is None:
    print("\nRecord fetch failed. If this is a 404 the deposit may be restricted or")
    print("unpublished; if 410 it was removed. Check in a browser at")
    print(f"  https://zenodo.org/records/{RECID}")
    sys.exit(1)

show(rec)

links = rec.get("links", {}) or {}
vurl = links.get("versions")
concept = rec.get("conceptrecid")

print("\nVersions")
data = get(vurl) if vurl else None
if data is None and concept:
    print("  versions link absent or failed; falling back to a concept query")
    data = get("https://zenodo.org/api/records"
               f"?q=conceptrecid:{concept}&all_versions=true&size=50&sort=oldest")

if data is None:
    print("  could not enumerate versions; use the Versions panel at")
    print(f"    https://zenodo.org/records/{RECID}")
else:
    hits = data.get("hits", {}).get("hits", []) if isinstance(data, dict) else []
    if not hits:
        print("  no versions returned")
    for h in sorted(hits, key=lambda d: str(d.get("created", ""))):
        show(h, indent="  ")

print("\nWhat to take from this")
print("  - the version whose created date is closest after the decade-10 scan is")
print("    the deposit to cite for reproducibility in the code-availability section")
print("  - a version created between 2026-07-23T20:43Z (preregistration pushed)")
print("    and the start of the scan would be immutable third-party custody of the")
print("    predictions, which would let section 7 drop its hedge")
