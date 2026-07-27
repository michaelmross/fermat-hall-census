#!/usr/bin/env python3
"""Perfect-power channel readout and cross-certification against the companion census.

    python3 verification/power_channel.py

Lists every k(x) that is a perfect m-th power, m >= 7, and checks each against
the companion {2,3,m} census.  The two programs record different orientations
of the same signature, so matching is numeric rather than by record equality:
a Hall event x^3 + a^m = y^2 is sought in the companion by its value set.
"""
import glob, json, os

def load(pat):
    out = []
    for p in sorted(glob.glob(pat)):
        try:
            for line in open(p):
                line = line.strip()
                if line:
                    out.append(json.loads(line))
        except OSError:
            pass
    return out

hits, seen = [], set()
for r in load("data/hall/w*/power_hits.jsonl") + load("data/hall/power_hits.jsonl"):
    key = json.dumps(r, sort_keys=True)
    if key not in seen:
        seen.add(key); hits.append(r)

print(f"perfect-power events (deduplicated): {len(hits)}")
by_decade = {}
for r in hits:
    d = len(str(r["x"])) - 1
    by_decade.setdefault(d, []).append(r)
for d in sorted(by_decade):
    for r in by_decade[d]:
        print(f"  decade {d:2d}  {r['equation']:48s} proper={r.get('proper')}")
top = max(r["x"] for r in hits)
print(f"\n  largest x = {top}; decades above {len(str(top))-1} contributed none")

# --- numeric cross-certification -------------------------------------------
comp, comp_seen, raw = [], set(), 0
for p in ("data/fc23m/hits.jsonl", "data/fc23m/run1/hits.jsonl",
          "data/fc23m/gapfill/hits.jsonl"):
    if os.path.exists(p):
        for line in open(p):
            line = line.strip()
            if not line:
                continue
            raw += 1
            rec = json.loads(line)
            key = json.dumps(rec, sort_keys=True)
            if key not in comp_seen:
                comp_seen.add(key)
                comp.append((p, rec))
print(f"\ncompanion {{2,3,m}} records: {len(comp)} distinct, from {raw} lines")
if raw != len(comp):
    print(f"  ({raw - len(comp)} duplicate lines: the merged hits.jsonl is the"
          f" union of run1 and gapfill, so all three overlap)")

def values(rec):
    """every integer appearing in a record, as a set of strings"""
    out = set()
    def walk(v):
        if isinstance(v, dict):
            for w in v.values(): walk(w)
        elif isinstance(v, (list, tuple)):
            for w in v: walk(w)
        elif isinstance(v, int):
            out.add(str(v))
        elif isinstance(v, str):
            for tok in "".join(c if c.isdigit() else " " for c in v).split():
                out.add(tok)
    walk(rec)
    return out

print("\ncross-certification (each Hall event must appear in the companion census)")
missing = 0
for r in hits:
    want = {str(r["x"]), str(r["y"]), str(r["a"]), str(r["m"])}
    match = [p for p, c in comp if want <= values(c)]
    if match:
        print(f"  FOUND  {r['equation']:48s} in {match[0]}")
    else:
        loose = [p for p, c in comp if {str(r['x']), str(r['y'])} <= values(c)]
        print(f"  ABSENT {r['equation']:48s}"
              + (f" (partial match in {loose[0]})" if loose else ""))
        missing += 1
print(f"\n{len(hits)-missing}/{len(hits)} cross-certified")
if missing:
    print("Inspect any ABSENT row before repeating the cross-certification claim:")
    print("it may be a schema difference rather than a genuine omission.")
