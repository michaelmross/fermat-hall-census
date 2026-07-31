#!/usr/bin/env python3
r"""Independent re-certification of the Mordell rank census by mwrank (eclib).
IMPORTANT: Special setup instructions below.
The census determines every rank by 2-descent in PARI/GP (\texttt{ellrank}).
This re-runs the same classes under mwrank, a separate implementation of
2-descent in a separate library, and diffs the two. The critical assertion is
one-sided: every class the census calls rank zero must be certified rank zero
here, since Proposition "rank zero" turns those into a proof that they carry
no proper solution. Agreement at positive rank is reported but is not load
bearing.

Input: the census's own rank table, one class per line,
    k0|rank_lower|rank_upper|torsion_order
(the format written by analysis/skeleton_ranks.gp; the fourth field is the
torsion order from elltors, not a multiplicity). No hand-built CSV is needed;
--pari-ranks also accepts plain "k0,rank" if that is what you have.

Expected outcome on this census: 647 of 670 determined, no disagreements, and
23 returning "0 <= rank <= selmer-rank = 2". Those 23 are not a failure of
either program -- they are the classes where ellrank's third component reports
Sha[2] = (Z/2)^2, so the 2-Selmer rank exceeds the rank by 2 and mwrank's
plain descent cannot attribute the excess to Sha. PARI finds 25 such classes;
the two missing from the 23 are k0 = -22^3 and -46^3, which have rational
2-torsion and so admit mwrank's second descent. See analysis/sha_census.gp.

Output: mwrank_results.json (per class: k0, PARI rank, mwrank rank or bounds,
runtime, raw tail of the mwrank output) and a summary to stdout. Exit status
is nonzero if any class disagrees or any rank-zero class fails to certify.

Setup
-----
Requires the mwrank binary (eclib) on PATH. Python needs nothing beyond the
standard library.

  Debian / Ubuntu / WSL:
      sudo apt update
      sudo apt install -y eclib-tools
      mwrank -V                     # confirm; expect /usr/bin/mwrank

  Windows: there is no native eclib build. Use WSL as above; the repository
  is reachable from WSL under /mnt/c, e.g.
      cd /mnt/c/Users/<you>/.../fermat-hall-census
      python3 verification/mwrank_check.py --pari-ranks analysis/skeleton_ranks.txt

  macOS:  brew install eclib   (or use Sage, below)

  Sage:   Sage ships eclib but not on the system PATH. Run this script from
          inside a Sage shell rather than trying to wrap the binary:
              sage -sh
              python3 verification/mwrank_check.py --pari-ranks ...
          Wrapping as `--mwrank "sage -sh -c mwrank"` does NOT work: the
          value of --mwrank is executed directly, not through a shell.

  Binary elsewhere:  --mwrank /full/path/to/mwrank

Expected outcome and exit status
-------------------------------
On this census: 647 of 670 determined, zero disagreements with PARI, 23
returning bounds, of which 22 are rank-zero classes. That is the correct
result, not a failure, but the script still exits nonzero because unpinned
rank-zero classes are a failed assertion by design (they carry the rank-zero
proposition). A CI wrapper should therefore compare against the expected
counts rather than the exit status alone. Runtime is about 30-45 seconds
single-threaded; --procs N parallelizes.

Usage:
  mwrank_check.py --pari-ranks skeleton_ranks.txt
  mwrank_check.py --pari-ranks skeleton_ranks.txt --only-zero
  mwrank_check.py --pari-ranks skeleton_ranks.txt --procs 4 --timeout 600
"""
__version__ = "2"
__features__ = "reads skeleton_ranks.txt directly; fixed eclib rank parsing"

import argparse, json, re, subprocess, sys, time
from multiprocessing import Pool

RANK_EXACT = re.compile(r"Rank\s*=\s*(\d+)")
# eclib reports an undetermined rank as e.g. "0 <= rank <= selmer-rank = 2"
RANK_BOUNDS = re.compile(r"(\d+)\s*<=\s*rank\s*<=\s*(?:selmer-rank\s*=\s*)?(\d+)",
                         re.I)
RANK_BOUNDS_ALT = re.compile(r"Rank\s*(?:is\s*)?(?:at least|>=)\s*(\d+).*?"
                             r"(?:at most|<=)\s*(\d+)", re.S | re.I)


def parse_ranks(path):
    """census rank table: 'k0|lower|upper|mult' or 'k0,rank'"""
    out = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "|" in line:
            p = line.split("|")
            if len(p) >= 3 and p[1] not in ("ERR", ""):
                lo, hi = int(p[1]), int(p[2])
                out[int(p[0])] = lo if lo == hi else None
        else:
            k, r = line.split(",")
            out[int(k)] = int(r)
    return out


def run_one(args):
    k0, timeout, exe = args
    t0 = time.time()
    try:
        r = subprocess.run([exe, "-q", "-v", "0"],
                           input=f"0 0 0 0 {k0}\n",
                           capture_output=True, text=True, timeout=timeout)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return {"k0": k0, "status": "timeout", "seconds": timeout}
    m = RANK_EXACT.search(out)
    if m:
        return {"k0": k0, "status": "ok", "rank": int(m.group(1)),
                "lower": int(m.group(1)), "upper": int(m.group(1)),
                "seconds": round(time.time() - t0, 2), "raw": out.strip()[-400:]}
    b = RANK_BOUNDS.search(out) or RANK_BOUNDS_ALT.search(out)
    if b:
        return {"k0": k0, "status": "bounds", "rank": None,
                "lower": int(b.group(1)), "upper": int(b.group(2)),
                "seconds": round(time.time() - t0, 2), "raw": out.strip()[-400:]}
    return {"k0": k0, "status": "unparsed", "rank": None,
            "seconds": round(time.time() - t0, 2), "raw": out.strip()[-400:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pari-ranks", required=True)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--only-zero", action="store_true",
                    help="check only the classes the census calls rank zero")
    ap.add_argument("--out", default="mwrank_results.json")
    ap.add_argument("--mwrank", default="mwrank",
                    help="path to the mwrank binary if it is not on PATH")
    ns = ap.parse_args()

    print(f"mwrank_check v{__version__} ({__features__})", flush=True)

    import shutil
    exe = shutil.which(ns.mwrank)
    if exe is None:
        raise SystemExit(
            f"'{ns.mwrank}' is not on PATH.\n"
            "  Debian/Ubuntu:  sudo apt install eclib-tools\n"
            "  check it landed: dpkg -L eclib-tools | grep bin\n"
            "  Sage:           --mwrank \"sage -sh -c mwrank\" will not work; "
            "instead run this script inside `sage -sh`\n"
            "  elsewhere:      --mwrank /full/path/to/mwrank\n"
            "There is no native Windows build of eclib; use WSL.")
    print(f"using {exe}", flush=True)

    pari = parse_ranks(ns.pari_ranks)
    ks = sorted(pari, key=abs)
    if ns.only_zero:
        ks = [k for k in ks if pari[k] == 0]
    print(f"classes to check: {len(ks)} of {len(pari)}"
          f"{' (rank zero only)' if ns.only_zero else ''}")

    jobs = [(k, ns.timeout, ns.mwrank) for k in ks]
    results = []
    if ns.procs > 1:
        with Pool(ns.procs) as pool:
            for i, res in enumerate(pool.imap(run_one, jobs)):
                results.append(res)
                if (i + 1) % 50 == 0:
                    print(f"  ...{i+1}/{len(ks)}", file=sys.stderr, flush=True)
    else:
        for i, j in enumerate(jobs):
            results.append(run_one(j))
            if (i + 1) % 50 == 0:
                print(f"  ...{i+1}/{len(ks)}", file=sys.stderr, flush=True)

    disagree, unresolved, zero_fail = [], [], []
    for r in results:
        k, pr = r["k0"], pari[r["k0"]]
        if r["status"] == "ok":
            if pr is not None and r["rank"] != pr:
                disagree.append((k, pr, r["rank"]))
        elif r["status"] == "bounds":
            lo, hi = r["lower"], r["upper"]
            if pr is not None and not (lo <= pr <= hi):
                disagree.append((k, pr, f"[{lo},{hi}]"))
            else:
                unresolved.append((k, f"bounds [{lo},{hi}], PARI {pr}"))
                if pr == 0:
                    zero_fail.append(k)
        else:
            unresolved.append((k, r["status"]))
            if pr == 0:
                zero_fail.append(k)

    json.dump({"pari_source": ns.pari_ranks, "classes": len(ks),
               "results": results}, open(ns.out, "w"), indent=1)

    n_ok = sum(1 for r in results if r["status"] == "ok")
    slow = sorted(results, key=lambda r: -r.get("seconds", 0))[:3]
    print(f"certified by mwrank: {n_ok}/{len(ks)}")
    print(f"disagreements with PARI: {len(disagree)}")
    for k, a, b in disagree[:20]:
        print(f"  DISAGREE k0={k}: PARI {a}, mwrank {b}")
    print(f"unresolved (timeout/bounds/unparsed): {len(unresolved)}")
    for k, st in unresolved[:20]:
        print(f"  {st}: k0={k}")
    if zero_fail:
        print(f"rank-zero classes mwrank does not pin down: {len(zero_fail)} "
              f"-- these carry the rank-zero proposition, so the independent "
              f"certification is partial exactly where it is load bearing")
        print("  " + ", ".join(str(k) for k in zero_fail))
    print(f"slowest: " + ", ".join(f"k0={r['k0']} {r.get('seconds',0)}s"
                                   for r in slow))
    print(f"archived {ns.out}")
    sys.exit(1 if (disagree or zero_fail) else 0)


if __name__ == "__main__":
    main()