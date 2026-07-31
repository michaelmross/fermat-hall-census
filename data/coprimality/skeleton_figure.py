#!/usr/bin/env python3
"""Regenerate Figure 1 of the coprimality note: the Mordell skeleton.

  M. M. Ross, "Coprimality Density and the Proper-Solution Deficit in the
  Generalized Fermat Family {2,3,m}"

Reads the pipeline's own outputs so the figure cannot drift from the text:

  skeleton.json               written by coprimality_pipeline.py skeleton
  coprimality_results.json    written by coprimality_pipeline.py all
                              (supplies the E(1000) mass per anchor-orientation
                               and the rank of every class)

Writes skeleton_ranks_note.pdf (vector, for the paper) and
skeleton_ranks_note.png (raster, for the repository README), plus
skeleton_ranks_figure_numbers.txt recording every number that appears in the
caption, so the caption can be checked against the figure mechanically.
A bare --out name is anchored to this script's own directory, mirroring the
input resolution below, so a rerun from any working directory replaces the
same canonical artifacts rather than writing fresh copies into the working
directory and leaving a stale PDF beside the paper.

Panel (a) positions one point per sixth-power-free class k0 by
sgn(k0) * log10|k0|, banded by Mordell-Weil rank with deterministic vertical
jitter. The horizontal axis is display spread only: the anchor-derived k0 are
highly non-uniform, so any apparent rank-versus-|k0| trend is confounded by
the selection and is not offered as a covariate.

Panel (b) contrasts, per rank, the share of classes against the share of the
E(1000) diagnostic mass.

Usage:
  python skeleton_figure.py
  python skeleton_figure.py --results coprimality_results.json --skeleton skeleton.json
"""
import argparse
import json
import math
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# the seven known coprime solutions and the classes that own them
KNOWN = {2: "$2^7$", 8: "$2^9$", 17: "$17^7$", 65: "$65^7$",
         113: "$113^7$", 1849: "$43^8$", -1089: "$33^8$"}

RANK_COLOR = {0: "#7f7f7f", 1: "#1f77b4", 2: "#d62728", 3: "#2ca02c"}
JITTER_SEED = 20260730


def resolve(path):
    """Look in the working directory, then beside this script, so the same
    command works from the repository root and from analysis/."""
    import os
    if os.path.exists(path):
        return path
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        os.path.basename(path))
    return here if os.path.exists(here) else None


def load(results_path, skeleton_path):
    import os
    r, s = resolve(results_path), resolve(skeleton_path)
    if r is None or s is None:
        missing = [p for p, v in ((results_path, r), (skeleton_path, s))
                   if v is None]
        raise SystemExit(
            "missing input: " + ", ".join(missing) + "\n"
            "Looked in the working directory and beside this script.\n"
            "These are written by the pipeline (about 30 s); run it from "
            "wherever coprimality_pipeline.py lives:\n"
            "    python coprimality_pipeline.py all --ranks skeleton_ranks.txt")
    results_path, skeleton_path = r, s
    res = json.load(open(results_path))
    skel = json.load(open(skeleton_path))
    ranks = {int(k): v for k, v in res["ranks"]["ranks"].items()}
    per_anchor = res["local"]["per_anchor"]

    def mass(s, sign):
        d = per_anchor[str(s)]["mass_1000"]
        return d[str(sign)] if str(sign) in d else d[sign]

    class_mass = defaultdict(float)
    class_orients = defaultdict(int)
    for k, members in skel.items():
        for s, a, m, sign in members:
            class_mass[int(k)] += mass(int(s), sign)
            class_orients[int(k)] += 1
    return ranks, dict(class_mass), dict(class_orients)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="coprimality_results.json")
    ap.add_argument("--skeleton", default="skeleton.json")
    ap.add_argument("--out", default="skeleton_ranks_note")
    ns = ap.parse_args()

    # Outputs mirror the input resolution of resolve(): a bare name is
    # anchored to this script's directory, so reruns from the repository
    # root and from analysis/ replace the same files. An explicit path
    # (anything containing a separator) is honored as given.
    if not os.path.dirname(ns.out):
        ns.out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              ns.out)
    out_dir = os.path.dirname(ns.out)
    numbers_path = os.path.join(
        out_dir, os.path.basename(ns.out).replace("_note", "")
        + "_figure_numbers.txt")

    ranks, class_mass, class_orients = load(ns.results, ns.skeleton)
    total_mass = sum(class_mass.values())
    n_class = len(ranks)

    by_rank_classes = defaultdict(int)
    by_rank_mass = defaultdict(float)
    for k, r in ranks.items():
        by_rank_classes[r] += 1
        by_rank_mass[r] += class_mass.get(k, 0.0)
    rank_list = sorted(by_rank_classes)

    rng = np.random.default_rng(JITTER_SEED)
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.0, 4.1), gridspec_kw={"width_ratios": [2.35, 1]})

    # ---- panel (a): the skeleton ------------------------------------------
    for r in rank_list:
        ks = [k for k in ranks if ranks[k] == r and abs(k) != 1]
        x = [math.copysign(math.log10(abs(k)), k) for k in ks]
        y = r + rng.uniform(-0.26, 0.26, size=len(ks))
        axA.scatter(x, y, s=13, c=RANK_COLOR[r], alpha=0.72,
                    linewidths=0, zorder=2)

    # k0 = +-1: the 6 | m layer, drawn as open diamonds at their true position
    for k in (1, -1):
        if k in ranks:
            axA.scatter([0.0], [ranks[k]], s=95, facecolors="none",
                        edgecolors="black", marker="D", linewidths=1.2,
                        zorder=4)
    axA.annotate("$k_0 = \\pm 1$\n($6 \\mid m$ layer)", xy=(0.0, 0.0),
                 xytext=(0.55, -0.62), fontsize=7.5, ha="left",
                 arrowprops=dict(arrowstyle="-", lw=0.6, color="0.35"))

    # the seven solution-bearing classes; labels alternate above/below when
    # two stars sit close together on the same rank band
    stars = sorted(((math.copysign(math.log10(abs(k)), k), ranks[k], lab)
                    for k, lab in KNOWN.items() if k in ranks))
    prev_x, prev_y, flip = -99.0, None, False
    for x, y, label in stars:
        axA.scatter([x], [y], s=95, marker="*", c="black", zorder=5)
        if y == prev_y and x - prev_x < 0.45:
            flip = not flip
        else:
            flip = False
        axA.annotate(label, xy=(x, y), xytext=(0, -16 if flip else 9),
                     textcoords="offset points", fontsize=7.5, ha="center",
                     zorder=6)
        prev_x, prev_y = x, y

    axA.axhspan(-0.42, 0.42, color="0.90", zorder=0)
    axA.text(-7.6, 0.0, "rank-0 floor", fontsize=7.5, color="0.35",
             va="center", zorder=1)
    axA.set_yticks(rank_list)
    axA.set_yticklabels([f"rank {r}" for r in rank_list])
    axA.set_ylim(-0.85, max(rank_list) + 0.7)
    axA.set_xlabel("$\\mathrm{sgn}(k_0)\\,\\log_{10}|k_0|$"
                   "   (display spread only, not a covariate)")
    axA.set_title(f"(a) the {n_class} sixth-power-free classes, "
                  "by Mordell\u2013Weil rank", fontsize=9.5, loc="left")
    axA.grid(axis="x", lw=0.3, color="0.85")
    axA.set_axisbelow(True)
    for side in ("top", "right"):
        axA.spines[side].set_visible(False)
    axA.legend(handles=[
        Line2D([], [], marker="*", color="black", linestyle="none",
               markersize=9, label="owns a known proper solution"),
        Line2D([], [], marker="D", color="black", linestyle="none",
               markerfacecolor="none", markersize=7,
               label="$k_0 = \\pm 1$")],
        fontsize=7.5, loc="upper left", frameon=False)

    # ---- panel (b): classes against mass ----------------------------------
    idx = np.arange(len(rank_list))
    w = 0.38
    cshare = [by_rank_classes[r] / n_class for r in rank_list]
    mshare = [by_rank_mass[r] / total_mass for r in rank_list]
    for i, r in enumerate(rank_list):
        axB.bar(i - w / 2, cshare[i], w, color=RANK_COLOR[r], alpha=0.32,
                edgecolor=RANK_COLOR[r], linewidth=0.8)
        axB.bar(i + w / 2, mshare[i], w, color=RANK_COLOR[r])
        axB.text(i - w / 2, cshare[i] + 0.012, f"{cshare[i]:.1%}",
                 ha="center", fontsize=7.5, color="0.25")
        axB.text(i + w / 2, mshare[i] + 0.012, f"{mshare[i]:.1%}",
                 ha="center", fontsize=7.5, color="0.25")
    axB.set_xticks(idx)
    axB.set_xticklabels([f"rank {r}" for r in rank_list])
    axB.set_ylim(0, max(max(cshare), max(mshare)) * 1.22)
    axB.set_ylabel("share")
    axB.set_title("(b) share of classes (pale) against share of\n"
                  "$E(1000)$ diagnostic mass (solid)", fontsize=9.5,
                  loc="left")
    axB.grid(axis="y", lw=0.3, color="0.85")
    axB.set_axisbelow(True)
    for side in ("top", "right"):
        axB.spines[side].set_visible(False)

    fig.tight_layout()
    fig.savefig(f"{ns.out}.pdf")
    fig.savefig(f"{ns.out}.png", dpi=200)

    # ---- the caption's numbers, recorded for checking ---------------------
    lines = [f"classes: {n_class}",
             f"E(1000) total mass: {total_mass:.4f}",
             f"solution-bearing classes plotted: "
             f"{sum(1 for k in KNOWN if k in ranks)} of {len(KNOWN)}",
             "all solution-bearing classes have positive rank: "
             f"{all(ranks[k] > 0 for k in KNOWN if k in ranks)}",
             f"max |k0|: {max(abs(k) for k in ranks)}"]
    for r in rank_list:
        lines.append(f"rank {r}: classes {by_rank_classes[r]} "
                     f"({by_rank_classes[r]/n_class:.1%}), "
                     f"orientations {sum(class_orients.get(k,0) for k in ranks if ranks[k]==r)}, "
                     f"mass {by_rank_mass[r]:.2f} "
                     f"({by_rank_mass[r]/total_mass:.1%})")
    open(numbers_path, "w").write(
        "\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {os.path.abspath(ns.out)}.pdf, "
          f"{os.path.abspath(ns.out)}.png,\n      {numbers_path}")


if __name__ == "__main__":
    main()
