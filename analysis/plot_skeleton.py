#!/usr/bin/env python3
"""Draw Figure 1 of the coprimality note: the Mordell skeleton.

Panel (a): the sixth-power-free curve classes by Mordell-Weil rank, with the
solution-bearing classes ringed and the 6 | m layer (k0 = +-1) marked.
Panel (b): share of classes against share of the E(1000) diagnostic mass,
by rank, in the same palette.

Inputs (both produced by the pipeline):
    skeleton_ranks.txt        from analysis/skeleton_ranks.gp
    coprimality_results.json  from `coprimality_pipeline.py all`

Usage:
    python3 plot_skeleton.py --out skeleton_ranks_note.pdf
"""
import argparse
import json
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.patches import Patch     # noqa: E402
from matplotlib.ticker import PercentFormatter  # noqa: E402

FS = 8.5
COLORS = {0: "#a5a9ad", 1: "#4878cf", 2: "#e08214", 3: "#c0392b"}

# Solution-bearing classes, with the label anchors placed in the empty
# corridors between rank bands; leader lines do the pointing.
SOLUTIONS = {
    2:     (r"$2^7$",    (-0.75, 0.62)),
    8:     (r"$2^9$",    (0.95, 0.55)),
    17:    (r"$17^7$",   (0.05, 1.55)),
    65:    (r"$65^7$",   (1.85, 1.50)),
    113:   (r"$113^7$",  (3.15, 3.32)),
    1849:  (r"$43^8$",   (3.40, 2.45)),
    -1089: (r"$33^8$",   (-3.70, 0.58)),
}


def load_ranks(path):
    ranks = {}
    for line in open(path):
        parts = line.strip().split("|")
        if len(parts) == 4 and parts[1] != "ERR":
            ranks[int(parts[0])] = int(parts[1])
    return ranks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranks", default="skeleton_ranks.txt")
    ap.add_argument("--results", default="coprimality_results.json")
    ap.add_argument("--out", default="skeleton_ranks_note.pdf")
    ap.add_argument("--seed", type=int, default=7, help="jitter seed")
    a = ap.parse_args()

    ranks = load_ranks(a.ranks)
    res = json.load(open(a.results))
    rk = res["ranks"]
    n_class = {int(k): v for k, v in rk["by_rank_classes"].items()}
    mass = {int(k): v for k, v in rk["by_rank_mass"].items()}
    tot_mass = sum(mass.values())
    tot_class = sum(n_class.values())

    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                         "font.size": FS})
    rng = np.random.default_rng(a.seed)
    xs, ys, cs = [], [], []
    ypos = {}
    for k, r in ranks.items():
        x = math.copysign(math.log10(abs(k)) if abs(k) > 1 else 0.0, k)
        y = r + rng.uniform(-0.15, 0.15)
        xs.append(x); ys.append(y); cs.append(COLORS[r]); ypos[k] = y

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 2.6),
                                   gridspec_kw={"width_ratios": [1.75, 1]})

    # ---- panel (a) ----
    ax1.scatter(xs, ys, c=cs, s=8, alpha=0.8, linewidths=0, zorder=2)
    lead = dict(arrowstyle="-", lw=0.55, color="0.35", shrinkA=1, shrinkB=5)
    for k, (label, tp) in SOLUTIONS.items():
        if k not in ypos:
            continue
        x = math.copysign(math.log10(abs(k)), k)
        ax1.scatter([x], [ypos[k]], marker="o", s=52, facecolors="none",
                    edgecolors="black", linewidths=0.85, zorder=4)
        ax1.annotate(label, xy=(x, ypos[k]), xytext=tp, textcoords="data",
                     fontsize=7.5, ha="center", va="center", zorder=5,
                     arrowprops=lead)
    for k, dx in ((1, 0.03), (-1, -0.03)):
        if k in ypos:
            ax1.scatter([dx], [ypos[k]], marker="D", s=26, facecolors="none",
                        edgecolors="black", linewidths=0.85, zorder=4)
    if -1 in ypos:
        ax1.annotate(r"$k_0=\pm1$", xy=(0, ypos[-1] - 0.10), xytext=(0.0, -0.35),
                     textcoords="data", fontsize=7.5, ha="center", va="center",
                     zorder=5, arrowprops=dict(arrowstyle="-", lw=0.55,
                                               color="0.35", shrinkA=2, shrinkB=3))
    ax1.set_xlabel(r"$\mathrm{sgn}(k_0)\cdot\log_{10}|k_0|$", fontsize=FS + 1.0)
    ax1.set_ylabel("rank (jittered)", fontsize=FS)
    ax1.set_yticks([0, 1, 2, 3])
    ax1.set_ylim(-0.58, 3.62)
    ax1.set_xlim(-7.9, 8.4)
    ax1.tick_params(labelsize=FS)
    for yl in (0.5, 1.5, 2.5):
        ax1.axhline(yl, color="0.90", lw=0.5, zorder=1)
    ax1.text(0.015, 0.965, "(a)", transform=ax1.transAxes, fontsize=FS + 0.5,
             va="top", weight="bold")

    # ---- panel (b) ----
    rs = sorted(n_class)
    w = 0.38
    fmt = lambda v: f"{v:.1f}" if v < 1 else f"{v:.0f}"
    for i, r in enumerate(rs):
        c_pct = 100 * n_class[r] / tot_class
        m_pct = 100 * mass[r] / tot_mass
        ax2.bar(i - w / 2, c_pct, w, color=COLORS[r], alpha=0.42)
        ax2.bar(i + w / 2, m_pct, w, color=COLORS[r])
        ax2.text(i - w / 2, c_pct + 1.3, fmt(c_pct), ha="center", fontsize=7.5)
        ax2.text(i + w / 2, m_pct + 1.3, fmt(m_pct), ha="center", fontsize=7.5)
    ax2.set_xticks(range(len(rs)))
    ax2.set_xticklabels([str(r) for r in rs])
    ax2.set_xlabel("rank", fontsize=FS)
    ax2.set_ylim(0, 62)
    ax2.tick_params(labelsize=FS)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax2.legend(handles=[Patch(facecolor="0.45", alpha=0.42, label="classes"),
                        Patch(facecolor="0.45", label=r"$E(1000)$ mass")],
               fontsize=6.5, frameon=False, loc="upper right",
               handlelength=1.1, handletextpad=0.45, borderaxespad=0.3,
               labelspacing=0.3)
    ax2.text(0.035, 0.965, "(b)", transform=ax2.transAxes, fontsize=FS + 0.5,
             va="top", weight="bold")

    plt.tight_layout(pad=0.4)
    plt.savefig(a.out)
    if a.out.endswith(".pdf"):
        plt.savefig(a.out[:-4] + ".png", dpi=200)
    print(f"wrote {a.out}")
    print(f"  classes by rank: {dict(sorted(n_class.items()))}")
    print(f"  mass shares:     "
          f"{ {r: round(100 * mass[r] / tot_mass, 1) for r in rs} }")


if __name__ == "__main__":
    main()
