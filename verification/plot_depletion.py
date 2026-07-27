#!/usr/bin/env python3
"""Figure: the depletion collapse, its local slope, and the dispersion calibration.

  python3 analysis/plot_depletion.py --out depletion_note.pdf

Panel (a) plots the family-subtracted depletion D against the scaling variable
s = x^(theta-3/2) for every cell eligible under the frozen rule (uniform-model
expectation >= 100), with the two registered decade-10 cells ringed.
Panel (b) plots the within-column local decay exponent, which tests curvature
without assuming the one-variable collapse.
Panel (c) plots the measured Fano factor of family-subtracted counts over
equal-model-mass blocks against the Poisson surrogate's assumed value of 1.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.linewidth": 0.7,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": False,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})
# uncomment for exact LaTeX type matching with the surrounding document
# plt.rcParams["text.usetex"] = True

GREY, BLUE, ORANGE, RED = "#8c8c8c", "#4878b8", "#e08020", "#c03028"

# (log10 s, D%, sigma_D%) -- theta = 0.9 column, decades 4..9 then decade 10
T9 = [(-2.70, 32.88, 7.66), (-3.30, 20.63, 5.02), (-3.90, 16.37, 3.26),
      (-4.50, 11.61, 2.08), (-5.10, 8.87, 1.31), (-5.70, 6.48, 0.83)]
T9_OUT = (-6.30, 3.68, 0.53)
# theta = 0.8 column, decades 5..9 then decade 10 (central value negative)
T8 = [(-3.85, 34.23, 9.29), (-4.55, 16.87, 7.12), (-5.25, 8.54, 5.18),
      (-5.95, 7.00, 3.66), (-6.65, 4.67, 2.72)]
T8_OUT = (-7.35, -0.92, 1.95)

# within-column local exponents, theta = 0.9, plotted at the step midpoint
LOCAL = [(-3.00, 0.337, 0.244), (-3.60, 0.167, 0.228), (-4.20, 0.249, 0.194),
         (-4.80, 0.195, 0.168), (-5.40, 0.227, 0.142), (-6.00, 0.411, 0.139)]
COL_FIT, COL_FIT_SD = 0.221, 0.007

# measured Fano factors over equal-model-mass blocks
FANO9 = [(-3.90, 0.526, 0.258), (-4.50, 0.663, 0.155), (-4.98, 0.887, 0.138)]
FANO10 = [(-3.25, 1.100, 0.114), (-3.75, 1.041, 0.071)]
FANO9_MEAN, FANO9_SD = 0.752, 0.096


def panel_a(ax):
    xs = np.linspace(-7.6, -2.4, 100)
    ax.plot(xs, 100 * 0.864 * 10 ** (0.191 * xs), ls="--", lw=0.9, color=GREY, zorder=1)
    ax.plot(xs, 100 * 1.425 * 10 ** (0.226 * xs), ls="-", lw=0.9, color="#c8c8c8", zorder=1)

    for data, col, mk in ((T9, BLUE, "o"), (T8, ORANGE, "s")):
        a = np.array(data)
        ax.errorbar(a[:, 0], a[:, 1], yerr=a[:, 2], fmt=mk, ms=3.4, color=col,
                    ecolor=col, elinewidth=0.7, capsize=1.6, mew=0, zorder=3)

    lx, d, sd = T9_OUT
    ax.errorbar([lx], [d], yerr=[sd], fmt="o", ms=3.4, color=BLUE, ecolor=BLUE,
                elinewidth=0.7, capsize=1.6, mew=0, zorder=4)
    ax.plot([lx], [d], "o", ms=9, mfc="none", mec="k", mew=0.7, zorder=5)

    lx, d, sd = T8_OUT
    cap = d + sd
    ax.plot([lx], [cap], "_", ms=7, color=ORANGE, mew=1.0, zorder=4)
    ax.annotate("", xy=(lx, 0.62), xytext=(lx, cap),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=ORANGE,
                                mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.plot([lx], [cap], "s", ms=9, mfc="none", mec="k", mew=0.7, zorder=5)

    dec = ("dec." + chr(92) + ",10") if plt.rcParams["text.usetex"] else "dec. 10"
    ax.annotate(dec, xy=(-6.30, 3.68), xytext=(-6.02, 2.15), fontsize=7)
    ax.annotate(dec, xy=(-7.35, cap), xytext=(-7.60, 1.85), fontsize=7)

    ax.set_yscale("log")
    ax.set_xlim(-7.75, -2.35); ax.set_ylim(0.5, 60)
    ax.set_yticks([1, 2, 5, 10, 20, 50])
    ax.set_yticklabels(["1", "2", "5", "10", "20", "50"])
    ax.set_xlabel(r"$\log_{10} s$, \ $s = x^{\theta-3/2}$" if plt.rcParams["text.usetex"]
                  else r"$\log_{10} s$,  $s = x^{\theta-3/2}$")
    ax.set_ylabel(r"non-family depletion $D$ (\%)" if plt.rcParams["text.usetex"]
                  else r"non-family depletion $D$ (%)")
    ax.grid(axis="y", lw=0.4, color="#e4e4e4", zorder=0)
    ax.set_axisbelow(True)

    h = [plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=3.4, label=r"$\theta=0.9$"),
         plt.Line2D([], [], marker="s", ls="", color=ORANGE, ms=3.4, label=r"$\theta=0.8$"),
         plt.Line2D([], [], ls="--", lw=0.9, color=GREY, label=r"registered, $\gamma=0.191$"),
         plt.Line2D([], [], ls="-", lw=0.9, color="#c8c8c8", label=r"11-cell, $\gamma=0.226$")]
    ax.legend(handles=h, loc="lower right", fontsize=6.5, frameon=False,
              handlelength=1.6, borderaxespad=0.6, labelspacing=0.3)
    ax.text(0.025, 0.945, r"$\mathbf{(a)}$", transform=ax.transAxes, fontsize=9)


def panel_b(ax):
    a = np.array(LOCAL)
    ax.axhspan(COL_FIT - COL_FIT_SD, COL_FIT + COL_FIT_SD, color="#e8e8e8", zorder=0)
    ax.axhline(COL_FIT, lw=0.8, color=GREY, zorder=1)
    ax.errorbar(a[:-1, 0], a[:-1, 1], yerr=a[:-1, 2], fmt="o", ms=3.4, color=BLUE,
                ecolor=BLUE, elinewidth=0.7, capsize=1.6, mew=0, zorder=3)
    ax.errorbar(a[-1:, 0], a[-1:, 1], yerr=a[-1:, 2], fmt="o", ms=3.4, color=BLUE,
                ecolor=BLUE, elinewidth=0.7, capsize=1.6, mew=0, zorder=3)
    ax.plot(a[-1, 0], a[-1, 1], "o", ms=9, mfc="none", mec="k", mew=0.7, zorder=5)
    ax.set_xlim(-6.45, -2.55); ax.set_ylim(-0.05, 0.62)
    ax.set_xlabel(r"$\log_{10} s$")
    ax.set_ylabel(r"local exponent $\gamma_{\mathrm{loc}}$")
    ax.grid(axis="y", lw=0.4, color="#e4e4e4", zorder=0)
    ax.set_axisbelow(True)
    ax.text(0.05, 0.945, r"$\mathbf{(b)}$", transform=ax.transAxes, fontsize=9)


def panel_c(ax):
    ax.axhline(1.0, lw=0.9, color="k", ls=":", zorder=1)
    ax.axhspan(FANO9_MEAN - FANO9_SD, FANO9_MEAN + FANO9_SD, color="#e8eef6", zorder=0)
    ax.axhline(FANO9_MEAN, lw=0.8, color=BLUE, zorder=1)
    for data, col, mk in ((FANO9, BLUE, "o"), (FANO10, RED, "^")):
        a = np.array(data)
        ax.errorbar(a[:, 0], a[:, 1], yerr=a[:, 2], fmt=mk, ms=3.4, color=col,
                    ecolor=col, elinewidth=0.7, capsize=1.6, mew=0, zorder=3)
    ax.set_xlim(-5.35, -2.90); ax.set_ylim(0.15, 1.45)
    ax.set_xlabel(r"$\log_{10} s$")
    ax.set_ylabel(r"Fano factor $F$")
    ax.grid(axis="y", lw=0.4, color="#e4e4e4", zorder=0)
    ax.set_axisbelow(True)
    h = [plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=3.4, label=r"$\theta=0.9$"),
         plt.Line2D([], [], marker="^", ls="", color=RED, ms=3.4, label=r"$\theta=1.0$"),
         plt.Line2D([], [], ls=":", lw=0.9, color="k", label="Poisson")]
    ax.legend(handles=h, loc="lower right", fontsize=6.5, frameon=False,
              handlelength=1.6, borderaxespad=0.5, labelspacing=0.3)
    ax.text(0.05, 0.945, r"$\mathbf{(c)}$", transform=ax.transAxes, fontsize=9)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="depletion_note.pdf")
    a = ap.parse_args()

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.55),
                             gridspec_kw=dict(width_ratios=[2.05, 1.0, 1.0], wspace=0.34))
    panel_a(axes[0]); panel_b(axes[1]); panel_c(axes[2])
    fig.subplots_adjust(left=0.075, right=0.985, top=0.965, bottom=0.185)
    fig.savefig(a.out)
    fig.savefig(a.out.replace(".pdf", ".png"), dpi=220)
    print("wrote", a.out)
