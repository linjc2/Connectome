"""Recreate Fig. 5 of Kovacs et al. 2020 (the network of significant innexin rules) from our Cook and Varshney runs.

Usage (from the repository root): python src/plot_fig5.py
Boxes are innexins; a line or loop is a positive rule with z >= 2 in Cook (alpha 0.215), labeled with its weight O and
drawn thicker for larger weights. Dark blue: also significant in Varshney (alpha 0.12; z >= 1.95, because the SI prints
1.96 as 2.0 and the paper counts it); light blue: Cook only. * marks the five rules the paper lists as supported by
earlier experiments. Box positions follow the paper's layout; which links exist, their colors and weights come from
our data. Writes results/figures/fig5_rule_network.png.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path
from plot_si_figures import GENES, INK, MUTED
from scm import load

DARK, LIGHT = "#1f4fe0", "#33c3ef"
# Box centres, read off the paper's Fig. 5 (y grows downward). unc-7 has no significant rule and is not drawn.
POS = {"inx-9": (435, 77), "inx-17": (617, 253), "inx-13": (812, 255), "inx-8": (1020, 257), "inx-12": (1243, 257),
       "che-7": (408, 479), "inx-2": (413, 646), "inx-10": (545, 756), "inx-3": (658, 648), "inx-5": (410, 860),
       "inx-11": (662, 860), "inx-6": (845, 650), "inx-18": (842, 860), "inx-1": (1038, 860), "inx-7": (1038, 648),
       "unc-9": (1221, 650), "inx-19": (1221, 860)}
LOOP_SIDE = {"inx-9": "left"}  # every other loop sits on top of its box
# Rules the paper reports as found in earlier experiments (main text, refs 46, 50)
KNOWN = {("inx-19", "inx-19"), ("unc-9", "unc-9"), ("inx-10", "inx-11"), ("inx-3", "inx-3"), ("inx-6", "inx-6")}
HALF = 40  # half the box size


def width(w):
    return 2 + 12 * w


def loop(ax, g, w, col, star):
    x, y = POS[g]
    if LOOP_SIDE.get(g) == "left":
        base_a, base_b, tip = (x - HALF, y - 12), (x - HALF, y + 12), (x - HALF - 95, y)
        c1, c2 = (tip[0] - 10, y - 70), (tip[0] - 10, y + 70)
        verts = [base_a, (base_a[0] - 40, base_a[1] - 10), c1, tip, c2, (base_b[0] - 40, base_b[1] + 10), base_b]
        label = (tip[0] + 38, y)
    else:
        base_a, base_b, tip = (x - 12, y - HALF), (x + 12, y - HALF), (x, y - HALF - 95)
        c1, c2 = (x - 70, tip[1] - 10), (x + 70, tip[1] - 10)
        verts = [base_a, (base_a[0] - 10, base_a[1] - 40), c1, tip, c2, (base_b[0] + 10, base_b[1] - 40), base_b]
        label = (x, tip[1] + 38)
    path = Path(verts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(PathPatch(path, fill=False, edgecolor=col, linewidth=width(w), capstyle="round", zorder=1))
    ax.text(*label, f"{w:.2f}", ha="center", va="center", fontsize=11, color=MUTED, zorder=3)
    if star:
        ax.text(x + 48, y - HALF - 82, "*", fontsize=22, ha="center", va="center", color=INK, fontweight="bold")


def link(ax, a, b, w, col, star):
    (x1, y1), (x2, y2) = POS[a], POS[b]
    ax.plot([x1, x2], [y1, y2], color=col, linewidth=width(w), solid_capstyle="butt", zorder=1)
    d = np.array([x2 - x1, y2 - y1], float)
    n = np.array([-d[1], d[0]]) / np.hypot(*d)
    if n[1] > 0:  # put labels above the line, as in the paper
        n = -n
    mx, my = (x1 + x2) / 2 + 26 * n[0], (y1 + y2) / 2 + 26 * n[1]
    ax.text(mx, my, f"{w:.2f}", ha="center", va="center", fontsize=11, color=MUTED, zorder=3)
    if star:
        ax.text(x2 + 26, y2 - HALF - 22, "*", fontsize=22, ha="center", va="center", color=INK, fontweight="bold")


def main():
    O = load("results/authors_inputs/rule_matrix_0.2150.csv")
    Zc = np.nan_to_num(load("results/authors_inputs/zscores_0.2150.csv"))
    Zv = np.nan_to_num(load("results/varshney/zscores_0.1200.csv"))
    fig, ax = plt.subplots(figsize=(11, 6.6))
    n_rules = 0
    for i in range(18):
        for j in range(i, 18):
            if Zc[i, j] < 2:
                continue
            a, b = GENES[i], GENES[j]
            n_rules += 1
            col = DARK if Zv[i, j] >= 1.95 else LIGHT
            star = (a, b) in KNOWN or (b, a) in KNOWN
            if a == b:
                loop(ax, a, O[i, j], col, star)
            else:
                # draw from the box on the left/top so the star lands at the far end, as in the paper
                a, b = sorted((a, b), key=lambda g: (POS[g][0], POS[g][1]))
                link(ax, a, b, O[i, j], col, star)
    for g, (x, y) in POS.items():
        ax.add_patch(FancyBboxPatch((x - HALF, y - HALF), 2 * HALF, 2 * HALF, boxstyle="square,pad=0",
                                    facecolor="white", edgecolor=INK, linewidth=1.8, zorder=2))
        ax.text(x, y, g, ha="center", va="center", fontsize=9.5, color=INK, zorder=3)
    ax.plot([], [], color=DARK, linewidth=5, label="significant in Cook and Varshney")
    ax.plot([], [], color=LIGHT, linewidth=5, label="significant in Cook only")
    ax.plot([], [], " ", label="* found in earlier experiments")
    ax.legend(loc="upper left", bbox_to_anchor=(0.57, 0.56), frameon=False, fontsize=10, labelcolor=INK)
    ax.set_title(f"Fig. 5 recreated: the {n_rules} positive innexin rules with z ≥ 2 (our Cook and Varshney runs)",
                 fontsize=12, color=INK, loc="left")
    ax.set_xlim(230, 1330)
    ax.set_ylim(930, -40)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig("results/figures/fig5_rule_network.png", dpi=150)
    plt.close(fig)
    print(f"{n_rules} rules drawn")


if __name__ == "__main__":
    main()
