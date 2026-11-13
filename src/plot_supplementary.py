"""Explanatory figures for supplementary/: the input matrices, the z-scores of O, the rule network, Cook vs Varshney,
a small Kronecker product example, and the ridge path (rule weights as alpha changes).

Usage (from the repository root): python src/plot_supplementary.py
Needs the SCM runs in results/authors_inputs (Cook) and results/varshney (alpha 0.12).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from plot_si_figures import GENES, INK, MUTED, GRID, NEG, style, table_heatmap
from scm import load

REF = "data/reference_authors"
OUT = "supplementary/figures"
ACCENT, OTHER = NEG, "#a8b1bb"
BINARY = ListedColormap(["#ffffff", INK])
# Hand-placed labels for points that sit close together in cook_vs_varshney.png
LABEL_OFFSET = {"inx-13/inx-13": (-7, 3), "unc-9/unc-9": (-6, 8), "inx-19/inx-19": (6, -2),
                "inx-10/inx-2": (-7, 4), "inx-18/inx-6": (6, 6), "inx-12/inx-8": (6, -9)}


def kronecker_example():
    """K = X ⊗ X for 2 neurons and 2 genes: one row per neuron pair, one column per gene pair."""
    X = np.array([[1, 0], [1, 1]])
    K = np.kron(X, X)
    pairs = ["(1,1)", "(1,2)", "(2,1)", "(2,2)"]
    genes = ["(g1,g1)", "(g1,g2)", "(g2,g1)", "(g2,g2)"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), gridspec_kw={"width_ratios": [1, 1.6]})
    for ax, M, rows, cols, title in [
        (axes[0], X, ["neuron 1", "neuron 2"], ["g1", "g2"], "X: which genes each neuron expresses"),
        (axes[1], K, [f"pair {p}" for p in pairs], genes, "K = X ⊗ X: entry = X[i,a] · X[j,b]")]:
        ax.imshow(M, cmap=ListedColormap(["#ffffff", "#c6dbef"]), vmin=0, vmax=1)
        for (r, c), v in np.ndenumerate(M):
            ax.text(c, r, str(v), ha="center", va="center", fontsize=12, color=INK)
        ax.set_xticks(range(M.shape[1]), cols, fontsize=9)
        ax.set_yticks(range(M.shape[0]), rows, fontsize=9)
        ax.xaxis.tick_top()
        ax.tick_params(length=0, labelcolor=INK)
        ax.set_xticks(np.arange(-0.5, M.shape[1]), minor=True)
        ax.set_yticks(np.arange(-0.5, M.shape[0]), minor=True)
        ax.grid(which="minor", color=GRID, linewidth=1)
        ax.tick_params(which="minor", length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(title, fontsize=10, color=INK, loc="left", pad=28)
    axes[1].set_xlabel("columns: gene pair (a, b), the rule O[a,b] it switches on", fontsize=9, color=MUTED)
    fig.tight_layout()
    fig.savefig(f"{OUT}/kronecker_example.png", dpi=150)
    plt.close(fig)


def ridge_path(C, B, X):
    """Rule weights as the ridge penalty alpha changes, for a few rules (Cook data)."""
    keep = C.flatten() == 1
    kp, bp = np.kron(X, X)[keep], B.flatten()[keep]
    kk, kb = kp.T @ kp, kp.T @ bp
    alphas = np.logspace(-4, 1.5, 60)
    O = np.array([np.linalg.solve(kk + a * np.eye(324), kb) for a in alphas])
    rules = [("inx-12", "inx-12"), ("inx-9", "inx-9"), ("inx-10", "inx-2"), ("inx-10", "inx-11"),
             ("inx-12", "inx-9"), ("unc-9", "unc-9")]
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"]  # Okabe-Ito, colorblind-safe
    fig, ax = plt.subplots(figsize=(8, 5))
    for (a, b), col in zip(rules, colors):
        w = O[:, GENES.index(a) * 18 + GENES.index(b)]
        label = f"{a}/{b}" + (" (= inx-10/inx-5)" if b == "inx-11" else "")
        ax.plot(alphas, w, color=col, linewidth=2, label=label)
    ax.axvline(0.215, color=MUTED, linestyle="--", linewidth=1)
    ax.annotate("α = 0.215\n(chosen)", (0.215, ax.get_ylim()[1]), xytext=(4, -24), textcoords="offset points",
                fontsize=9, color=INK)
    ax.axhline(0, color=MUTED, linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlim(alphas[0], alphas[-1])
    ax.set_xlabel("ridge penalty α (log scale)", color=INK)
    ax.set_ylabel("rule weight O[a,b]", color=INK)
    ax.set_title("Ridge path: weights shrink toward 0 as α grows (Cook data)", fontsize=11, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower left")
    style(ax)
    fig.tight_layout()
    fig.savefig(f"{OUT}/ridge_path.png", dpi=150)
    plt.close(fig)


def matrices(C, B, X):
    """The three inputs of B = X O X^T on the 185 anterior neurons (Cook data)."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1, 1, 0.45]})
    panels = [
        (C, "C: contact (185 × 185)", f"{int(np.triu(C).sum()):,} touching pairs"),
        (B, "B: gap junctions (185 × 185)", f"{int(np.triu(B).sum())} links, all inside C"),
        (X, "X: expression\n(185 × 18)", f"{int((X.sum(1) == 0).sum())} neurons express none"),
    ]
    for ax, (M, title, sub) in zip(axes, panels):
        ax.imshow(M, cmap=BINARY, aspect="auto", interpolation="nearest")
        ax.set_title(f"{title}\n{sub}", fontsize=10, color=INK, loc="left")
        ax.tick_params(length=0, labelsize=7, labelcolor=MUTED)
        for s in ax.spines.values():
            s.set_color(GRID)
    axes[0].set_ylabel("neuron", color=INK)
    for ax in axes[:2]:
        ax.set_xlabel("neuron", color=INK)
    axes[2].set_xticks(range(18), GENES, rotation=90, fontsize=6)
    axes[2].set_xlabel("innexin", color=INK)
    fig.suptitle("The model's inputs: B ≈ X O Xᵀ, fitted only on touching pairs (C = 1)",
                 fontsize=12, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/model_matrices.png", dpi=150)
    plt.close(fig)


def o_zscores(Z):
    """The output: z-score of every rule, as in the bottom table of SI Fig. S3. Bold = significant (z >= 2)."""
    fig, ax = plt.subplots(figsize=(9, 8))
    table_heatmap(ax, Z, ".1f", "", nan_text="-")
    for k, t in enumerate(ax.texts):  # one label per cell, row by row
        t.set_fontsize(8)
        if Z.flat[k] >= 2:  # unrounded: a few cells print as 2.0 but are 1.96-1.99
            t.set_fontweight("bold")
    ax.set_title("z-scores of the rule matrix O (Cook data, α = 0.215)\n"
                 "red = connect more than chance, blue = less; bold = significant (z ≥ 2); "
                 "- = undefined (inx-5/inx-11 duplicate)", fontsize=10, color=INK, loc="left", pad=48)
    fig.tight_layout()
    fig.savefig(f"{OUT}/o_zscores.png", dpi=150)
    plt.close(fig)


def rule_network(Zc, Zv):
    """Genes on a circle; a line for each rule with z >= 2 in Cook, a ring for same-gene rules."""
    ang = np.pi / 2 - 2 * np.pi * np.arange(18) / 18
    pos = np.c_[np.cos(ang), np.sin(ang)]
    fig, ax = plt.subplots(figsize=(7, 7))
    n_rules = 0
    for i in range(18):
        for j in range(i, 18):
            if not Zc[i, j] >= 2:  # NaN (no variance) is not significant
                continue
            n_rules += 1
            shared = Zv[i, j] >= 1.95  # inx-8/inx-12 is 1.96 in Varshney, printed 2.0 in the SI
            col, lw = (ACCENT, 3) if shared else (OTHER, 1.8)
            if i == j:
                ax.add_patch(plt.Circle(pos[i], 0.085, fill=False, color=col, linewidth=lw, zorder=2))
            else:
                ax.plot(*pos[[i, j]].T, color=col, linewidth=lw, zorder=1, solid_capstyle="round")
    ax.scatter(*pos.T, s=260, color="white", edgecolor=MUTED, linewidth=1, zorder=3)
    for k, g in enumerate(GENES):
        ax.annotate(g, pos[k] * 1.17, ha="center", va="center", fontsize=10, color=INK)
    ax.legend(handles=[Line2D([], [], color=ACCENT, lw=3, label="significant in Cook and Varshney (12)"),
                       Line2D([], [], color=OTHER, lw=1.8, label="significant in Cook only (7)"),
                       Line2D([], [], color=MUTED, lw=0, marker="o", markersize=12, markerfacecolor="white",
                              markeredgewidth=2, label="ring = same-gene (homotypic) rule")],
              loc="lower center", bbox_to_anchor=(0.5, -0.12), frameon=False, labelcolor=INK, fontsize=9)
    ax.set_title(f"The {n_rules} significant rules (z ≥ 2, Cook data)", fontsize=12, color=INK, loc="left")
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(f"{OUT}/rule_network.png", dpi=150)
    plt.close(fig)


def cook_vs_varshney(Zc, Zv):
    """z-score of each of the 171 gene pairs in the two connectomes."""
    iu = np.triu_indices(18)
    zc, zv = Zc[iu], Zv[iu]
    ok = np.isfinite(zc) & np.isfinite(zv)
    shared = ok & (zc >= 2) & (zv >= 1.95)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.add_patch(plt.Rectangle((2, 2), 6, 6, color=GRID, alpha=0.45, zorder=0, linewidth=0))
    ax.scatter(zc[ok & ~shared], zv[ok & ~shared], s=36, color=OTHER, edgecolor="white", linewidth=1, zorder=2)
    ax.scatter(zc[shared], zv[shared], s=60, color=ACCENT, edgecolor="white", linewidth=1.5, zorder=3)
    for (i, j), x, y in zip(np.array(iu).T[shared], zc[shared], zv[shared]):
        name = f"{GENES[i]}/{GENES[j]}"
        ax.annotate(name, (x, y), xytext=LABEL_OFFSET.get(name, (6, 3)), textcoords="offset points",
                    fontsize=7.5, color=INK, ha="right" if LABEL_OFFSET.get(name, (6, 3))[0] < 0 else "left")
    r = np.corrcoef(zc[ok], zv[ok])[0, 1]
    ax.set_xlabel("z-score, Cook connectome", color=INK)
    ax.set_ylabel("z-score, Varshney connectome", color=INK)
    ax.set_title(f"Same rules from two connectomes (r = {r:.2f}); shaded: z ≥ 2 in both", fontsize=11,
                 color=INK, loc="left")
    ax.set_xlim(-5.5, 8)
    ax.set_ylim(-5.5, 8)
    style(ax)
    fig.tight_layout()
    fig.savefig(f"{OUT}/cook_vs_varshney.png", dpi=150)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    C, X = load(f"{REF}/ContactSubgraphMatrix.csv", int), load(f"{REF}/INXExpressionJustContact.csv", int)
    B = load(f"{REF}/GapJunctContact.csv", int)

    Zc, Zv = load("results/authors_inputs/zscores_0.2150.csv"), load("results/varshney/zscores_0.1200.csv")
    matrices(C, B, X)
    kronecker_example()
    ridge_path(C, B, X)
    o_zscores(Zc)
    rule_network(Zc, Zv)
    cook_vs_varshney(Zc, Zv)


if __name__ == "__main__":
    main()
