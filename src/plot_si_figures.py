"""Recreate SI Figs. S2-S6 of Kovacs et al. 2020 from our data.

Usage (from the repository root): python src/plot_si_figures.py
Needs the SCM runs in results/authors_inputs (Cook, alpha 0.215) and results/varshney (alpha 0.12).
Writes PNGs and the computed numbers (alpha scans, Moore-Penrose solutions) to results/figures/.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from scm import load, ridge_operator

GENES = ["inx-1", "inx-10", "inx-11", "inx-12", "inx-13", "inx-17", "inx-18", "inx-19", "inx-2",
         "inx-3", "inx-5", "inx-6", "inx-7", "inx-8", "inx-9", "che-7", "unc-7", "unc-9"]
REF = "data/reference_authors"
OUT = "results/figures"
C = load(f"{REF}/ContactSubgraphMatrix.csv", int)
X = load(f"{REF}/INXExpressionJustContact.csv", int)
DATA = {"Cook": (load(f"{REF}/GapJunctContact.csv", int), 0.215, "results/authors_inputs"),
        "Varshney": (load("data/processed/B_varshney.csv", int), 0.12, "results/varshney")}
C_BG = 8.55e-6  # constant background subtracted on the y-axis of SI Fig. S2
INK, MUTED, GRID = "#1f2328", "#59636e", "#d8dee4"
POS, NEG = "#d6604d", "#4393c3"
# Pale blue-white-red, like the SI tables, so dark text stays readable on every cell
_rdbu = plt.get_cmap("RdBu_r")(np.linspace(0, 1, 256))
PALE_RDBU = ListedColormap(0.45 * _rdbu + 0.55 * np.array([1, 1, 1, 1]))


def contact_design(B):
    """Rows of kron(X, X) and entries of B for contacting pairs only (as in SCM_alpha.py)."""
    k = np.kron(X, X)
    keep = C.flatten() == 1
    return k[keep], B.flatten()[keep]


def alpha_scan(B, alphas):
    """Generalized cross-validation score r^2 / tau^2 for each alpha (SCM_alpha.py)."""
    kp, bp = contact_design(B)
    kk = kp.T @ kp
    scores = []
    for a in alphas:
        ka = np.linalg.solve(kk + a * np.eye(kk.shape[0]), kp.T)
        tau = kp.shape[0] - np.trace(ka @ kp)  # trace(I - kp ka) without the n x n matrix
        r = np.linalg.norm(bp - kp @ (ka @ bp))
        scores.append(r * r / tau / tau)
    return np.array(scores)


def moore_penrose(B):
    """alpha -> 0 solution with the pseudo-inverse (rule_matrix_MP.csv in the authors' code)."""
    kp, bp = contact_design(B)
    return (np.linalg.pinv(kp) @ bp).reshape(X.shape[1], X.shape[1])


def style(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def table_heatmap(ax, M, fmt, title, nan_text=""):
    """Annotated gene-by-gene table, red for positive and blue for negative, white at 0."""
    lim = np.nanmax(np.abs(M))
    ax.imshow(np.ma.masked_invalid(M), cmap=PALE_RDBU, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for i in range(18):
        for j in range(18):
            txt = nan_text if np.isnan(M[i, j]) else format(M[i, j], fmt)
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.5, color=INK)
    ax.set_xticks(range(18), GENES, rotation=90, fontsize=7)
    ax.set_yticks(range(18), GENES, fontsize=7)
    ax.xaxis.tick_top()
    ax.tick_params(length=0, labelcolor=INK)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, fontsize=9, color=INK, loc="left", pad=40)


def fig_tables(top, bottom, titles, fmts, path, suptitle):
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 15))
    for ax, M, t, f in zip(axes, [top, bottom], titles, fmts):
        table_heatmap(ax, M, f, t, nan_text="-")
    fig.suptitle(suptitle, fontsize=10, color=INK, x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    runs = {}
    for name, (B, a, d) in DATA.items():
        tag = f"{a:.4f}"
        runs[name] = (load(f"{d}/rule_matrix_{tag}.csv"), load(f"{d}/zscores_{tag}.csv"), a)

    # Fig. S2: alpha selection (Cook as in the paper; Varshney added to check its alpha = 0.12)
    alphas = np.unique(np.round(np.concatenate([[0.001, 0.01], np.arange(0.05, 1.0001, 0.05),
                                                np.arange(0.10, 0.30001, 0.005)]), 4))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, name in zip(axes, DATA):
        s = alpha_scan(DATA[name][0], alphas)
        np.savetxt(f"{OUT}/alpha_scan_{name.lower()}.csv", np.c_[alphas, s], delimiter=",",
                   header="alpha,r2_over_tau2", comments="")
        best = alphas[np.argmin(s)]
        y = s - C_BG if name == "Cook" else s - s.min()
        ax.plot(alphas, y, color=POS, linewidth=2, marker="+", markersize=6)
        ax.axvline(best, color=MUTED, linewidth=1, linestyle="--")
        ax.annotate(f"minimum at α = {best:.3f}", (best, y.max() * 0.9), xytext=(8, 0),
                    textcoords="offset points", color=INK, fontsize=9)
        ax.set_xlabel("α", color=INK)
        ax.set_ylabel("r²/τ² − c" if name == "Cook" else "r²/τ² − min", color=INK)
        ax.set_title(f"{name}" + (" (SI Fig. S2, c = 8.55e-6)" if name == "Cook" else " (not in SI)"),
                     fontsize=10, color=INK, loc="left")
        style(ax)
        print(f"S2 {name}: GCV minimum at alpha = {best:.3f}, r2/tau2 = {s.min():.6e}")
    fig.tight_layout()
    fig.savefig(f"{OUT}/figS2_alpha_selection.png", dpi=150)
    plt.close(fig)

    # Figs. S3 and S4: rule weights and z-scores
    for fig_id, name in [("S3", "Cook"), ("S4", "Varshney")]:
        O, Z, a = runs[name]
        fig_tables(O, Z, [f"Rule weights O (α = {a})", "z-scores (- = no variance in the randomized ensemble)"],
                   [".2f", ".1f"], f"{OUT}/fig{fig_id}_{name.lower()}_rules.png",
                   f"Fig. {fig_id}. Solution of the SCM for {name} et al. data (our run)")

    # Fig. S5: cumulative eigenvalue distribution of the Cook O
    O, Z, _ = runs["Cook"]
    sig = np.abs(np.nan_to_num(Z)) >= 2
    variants = [("full", O, "+", POS), ("sig. (|z| ≥ 2)", np.where(sig, O, 0), "x", NEG),
                ("sig. pos. (z ≥ 2)", np.where(np.nan_to_num(Z) >= 2, O, 0), "s", MUTED)]
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    for label, M, mk, col in variants:
        ev = np.sort(np.linalg.eigvalsh((M + M.T) / 2))
        ax.plot(ev, np.arange(1, 19), linestyle="none", marker=mk, markersize=8, color=col,
                markerfacecolor="none" if mk == "s" else col, label=f"{label}: {int((ev < 0).sum())} negative")
        print(f"S5 {label}: {int((ev < 0).sum())} negative eigenvalues, range {ev[0]:.2f} to {ev[-1]:.2f}")
    ax.set_xlabel("Eigenvalue", color=INK)
    ax.set_ylabel("Cumulative number", color=INK)
    ax.set_title("Fig. S5. Eigenvalues of O (Cook, our run)", fontsize=10, color=INK, loc="left")
    ax.legend(frameon=False, loc="upper left", labelcolor=INK)
    style(ax)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figS5_eigenvalues.png", dpi=150)
    plt.close(fig)

    # Fig. S6: Moore-Penrose solutions
    mp = {}
    for name in DATA:
        mp[name] = moore_penrose(DATA[name][0])
        np.savetxt(f"{OUT}/rule_matrix_MP_{name.lower()}.csv", mp[name], delimiter=",")
    fig_tables(mp["Cook"], mp["Varshney"], ["Cook et al.", "Varshney et al."], [".2f", ".2f"],
               f"{OUT}/figS6_moore_penrose.png",
               "Fig. S6. Solution of the SCM with the Moore-Penrose pseudoinverse (α → 0), our run")


if __name__ == "__main__":
    main()
