"""Compare the SCM rules across four connectomes: Cook, Varshney and the two Witvliet et al. 2021 adults.

Usage (from the repository root): python src/compare_connectomes.py [gcv|fixed]
  gcv   (default) each connectome at its own cross-validated alpha (0.215, 0.12, 15.849, 6.3096)
  fixed all four at alpha = 0.215
Needs the SCM runs in results/ (see README). Writes results/connectome_comparison_<mode>.csv and
results/figures/rules_across_connectomes_<mode>.png, and prints a summary.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from plot_si_figures import GENES, PALE_RDBU, INK
from scm import load

REF = "data/reference_authors"
MODE = sys.argv[1] if len(sys.argv) > 1 else "gcv"
ALPHA = {"gcv": ["0.2150", "0.1200", "15.8490", "6.3096"], "fixed": ["0.2150"] * 4}[MODE]
RUNS = {  # name: (results folder, B)
    "Cook": ("results/authors_inputs", f"{REF}/GapJunctContact.csv"),
    "Varshney": ("results/varshney", "data/processed/B_varshney.csv"),
    "Witvliet adult 1": ("results/witvliet_adult1", "data/processed/B_witvliet_adult1.csv"),
    "Witvliet adult 2": ("results/witvliet_adult2", "data/processed/B_witvliet_adult2.csv"),
}
IU = np.triu_indices(18)
PAIRS = [f"{GENES[i]}/{GENES[j]}" for i, j in zip(*IU)]


def direct_evidence(B, C, X):
    """For each gene pair (a, b): gap junctions among touching neuron pairs where the rule applies."""
    n = np.zeros(len(PAIRS), dtype=int)
    for k, (a, b) in enumerate(zip(*IU)):
        applies = np.outer(X[:, a], X[:, b]) | np.outer(X[:, b], X[:, a])
        n[k] = int(np.triu(B * C * applies).sum())
    return n


def main():
    C, X = load(f"{REF}/ContactSubgraphMatrix.csv", int), load(f"{REF}/INXExpressionJustContact.csv", int)
    names = list(RUNS)
    zs = np.array([load(f"{d}/zscores_{a}.csv")[IU] for (d, _), a in zip(RUNS.values(), ALPHA)])
    gj = np.array([direct_evidence(load(b, int), C, X) for _, b in RUNS.values()])
    with open(f"results/connectome_comparison_{MODE}.csv", "w") as f:
        f.write("gene_pair," + ",".join(f"z {n},GJs {n}" for n in names) + "\n")
        for k, p in enumerate(PAIRS):
            f.write(p + "," + ",".join(("" if np.isnan(zs[d, k]) else f"{zs[d, k]:.3f}") + f",{gj[d, k]}"
                                       for d in range(len(names))) + "\n")

    pos = np.nan_to_num(zs) >= 2
    real = pos & (gj > 0)  # significant and backed by at least one gap junction where the rule applies
    neg = np.nan_to_num(zs) <= -2
    cook_rules, both_rules = pos[0], pos[0] & (np.nan_to_num(zs[1]) >= 1.95)  # the paper's 19 and 12
    print(f"alpha: {dict(zip(names, ALPHA))}")
    print("Positive rules (z >= 2):", {n: int(pos[k].sum()) for k, n in enumerate(names)})
    print("  of which with no gap junction where the rule applies:",
          {n: [PAIRS[i] for i in np.where(pos[k] & (gj[k] == 0))[0]] for k, n in enumerate(names)})
    print("Negative rules (z <= -2):", {n: int(neg[k].sum()) for k, n in enumerate(names)})
    for k in (2, 3):
        print(f"{names[k]}: {int((real[k] & cook_rules).sum())} of Cook's 19 rules, "
              f"{int((real[k] & both_rules).sum())} of the 12 shared Cook/Varshney rules")
    print("Witvliet adult 1 vs adult 2 rules shared:", int((real[2] & real[3]).sum()))
    print("Rules in all four:", [PAIRS[i] for i in np.where(real.all(0))[0]])
    print("Rules in 3 of 4:", [PAIRS[i] for i in np.where(real.sum(0) == 3)[0]])
    print("Negative rules in 3 or 4:", [PAIRS[i] for i in np.where(neg.sum(0) >= 3)[0]])
    ok = np.isfinite(zs).all(0)
    print("Correlation of z-scores (Cook, Varshney, Witvliet 1, Witvliet 2):")
    print(np.round(np.corrcoef(zs[:, ok]), 2))

    # Figure: every rule significant (z >= 2) in at least one connectome, sorted by how many agree
    sel = np.where(pos.any(0))[0]
    order = sel[np.lexsort((-np.nanmean(zs[:, sel], 0), -real[:, sel].sum(0)))]
    M = zs[:, order].T
    fig, ax = plt.subplots(figsize=(6.4, 0.3 * len(order) + 1.9))
    ax.imshow(M, cmap=PALE_RDBU, norm=TwoSlopeNorm(0, -6, 6), aspect="auto")
    for (r, c), v in np.ndenumerate(M):
        if np.isfinite(v):
            flag = "*" if v >= 2 and gj[c, order[r]] == 0 else ""
            ax.text(c, r, f"{v:.1f}{flag}", ha="center", va="center", fontsize=7.5, color=INK,
                    fontweight="bold" if v >= 2 and not flag else "normal")
    ax.set_xticks(range(len(names)), [f"{n.replace('Witvliet ', 'Witvliet' + chr(10))}\nα = {float(a):g}"
                                      for n, a in zip(names, ALPHA)], fontsize=8)
    ax.set_yticks(range(len(order)), [f"{PAIRS[i]}  ({int(real[:, i].sum())}/4)" for i in order], fontsize=7.5)
    ax.xaxis.tick_top()
    ax.tick_params(length=0, labelcolor=INK)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("z-scores of every rule significant in at least one connectome\n"
                 "bold = z ≥ 2; * = z ≥ 2 but no gap junction where the rule applies;\n"
                 "(n/4) = connectomes where it is significant with gap junctions behind it",
                 fontsize=9, color=INK, loc="left", pad=46)
    fig.tight_layout()
    fig.savefig(f"results/figures/rules_across_connectomes_{MODE}.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
