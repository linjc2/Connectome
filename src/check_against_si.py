"""Compare our rule matrices and z-scores with the tables printed in the SI (Figs. S3, S4, S6).

Usage (from the repository root): python src/check_against_si.py [path/to/pnas.2009093117.sapp.pdf]
The SI PDF is not in the repository (copyright); download it from the paper's PNAS page.
Run src/plot_si_figures.py first: it writes the Moore-Penrose solutions this compares with Fig. S6.
"""
import sys
import numpy as np
import pymupdf
from plot_si_figures import GENES

PDF = sys.argv[1] if len(sys.argv) > 1 else "pnas.2009093117.sapp.pdf"


def page_tables(page):
    """The 18 x 18 tables on one SI page, in order. '-' (undefined z-score) becomes NaN."""
    toks = page.get_text().split()
    tables = []
    for start in [i for i in range(len(toks) - 18) if toks[i:i + 18] == GENES]:
        pos, M = start + 18, np.full((18, 18), np.nan)
        for i, g in enumerate(GENES):
            assert toks[pos] == g, f"expected row {g}, found {toks[pos]}"
            pos += 1
            for j in range(18):
                M[i, j] = np.nan if toks[pos] == "-" else float(toks[pos])
                pos += 1
        tables.append(M)
    return tables


def compare(label, ours, si, decimals):
    ok = np.isfinite(si) & np.isfinite(ours)
    diff = np.abs(ours - si)[ok].max()
    tol = 0.5 * 10 ** -decimals + 1e-9
    print(f"{label:40s} max |ours - SI| = {diff:.4f}  (SI rounds to {10 ** -decimals:g}): "
          f"{'MATCH' if diff <= tol else 'DIFFERENT'}")
    return diff <= tol


def main():
    doc = pymupdf.open(PDF)
    s3, s4, s6 = page_tables(doc[3]), page_tables(doc[4]), page_tables(doc[6])
    load = lambda p: np.loadtxt(p, delimiter=",")
    checks = [
        ("Fig. S3 Cook rule weights (alpha 0.215)", load("results/authors_inputs/rule_matrix_0.2150.csv"), s3[0], 2),
        ("Fig. S3 Cook z-scores", load("results/authors_inputs/zscores_0.2150.csv"), s3[1], 1),
        ("Fig. S4 Varshney rule weights (alpha 0.12)", load("results/varshney/rule_matrix_0.1200.csv"), s4[0], 2),
        ("Fig. S4 Varshney z-scores", load("results/varshney/zscores_0.1200.csv"), s4[1], 1),
        ("Fig. S6 Cook Moore-Penrose", load("results/figures/rule_matrix_MP_cook.csv"), s6[0], 2),
        ("Fig. S6 Varshney Moore-Penrose", load("results/figures/rule_matrix_MP_varshney.csv"), s6[1], 2),
    ]
    results = [compare(*c) for c in checks]
    print(f"{sum(results)} of {len(results)} tables match the SI")


if __name__ == "__main__":
    main()
