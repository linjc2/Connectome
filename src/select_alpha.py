"""Choose the ridge alpha for a gap-junction matrix B by generalized cross-validation (minimum of r^2/tau^2).

Usage (from the repository root): python src/select_alpha.py data/processed/B_witvliet_adult1.csv [out.csv]
Uses the authors' C and X, as in SI Fig. S2.
"""
import sys
import numpy as np
from scm import load
from plot_si_figures import alpha_scan

B = load(sys.argv[1], int)
alphas = np.unique(np.round(np.concatenate([[0.001, 0.01], np.arange(0.05, 1.0001, 0.05),
                                            np.arange(0.02, 0.40001, 0.005), np.logspace(0, 3, 31)]), 4))
scores = alpha_scan(B, alphas)
if len(sys.argv) > 2:
    np.savetxt(sys.argv[2], np.c_[alphas, scores], delimiter=",", header="alpha,r2_over_tau2", comments="")
print(f"{sys.argv[1]}: GCV minimum at alpha = {alphas[np.argmin(scores)]:.3f}")
