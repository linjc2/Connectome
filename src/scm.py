"""Spatial Connectome Model (Kovacs et al. 2020), Python 3 port of the authors' SCM_alpha*.py.

Usage: python scm.py C.csv B.csv X.csv alpha outdir
Needs alpha_matrix.out (compiled from alpha_matrix.cpp) in results/bin.
"""
import subprocess, sys, os, tempfile
import numpy as np

def load(p, dtype=float):
    return np.loadtxt(p, delimiter=",", dtype=dtype)

def ridge_operator(C, X, alpha):
    k = np.kron(X, X)
    kp = k[C.flatten() == 1]
    ka = np.linalg.solve(kp.T @ kp + alpha * np.eye(kp.shape[1]), kp.T)
    return kp, ka

def solve(C, B, X, alpha):
    kp, ka = ridge_operator(C, X, alpha)
    bp = B.flatten()[C.flatten() == 1]
    return (ka @ bp).reshape(X.shape[1], X.shape[1])

def randomized_moments(C, B, X, alpha, binpath, workdir):
    """Mean and variance of O over the degree-preserving ensemble restricted to contacts."""
    n = C.shape[0]
    bt, ct = os.path.join(workdir, "B.txt"), os.path.join(workdir, "C.txt")
    np.savetxt(bt, B.astype(int), fmt="%d", delimiter="\t", newline="\t")
    np.savetxt(ct, C.astype(int), fmt="%d", delimiter="\t", newline="\t")
    subprocess.run([binpath, "B.txt", "C.txt", str(n), "100"], cwd=workdir, check=True)
    P = load(os.path.join(workdir, "alpha_exp_B.txt_C.txt_100.csv"))
    kp, ka = ridge_operator(C, X, alpha)
    pp = P.flatten()[C.flatten() == 1]
    mean = (ka @ pp).reshape(X.shape[1], X.shape[1])
    var = ((ka * ka) @ (pp * (1.0 - pp))).reshape(X.shape[1], X.shape[1])
    return mean, var

if __name__ == "__main__":
    Cp, Bp, Xp, alpha, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
    os.makedirs(out, exist_ok=True)
    C, B, X = load(Cp, int), load(Bp, int), load(Xp, int)
    O = solve(C, B, X, alpha)
    np.savetxt(os.path.join(out, f"rule_matrix_{alpha:.4f}.csv"), O, delimiter=",")
    binpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "bin", "alpha_matrix.out")
    with tempfile.TemporaryDirectory() as wd:
        mean, var = randomized_moments(C, B, X, alpha, binpath, wd)
    z = (O - mean) / np.sqrt(var)
    np.savetxt(os.path.join(out, f"zscores_{alpha:.4f}.csv"), z, delimiter=",")
