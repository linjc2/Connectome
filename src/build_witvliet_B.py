"""Build 185-neuron GJ matrices B from the two Witvliet et al. 2021 adult connectomes (binarized, symmetrized, masked by C).

Witvliet et al. reconstructed new animals, independent of the White et al. 1986 series behind Cook and Varshney.
"""
import numpy as np, pandas as pd
names = open("data/reference_authors/NeurNamesContact.csv").read().split()
idx = {n: i for i, n in enumerate(names)}
C = np.loadtxt("data/reference_authors/ContactSubgraphMatrix.csv", delimiter=",", dtype=int)
for k in (1, 2):
    df = pd.read_excel(f"data/raw/connectome/witvliet2021_adult{k}.xlsx")
    present = set(df["pre"].astype(str)) | set(df["post"].astype(str))
    W = np.zeros((len(names), len(names)), dtype=int)
    for s, t in df.loc[df["type"] == "electrical", ["pre", "post"]].astype(str).itertuples(index=False):
        if s in idx and t in idx:
            W[idx[s], idx[t]] = W[idx[t], idx[s]] = 1
    B = W * C
    np.savetxt(f"data/processed/B_witvliet_adult{k}.csv", B, fmt="%d", delimiter=",")
    print(f"adult {k}: links {int(np.triu(B).sum())}, self {int(np.trace(B))}, "
          f"removed outside C {int(np.triu(W).sum() - np.triu(B).sum())}, "
          f"missing neurons {sorted(n for n in names if n not in present)}")
