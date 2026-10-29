"""Build the 185-neuron Cook GJ matrix B from the WormWiring hermaphrodite edge list (binarized, symmetrized).

This version matches the authors' GapJunctContact.csv exactly (601 links, 7 self-loops).
"""
import numpy as np, pandas as pd
names = open("data/reference_authors/NeurNamesContact.csv").read().split()
idx = {n: i for i, n in enumerate(names)}
df = pd.read_csv("data/raw/connectome/herm_full_edgelist.csv", skipinitialspace=True)
for c in ["Source", "Target", "Type"]:
    df[c] = df[c].astype(str).str.strip()
B = np.zeros((len(names), len(names)), dtype=int)
for s, t in df.loc[df["Type"] == "electrical", ["Source", "Target"]].itertuples(index=False):
    if s in idx and t in idx:
        B[idx[s], idx[t]] = B[idx[t], idx[s]] = 1
np.savetxt("data/processed/B_cook_wormwiring.csv", B, fmt="%d", delimiter=",")
A = np.loadtxt("data/reference_authors/GapJunctContact.csv", delimiter=",", dtype=int)
print("links", int(np.triu(B).sum()), "self", int(np.trace(B)), "differences from authors' B", int(np.triu(B != A).sum()))
