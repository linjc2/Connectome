"""Build the 185-neuron Varshney 2011 GJ matrix B from ConnOrdered_040903.mat (Ag_t_ordered, binarized, symmetrized).

This version matches the authors' CEO_Gap_OnlyContact.csv exactly (301 links, 2 self-loops).
"""
import numpy as np, scipy.io as sio
names = open("data/reference_authors/NeurNamesContact.csv").read().split()
m = sio.loadmat("data/raw/connectome/ConnOrdered_040903.mat")
order = [str(n[0][0]).strip() for n in m["Neuron_ordered"]]
G = (m["Ag_t_ordered"].toarray() > 0).astype(int)
G = np.maximum(G, G.T)
pos = [order.index(n) for n in names]
B = G[np.ix_(pos, pos)]
np.savetxt("data/processed/B_varshney.csv", B, fmt="%d", delimiter=",")
A = np.loadtxt("data/reference_authors/CEO_Gap_OnlyContact.csv", delimiter=",", dtype=int)
print("links", int(np.triu(B).sum()), "self", int(np.trace(B)), "differences from authors' B", int(np.triu(B != A).sum()))
