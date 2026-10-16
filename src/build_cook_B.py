"""Build the 185-neuron Cook 2019 GJ matrix B from SI 5 (binarized, masked by contact matrix C)."""
import numpy as np, pandas as pd
names = open("data/reference_authors/NeurNamesContact.csv").read().split()
C = np.loadtxt("data/reference_authors/ContactSubgraphMatrix.csv", delimiter=",", dtype=int)
df = pd.read_excel("data/raw/connectome/cook2019_SI5_adjacency_matrices_original.xlsx",
                   sheet_name="herm gap jn symmetric", header=None)
M = pd.DataFrame(df.iloc[3:, 3:].apply(pd.to_numeric, errors="coerce").fillna(0).values,
                 index=df.iloc[3:, 2].tolist(), columns=df.iloc[2, 3:].tolist())
B = (M.loc[names, names].values > 0).astype(int) * C
np.savetxt("data/processed/B_cook_rebuilt.csv", B, fmt="%d", delimiter=",")
print("links", int(np.triu(B).sum()), "self", int(np.trace(B)))
