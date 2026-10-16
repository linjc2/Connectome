# Raw data manifest (downloaded 2026-10-16)

| File | Source | Used for |
|------|--------|----------|
| connectome/ConnOrdered_040903.mat | https://github.com/lrvarshney/elegans (Varshney 2011) | GJ connectome B, robustness run |
| connectome/NeuronTypeOrdered_040903.mat | same | neuron class labels |
| connectome/cook2019_SI5_adjacency_matrices_corrected_2020-07.xlsx | https://github.com/openworm/ConnectomeToolbox (Cook 2019 SI 5, corrected July 2020) | GJ connectome B, main run |
| expression/Altun2009_DevDyn_innexin_map.pdf | https://www.wormatlas.org/papers/Altun%202009.pdf | innexin expression X (figures/tables, needs manual curation) |

## Still missing
- contact/: Brittin et al. 2018 contactome (anterior 185 neurons). No direct download found. Fallback: authors' `ContactSubgraphMatrix.csv` and `NeurNamesContact.csv` from the Zenodo package.
- expression/: curated binary X (Hobert lab, ref 5). Fallback: authors' `INXExpressionJustContact.csv`.
