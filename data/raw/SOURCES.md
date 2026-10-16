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

## Validation against authors' processed data (data/reference_authors/, from Zenodo 4027588)

Authors' matrices cover 185 anterior neurons: B (Cook) 601 links / 7 self-loops, B (Varshney) 301 / 2, C 5,557 contacting pairs, X is 185 x 18.
All 594 off-diagonal GJ pairs fall inside C (none outside). 43 of 185 neurons (23%) have no innexin in X.

- **Varshney: exact match.** `ConnOrdered_040903.mat`, matrix `Ag_t_ordered`, symmetrized and subset to the 185 neurons gives 301 links, 2 self-loops, 0 mismatches. Use `Ag_t_ordered` (not `A_init_t_ordered`, which is the chemical network).
- **Cook: close but not identical.** SI 5 "hermaphrodite gap jn symmetric" (binarized, 185-neuron subset) gives 611 links / 7 self-loops vs 601 / 7. 78 edges differ (44 only in the download, 34 only in the authors' file). Likely cause: the file is the July 2020 correction, newer than the version the authors used (preprint May 2020). Not yet confirmed. Need to check whether the SI 5 weights or a threshold explain it, and try the original 2019 supplement.
- **Contact C and expression X:** taken from the authors' files for now. Note `inx-5` and `inx-11` columns of X are identical (columns 3 and 11), matching the identical rows in SI Fig. S3.
