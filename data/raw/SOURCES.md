# Raw data manifest (downloaded 2026-10-16)

| File | Source | Used for |
|------|--------|----------|
| connectome/ConnOrdered_040903.mat | https://github.com/lrvarshney/elegans (Varshney 2011) | GJ connectome B, robustness run |
| connectome/NeuronTypeOrdered_040903.mat | same | neuron class labels |
| connectome/cook2019_SI5_adjacency_matrices_original.xlsx | https://github.com/openworm/ConnectomeToolbox (Cook 2019 SI 5, original) | GJ connectome B, main run (`src/build_cook_B.py`) |
| connectome/cook2019_SI5_adjacency_matrices_corrected_2020-07.xlsx | same (Cook 2019 SI 5, corrected July 2020) | check only; same gap-junction result as the original |

All files above are fetched and checksum-verified by `src/download_data.sh`.
| expression/Altun2009_DevDyn_innexin_map.pdf | https://www.wormatlas.org/papers/Altun%202009.pdf | innexin expression X (figures/tables, needs manual curation) |

## Still missing
- contact/: Brittin et al. 2018 contactome (anterior 185 neurons). No direct download found. Fallback: authors' `ContactSubgraphMatrix.csv` and `NeurNamesContact.csv` from the Zenodo package.
- expression/: curated binary X (Hobert lab, ref 5). Fallback: authors' `INXExpressionJustContact.csv`.

## Validation against authors' processed data (data/reference_authors/, from Zenodo 4027588)

Authors' matrices cover 185 anterior neurons: B (Cook) 601 links / 7 self-loops, B (Varshney) 301 / 2, C 5,557 contacting pairs, X is 185 x 18.
All 594 off-diagonal GJ pairs fall inside C (none outside). 43 of 185 neurons (23%) have no innexin in X.

- **Varshney: exact match.** `ConnOrdered_040903.mat`, matrix `Ag_t_ordered`, symmetrized and subset to the 185 neurons gives 301 links, 2 self-loops, 0 mismatches. Use `Ag_t_ordered` (not `A_init_t_ordered`, which is the chemical network).
- **Cook: close but not identical.** SI 5 "hermaphrodite gap jn symmetric" (binarized, 185-neuron subset) gives 611 links / 7 self-loops vs 601 / 7. 78 edges differ (44 only in the download, 34 only in the authors' file). Investigated (see below); the cause is still open.
- **Contact C and expression X:** taken from the authors' files for now. Note `inx-5` and `inx-11` columns of X are identical (columns 3 and 11), matching the identical rows in SI Fig. S3.

## Cook mismatch investigation (611 vs 601 links; 78 differing edges)

Ruled out:
- **Weight threshold.** Binarizing at weight >= 1 is the best fit; >= 2 gives 385 links and 249 mismatches.
- **July 2020 correction.** `cook2019_SI5_adjacency_matrices_original.xlsx` (the pre-correction SI 5, also from the openworm ConnectomeToolbox repo, 4.4 MB) gives the identical 611 links and the same 78 differences.
- **Symmetric vs asymmetric sheet.** Both give the same result.

Partly explained:
- 9 of the 44 Cook-only edges lie outside the contact matrix C. The authors' B has no GJs outside C, so masking B with C removes these. That leaves about 69 differing edges.

Still unexplained:
- 35 Cook-only edges inside C (mostly weight 1-2) and 34 authors-only edges that are absent from SI 5 entirely. 14 of the 44 Cook-only edges appear in Varshney; only 1 of the 34 authors-only edges does.
- Differences concentrate on SIB* neurons (SIBVL 9, SIBDL 8, SIBDR 7, SIBVR 7), then PVNR, AVBR, RIVL, ADEL.
- Agreement is Jaccard 0.878. The authors probably used a Cook-lab version of the data (e.g. an earlier or WormWiring-hosted release) that differs from SI 5 for these neurons. Next step: check the WormWiring Cook adjacency files and, if needed, ask the authors.
