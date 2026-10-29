# Data collection outline: reproducing Kovács, Barabási & Barabási (PNAS 2020)

Goal: rebuild the three inputs to the Spatial Connectome Model (SCM) from primary
sources, then run it and compare against the authors' published outputs.

Model: `B ≈ X O Xᵀ`, restricted to neuron pairs in physical contact (`C`).
- **B**: gap-junction (GJ) connectome (edge list over contacting pairs)
- **X**: neurons × innexin expression matrix (binary, 18 genes)
- **C**: contact matrix (which neuron pairs touch)
- **O**: innexin-by-innexin rule matrix (ridge regression, α = 0.215)

## Inputs to collect (bottom layer)

| # | Input | Paper's source | Where to get it | Output file |
|---|-------|----------------|-----------------|-------------|
| 1 | GJ connectome B (main) | Cook et al. 2019, *Nature* 571:63 (ref 33). 1,051 undirected links, 11 self-loops | Supplementary data of Cook et al. 2019; also wormwiring.org | `data/raw/connectome/cook2019_*` |
| 2 | GJ connectome B (robustness) | Varshney et al. 2011, *PLoS Comput Biol* 7(2) (SI ref 1). 517 links, 3 self-loops | Paper's supplementary data; wormatlas.org | `data/raw/connectome/varshney2011_*` |
| 3 | Contact matrix C | Brittin et al. 2018 bioRxiv:485771 (ref 31) via wormwiring.org; WormBase (ref 34). Anterior 185 neurons; 5,592 contacting pairs (~33% of pairs), 601 form GJs | wormwiring.org; WormBase | `data/raw/contact/` |
| 4 | Innexin expression X | Altun et al. 2009, *Dev Dyn* 238:1936 (ref 45): 18 of 25 innexins. Binarized and curated per Hobert, Glenwinkel & White 2016, *Curr Biol* 26:R1197 (ref 5) | Altun 2009 + WormAtlas gap-junction maps; Hobert lab (authors thank Oliver Hobert for sharing data) | `data/raw/expression/` |
| 5 | Neuron name list / class map | White et al. 1986; WormAtlas | wormatlas.org | `data/raw/neuron_names.csv` |

The 18 innexins in X: inx-1, -2, -3, -5, -6, -7, -8, -9, -10, -11, -12, -13, -17, -18, -19, che-7, unc-7, unc-9.
Expression is only resolved to neuron class, not individual neurons.

## Steps

1. **Harmonize neuron identifiers.** Map every source to one naming scheme
   (302 neurons; 279 somatic connected; 185 anterior in the contact set).
   Handle L/R and class-level names (expression is per class).
2. **Build B.** Parse Cook 2019 GJ edges, symmetrize, keep self-loops, restrict to the 185 anterior neurons. Repeat for Varshney.
3. **Build C.** Extract the contact list, restrict to the same 185 neurons, check count (5,592 pairs).
4. **Build X.** Encode Altun/Hobert expression as 0/1 per neuron class, expand to neurons.
5. **Validate against the authors' processed data** (see below): row/column counts, 601 GJ pairs among contacting pairs, 18 innexin columns.
6. **Run the SCM** (ridge regression, α from generalized cross-validation, subgraph randomization for z-scores, 100 iterations). Expect 19 rules with z ≥ 2.
7. **Robustness run** on Varshney. Expect 12 shared rules (inx-3–inx-3, inx-9–inx-17, inx-19–inx-19, unc-9–unc-9, inx-9–inx-9, inx-13–inx-13, inx-12–inx-12, inx-2–inx-10, inx-6–inx-18, inx-7–inx-7, inx-8–inx-12, che-7–inx-17).
8. **Compare** to SI Figs S3, S4, S6; check eigenvalues of O (at least 4 negative).

## Ground truth already on this machine

- `4027588.zip` is the Zenodo record (DOI 10.5281/zenodo.4027588) for the paper. It holds `SCM-v1.0.zip`.
- `kpisti/kpisti-SCM-6cc8731/` has the authors' code (`SCM_alpha.py`, `SCM_alpha_rnd.py`, `SCM_alpha_rndpinv.py`, `alpha_matrix.cpp`) and processed inputs: `GapJunctContact.csv`, `ContactSubgraphMatrix.csv`, `INXExpressionJustContact.csv`, `NeurNamesContact.csv`, plus their outputs (`rule_matrix_0.2150*.csv`).
- Use these to check your rebuilt B, C, X. Build from raw first, then diff.

## Open questions

- Exact form of the curated expression table (Hobert lab) may need to be requested; Altun 2009 alone is images and text.
- Cook 2019 distinguishes hermaphrodite and male; use hermaphrodite.
- ~~Confirm which Cook 2019 supplement gives the GJ list.~~ Resolved: the WormWiring edge list `herm_full_edgelist.csv`, not SI 5 (see `data/raw/SOURCES.md`).
- The unusual `inx-5` / `inx-11` rows in SI Fig. S3 are identical, which suggests inx-5 is a duplicate of inx-11 in X. Check the data.
- Heteromeric rules and chemical synapses are out of scope for the paper (SI text).

## Reference papers

Kept locally, not committed: `kovács-et-al-2020-...pdf` and `pnas.2009093117.sapp.pdf` (both gitignored).
