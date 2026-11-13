# Connectome

Reproduction of Kovács, Barabási & Barabási (2020), "Uncovering the genetic blueprint
of the C. elegans nervous system", PNAS 117(52). Built from the raw data up.

- Plan and data sources: [docs/DATA_COLLECTION_OUTLINE.md](docs/DATA_COLLECTION_OUTLINE.md)
- `data/raw/` downloaded source data (connectomes, contact, expression)
- `data/processed/` rebuilt B, C, X matrices
- `src/` pipeline code, `results/` outputs and recreated SI figures
- [`supplementary/`](supplementary/README.md) explanation of the model, procedure, data and terms, with graphics

## How to run

Needs Python 3, a C++ compiler (`g++`) and `curl`. Run from the repository root.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
bash src/download_data.sh
g++ -O3 -o results/bin/alpha_matrix.out src/alpha_matrix.cpp
```

Reproduce the authors' rules from their processed inputs:

```bash
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/reference_authors/GapJunctContact.csv data/reference_authors/INXExpressionJustContact.csv 0.215 results/authors_inputs
```

Rebuild the authors' B exactly from the WormWiring edge list (prints 0 differences):

```bash
python src/build_cook_B_wormwiring.py
```

Robustness run on the Varshney 2011 connectome (the paper uses alpha = 0.12 here):

```bash
python src/build_varshney_B.py
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_varshney.csv data/reference_authors/INXExpressionJustContact.csv 0.12 results/varshney
```

Independent replication on the two Witvliet et al. 2021 adult connectomes (new animals), at each one's
cross-validated alpha and at alpha = 0.215, then compare the rules across all four connectomes:

```bash
python src/build_witvliet_B.py
python src/select_alpha.py data/processed/B_witvliet_adult1.csv results/witvliet_adult1/alpha_scan.csv
python src/select_alpha.py data/processed/B_witvliet_adult2.csv results/witvliet_adult2/alpha_scan.csv
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_witvliet_adult1.csv data/reference_authors/INXExpressionJustContact.csv 15.849 results/witvliet_adult1
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_witvliet_adult2.csv data/reference_authors/INXExpressionJustContact.csv 6.3096 results/witvliet_adult2
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_witvliet_adult1.csv data/reference_authors/INXExpressionJustContact.csv 0.215 results/witvliet_adult1
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_witvliet_adult2.csv data/reference_authors/INXExpressionJustContact.csv 0.215 results/witvliet_adult2
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_varshney.csv data/reference_authors/INXExpressionJustContact.csv 0.215 results/varshney
python src/compare_connectomes.py gcv
python src/compare_connectomes.py fixed
```

Rebuild B from Cook 2019 SI 5 instead and rerun with it:

```bash
python src/build_cook_B.py
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_cook_rebuilt.csv data/reference_authors/INXExpressionJustContact.csv 0.215 results/cook_rebuilt
```

Recreate SI Figs. S2-S6 in `results/figures/` (needs the Cook and Varshney runs above; takes a few minutes), then
compare our numbers with the tables printed in the SI PDF (not in the repository; download it from PNAS):

```bash
python src/plot_si_figures.py
python src/check_against_si.py pnas.2009093117.sapp.pdf
python src/plot_supplementary.py
```

See [results/README.md](results/README.md) for what the outputs mean.
