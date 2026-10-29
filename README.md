# Connectome

Reproduction of Kovács, Barabási & Barabási (2020), "Uncovering the genetic blueprint
of the C. elegans nervous system", PNAS 117(52). Built from the raw data up.

- Plan and data sources: [docs/DATA_COLLECTION_OUTLINE.md](docs/DATA_COLLECTION_OUTLINE.md)
- `data/raw/` downloaded source data (connectomes, contact, expression)
- `data/processed/` rebuilt B, C, X matrices
- `src/` pipeline code, `results/` outputs

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

Rebuild B from Cook 2019 SI 5 instead and rerun with it:

```bash
python src/build_cook_B.py
python src/scm.py data/reference_authors/ContactSubgraphMatrix.csv data/processed/B_cook_rebuilt.csv data/reference_authors/INXExpressionJustContact.csv 0.215 results/cook_rebuilt
```

See [results/README.md](results/README.md) for what the outputs mean.
