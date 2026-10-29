# SCM results (alpha = 0.215, Cook GJ connectome, 185 anterior neurons)

Code: `src/scm.py` (Python 3 port of the authors' SCM_alpha*.py), `src/alpha_matrix.cpp` (authors' randomization, unchanged),
`src/build_cook_B.py` (rebuilds B from Cook 2019 SI 5). Compile: `g++ -O3 -o results/bin/alpha_matrix.out src/alpha_matrix.cpp`.

| Run | Folder | Result |
|-----|--------|--------|
| Authors' B, C, X | `authors_inputs/` | Matches the authors' `rule_matrix_0.2150.csv` and z-scores to ~1e-12. 19 positive rules with z >= 2, as in the paper. 7 negative eigenvalues of the symmetrized O. |
| Rebuilt Cook B (SI 5, binarized, masked by C; 602 links) with authors' C and X | `cook_rebuilt/` | 19 rules, 14 shared with the authors' set. Correlation of O with the authors' O: 0.973 (mean abs diff 0.031, max 0.374). |

The authors' B is the WormWiring Cook edge list, not SI 5: `src/build_cook_B_wormwiring.py` rebuilds it with 0 differences
(see `data/raw/SOURCES.md`), so the `authors_inputs/` run is also the run on B rebuilt from raw data. The `cook_rebuilt/` run
shows how much the rules change between the two Cook releases.

Rebuilt-B (SI 5) differences (all near the z = 2 cutoff):
- Lost: inx-1/inx-18 (z 1.6), inx-10/inx-2 (1.9), inx-13/inx-13 (1.4), inx-6/inx-6 (1.7), inx-7/inx-7 (1.9)
- Gained: inx-11/inx-12 (2.1), inx-12/inx-5 (2.1), inx-13/inx-7 (2.0), inx-17/inx-19 (2.3), inx-7/unc-7 (2.2)
- Strong rules are stable: inx-12/inx-12, inx-9/inx-9, inx-10 with inx-11, inx-3, inx-5, inx-19/inx-19, unc-9/unc-9.
