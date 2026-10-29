# SCM results (alpha = 0.215, Cook GJ connectome, 185 anterior neurons)

Code: `src/scm.py` (Python 3 port of the authors' SCM_alpha*.py), `src/alpha_matrix.cpp` (authors' randomization, unchanged),
`src/build_cook_B.py` (rebuilds B from Cook 2019 SI 5). Compile: `g++ -O3 -o results/bin/alpha_matrix.out src/alpha_matrix.cpp`.

| Run | Folder | Result |
|-----|--------|--------|
| Authors' B, C, X | `authors_inputs/` | Matches the authors' `rule_matrix_0.2150.csv` and z-scores to ~1e-12. 19 positive rules with z >= 2, as in the paper. 7 negative eigenvalues of the symmetrized O. |
| Varshney B (rebuilt from `ConnOrdered_040903.mat`, identical to the authors' `CEO_Gap_OnlyContact.csv`; 301 links), alpha = 0.12 as in SI Fig. S4 | `varshney/` | Matches SI Fig. S4 to its printed precision (O within 0.005, z within 0.05). 21 rules with z >= 2, 11 shared with Cook. 6 negative eigenvalues. |
| Rebuilt Cook B (SI 5, binarized, masked by C; 602 links) with authors' C and X | `cook_rebuilt/` | 19 rules, 14 shared with the authors' set. Correlation of O with the authors' O: 0.973 (mean abs diff 0.031, max 0.374). |

The authors' B is the WormWiring Cook edge list, not SI 5: `src/build_cook_B_wormwiring.py` rebuilds it with 0 differences
(see `data/raw/SOURCES.md`), so the `authors_inputs/` run is also the run on B rebuilt from raw data. The `cook_rebuilt/` run
shows how much the rules change between the two Cook releases.

Rebuilt-B (SI 5) differences (all near the z = 2 cutoff):
- Lost: inx-1/inx-18 (z 1.6), inx-10/inx-2 (1.9), inx-13/inx-13 (1.4), inx-6/inx-6 (1.7), inx-7/inx-7 (1.9)
- Gained: inx-11/inx-12 (2.1), inx-12/inx-5 (2.1), inx-13/inx-7 (2.0), inx-17/inx-19 (2.3), inx-7/unc-7 (2.2)
- Strong rules are stable: inx-12/inx-12, inx-9/inx-9, inx-10 with inx-11, inx-3, inx-5, inx-19/inx-19, unc-9/unc-9.

## Gene order

Rows and columns of every `rule_matrix_*.csv` and `zscores_*.csv` (and the columns of X) follow SI Figs. S3/S4:
inx-1, inx-10, inx-11, inx-12, inx-13, inx-17, inx-18, inx-19, inx-2, inx-3, inx-5, inx-6, inx-7, inx-8, inx-9, che-7, unc-7, unc-9.
(inx-5 and inx-11 have identical expression in X, so their rows are identical.)

The `authors_inputs/` run matches SI Fig. S3 to its printed precision (O within 0.005, z within 0.05).

## Varshney robustness run (SI text, Fig. S4)

The paper reports 12 rules shared between Cook and Varshney. We reproduce all 12; 11 have z >= 2 in both runs, and the 12th,
inx-8/inx-12, has Varshney z = 1.96, which SI Fig. S4 prints as 2.0. No other rules are shared.

| Shared rule | Cook z | Varshney z |
|-------------|--------|------------|
| inx-12/inx-12 | 5.61 | 7.04 |
| inx-9/inx-9 | 5.00 | 5.53 |
| inx-3/inx-3 | 2.83 | 6.51 |
| unc-9/unc-9 | 3.15 | 4.55 |
| inx-19/inx-19 | 3.44 | 4.49 |
| inx-13/inx-13 | 2.40 | 4.44 |
| inx-7/inx-7 | 2.36 | 3.85 |
| che-7/inx-17 | 2.84 | 3.19 |
| inx-9/inx-17 | 2.45 | 2.65 |
| inx-6/inx-18 | 2.42 | 2.16 |
| inx-2/inx-10 | 2.09 | 2.16 |
| inx-8/inx-12 | 2.57 | 1.96 |

- Cook only: inx-1/inx-18, inx-10/inx-11, inx-10/inx-3, inx-10/inx-5, inx-13/inx-17, inx-13/inx-8, inx-6/inx-6 (plus inx-8/inx-12 above).
- Varshney only: inx-1/inx-11, inx-1/inx-17, inx-1/inx-5, inx-11/inx-12, inx-12/inx-5, inx-13/inx-18, inx-18/inx-18, inx-18/unc-7, inx-2/inx-2, inx-6/che-7.
- Correlation of the two O matrices: 0.82 (z-scores: 0.83). Varshney has half as many links (301 vs 601), so it is the noisier fit.
- `varshney/*_0.2150.*` is the same run at the Cook alpha for comparison: 20 rules, the same 11 shared.
