# Runs N' Poses covalent entries: Chai-1, and the warhead distance across two models

Companion to `../rnp-boltz2/`. The same six Runs N' Poses entries whose benchmark ligand is covalently bonded in the deposited coordinates with no `_struct_conn` record, predicted again with Chai-1 under the same two conditions, and a single measurement of the warhead-to-nucleophile distance applied to both models' predictions and to the crystals.

The question this answers was put by Fernando M. Delgado-Chaves: given the same SMILES-only input, does a second co-folding model also place the reactive atom at van der Waals contact rather than bonded? If every model does, the behaviour belongs to the input format. If only one does, it belongs to that model.

**It is not a property of the input format.** Given identical SMILES-only input, Chai-1 closes the warhead to covalent bond distance in four of the six entries. Boltz-2 does so in none of them.

## Result

Warhead-to-nucleophile distance (Å), SMILES-only arm — the input the benchmark actually uses. Ranges across all 5 samples. `xtal` is the deposited bond length. Pose is symmetry-corrected ligand heavy-atom RMSD after pocket-Cα superposition.

| Entry | Ligand | Bond | xtal | Boltz-2 pose | Boltz-2 warhead | Chai-1 pose | Chai-1 warhead |
|---|---|---|---|---|---|---|---|
| 7EI4 | J3U | Cys90 SG–CAT | 1.76 | 0.97–1.22 | 2.74–3.22 | 3.86–11.32 | 2.16–19.41 |
| 7F1V | 7XF | Lys211 NZ–CAL | 1.32 | 0.55–0.78 | 3.87–3.91 | 0.37–0.48 | 2.38–2.67 |
| 7YCC | IQC | Cys12 SG–C36 | 1.82 | 1.02–1.12 | 2.75–2.78 | 0.80–1.59 | **1.65–1.82** |
| 7YCE | IQN | Cys12 SG–C31 | 1.84 | 0.77–0.85 | 2.74–2.78 | 0.74–0.80 | **1.65–1.79** |
| 8DVH | BO2 | Ser445 OG–B26 | 1.47 | 0.57–1.11 | 2.64–2.69 | 2.14–2.55 | **1.55–1.62** |
| 8X6R | Y9D | Cys12 SG–C46 | 1.80 | 0.42–0.61 | 2.74–2.77 | 0.46–0.67 | **1.64–1.87** |

Samples at bond distance (< 2.0 Å), out of 30:

| | Boltz-2 | Chai-1 |
|---|---|---|
| SMILES-only | **0 / 30**, median 2.75 Å | **20 / 30**, median 1.79 Å |
| bond declared | 30 / 30, median 1.62 Å | 30 / 30, median 1.65 Å |

Three observations:

1. Boltz-2's SMILES-only distances are a hard stop rather than a near miss. Across four different proteins and two warhead chemistries (acrylamide, boronic acid) they fall between 2.64 and 2.78 Å — a 0.14 Å window across unrelated systems.
2. On those four entries, declaring the bond barely moves Chai-1 (8X6R 1.78 → 1.72 Å). The same declaration moves Boltz-2 a full ångström (2.75 → 1.61 Å).
3. **8DVH separates pose from chemistry in both directions.** Chai-1 fails a 2.0 Å pose cut (2.14–2.55 Å) while placing the warhead at bond distance (1.55–1.62 Å against a deposited 1.47 Å). Boltz-2 passes the pose cut (0.57–1.11 Å) with the warhead a full ångström off. On that entry a pose-RMSD cut ranks the two models in the opposite order from the geometry the missing annotation is about.

## What this does not show

Chai-1 is not uniformly better. Its pose is better on 7F1V, comparable on 7YCC / 7YCE / 8X6R, worse on 8DVH, and it fails 7EI4 outright (RMSD 3.86–11.32 Å; the ligand flips in 3 of 5 samples, warhead 18–19 Å away). Chai-1 reports that itself: mean aggregate score 0.444 on 7EI4 against 0.93–0.96 on the other five. 7F1V is the one entry where neither model reaches bond distance.

The claim is confined to the warhead distance: stopping at van der Waals contact under SMILES-only input is Boltz-2 behaviour, not a consequence of the input format.

## How the warhead was located

Neither model names predicted ligand atoms as the chemical component dictionary does, and the two disagree with each other — for 8X6R's warhead Boltz-2 writes `C77`, Chai-1 writes `C15_1`, the CCD says `C46`. `warhead_compare.py` therefore uses no names: it infers each predicted ligand's bond graph from its own coordinates (covalent radii, 0.45 Å tolerance), graph-matches that to the CCD definition, and takes the warhead to be whichever predicted atom maps to the CCD warhead atom. The crystals are measured the same way.

Checks on the measurement:

- Applied to the deposited coordinates it returns the bond lengths: 1.76 / 1.32 / 1.82 / 1.84 / 1.47 / 1.80 Å.
- For 8X6R's Boltz-2 prediction it selects atom `C77` at 2.75 Å, matching the atom and distance used in the preprint figure, which were obtained by a different route.
- The Chai-1 distances agree with those produced by `chai_rnp.py score`, a separate code path.
- The warhead maps to a single candidate atom in all 126 measurements; no symmetry ambiguity arises.

Independently, `chai_rnp.py score` verifies that every SMILES bond in each prediction spans 0.9–2.1 Å. All 60 Chai-1 predictions pass (`bad_bonds` = 0), confirming the atom-order mapping.

Note on `../rnp-boltz2/README.md`: the "2.6–3.5 Å" recorded there is the nucleophile to *nearest ligand atom* distance, which here is 2.64–3.55 Å. The warhead specifically is 2.64–3.91 Å. The two differ only on 7F1V, where Boltz-2's warhead is not the closest atom to the lysine. Both are reported per sample in `warhead_compare.csv`.

## Setup

Chai-1 0.6.1, Python 3.12 venv, NVIDIA A100 40 GB. 5 diffusion samples, seed 42, MSAs from the ColabFold server. One copy of the chain carrying the site; GDP and Mg for the three KRAS entries. The two conditions share identical protein input. Entries, SMILES, nucleophile and warhead atoms are the same as `../rnp-boltz2/manifest.json`.

- `smiles` arm: ligand as SMILES, no bond — the benchmark input.
- `bond` arm: the same, plus a covalent restraint from the nucleophile to the warhead atom.

## Caveats

- SMILES inputs were generated with RDKit from the CCD definition files, without stereochemistry flags.
- Two models, six entries. Chai-1's training cutoff is not stated as precisely as Boltz-2's; the within-entry and between-arm comparisons do not depend on it.
- AF3, the other half of the original question, is untested — the weights are gated.
- Predicted structure files (60 CIFs, 49 MB) are not included here.

## Files

| File | Contents |
|---|---|
| `warhead_compare.csv` | 126 rows: every prediction from both models and the six crystals — warhead distance, nearest-ligand-atom distance, whether they are the same atom, warhead candidate count |
| `results.csv` | 60 rows: entry, ligand, condition, sample, Chai-1 aggregate score, ligand RMSD, warhead and nearest distances, centroid offset, bad-bond count, symmetry matches, pocket Cα count |
| `summary.md` | RMSD table, top-ranked / best of 5 |
| `manifest.json` | Entry → chain, sequence index of the nucleophile, SMILES used, cofactors |
| `warhead_compare.py` | The cross-model warhead measurement, as run |
| `chai_rnp.py` | Preparation, prediction and scoring script, as run |
| `chai_rnp_covalent.ipynb` | The notebook the run was executed from |
