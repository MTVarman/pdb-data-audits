# Runs N' Poses covalent entries: Boltz-2, SMILES-only vs bond-declared

Supplementary Table S12 of the preprint (v1.8). The six Runs N' Poses entries whose benchmark ligand is covalently bonded in the deposited coordinates with no `_struct_conn` record (see `../strict-cut-491/` and `../benchmark-overlap/strict_subset_in_runs_n_poses.csv`), each predicted twice with Boltz-2 2.2.1: ligand as SMILES with no bond (the benchmark input), and ligand as chemical-component code with a bond constraint.

## Result

Ligand heavy-atom RMSD (Å) to the deposited pose, pocket-Cα superposition (Cα within 10 Å of the ligand), symmetry-corrected. Top-ranked sample / best of 5.

| Entry | Ligand | Bond | SMILES-only | Bond declared |
|---|---|---|---|---|
| 7EI4 | J3U | Cys90 SG–CAT | 1.08 / 0.98 | 0.91 / 0.75 |
| 7F1V | 7XF | Lys211 NZ–CAL | 0.78 / 0.55 | 0.42 / 0.33 |
| 7YCC | IQC | Cys12 SG–C36 | 1.08 / 1.02 | 0.96 / 0.89 |
| 7YCE | IQN | Cys12 SG–C31 | 0.77 / 0.77 | 0.77 / 0.67 |
| 8DVH | BO2 | Ser445 OG–B26 | 0.59 / 0.57 | 0.48 / 0.43 |
| 8X6R | Y9D | Cys12 SG–C46 | 0.43 / 0.42 | 0.53 / 0.53 |

- Bond declared: nucleophile–ligand distance 1.55–1.93 Å in all 30 samples.
- SMILES only: ligand centroid within 0.1–0.75 Å of the deposited centroid in all 30 samples; nucleophile to nearest ligand atom 2.6–3.5 Å (contact distance, not bonded).
- All 60 samples: complex pLDDT and pTM 0.94–0.98.

## Setup

Boltz-2 2.2.1, Python 3.12 venv, NVIDIA A100. `--use_msa_server --use_potentials --no_kernels --diffusion_samples 5 --seed 42`. One copy of the chain carrying the site; GDP and Mg for the three KRAS entries. Residue indices mapped from author numbering through `_pdbx_poly_seq_scheme`. The two conditions share identical protein input.

## Caveats

- SMILES inputs were generated with RDKit from the CCD definition files, without stereochemistry flags.
- Reported Boltz-2 training cutoffs differ (2021-09-30, 2023-06-01). Under the later date five of the six entries pre-date it; only 8X6R (released 2024) post-dates both. The within-entry comparison does not depend on this.
- One method, six entries.

## Files

| File | Contents |
|---|---|
| `summary.md` | The table above, three decimals |
| `results.csv` | Every sample: entry, ligand, condition, model rank, confidence, ipTM, ligand RMSD, symmetry matches, pocket Cα count |
| `manifest.json` | Entry → chain, sequence index of the nucleophile, SMILES used, cofactors |
| `inputs_yaml.zip` | The 12 Boltz input YAMLs |
| `run_rnp_covalent_as_run.py` | Preparation and scoring script, as run |
