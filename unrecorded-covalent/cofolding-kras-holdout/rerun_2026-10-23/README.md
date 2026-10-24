# Re-run of the co-folding arms lost on 22 October 2026

Working files for the covalent KRAS temporal-holdout benchmark (Chai-1, Boltz-2, Protenix v1.0.0 and 20250630). The per-sample rows for the Protenix-20250630 arm, the Boltz-2 re-run on all 21 entries, and both models' 40-entry control arms were written on 22 October to a Google Drive mount on an account that was over quota; the mount listed the files while the upload failed, and the runtime was reclaimed. This folder re-creates them.

`kras_cofold_rerun_2026-10-23.ipynb` is the Colab notebook (A100). It embeds everything else here and runs `queue.sh` under `nohup`; each stage ends with `persist.sh`, which scores, computes symmetry-corrected pose RMSD, and copies results plus a tar of the predicted CIFs to Drive before the next stage starts.

| stage | model | set | jobs |
|---|---|---|---|
| 1a | Protenix-20250630 | Tier B + C, both arms | 30 |
| 1b | Protenix-20250630 | Tier A, both arms | 12 |
| 2 | Boltz-2 | all 21 entries, both arms (15 bond arms fail on the frozen CCD by design) | 42 |
| 3 | Boltz-2 | 40 pre-cutoff controls, SMILES only | 40 |
| 4 | Protenix v1.0.0 | 40 pre-cutoff controls, SMILES only | 40 |

`kras_holdout.py` is the audited runner (sha256 d17a5e53…, residue-index fix applied). `poserms.py` is the symmetry-corrected pose RMSD, extended here with Boltz-2's atom-naming scheme so Boltz-2 outputs can be pose-scored too. `collect.py` flattens raw model output into `out/<model>/<entry>_<arm>/sample_i.cif` from a directory outside `out/`, which is what prevents the double count that inflated the 21 October Boltz-2 Tier C figures.

Settings unchanged from the original runs: seed 42, 5 diffusion samples, 200 steps, ColabFold MSAs, nucleotide and Mg present.
