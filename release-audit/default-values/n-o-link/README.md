# cctbx links a lysine nitrogen to a ligand oxygen at 1.297 A
Date: 2026-11-11

## Question
Thirteen unrecorded protein-ligand pairs in the 4 Nov and 11 Nov PDB releases sit at exact repeated distances, eleven of them N-O at about 1.297 A, all in entries that declare PHENIX. cctbx's generic table has ("N","O") = 1.297 (mmtbx/monomer_library/bondlength_defaults.py line 592, four lines above the Mg-N 1.895 value). Does the restraint builder really create that bond between a side chain and a ligand?

## Test
Fragment from 36YY (X-ray, 1.9 A): Lys68 of chain D and ATP 201, deposited NZ-O1B 1.298 A. The ATP was moved along the NZ-O1B line to set other separations. cctbx-base 2025.11, default pdb_interpretation parameters, dictionaries geostd f4d21cc + mon_lib c382734.

## Result
| NZ-O1B in the input (A) | Bond restraint made | Ideal | Sigma |
|---|---|---|---|
| 1.298 (as deposited) | yes | 1.297 | 0.01 |
| 2.00 | yes | 1.297 | 0.01 |
| 2.10 | yes | 1.297 | 0.01 |
| 2.15 | yes | 1.297 | 0.01 |
| 2.20 | yes | 1.297 | 0.01 |
| 2.25 | no | | |
| 2.29 | no | | |
| 2.40 | no | | |
| 2.80 | no | | |
| 3.10 | no | | |

At default settings, a lysine NZ within 2.20 A of an ATP phosphate oxygen is given a covalent bond restraint with ideal 1.297 A and sigma 0.01 A. At 2.25 A and beyond there is none. A normal salt bridge (2.7-3.0 A) is left alone; a pair that has drifted to 2.2 A or closer is pulled the rest of the way to 1.297 A and held there.

## Not tested
- Other residue-ligand pairs (Arg NH1, His, Asn; phosphate, carboxylate and heme oxygens). The other twelve pairs are value matches only.
- Older builds, and real-space refinement specifically.
- The C-O 1.374 value on the same code path.
- How many archive entries carry an N-O pair at 1.297 A.

## Files
- `n_o_link_output.txt`: the restraint builder's bond lines for each input.
- `mk.py`, `run.py`: fragment builder and runner (`mk.py <36YY.cif.gz>`; set `MMTBX_CCP4_MONOMER_LIB`).
- `frag_*.pdb`: the ten inputs.
