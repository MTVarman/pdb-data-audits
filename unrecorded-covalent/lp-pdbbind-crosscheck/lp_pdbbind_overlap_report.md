# Unrecorded covalent scan vs LP-PDBBind

Source: https://raw.githubusercontent.com/THGLab/LP-PDBBind/master/dataset/LP_PDBBind.csv
SHA-256 of downloaded file: `7bb4e54e66c0c58ab74263edb8c2a2d677c2474cfc4ae4dec1eb31a231b10f8a`
Comparison is at **PDB entry level only**: LP-PDBBind has no ligand column, so a hit means the
entry is in both tables, not that the same ligand is involved.

## 1. Headline contradiction

- Headline entries in the scan: 2568; present in LP-PDBBind: **169**
- Of those, labelled covalent=False: **139** (covalent=True: 30)

Covalent=False entries by triage class (an entry is counted in each class it has sites in):

- known cofactor or PTM: 8 entries (15 sites)
- bond-length consistent adduct: 39 entries (64 sites)
- likely artefact: 98 entries (173 sites)
- Entries whose only sites are likely artefact (discard): 92

## 2. Strict subset (class 2, X-ray, <= 2.5 A)

- Strict headline entries: 289; in LP-PDBBind: 34
- Labelled covalent=False: **27 entries (38 sites)**; by split: test 6, train 19, val 2
- Excluding DTT (a buffer reducing agent, not a binding ligand): **25 entries (34 sites)**;
  by split: test 4, train 19, val 2
- Bond types across the drug-like hits: C-O 11, C-S 10, C-N 9, S-S 2, B-O 2

Full per-site table in `lp_pdbbind_unrecorded_covalent_hits.csv`.

Catalytic residues represented: beta-lactamase Ser70 (boronate adducts, 1NXY/1NYM),
proteasome Thr1 (38N, four chains of 4QLQ, both OG1 and backbone N), 3C-like protease
Cys145 (2VJ1), papain-family Cys25 (3O0U, 2R9M, 3N3G), glycosidase Glu526 (6NZG),
serine hydrolase Ser81 (3FV7/3FYZ/3FZC), His212 (5YKP/5YR4/5YR5/5YR6).

## 3. Sensitivity control (detector vs their covalent=True)

- LP-PDBBind covalent=True: 893; also a headline site in the scan: **30 / 893 = 3.4%**
- Those by class: known cofactor or PTM 0, bond-length consistent adduct 11, likely artefact 20
- covalent=True entries with any contact in the verified-sites table: 49

The scan only reports bonds with no connectivity record, so recorded covalent bonds are
excluded by design; a low share is therefore expected and does not by itself test the detector.
The 30 entries that LP-PDBBind marks covalent AND that still carry an unrecorded contact are a
separate observation: the dataset curators knew the complex was covalent while the coordinate
file still lacks the bond record.

## Known limits

- Entry-level matching only. The LP-PDBBind binding ligand is not proven to be the same
  component carrying the unrecorded bond; per-entry confirmation is required before the
  claim is made in public. The DTT entries are a worked example of where that inference fails.
- Class 2 means the contact distance falls inside a bond-length window for that element pair.
  It is evidence of a covalent bond, not proof of one.
