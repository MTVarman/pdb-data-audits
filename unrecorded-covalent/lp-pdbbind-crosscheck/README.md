# LP-PDBBind × the strict cut: 20 covalent complexes labelled non-covalent, and why

**Date:** 2026-10-17
**LP-PDBBind source:** `THGLab/LP-PDBBind`, `dataset/LP_PDBBind.csv`, SHA-256 `7bb4e54e66c0c58ab74263edb8c2a2d677c2474cfc4ae4dec1eb31a231b10f8a`; curation code in `dataset_creation/CL_filters.ipynb` and `dataset_creation/prepare_info.ipynb`, same commit.
**Scan source:** the 491-site strict cut in `../strict-cut-491/`.
**Reproduce:** `python3 -I verify_lp.py LP_PDBBind.csv unrecorded_covalent_tight_cut_491_sites.csv github_issue_8_as_filed.md`

## Numbers

| | Entries | Sites |
|---|---:|---:|
| Strict cut (class 2, X-ray ≤ 2.5 Å) | 289 | 491 |
| …present in LP-PDBBind | 34 | — |
| …labelled `covalent=False` | 27 | 38 |
| …after per-entry ligand verification: **confirmed or plausible covalent inhibitor of the LP-PDBBind ligand** | **20** | **29** |
| …doubtful (clash or synthetic construct) | 5 | 5 |
| …excluded (bond is to buffer DTT, not the tracked ligand) | 2 | 4 |

The 7 strict-cut entries LP-PDBBind *does* label covalent: 1K2I, 2OP9, 4WX4, 5D6F, 5OM9 (train); 6H0U, 6S9W (val).

Splits for the 20: train 16, val 2 (2R9M, 3N3G), test 2 (3O0U, 5Z1D).

Per-site table with LP-PDBBind label, split, Clean Level flags and SMILES: `strict_vs_lp_27_entries.csv`. Entry-level overlap report from the scan pipeline: `lp_pdbbind_overlap_report.md`, `lp_pdbbind_unrecorded_covalent_hits.csv`.

## How LP-PDBBind assigns `covalent` (from the notebook)

`CL_filters.ipynb`, cells 22–28:

```
cov = read CovBinderInPDB_2022Q4_AllRecords.csv
covalent_in_pdbbind = [entries whose PDB ID is in cov]
false_covalent = [hand-curated list of 26 IDs]
covalent = True for (covalent_in_pdbbind − false_covalent), else False
```

The flag is CovBinderInPDB membership minus a hand list. No geometry is consulted. An entry whose covalent bond has no `_struct_conn` / LINK record is invisible to CovBinderInPDB and therefore to LP-PDBBind.

## Three findings

### 1. 3N3G was flagged covalent and manually overridden to non-covalent

3N3G is the first entry in the `false_covalent` hand list. CovBinderInPDB recorded it as covalent; the curators removed it. In the deposited structure, Cys25 SG sits **1.721 Å** from the nitrile carbon of ligand 935 — a thioimidate bond length. The LP-PDBBind SMILES for 3N3G carries an intact nitrile (`C#N`). The entry passes **CL1, CL2 and CL3** and sits in the `val` split as `covalent=False`.

Our density validation marks the 3N3G ligand `ligand_uncertain` (RSCC 0.685, occupancy 0.75), so a density-based curator judgement is possible. 3N3G is the only `false_covalent` entry that falls in the strict cut.

The cleaner companion is **2R9M** (val, `covalent=False`, CL1=True): the same Cys25 thioimidate chemistry, both chains density-supported (RSCC 0.923 / 0.930), C–S at 1.878 / 1.892 Å, and a SMILES that encodes the *reacted* `C=N`. Same warhead class, same catalytic residue, same split, same label — one missed by CovBinderInPDB, one overruled by hand, with opposite SMILES reaction states.

### 2. The CL1 1.75 Å clash filter is an accidental, chemically biased covalent filter

Cell 4 sets `CL1 = False` when the minimum heavy-atom distance between PDBBind's `_pocket.pdb` and `_ligand.mol2` is below 1.75 Å. Covalent bond lengths straddle that number: C–O, C–N and B–O adducts sit at 1.3–1.55 Å and are removed; C–S and S–S adducts sit at 1.7–2.0 Å and pass.

| CL1 | Entries | Chemistry |
|---|---|---|
| `False` (11) | 1NYM, 2A4R, 3FV7, 3FYZ, 3FZC, 4QLQ, 5YKP, 5YR4, 5YR5, 5YR6, 6NZG | all Ser/Thr/His/Glu adducts (C–O, C–N, B–O) |
| `True` (9) | 3O0U, 5Z1D, 2R9M, 3N3G, 2VJ1, 4YV2, 5YJO, 6IJL, 1NXY | 8 cysteine adducts (C–S, S–S) + one boronate |

LP-PDBBind's *clean* tiers are therefore enriched for cysteine covalent adducts labelled non-covalent. Both test-split cases and both val-split cases are cysteine and CL1=True.

Open: 1NXY (1.525 Å), 3N3G (1.721), 4KB7 (1.728) and 4U6W (1.749) have deposited contacts under 1.75 Å yet pass CL1; 1NYM, with the same chemistry and distance as 1NXY, does not. PDBBind's pocket/ligand files likely differ from the deposited coordinates. Not verified.

### 3. SMILES are coordinate-derived and encode the post-reaction adduct

`prepare_info.ipynb` cell 7: `Chem.MolFromMol2File(…_ligand.mol2)` → `RemoveHs` → `MolToSmiles`. The SMILES is whatever bond orders PDBBind's mol2 carried, so reacted adducts appear as products:

| Entry (split) | Warhead | LP-PDBBind SMILES shows |
|---|---|---|
| 5Z1D (test) | acrylamide | saturated propanamide `CCC(=O)N…` |
| 3O0U (test), 2R9M (val) | nitrile | thioimidate `C=N` |
| 3N3G (val) | nitrile | intact `C#N` |
| 2A4R (train) | α-ketoamide | α-hydroxy amide |
| 3FV7, 3FYZ, 3FZC (train) | 6-alkylidene penam sulfone | ring-opened acyl-enzyme product |
| 4QLQ (train) | epoxyketone | opened diol |
| 5YKP, 5YR4–5YR6 (train) | spiroepoxide | ring-opened tertiary alcohol |
| 6NZG (train) | cyclitol epoxide | ring-opened |
| 1NXY, 1NYM (train) | boronic acid | `B(O)O` |

## Density and clash cross-check (from the strict-cut table)

Of the 29 confirmed/plausible sites, **24 are `density_supported`**; five are `ligand_uncertain`: 5Z1D (RSCC 0.752), 3N3G (0.685), 3FV7 (0.751), 3FZC (0.756), 6IJL (0.757). Fifteen of the 20 entries are fully supported. None of the 38 sites is recorded as a steric-clash outlier at the flagged pair in the wwPDB validation report.

## Doubtful and excluded, with reasons

- **3AAV** (test): bis-benzamidine Cu(II) chelate; 1.550 Å sits at the C–O window edge; reversible S1-pocket binder. Clash.
- **3KIG** (test): sulfonamide S(VI) cannot form a disulfide with Cys64. Valence.
- **4KB7, 4U6W** (train): Met SD to an unactivated carbon at 1.73–1.75 Å would require a sulfonium. Clash / alternate conformation.
- **4R1E** (train): Gly802 N to fragment 3EC at 1.316 Å is an amide bond *within* a synthetic peptide–fragment construct that LP-PDBBind extracted as the ligand. Not a protein–ligand reaction.
- **3BMO, 3JQF** (test): Cys168–DTT disulfides; DTT is a buffer reagent (`reagent_not_drug=yes`), not the tracked inhibitor. The pocket carries an unrecorded buffer adduct; the `covalent=False` label for the inhibitor is correct.

## Relation to issue #8

Seven of the 20 (1NXY, 2A4R, 2R9M, 3FYZ, 3O0U, 4QLQ, 5Z1D) are already named in `THGLab/LP-PDBBind#8`. Thirteen are new. The `false_covalent` override of 3N3G and the CL1 chemistry bias are new observations.

## Not verified here

- That CovBinderInPDB harvests only from `_struct_conn` / LINK / CONECT records.
- The cause of the four CL1 pass-throughs above.
