# SARS-CoV-2: unrecorded covalent bonds

`sars_cov2_unrecorded_covalent_appendix.csv` lists every SARS-CoV-2 site in the strict cut: class 2, X-ray, ≤ 2.5 Å, a protein nucleophile within covalent bond length of a ligand atom, and no `_struct_conn` record. There are 73 sites in 46 entries.

It is a subset of `../strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv` (SHA-256 `05be1431…4edb7d0`), filtered on organism. `python3 -I build.py ../strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv` rebuilds it byte for byte. Density and clash values are taken from the wwPDB validation reports.

## Numbers

| Protein | Sites | Entries | Density-supported | Reading |
|---|---:|---:|---:|---|
| Mpro (nsp5), Cys145 | 66 | 41 | 66 / 66 | Covalent inhibitors on the catalytic cysteine; C–S 1.70–1.93 Å |
| PLpro (nsp3), Cys111 | 3 | 1 (9Z0D) | 3 / 3 | Irreversible inhibitor on the catalytic cysteine, three chains |
| PLpro C111S, Glu203 | 3 | 3 (13MX, 13NJ, 13NL) | 0 / 3 | Fragment screen; likely modelling artifact |
| PLpro C111S, Asp164 | 1 | 1 (13QF) | 1 / 1 | Fragment screen; chemically unexpected |
| **Total** | **73** | **46** | **70 sites / 43 entries** | |

Density-supported means both partners have RSCC ≥ 0.8, RSR ≤ 0.3 and occupancy > 0. The thresholds were fixed before any values were read.

Every Mpro site passes. These are not poorly fitted models. They are covalent inhibitors, modelled bonded to the catalytic cysteine, with no record in the file that the two atoms are bonded. Any pipeline that takes bonds from `_struct_conn` rather than from geometry will read them as non-covalent.

Four ligand series account for 49 of the 66 Mpro sites: 4WI (18), HUR (11), A1EN0 (10) and 80I (10).

## Restraint signature

Seventeen of the eighteen 4WI chains, across nine separate entries, refine to C–S distances within 1.761–1.770 Å. The remaining chain, 8IGY, is at 1.846 Å. Independent crystals do not agree to 0.009 Å, so this spread reflects the refinement restraint target rather than a measurement. The same applies across the table: a deposited bond length mostly records the link restraint the depositor supplied. These distances are useful for finding bonds and should not be used to parameterize them.

## The four PLpro C111S rows

All four come from one PanDDA fragment screen on catalytically inactive PLpro-C111S (Wang et al. 2026, *Int J Biol Macromol* 347:150689).

- **Glu203 (13MX, 13NJ, 13NL):** three unrelated fragments (A1CUM, a thiourea; XZT, a succinimide; A1CUW) each sit 1.375–1.379 Å from Glu203 OE2. None of them has an electrophile that could form an ester with a carboxylate. Ligand RSCC is 0.75, 0.64 and 0.65, and occupancy is below 1. The near-identical distance across three different fragments suggests a link restraint was active during refinement even though no bond is recorded in the deposited file.
- **Asp164 (13QF):** a bromophenol fragment (A1AW7) sits 1.378 Å from Asp164 OD1 with good density (ligand RSCC 0.93). The chemistry is still unexpected and needs inspection against the map.

For a covalent-inhibitor claim, the defensible SARS-CoV-2 set is the 42 entries with a catalytic-cysteine site: 41 Mpro entries plus 9Z0D.

## Columns

`pdb_id, protein, chain, residue, protein_atom, ligand, ligand_residue_number, ligand_atom, distance_A, resolution_A, release_date, protein_rscc, protein_rsr, ligand_rscc, ligand_rsr, ligand_occupancy, density_class, clash_class, interpretation`

`ligand` is the three- or five-character CCD code. `clash_class` is `clash_between_same_residue_and_ligand_other_atoms` for three rows (7WYM A, 8HQG B, 8HVX A). Those clashes are between other atoms, not the bonded pair.

## Not included

The strict cut contains 13 further sites in other coronaviruses: SARS-CoV and SARS-related entries (2OP9, 2VJ1, 6LNY, 6LO0, 7EO8, 7VLO, 7WQI), MERS (7XRY) and HCoV-229E (7YRZ, 8IM6).
