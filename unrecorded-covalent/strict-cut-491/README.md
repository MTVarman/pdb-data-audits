# Unrecorded covalent contacts: tight cut, 491 sites in 289 PDB entries

`unrecorded_covalent_tight_cut_491_sites.csv` has one row per flagged contact between a protein residue and a ligand, with density, clash and protein-identity information added. UTF-8 CSV, 28 columns, 491 rows.

## What the cut is

A scan of 199,603 experimental PDB entries (snapshot 2026-10-13) looked for ligand atoms lying within 2.0 Å of a protein nucleophile (the scan's list of nucleophilic residue atoms) where the entry's own bond records (`_struct_conn`) list no covalent bond. The ligand had to be organic and the contacted ligand atom could not be a metal. Each site was then triaged by rule:

- Class 1: known cofactor or post-translational attachment.
- Class 2: the distance falls in the bond-length window for that element pair (for example C–O 1.3–1.55 Å, C–S 1.7–1.95 Å) and the ligand does not match the known-cofactor rules (class 1) or the nucleotide-like and other artefact rules (class 3).
- Class 3: likely artefact.

**This file holds class 2 only, restricted to X-ray structures at 2.5 Å resolution or better: 491 sites in 289 entries.** Class 2 means only that the distance fits a bond. It does not show that a bond exists.

## Density validation

Fit values are the wwPDB validation report values for each entry, computed by the wwPDB against the deposited map. Nothing was recomputed. Thresholds were fixed before any results were read. A residue or ligand is *well fitted* when real-space correlation (RSCC) ≥ 0.8, real-space R (RSR) ≤ 0.3 and occupancy > 0. If a partner has several report rows (alternate conformers), the deposited conformer is used, otherwise the row with the lowest RSCC.

- `density_supported`: both partners well fitted. **403 of 491 sites, in 233 of 289 entries.**
- `ligand_uncertain` (protein well fitted, ligand not): 56 sites.
- `protein_uncertain` (ligand well fitted, protein not): 5 sites.
- `unsupported` (neither well fitted): 6 sites.
- `unknown` (report missing, or no density values for a partner): 21 sites.

## Clash check

The wwPDB reports were searched for a steric-clash outlier at the exact protein-atom/ligand-atom pair of each site. **None of the 491 flagged atom pairs is recorded as a clash outlier.** Qualifications:

- At 45 sites the report records a clash between the same residue and ligand at *other* atoms.
- At 6 sites the report holds no clash data.

`clash_class` gives the result per site.

## Worked example: 6T5V

Site `6T5V|A|CYS12|SG|A|MKZ203|C`: Cys12 SG to ligand MKZ203 C, 1.816 Å, released 2020-02-19, X-ray at 1.31 Å. The distance fits a C–S bond.

| | RSCC | RSR | occupancy |
|---|---|---|---|
| protein Cys12 | 0.980 | 0.049 | 1.000 |
| ligand MKZ203 | 0.951 | 0.068 | 1.000 |

Both partners are well fitted, so the class is `density_supported`. The validation report for this entry (clashscore 0.72) records no clash at the Cys12 SG–MKZ C pair, so `clash_class` is `no_clash_recorded_for_pair`. Protein: GTPase KRas (P01116), *Homo sapiens*.

## What good fit does and does not establish

Good real-space fit shows that the atoms are where the deposited density says they are. It does not by itself prove a covalent bond. A well-fitted residue next to a well-fitted ligand at 1.3–2.0 Å is consistent with a bond, but density at this resolution cannot always tell a bond from a close non-bonded contact, and refinement choices can fit a model to density. A poor fit does not prove there is no bond. It means the coordinates do not support the claim strongly enough to rely on. Fit values are averaged over a residue or ligand, so a well-fitted ligand with one badly fitted atom at the contact still counts as well fitted. `ligand_category` is a heuristic based on the component's identity (pharmacological, cofactor_or_metabolite, additive_or_reagent), not a chemical assessment.

## Provenance (the per-row source columns were left out)

- **Scan:** RCSB PDB mmCIF files and metadata, snapshot 2026-10-13 (ID 20261013T233952Z). Every download was logged with URL, retrieval time and SHA-256.
- **Validation:** the wwPDB validation report XML for each entry, as served by PDBe, retrieved 2026-10-17 between 09:50 and 09:52 UTC (report schema v6.06 or v6.07). One layer holds the real-space fit values and the other holds the clash records.
- **Protein identity:** the UniProt accession comes from the entry's own mmCIF reference records (the RCSB/SIFTS alignment is used where the mmCIF accession is inactive in UniProt). Names come from the UniProt REST service, retrieved 2026-10-17.
- Per-row URLs, hashes and retrieval times are kept in the full project tables and can be supplied on request.

## Notes on reading the file

- `site_id` joins the entry, chain, residue, atoms and ligand with `|`. It is not unique for three pairs of rows (5BPK chains C and D, 6G3Z chains A and B, 6P4E chains A and B). Within each pair the two rows are identical in every source column. They are kept so that the totals stay at 491.
- 12 sites in 7 entries have a blank `uniprot_accession`. UniProt holds no sequence for that chain: peptides, designed proteins or unreferenced sequences.
- Blank fit values occur only at the 21 `unknown` sites.
- The SARS-CoV-2 replicase polyprotein accessions account for 46 of the 289 entries (about 16%), mostly the main protease. Otherwise the proteins are widely spread (159 distinct accessions).
