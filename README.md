# PDB data audits

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23441847.svg)](https://doi.org/10.5281/zenodo.23441847)

Audits of public protein-structure data for errors that carry into the models trained on it: protein–ligand covalent bonds present in the coordinates but absent from the connectivity records, benchmark labels that contradict the deposited structure, and mutations recorded as wildtype because the mutated residue was not modeled. Where the error is in a maintained dataset, it is reported to the maintainers with a reproducer.

Michael Varman · independent · Wasilla, Alaska · [ORCID 0009-0001-9453-0226](https://orcid.org/0009-0001-9453-0226) · [Morphogenesis Analytics](https://morphalytica.com)

## The three numbers

A scan of 199,603 experimental PDB entries (RCSB snapshot `20261013T233952Z`) for ligand atoms within 2.0 Å of a protein nucleophile with no `_struct_conn` record:

| | Sites | Entries | Definition |
|---|---:|---:|---|
| Scan floor | 16,827 | 2,568 | Organic ligand, protein nucleophile, ≤ 2.0 Å, no bond record. All classes, methods, resolutions. |
| Chemistry-consistent | 1,705 | 477 | Triage class 2: distance inside the element-pair bond-length window; not a known cofactor, PTM or artefact. |
| **Strict** | **485** | **289** | Class 2, X-ray, ≤ 2.5 Å. |
| …density-supported | 397 | 233 | Both partners RSCC ≥ 0.8, RSR ≤ 0.3, occupancy > 0 in the wwPDB validation report. Thresholds fixed before reading. |

The floor is honest as a scan result and should not be read as a drug-discovery claim: half of it is one bilin cofactor on cysteine, and four-fifths of it is cryo-EM at ~3 Å. The defensible number is 233 entries. More than half of the strict 289 were released since 2020.

None of the 485 strict sites is reported as a steric-clash outlier in its wwPDB validation report. The reason is in `unrecorded-covalent/clash-misclassification/`: the validation software perceives the bond from geometry, suppresses the clash accordingly, and writes nothing down. That perception stops at 2.01 Å for C–S. The bond record itself follows a similar ruler near 2.0 Å, so the same drug in the same entry is recorded at 1.981 Å and not at 2.000 Å. Past the line, an unrecorded bond is reported as a clash. Details in `unrecorded-covalent/two-rulers/`.

135 of the 485 sites (88 entries) did have the bond recorded in an earlier release. A wwPDB revision removed it, mostly on 2024-12-25, 2025-09-17 or 2025-10-01, and nothing re-filed it in any linkage category; 98 of those 135 are density-supported. Details in `unrecorded-covalent/remediation-history/`.

Removed bond records are not only cysteine. Across all 811 protein–ligand covalent sites removed by wwPDB revisions since May 2024 (archive-wide, from the Protenix training-set diff), 23 density-supported non-cysteine covalent sites on real drug targets lost their record, including a ligand bonded to the mutant Asp12 of KRAS G12D (8JHL), HIV-1 reverse transcriptase Lys102 and Tyr181, GSTP1 Tyr109 and 14-3-3σ Lys122. Details in `unrecorded-covalent/reactive-residues/`.

## Layout

```
unrecorded-covalent/
  strict-cut-491/            the 485-site table (28 columns: geometry, density, clash, UniProt) + README
  sars-cov-2/                73 SARS-CoV-2 sites / 46 entries; Mpro Cys145 66 of 66 density-supported
  remediation-history/       135 of the 485 were recorded, then removed by a 2024–25 revision and not re-filed
  lp-pdbbind-crosscheck/     20 LP-PDBBind entries that are covalent and labelled covalent=False; the mechanism
  benchmark-overlap/         PoseBusters (0 of 308), Runs N' Poses (7 entries / 23 sites)
  rnp-boltz2/                Boltz-2 on the 6 Runs N' Poses covalent entries, SMILES-only vs bond-declared
  rnp-chai1/                 Chai-1 on the same six; warhead-to-nucleophile distance measured across both models
  reactive-residues/         811 archive-removed covalent sites classified; 23 non-Cys covalent sites on drug targets
  clash-misclassification/   why zero clash outliers: probe's distance-based bond perception
  two-rulers/                the 2.01 Å validator line and the 2.0 Å record line, archive-wide; 132 density-supported Cys adducts past 2 Å
  validator-merge/           decision spec for merging resolution / occupancy / altloc gates into the triage
lp-pdbbind-label-conflicts/  11 identical-ligand twins labelled both covalent and non-covalent (issue #8, as filed)
rascore/                     KRAS switch-state table after correcting silent wildtype mislabels (issue #2)
release-audit/               weekly audit of each PDB release: heme c parity, overlapping atoms, unrecorded contacts, software defaults
priority/                    SHA-256 fingerprints of unpublished bundles, by date
```

Every table carries its source URL, snapshot or commit, and SHA-256 where a file was downloaded. Nothing here was recomputed from maps; density and clash values are the wwPDB's own.

## Public record

| Date | What | Where |
|---|---|---|
| Oct 2026 | These audit tables (dataset; v1.1.0 adds the SARS-CoV-2 table) | [Zenodo 10.5281/zenodo.23441847](https://doi.org/10.5281/zenodo.23441847) |
| Oct 2026 | Identical ligand on the same protein labelled both covalent and non-covalent, 11 entries | [THGLab/LP-PDBBind #8](https://github.com/THGLab/LP-PDBBind/issues/8) |
| Oct 2026 | Mutation status defaults to wildtype when the mutated residue isn't modeled (38 chains, mostly KRAS Q61H) | [mitch-parker/rascore #2](https://github.com/mitch-parker/rascore/issues/2) |
| Oct 2026 | *Silent Mutation Misclassification and Unrecorded Covalent Bonds in the Protein Data Bank* (preprint) | [Zenodo 10.5281/zenodo.23380304](https://doi.org/10.5281/zenodo.23380304) |

## Cite

Dataset: [10.5281/zenodo.23441847](https://doi.org/10.5281/zenodo.23441847) (all versions; resolves to the latest). Each release also has its own version DOI, listed on the Zenodo record; v1.0.0 is [10.5281/zenodo.23441848](https://doi.org/10.5281/zenodo.23441848). Machine-readable citation in `CITATION.cff`.

## License

Data and documents: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Scripts: MIT.
