# Reactive non-cysteine residues hidden by deleted bond records (25 Oct 2026)

Source: archive_removed_811_sites.csv (21 Oct: 811 protein–ligand covalent sites in 524 entries whose
_struct_conn record was removed by wwPDB revisions after 2024-05-22). Re-read today from current RCSB
mmCIF + validation XML for all 524 entries.

## Classification of the 811 (804 located; 7 atoms not found)
- Chemistry: 271 plausible bonds (protein S/O/N/Se to ligand C/B/P/S within covalent radii + 0.4 Å),
  244 heteroatom pairs that are not bonds (removal correct), 113 C–C to inspect, 176 other.
- Ligand type: 439 drug-like/organic ligands, 341 nucleotides/cofactors, 16 additives, 8 sugars.
- Density (both partners RSCC >= 0.8, RSR <= 0.3): 298 supported, 122 not, 384 unknown
  (cryo-EM / NMR / old X-ray without structure factors).

Plausible bond AND organic ligand: 198 sites. By residue: Cys 91, Thr 19, Ser 15, Lys 15, Glu 12,
His 12, Tyr 9, Asp 6, Pro 4, Gly 4, other 11. Density-supported: Cys 60, Lys 14, Ser 9, Tyr 8,
Thr 6, His 4, Pro 4, Asp 3, Glu 3.

## Hand-picked: 23 density-supported non-cysteine covalent sites on real targets
(non_cys_covalent_targets.csv). Highlights:
- KRAS G12D 8JHL — ligand bonded to the mutant Asp12 itself, 1.44 Å, density-supported; record removed 2025-09-24.
- HIV-1 RT Lys102 (8U6G/8U6H/8U6T, acrylamide NNRTIs) and Tyr181 (7KRC/7KRD/7KRF, sulfur(VI) warheads).
- GSTP1 Tyr109 (5X79), 14-3-3σ Lys122 aldehyde fragments (8B2I, 8B5P), MetAP2 His231 (5D6F),
  OAT Lys292/Thr322 (6V8C, 7LNM), MIF N-terminal Pro (4OYQ, 4PKK), SARM1 Glu642 (7NAJ).
- Classic enzymology also present (serine proteases, beta-lactamases, proteasome Thr1, HaloTag-family Asp).

## Honest framing
- The chemistry in these papers is mostly published; what is lost is the archive record. Anyone searching
  the PDB for covalent Lys/Tyr/Asp ligands, and every model trained on current files, no longer sees them.
- Excluded on purpose: cryo-EM lipid/detergent "bonds" at ~1.37 Å (modelling artifacts), additives
  (LAC, DMF, CO2, CO3), lipidated peptides, nucleotides.
- Selection of the 23 is by hand from the 198; the full classified table is removed_811_classified.csv.

## Open
- Q-score check for the 281 cryo-EM sites; inspect the 113 C–C pairs.
- Extend from removed records to never-recorded sites (the 485 strict cut) for the same non-Cys map.
