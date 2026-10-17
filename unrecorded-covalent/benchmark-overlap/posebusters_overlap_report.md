# PoseBusters (degrado-lab list) vs unrecorded covalent scan

## Exact pair hits (benchmark's own ligand carries the unrecorded bond): **0**

- known cofactor or PTM: 0
- bond-length consistent adduct: 0
- likely artefact: 0

No benchmark (PDB, ligand) pair is a headline unrecorded covalent site.

## Same entry, different ligand (weaker): **2**

- known cofactor or PTM: 0
- bond-length consistent adduct: 0
- likely artefact: 2

| entry | ligand | class | chain/residue/atom | ligand atom | distance Å | method | resolution Å | year |
|---|---|---|---|---|---|---|---|---|
| 6M73 | LAC | likely artefact | A TYR211 OH | CA | 1.389 | X-ray | 1.7 | 2021 |
| 7C0U | EDO | likely artefact | B ASP336 OD1 | O1 | 1.857 | X-ray | 1.8 | 2021 |

Benchmark ligand for 6M73 is FNR; for 7C0U it is FGO. Both unrecorded contacts are
crystallisation additives (lactic acid, ethylene glycol), flagged by triage rule 3B.

## Method

Parsed 308 unique pairs (308 entries), uppercased, from the repository file
(SHA-256 `19d8b77b9f0ab9edf1703c6544bb88ca22c216a0c31c8b2cd015e574f6e26f16`,
source https://raw.githubusercontent.com/degrado-lab/posebusters-benchmark/main/posebuster_benchmark.txt).
Matched on entry + ligand component code against the 16,827 headline sites
(organic ligand, non-metal ligand atom, protein nucleophile, <=2.0 Å, no struct_conn row)
in the triage table; class is from the same table. Matching is by component code, so a hit
is not proof of the same ligand instance. The scan covers headline sites only; sub-threshold
or non-nucleophile contacts are not included.
