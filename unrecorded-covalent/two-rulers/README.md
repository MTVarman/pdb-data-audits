# Two rulers at 2.0 Å: why covalent bonds modelled a hair long vanish from the record

**Result.** Whether the PDB records a cysteine covalent bond, and whether its validation report calls it a bond or a clash, both depend on a distance line near 2.0 Å. Chemistry plays no part. The same drug in the same entry is recorded in one copy at 1.981 Å and left unrecorded in the next at 2.000 Å.

**Data.** Every entry in the PDB, 261,583 of 261,584, downloaded once from RCSB on 31 Oct 2026 (holdings retrieved 08:12 UTC). The one entry left out is 36ZA (11.2 million atoms), which exceeded available memory. Each download's URL, time and SHA-256 are logged. Contacts are ligand atoms within 2.5 Å of a protein nucleophile, read from model 1 of the asymmetric unit, using the same nucleophile list as `../strict-cut-491/`.

---

## Ruler 1: the validator's bond line (2.01 Å)

The wwPDB clash list comes from MolProbity's `probe`. Probe treats a pair as bonded when its distance is within the two covalent radii plus a 0.2 Å allowance (`probe.c`, `markBonds()`, `COVRADFUDGE` 0.2). For carbon–sulfur that is 0.77 + 1.04 + 0.2 = **2.01 Å**.

**Pre-registered prediction.** Unrecorded Cys SG–ligand C pairs at or under 2.01 Å would almost never appear in the clash list; pairs over 2.01 Å almost always would.

**Test.** All such pairs in X-ray entries, both partners at full occupancy, closest ligand atom per cysteine, read against each entry's own wwPDB validation report.

| | ≤ 2.01 Å | > 2.01 Å |
|---|---:|---:|
| **Unrecorded; no bond recorded anywhere between that Cys and that ligand** | **1 / 496** clash-listed | **1,606 / 1,686** (95%) |
| Unrecorded, but the ligand is recorded as bonded to that Cys through another atom | 0 / 32 | 9 / 1,393 |
| Recorded (`_struct_conn` covale), sampled at 1.95–2.50 Å | 0 / 506 | 3 / 101 |

Below the line, every 0.05 Å bin is at 0%. From 2.05 to 2.50 Å, every bin is at 91–98%. The 2.00–2.05 bin straddles the line and reads 69% (`probe_line_bins.csv`).

**What this settles.** This closes the open items in `../clash-misclassification/`:
- The validator treats a pair as bonded if it is **recorded**, or if it falls **within the covalent radii + 0.2 Å**.
- The production pipeline does read bond records: recorded long bonds are almost never flagged.
- Its geometric perception covers unrecorded C–S bonds only up to 2.01 Å.

Cofactors behave the same way. Unrecorded heme-c and bilin attachments past 2.01 Å are clash-listed 1,414 of 1,482 times.

## Ruler 2: the record line (≈ 2.0 Å)

**Question.** When one entry holds the same covalent bond in several copies (the same ligand atom to the same Cys in different chains), is the bond recorded in every copy?

**Test.** Groups are the same ligand component and ligand atom to Cys SG within one entry, full occupancy, no crossing alternate conformers, 1.6–2.5 Å. A group counts as mixed when at least one copy has a `covale` record and at least one has none, with no record on that residue pair either.

| | entries | recorded copies < 2.0 Å | unrecorded copies ≥ 2.0 Å |
|---|---:|---:|---:|
| Ligands (inhibitors, drugs, substrates) | 27 | 47 / 51 | **40 / 40** (closest at 2.000 Å) |
| Cofactors (heme c, bilins, flavins…) | 186 | 8,264 / 8,285 | 1,709 / 1,712 |

Within one entry, the chemistry, the depositor and the experiment are all the same. The copy modelled under 2.0 Å gets the bond record and the copy at 2.0 Å or over does not.

| Entry | Target | Recorded | Not recorded |
|---|---|---|---|
| 7VJL | FGFR4 | B 1.981 Å | **A 2.000 Å** |
| 9PZF | KRAS G12C | A 1.943 | B 2.003 |
| 9Q62 | PREP | C 1.866, B 1.957 | A 2.003 |
| 9EW5 | 14-3-3σ | A 1.942 | D 2.011 |
| 9Q65 | PREP | B 1.953, F 1.957 | D 2.014 |
| 7JXH | HER2 | six chains, 1.772–1.977 | G 2.083 |
| 7MAT | SARS-CoV-2 Mpro | B 1.766 | A 2.006 |

All 275 groups in 211 entries are in `same_bond_recorded_and_unrecorded.csv`.

## Together

A cysteine adduct modelled at 2.0–2.5 Å fails both rulers:
- It is not written into `_struct_conn`.
- The validation report files it as a severe steric clash, with overlaps above 1 Å.

Every downstream consumer reads an unbonded ligand. A dataset that filters clashes by distance either drops the entry or keeps it as non-covalent.

Archive-wide there are 353 such cysteine–ligand pairs in 216 entries, after excluding cofactors, reducing agents and pairs already bonded through another atom. **132 of those entries** are X-ray ≤ 2.5 Å with density support on both partners. 126 of them are not among the 233 density-supported entries in `../strict-cut-491/`, which stops at 2.0 Å. The list mixes covalent inhibitors with enzyme–substrate intermediates (IMP in IMPDH, CoA thioesters); the split by ligand category is still open.

## Fix

**Ruler 2 has a one-line fix.** If a bond is recorded in one copy of an entry, record it in every equivalent copy. No chemistry judgement is needed. **For Ruler 1,** a clash between a nucleophile and a reactive ligand atom at 2.0–2.5 Å is a candidate missing bond and should be reported as one.

## Established vs. not

**Established:**
- probe's distance rule, read from source
- the 2.01 Å split in clash listing across the archive
- record-honouring behaviour (recorded long bonds 3 of 101 flagged)
- the within-entry record split at 2.0 Å

**Not established:**
- the annotation step's exact rule. Ruler 2 is inferred from the records, not read from annotation software.
- why bonds end up modelled past 2.0 Å. One plausible loop: with no link restraint, refinement treats the pair as non-bonded and pushes it apart, which then keeps it from being recorded. Not tested.
- other element pairs, where probe's lines are C–O 1.63 Å and C–N 1.67 Å.

**Exceptions:** 4 recorded ligand copies sit over 2.0 Å (2C6Q 2.080; 2VK8 2.026 and 2.052; 6MD2 2.144). Not investigated.

## Files

| File | Rows |
|---|---:|
| `probe_line_bins.csv` | 0.05 Å bins, three arms |
| `probe_line_unrecorded_pairs.csv` | 3,612 unrecorded Cys–C pairs with clash, density, family |
| `probe_line_recorded_pairs.csv` | 607 recorded pairs, 1.95–2.50 Å |
| `same_bond_recorded_and_unrecorded.csv` | 275 groups, 211 entries |
| `scripts/` | `archive_index.py` (whole-PDB index), `archive_triage.py`, `archive_evidence.py`, `record_cliff.py` |

The full index is 32 shards (506 MB compressed). It is too large for this repository and is available on request.
