# Remediation history of the 491 strict sites

What happened to the bond record at each of the 491 strict-cut sites (289 entries, see `../strict-cut-491/`). In the current files none of them has a `_struct_conn` row for the bond. This folder checks whether the bond was ever recorded in an earlier release of the same entry and, if it was, what removed it.

## Sources

- **Current files:** RCSB snapshot `20261013T233952Z`, the same files the scan used.
- **wwPDB versioned archive** (`files-versioned.wwpdb.org`). It keeps only the latest minor version of each major version. 287 requests, 160 mmCIF files, every one with URL, timestamp, bytes and SHA-256: `downloads_versioned.jsonl`.
- **PDBj Snapshot Archive** (`snapshots.pdbj.org`, 34 snapshots 2005–2026). For each entry, the latest snapshot dated after initial release and before the first revision that lists `struct_conn`. 107 requests, 106 returned 200 (one 404, 2IWC at 20180321; the 20180101 snapshot was used instead): `downloads_v2.jsonl`.
- **Revision history:** each file's own `pdbx_audit_revision_history`, `_group`, `_category`, `_item` and `_details`.
- **Density:** the wwPDB validation values already used for `../strict-cut-491/` (both partners RSCC ≥ 0.8, RSR ≤ 0.3, occupancy > 0). Nothing was recomputed.

Matching is strict and on covale rows only: author chain, label asym, residue number, residue name, protein atom, ligand component, ligand atom and identity symmetry must all agree. A *sister* match is the same bond on another copy of the same entity.

## Result

| Class | Sites | Entries | Density-supported sites | Entries with a density-supported site |
|---|---:|---:|---:|---:|
| Never recorded; `struct_conn` never revised since release | 286 | 162 | 246 | 136 |
| **Recorded, then removed by a revision** | **141** | **88** | **104** | **67** |
| Revised; bond absent from every pre-revision file retrieved | 61 | 38 | 51 | 30 |
| Undeterminable (no pre-revision file anywhere) | 3 | 3 | 2 | 2 |
| Recorded now | 0 | 0 | | |
| All | 491 | 289 | 403 | 233 |

An entry with sites in more than one class is counted in each row, so the entry column does not sum to 289. Sister classes are empty: no site was recorded only on another copy of its entity, in any file.

"Removed" means the row is present in an earlier retrieved file, absent from every later one, and a revision listing `struct_conn` lies between the two. A revision that touches `struct_conn` is not by itself proof that a row was removed; the removal is read from the files on either side of it.

## The 141 removed bonds were not re-filed

In each current file, the removed bond was looked for in `pdbx_modification_feature`, in any `struct_conn` row of any type between the same two residues, in `pdbx_entity_branch_link`, in `pdbx_struct_mod_residue`, and as a residue merge or component rename.

| Result | Sites |
|---|---:|
| Re-filed in any linkage category | 0 |
| Removed, not represented elsewhere | 141 |

In all 141, the protein residue and the ligand are unchanged (same residue name at the same number; same ligand component). In all 141, the same two residues now appear as a `pdbx_validate_close_contact` pair, which is a validation listing, not a linkage record. Per site: `refiled_check.csv`.

## Which links the revisions removed

Every covale row in the 88 pre-revision files, kept or removed in the current file:

| Partner type | Kept | Removed | Removed that are flagged bonds |
|---|---:|---:|---:|
| Protein–small-molecule ligand | 10 | 166 | 141 |
| Protein–protein non-standard linkage | 47 | 6 | 0 |
| Protein–glycan | 6 | 0 | 0 |
| Glycan–glycan | 3 | 0 | 0 |
| Ligand–ligand | 3 | 3 | 0 |
| All | 69 | 175 | 141 |

No removed row reappears under another connection type. These 88 entries were selected because they lost a protein–ligand record, so 166 of 176 describes this set, not the archive. Within the same files, glycan links were kept 9 of 9 and protein–protein linkages 47 of 53. Row level: `covale_rows_detail.csv`.

## When

For 80 of the 88 entries, the removal window contains a revision on one of three dates: **2024-12-25 (34 entries), 2025-10-01 (28), 2025-09-17 (18)**. 72 of the 88 windows contain only minor-version revisions; 16 contain a major-version (x-0) release. On all three dates the revisions list `struct_conn`, `pdbx_validate_close_contact` and `pdbx_entry_details` (most often `has_protein_modification`), and some also list `pdbx_modification_feature`. 27 of the 80 revisions on those dates list no items at all. Across the 88 entries, `pdbx_audit_revision_details` carries one description beyond "Initial release": "Carbohydrate remediation" (1A65, 2020-07-29). None of the three dates has a free-text description.

The dates line up with the wwPDB protein-modification remediation. The [wwPDB remediation page](https://www.wwpdb.org/documentation/remediation) says 76,821 entries were updated with new PTM/PCM annotation in October 2024, and that about 12,900 entries would be updated with "link correction and residue split or merge" later in 2025.

Examples in the removed set: 6T5V (KRAS Cys12–MKZ; recorded in the 2024-01-01 snapshot, removed by the revision of 2024-05-01 or of 2025-10-01; both list `struct_conn` and no file between them survives), 7RP3 and 7RP4 (same ligand), and the SARS-CoV-2 Mpro Cys145 entries 7B3E, 7DPP and 7WYM.

## Files

| File | Rows | Contents |
|---|---:|---|
| `density_crosstab.csv` | 5 | The result table above, with the full density breakdown per class |
| `refiled_check.csv` | 141 | Each removed site, checked in the current file for every linkage category. Six `site_id` values appear twice, as they do in the strict-cut table |
| `covale_kept_removed_by_type.csv` | 6 | Kept and removed covale rows by partner type, 88 pre-revision files |
| `covale_rows_detail.csv` | 244 | The same, one row per covale row |
| `deleted_rows_strict.csv` | 221 | Every covale row present in an earlier retrieved file and absent from the current file, across all entries checked, with the revisions in the removal window and the strict match to a flagged site (exact 135, sister 9, none 77) |
| `entries_snapshot_check.csv` | 109 | Per entry that needed the PDBj snapshot archive: snapshots in the window, the snapshot used, its revision ordinal, `struct_conn` row counts and the exact/sister/none result |
| `downloads_versioned.jsonl` | 287 | wwPDB versioned-archive requests with SHA-256 |
| `downloads_v2.jsonl` | 107 | PDBj snapshot requests with SHA-256 |
| `split/sites_strict_match.part01–02.csv` | 491 | One row per strict-cut site: the class above, the files checked, the strict match in each, and other records on the same residue and ligand. Split in two for size; the header is repeated in each part |
| `split/revision_items.part01–03.csv` | 828 | Every revision item (`pdbx_audit_revision_item`), category and group listed after release for the 88 removed entries, with whether the revision falls in the removal window. Split in three, header repeated |
| `SHA256SUMS` | | SHA-256 of every table here. Each matches the hash recorded for the same file in the run directory that produced it |

## Notes

- The first pass of this audit misread quoted item names in `pdbx_audit_revision_item` (for example `'_struct_conn.pdbx_dist_value'`) as a one-row `_struct_conn`, which made 40 entries' row counts one too high. Revision dates and the set of files fetched were not affected. These tables come from the corrected parser.
- 6WVQ chain A Ser203 OG–UCY C10 is unrecorded. Chain B records Ser203 OG–UCY **P09**, a different ligand atom, so under the atom-strict rule it is not counted as the same bond.

## Open

- An archive-wide count: the same snapshot-against-current comparison over all entries, counting protein–ligand covale rows removed by the 2024-12-25, 2025-09-17 and 2025-10-01 revisions.
