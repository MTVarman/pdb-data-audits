# verify/

One script that re-derives the headline numbers in the top-level `README.md`
from the tables in this repository.

```
python3 verify/verify_claims.py
```

Standard library only. No network, no arguments, no install. Exit status is 1
if a number that should re-derive does not, and 0 otherwise.

## What it checks

| Claim | Source |
|---|---|
| strict cut, 485 sites / 289 entries | `strict-cut-491/…_491_sites.csv` |
| strict cut matches its own definition (class 2, X-ray, ≤ 2.5 Å) | same |
| density-supported, 397 sites / 233 entries | same, recomputed from RSCC / RSR / occupancy |
| `density_class` column agrees with that rule | same |
| no strict site is a clash outlier for its own bonded pair | same |
| more than half the strict entries released since 2020 | same |
| 135 sites / 88 entries had the bond in an earlier release | `remediation-history/deleted_rows_strict.csv` |
| 98 of those are density-supported | both of the above |
| 804 of the 811 archive-removed sites located | `reactive-residues/removed_811_classified.csv` |
| 23 density-supported non-cysteine sites, KRAS G12D 8JHL among them | `reactive-residues/non_cys_covalent_targets.csv` |

Density support is recomputed from the raw wwPDB validation columns rather
than read from the precomputed `density_class` column, and the column is then
checked against the recomputation.

## What it cannot check, and why

Two rows of the README table are counts over the whole scan, and the scan
output is not in this repository:

- **scan floor, 16,827 sites / 2,568 entries** — needs the full output of the
  scan over 199,603 entries.
- **chemistry-consistent, 1,705 sites / 477 entries** — needs the triage
  class 2 table as it stood before the X-ray and ≤ 2.5 Å cuts.

The script reports both by name instead of failing on them, and reports the
same way for any other input that is absent. A missing file is a stated gap,
never a crash: every remaining check still runs and the report is always
complete.

Two wording notes it surfaces rather than hides:

- The strict table is still named `…491…` from an earlier count; it holds 485
  rows, which is the number the README now uses.
- The 811 archive-wide total is itself one level up from this repository:
  `archive_removed_811_sites.csv` is not shipped, and `reactive-residues/README.md`
  records that 7 of the 811 atoms were not found on re-read, leaving the 804
  rows that are here.
