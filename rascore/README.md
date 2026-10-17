# Rascore within-class comparison, corrected labels

`TableS_rascore_corrected_within_class.csv` — mutant-vs-wildtype switch-state frequencies for KRAS after correcting the mutation labels that Rascore assigns as wildtype when the mutated residue is not modeled (reported as `mitch-parker/rascore#2`).

## Columns

| column | meaning |
|---|---|
| `nucleotide_class` | `2P` (GDP-bound), `3P` (GTP-analogue-bound), or `all` |
| `variant`, `switch`, `state` | KRAS mutant, switch I or II, Rascore conformational state |
| `mutant_n`, `mutant_pct`, `wt_n`, `wt_pct`, `delta` | counts and percentages after label correction |
| `p_corrected` | Fisher exact p, corrected labels |
| `p_labeled`, `delta_labeled` | the same test using Rascore's original labels |
| `bh_rank`, `bh_n_tests`, `bh_threshold`, `bh_pass` | Benjamini–Hochberg at FDR 0.05 *within each nucleotide class* (34, 29 and 60 tests) |
| `admissible` | `yes` only for within-class rows; `all`-class rows are reported but not admissible because they mix nucleotide states |
| `exclusion_reason` | e.g. `no_2P_structures` where a variant has no structures in that class |
| `ligand_defined_state` | the state is defined by a bound ligand, not by the switch alone |
| `significance_change` | `gained` / `lost` relative to the original labels |
| `sign_flip` | the direction of `delta` reversed after correction |

## What survives correction

Within-class BH pass: three 2P rows (G12C SW2 Y71in.2P-SP12, G12D SW2 Y71out.2P-SP2-A, G12C SW2 Y71out.2P-SP2-A) and two 3P rows (Q61H SW1 Outlier, Q61H SW2 Y71in.3P-R). Everything else fails the within-class threshold, including several comparisons that were nominally significant under the original labels (`significance_change = lost`).

Rows in the `all` class are shown for completeness and marked `admissible = no`.
