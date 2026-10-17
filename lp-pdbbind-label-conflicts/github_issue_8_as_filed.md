**Filed:** https://github.com/THGLab/LP-PDBBind/issues/8 (16 Oct 2026, from MTVarman)

**Title:** Identical ligand on the same protein labelled both covalent=True and covalent=False (11 entries)

---

Thanks for LP-PDBBind. The leakage-aware splits are the reason we use it as a reference set.

While auditing the `covalent` column I found label conflicts that are visible from `dataset/LP_PDBBind.csv` alone.

**Same ligand, same protein, opposite labels**

Eleven `covalent=False` entries have a twin with the identical `smiles` string and a ≥95% similar `seq` that is labelled `covalent=True`. Ten are in `train` and one is in `test`.

Six of the pairs have no amino-acid substitution between them. Three have identical sequences. The other three differ only by a few inserted or missing residues:

| covalent=False | split | covalent=True twin | protein difference | ligand |
|---|---|---|---|---|
| 3w07 | train | 2zz6 | identical | UMP |
| 3fyz | train | 3mbz | identical | ring-opened β-lactamase inhibitor |
| 4ob1 | train | 4ob2 | identical | butylboronic acid |
| 6gcx | train | 6gcr | indels only | anilinopyrimidine–squaramide |
| 6mzq | test | 6mzw | indels only | acrylamide kinase inhibitor |
| 2za3 | train | 3bar | indels only | UMP |

For these six, either the bond genuinely differs between two crystals of the same complex, which would be worth recording, or one of the two labels is wrong.

Five more pairs differ by one to three substitutions. A substitution at the nucleophile could legitimately stop the bond from forming, so these need a check rather than a fix:

| covalent=False | split | covalent=True twin | substitutions | ligand |
|---|---|---|---|---|
| 3ixg | train | 1c3b | 1 | benzothiophene-2-boronic acid |
| 2qcf | train | 3g3m | 2 | 5-fluoro-UMP |
| 3sv9 | train | 3sv6 | 1 | α-ketoamide |
| 1loq | train | 2zz6 | 1 | UMP |
| 3qzq | train | 4ght | 3 | Michael-acceptor ester, drawn saturated |

The exact-string version takes three lines and returns the first three pairs:

```python
import pandas as pd
df = pd.read_csv("dataset/LP_PDBBind.csv").dropna(subset=["smiles", "seq"])
n = df.groupby(["smiles", "seq"])["covalent"].nunique()
print(n[n > 1])
```

The remaining eight come from comparing sequences within each SMILES group (difflib, `autojunk=False`).

**Wider tier, needs case-by-case review.** 21 distinct SMILES carry both labels across 86 entries (37 False, 49 True). The 26 not listed above are the same ligand on a different protein, where differing labels can be legitimate. The full list is below. One pattern stands out: of 116 entries whose ligand contains a boronic acid, 104 are `covalent=True`.

<details>
<summary>All 37 False-labelled entries that share a SMILES with a True-labelled entry</summary>

| covalent=False | split | closest covalent=True twin | seq similarity | substitutions |
|---|---|---|---|---|
| 3w07 | train | 2zz6 | 1.000 | 0 |
| 3fyz | train | 3mbz | 1.000 | 0 |
| 4ob1 | train | 4ob2 | 1.000 | 0 |
| 3ixg | train | 1c3b | 0.997 | 1 |
| 2qcf | train | 3g3m | 0.992 | 2 |
| 6gcx | train | 6gcr | 0.992 | 0 |
| 6mzq | test | 6mzw | 0.991 | 0 |
| 2za3 | train | 3bar | 0.989 | 0 |
| 3sv9 | train | 3sv6 | 0.988 | 1 |
| 1loq | train | 2zz6 | 0.981 | 1 |
| 3qzq | train | 4ght | 0.969 | 3 |
| 3mg0 | train | 4qvm | 0.749 | — |
| 2qcd | train | 2qcn | 0.665 | — |
| 3g1d | train | 2zz6 | 0.660 | — |
| 6gxg | train | 6gxy | 0.656 | — |
| 3bcn | val | 3c9e | 0.501 | — |
| 2pre | test | 3c9e | 0.482 | — |
| 3ibc | train | 3gjs | 0.386 | — |
| 1nxy | train | 1ym1 | 0.376 | — |
| 3mke | train | 1yly | 0.333 | — |
| 3g1v | train | 3g3m | 0.264 | — |
| 2reg | train | 6nd3 | 0.175 | — |
| 1oba | train | 6nd3 | 0.170 | — |
| 3ppq | train | 6nd3 | 0.160 | — |
| 3r6u | train | 6nd3 | 0.158 | — |
| 5km3 | train | 2qcn | 0.153 | — |
| 1xz8 | train | 2zz6 | 0.145 | — |
| 3dxg | train | 3bar | 0.138 | — |
| 6b22 | train | 6b1h | 0.121 | — |
| 3gjt | train | 1iau | 0.113 | — |
| 5zbz | val | 3lle | 0.091 | — |
| 3as0 | val | 3lle | 0.082 | — |
| 2j4k | train | 3bar | 0.082 | — |
| 3arv | val | 3lle | 0.081 | — |
| 2clk | train | 1uod | 0.080 | — |
| 2ha3 | train | 6nd3 | 0.074 | — |
| 2prj | test | 5d9p | 0.058 | — |

</details>

**Related to #7.** @roccomoretti's note that 2jg8's SMILES shows an aldehyde where a carboxylate is expected looks like the same failure mode at single-entry level. Separately, an archive-wide scan I ran for protein–ligand bonds present in the coordinates but absent from `_struct_conn` flags 18 `covalent=False` entries here. In several, the `smiles` field already shows the warhead in reacted form: nitrile as C=N (2r9m, 3o0u), acrylamide saturated (5z1d), epoxyketone opened to a diol (4qlq), ketoamide as the hemiketal (2a4r). I'm happy to open that separately with residues, atoms and distances.

Method and full scan: https://doi.org/10.5281/zenodo.23380304
