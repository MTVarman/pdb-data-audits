# LP-PDBBind covalent-label audit — 16 Oct 2026

**Filed as GitHub issue #8:** https://github.com/THGLab/LP-PDBBind/issues/8 (from MTVarman, 16 Oct 2026, tags Rocco Moretti and references #7)

Checked the `covalent` column of LP-PDBBind (THGLab; 19,443 entries, 18,550 labelled False and 893 True) two ways: against itself, and against our archive-wide scan for unrecorded protein–ligand covalent bonds.

## Finding 1 — same ligand, same protein, opposite labels

Eleven entries labelled `covalent=False` have a twin with the identical SMILES and a ≥95% similar protein sequence labelled `covalent=True`. Ten are in **train** and one is in **test**.

**Six pairs have no amino-acid substitution at all.** Three are byte-identical; three differ only by a few inserted or missing residues:

| covalent=False | split | covalent=True twin | protein difference | ligand |
|---|---|---|---|---|
| 3w07 | train | 2zz6 | identical | UMP |
| 3fyz | train | 3mbz | identical | ring-opened β-lactamase inhibitor |
| 4ob1 | train | 4ob2 | identical | butylboronic acid |
| 6gcx | train | 6gcr | indels only | anilinopyrimidine–squaramide |
| 6mzq | test | 6mzw | indels only | acrylamide kinase inhibitor |
| 2za3 | train | 3bar | indels only | UMP |

For these, either the bond genuinely differs between two crystals of the same complex, or one label is wrong.

**Five pairs differ by 1–3 substitutions.** A substitution at the nucleophile could legitimately stop the bond forming, so these were filed as "check," not "wrong": 3ixg/1c3b (1, benzothiophene-2-boronic acid), 2qcf/3g3m (2, 5-fluoro-UMP), 3sv9/3sv6 (1, α-ketoamide), 1loq/2zz6 (1, UMP), 3qzq/4ght (3, Michael-acceptor ester drawn saturated).

3fyz is also one of the 18 flagged by the structural scan (Finding 3), found independently.

Method note: sequence similarity uses Python difflib with `autojunk=False`. With the default setting, difflib treats common amino acids as junk on sequences over 200 residues and badly underestimates similarity; the first pass used the default and found only 6.

Exact-string reproducer (returns the three byte-identical pairs):

```python
import pandas as pd
df = pd.read_csv("dataset/LP_PDBBind.csv").dropna(subset=["smiles", "seq"])
n = df.groupby(["smiles", "seq"])["covalent"].nunique()
print(n[n > 1])
```

## Finding 2 — same ligand, different protein (needs case-by-case review)

21 distinct SMILES carry both labels across 86 entries (37 False, 49 True). The 26 False members not in Finding 1 are the same ligand on a different protein, where differing labels can be legitimate.

Boronic acids stand out: 116 entries contain one, 104 labelled True and 12 False. `1nxy` has the identical SMILES of `1mxo`, `5w13`, `1ym1` and `1pi5` (all True), but those are different β-lactamases (≤38% similar), so it belongs in this tier, not Finding 1.

All 37 are in `lp_pdbbind_label_conflicts.csv` with closest twin, similarity, substitution count and `tier` (`no_substitution`, `few_substitutions`, `same_ligand_other_protein`).

## Finding 3 — structural scan (unrecorded bonds in the coordinates)

18 `covalent=False` entries have a protein–ligand bond to a catalytic residue that isn't in the deposited connectivity records. Split: train 14, val 2, test 2.

- train: 1nxy, 1nym, 2a4r, 2vj1, 3fv7, 3fyz, 3fzc, 4qlq, 4yv2, 5ykp, 5yr4, 5yr5, 5yr6, 6nzg
- val: 2r9m, 3n3g
- test: 3o0u, 5z1d

In several, LP-PDBBind's own SMILES already draws the warhead after it reacted: nitrile as C=N (2r9m, 3o0u), acrylamide saturated (5z1d), epoxyketone opened to a diol (4qlq), ketoamide as the hemiketal (2a4r).

3n3g is a cyanopyridine (a real cysteine-protease warhead) whose SMILES shows the nitrile intact. It stays in the 18 but isn't one of the reacted-SMILES examples.

Set aside on chemistry review: 5yjo and 6ijl (plausible, not confirmed); 3kig, 4u6w, 3aav, 4kb7 and 4r1e (doubtful); 3bmo and 3jqf (excluded, since LP's ligand of record there isn't the DTT).

Issue #8 offers to open Finding 3 as its own issue with residues, atoms and distances. Not filed yet.

## Related

Issue #7 on THGLab/LP-PDBBind, "Curation of SMILES strings?", was opened by Rocco Moretti (Rosetta, Meiler lab) in Aug 2025 and had no replies as of 16 Oct 2026. He flagged 2jg8's SMILES as an aldehyde where a carboxylate is expected, which is the same failure mode at the single-entry level.

## Status

- [x] **GitHub issue** filed: #8, 16 Oct 2026, 22:53 AKDT.
- [ ] **LinkedIn post**: drafted in `linkedin_post.md`. Post it Tuesday morning.
- [ ] Watch #8 for a maintainer reply or a reply from Moretti.
- [ ] If they want it: open Finding 3 as a separate issue with residues, atoms and distances.
- Don't quote the 16.7% unrecorded-covalent headline anywhere new until it's recomputed on class 1+2 entries for v1.8.

## Files

- `lp_pdbbind_label_conflicts.csv`: the 37 conflict rows
- `github_issue.md`: the issue exactly as filed
- `linkedin_post.md`: post draft for Tuesday

Scan method and archive results: Varman, *Silent Mutation Misclassification and Unrecorded Covalent Bonds in the Protein Data Bank*, v1.7 — https://doi.org/10.5281/zenodo.23380304
