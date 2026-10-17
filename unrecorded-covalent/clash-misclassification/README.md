# The archive does not file these as steric clashes. It files them as nothing at all.

**Result: hypothesis falsified, and the null is more interesting than the claim was. Run 17 Oct 2026.**

---

## What was claimed

That wwPDB validation, lacking a bond record, evaluates these contacts against van der Waals radii, reads a 1.8 Å bond as a >1.6 Å overlap, and reports it as a steric clash outlier — so the archive had been detecting every one of these and filing them under the wrong error type.

## What the data says

**0 of 491 sites. 0 of 289 entries.** Not one flagged contact appears as a clash outlier at the flagged atom pair.

It fails at the control. 6T5V's validation report contains four clash records in total — all hydrogens, on Gln22, Ile24, Ile55 and Arg149. None involve Cys12 or the ligand. The Cys12 and ligand rows carry no clash children at all.

Not a threshold effect. Expected overlap at that contact is about 1.68 Å against a 0.4 Å cutoff — four times over. Across all 491 sites the predicted overlap runs 1.55–1.98 Å, median 1.73 Å, and every single one would clear the cutoff. The largest overlap actually reported anywhere in the residue-level matches is 1.29 Å, between hydrogen-bearing atoms.

Loose matches exist and are not the contact: 45 sites have a clash somewhere between the same residue and ligand, but 22 of those involve only hydrogens and the rest are heavy atoms at 2.2–3.0 Å — neighbouring atoms, not the bond.

Cross-tabbed against density: no density-supported site is reported as a clash at the flagged pair. Neither is any other class.

Same behaviour under both report schema versions in play. Zero network retrievals — every report was re-parsed from the cache the density run left behind.

## Why the zero is the finding

**All 491 contacts are in the mmCIF close-contact table.** That is where the scan found them. So the pipeline measured every one of these pairs, wrote each one down as a close contact — and then left all 491 out of the clash list.

A 1.8 Å sulfur-to-carbon contact at four times the clash cutoff does not quietly fail to be flagged. Something in the pipeline removed it.

The live hypothesis, not yet established: the validation software performs its own distance-based bond perception and excludes pairs it perceives as bonded from clash analysis. If that holds, the pipeline is independently concluding these are covalent bonds — concluding it firmly enough to suppress a clash it would otherwise report — while `_struct_conn` stays empty and every downstream consumer reads the ligand as unbonded.

That is a stronger claim than the one that just died. Not "nobody noticed." **The software noticed, acted on it, and wrote nothing down.**

Supporting detail: 6T5V's reported clashscore is 0.72. A 1.68 Å overlap counted as a clash would not leave that number where it is.

## The mechanism, read from the source

Clash detection in the wwPDB report comes from MolProbity's `clashlist`, which drives the `probe` program. Two facts from probe's own source and MolProbity's own script explain the zero.

**1. probe excludes bonded pairs, and their whole 1–4 neighbourhood, from dot generation.** Bond connectivity is built by `markBonds()` into each atom's `bondedto` list. The `Maxbonded` flag defaults to **4** — exclusion bond depth, covering 1–4 interactions. Pairs inside that neighbourhood generate no contact dots. No dots means no overlap, which means no clash record can exist for them.

**2. probe's bond detection is distance-based by default.** It supports both distance-based and table-based bonding. The `UseStdBond` flag, which switches on standard residue bonding tables, defaults to **FALSE** — and MolProbity's `clashlist` has `-stdbonds` present but **commented out**, so it never turns it on. The wwPDB pipeline's clash analysis therefore runs with geometric bond perception.

So: probe measures a 1.816 Å sulfur-to-carbon separation, concludes from the distance alone that it is a bond, marks the pair bonded, and suppresses contact analysis across the pair and everything within four bonds of it. The contact never reaches the clash list. Not misfiled — deliberately and correctly excluded, by a program that decided it was looking at a covalent bond.

**This is the finding.** The validation pipeline independently perceives every one of these 491 bonds, from geometry, with the annotation-driven path switched off. It acts on that perception. And the perception exists only inside probe's memory for the length of the run. Nothing writes it down. `_struct_conn` stays empty, the mmCIF close-contact table keeps listing the pair as a non-bonded contact, and every downstream consumer reads an unbonded ligand sitting at bonding distance.

The recommendation is no longer "add a check." It is **write down the answer the pipeline already computed.**

### Established vs. not

Established: probe's 1–4 exclusion depth and its default-off standard bond tables; MolProbity's clashlist leaving `-stdbonds` inactive; zero clash records at all 491 flagged pairs; all 491 present in the close-contact table.

Not established: the exact distance threshold `markBonds()` applies; whether the wwPDB production pipeline runs the same clashlist flags as the public MolProbity distribution; whether probe reads `_struct_conn` or LINK records at all (if it does not, the geometric perception is wholly independent of the annotation, which is the strongest reading).

Open, and now cheap to answer: read `markBonds()` in `probe.c` and `stdconntable.c` directly for the distance rule, and confirm the production flags against a wwPDB pipeline description.

## Status of the related layers

- Density validation: 233 of 289 entries density-supported (80.6%), 403 of 491 sites. Unaffected by this.
- The scan, triage, classification rules and every published count: unchanged.

## Pre-registration

Before the reports were read, the prediction was that more than 70% of the 289 entries would carry a clash outlier at the flagged pair, and that under 40% would count as falsifying the mechanism. Observed: 0%. The hypothesis was falsified outright; the mechanism section above is what the null led to.
