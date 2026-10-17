# Merging the warhead validator into the covalent triage

**Status: not built. This is a decision document — read it before anything runs.**

The validator (`covalent_warhead_validator.py`, v2) is better than the triage chemistry in three specific places and worse in none that matter. Merging it is the right call. But it changes which number comes out of the archive scan, and that number is in a published preprint, so the decision is about what you're claiming, not about code.

---

## 1. The decision in one line

The validator gates on **resolution, occupancy and altloc**. The current scan does not. Applying those gates moves the headline from 2,568 entries toward roughly 289.

## 2. The three numbers

All from `analysis/unrecorded_covalent/` at snapshot 20261013T233952Z.

| | Sites | Entries | What it is |
|---|---:|---:|---|
| **Scan floor** | 16,827 | 2,568 | Organic ligand, protein nucleophile, ≤2.0 Å, no `_struct_conn` row. All classes, all methods, all resolutions. **This is the 16.7%.** |
| **Chemistry-consistent** | 1,705 | 477 | Triage class 2: bond-length-consistent adduct. Excludes known cofactors and artefacts. |
| **Defensible** | 491 | 289 | Class 2, X-ray only, ≤2.5 Å. |

Why the floor is hard to defend as a drug-discovery claim:

- **8,445 of the 16,827 sites are PEB** — phycoerythrobilin, a light-harvesting cofactor on cysteine. Half the headline is one bilin. It's a real covalent bond and nobody's error.
- **13,669 sites are cryo-EM against 3,143 X-ray**, median resolution 3.0 Å. At 3 Å a 2.0 Å contact is a much softer assertion than at 1.8 Å.
- Class 3 (likely artefact) is **9,863 sites in 1,903 entries** — more than half the entries in the floor.

The floor is honest as a scan result. It is not honest as "unrecorded covalent bonds that matter," and that's how a hostile reader will take it.

### The 289 is not a legacy problem — lead with this

Per release year in the strict subset: 2020: 20, 2021: 14, 2022: 24, **2023: 39**, 2024: 20, 2025: 17, 2026: 26. That is **160 of the 289 deposited since 2020**, more than half, with 2023 the peak year.

So these are not old files from before the conventions settled. They are well-resolved crystal structures, X-ray at ≤2.5 Å, arriving at a steady rate right now.

Better headline than any percentage, because it can't be waved off as ancient history and it doesn't invite an argument about denominators:

> 289 well-resolved crystal structures carry a covalent bond in the coordinates with no record of it, and more than half of them arrived since 2020.

## 3. What the validator does better

**Resolution gate.** Default ≤2.5 Å. The scan has no resolution gate at all.

**Occupancy gate.** Rejects zero-occupancy atoms on either side. Not currently checked.

**Altloc gate.** Rejects discordant alternate conformations. Your own sensitivity run dropped **157 pairs in a 1,000-entry sample** for altloc mismatches — so this is a real and frequent failure mode that currently isn't handled.

**Mechanism naming.** Class 2 currently says "C-S" or "B-O". The validator says "acrylamide Michael acceptor on Cys SG" or "tetrahedral boronate adduct." Same evidence, far more legible to a chemist.

**The mmCIF row.** Nothing in the current pipeline emits a proposed `_struct_conn` record. This is the genuinely new capability and the reason to do the merge at all — it turns a finding into a patch someone can apply.

## 4. What the triage does better — keep these

**Class 1 cofactor and PTM recognition.** Eleven curated families with their attachment atoms: bilins on Cys, heme c links, covalent flavin, PLP Schiff base on Lys, biotin and lipoate, retinal, fatty acyl, prenyl, phosphopantetheine, glutathionylation, N-terminal blocking groups. The validator has none of this and would misfile every one of them.

**Class 3A nucleotide rule.** 52 nucleotide-like components against Asp/Glu/Lys/Ser/Thr. GTP alone accounts for 1,618 sites. The validator would pass many of these as geometric candidates.

**Class 3B additive list.** 56 crystallisation components. The validator covers only BME and DTT.

**Class 3C element-pair rule.** Pairs with no covalent window (N-O, O-N, C-C) are the single largest artefact category — over 4,000 sites. The validator reaches the same outcome by a different route but the triage rule is explicit and auditable.

**Element resolution from the chemical-component formula.** The triage resolves two-letter symbols properly. The validator's atom-name heuristic does not — see bugs below.

## 5. The merge

**One chemistry engine.** The triage pipeline keeps ownership. The validator contributes gates and output, not a parallel classifier. Two engines that can disagree is worse than one engine that's wrong.

Order of operations:

1. **New pre-gates, before classification:** occupancy, altloc, resolution. Record each exclusion with its reason and keep the counts — these numbers are themselves a finding.
2. **Classification unchanged:** existing class 1 / 3A / 3B / 2 / 3C rules, in the existing order.
3. **Mechanism labelling, new:** for class 2 sites, attach a mechanism name from the warhead catalogue where the component matches, otherwise from the element pair and residue.
4. **mmCIF row, new:** emit for class 2 sites only. Never for class 1 (already known chemistry), never for class 3, never for geometric-only candidates.

**Resolution must be a parameter, not a constant.** Run the scan at no limit, at 3.0 Å and at 2.5 Å, and report all three. Burying a 2.5 Å default inside the code is how the headline silently changes without anyone deciding.

## 6. Bugs to fix before it runs

**`CA` and `CD` collide with metal symbols.** `STANDARD_PROTEIN_ATOMS` covers eight residues. Any other residue falls through to the heuristic, where `stripped[:2]` matches "CA" against calcium and "CD" against cadmium. So an alpha carbon on ALA or GLY, or a CD on ARG or PRO, is excluded as metal coordination. Worse on ligands, which never hit the lookup at all: a peptidomimetic with an atom named `CA1` reads as calcium and the site disappears. That hits 80I, HUR and 4WI — all three in your strict subset. **This silently deletes real findings and leaves no trace in the output.** Use the component formula for ligand atoms; never apply two-letter metal matching to a polymer atom.

**Warheads need a set of acceptable element pairs, not one.** `EVS` targets `{CYS, SER}` and `{SG, OG}` with a single pair `('S','C')`, so a serine hit gives O–C, fails the check and loses its mechanism. Same for `DON` on GLU and `38N` on a backbone N.

**Geometric candidates emit a record reading "Validated covalent bond."** On a result whose own confidence is `PLAUSIBLE`. Either carry the tier in the text or emit nothing below PROBABLE.

**The warhead descriptions are unverified.** The codes are real — 4WI, 80I, HUR, A1EN0 and BO2 all appear in your own strict-subset ligand counts. But "nirmatrelvir-derived nitrile" and "afatinib-derived acrylamide" are claims about *which drug*, asserted from memory by a language model. Check each against the component dictionary before any of it goes public. A wrong attribution inside a correction you're handing a depositor costs more than the correction is worth.

## 7. Verification

- 6T5V must survive every gate and produce a correct row for Cys12 SG to the inhibitor at 1.816 Å.
- 6OIM and 8UW2 must stay unflagged — their bonds are already recorded.
- The eleven LP-PDBBind conflict pairs should be run through and their outcomes recorded.
- Re-run the 1,000-entry sensitivity sample and confirm recall hasn't dropped from 0.79 for any reason other than the new gates.
- Report how many sites each new gate removed, separately.

## 8. What this means for the preprint

v1.7 leads with 16.7%. After this, v1.8 should give three numbers with their definitions, and lead with the defensible one. The floor stays in the paper as a scan result with its composition stated — half of it is one bilin, four fifths of it is cryo-EM at 3 Å.

That is a smaller headline and a much harder one to attack. It is also the version that survives someone checking.
