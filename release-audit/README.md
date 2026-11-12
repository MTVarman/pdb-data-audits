# Weekly PDB release audit

The PDB releases new entries every Wednesday. This runs four checks against
that week's release and commits the result, dated.

Each week lands in `weeks/<YYYY-MM-DD>_findings.csv` with a short
`_brief.md` beside it. Nothing is edited afterwards. A wrong call stays in the
record with the week it was made.

## The checks

**PARITY** — heme c (`HEC`) is defined by its covalent attachment: two
cysteine thioethers from the C-x-x-C-H motif. The check finds CXXCH in the
polymer sequence first, then looks for the heme it belongs to. Sequence is the
one thing in the file refinement cannot move, so the attachment site is fixed
independently of where the atoms ended up — and where the atoms ended up is
what is being audited.

Per site, one of four states:

| state | meaning |
|---|---|
| `OK` | both thioethers within bonding range, both carry a `covale` record |
| `UNRECORDED` | both within bonding range, at least one has no record |
| `UNBONDED_HIS_ON_FE` | the motif's histidine coordinates the iron; its cysteines are modelled beyond bonding range — the site is half built |
| `UNBONDED_DISPLACED` | motif and heme both present, but neither the histidine nor the cysteines are in contact |
| `NO_HEME` | motif present, no heme c near it |

`UNBONDED_HIS_ON_FE` is the strong one. The heme cannot be "somewhere else"
when the same motif's own histidine is holding its iron.

This is the check worth having, because no distance threshold decides it. The
motif decides it. Either the site has its two thioethers or it does not.

**GEOMETRY** — a contact shorter than the shortest real bond for that pair of
elements, double and triple bonds included (1.20 Å for N–C, 1.15 Å for O–C,
1.55 Å for S–C, 1.90 Å for S–S; the full table is `BOND_FLOOR` in the script).
These are atoms modelled on top of each other, not chemistry. A retinal Schiff
base at 1.3 Å or a glycosidic bond at 1.4 Å is a bond and is not listed.

**CONTACT** — a nucleophilic protein atom within 2.3 Å of a ligand heavy atom
carrying no `covale` or `disulf` record, reported in two bands (`<2.0` and
`2.0–2.3`). The band is stated, not judged, because 2.0 Å is where standard
bond perception stops looking and that is the thing under audit.

Hydrogens are stripped, zero-occupancy atoms skipped, altlocs must match, and
metals are excluded — a metal coordination bond is not a covalent bond and is
not counted as one.

**DEFAULT** — a distance that equals a refinement-software default to the
third decimal. Two metal defaults in cctbx, the library under PHENIX, are
traced to the code and reproduced in
[doi:10.5281/zenodo.24150338](https://doi.org/10.5281/zenodo.24150338). The
N–O value is reproduced in [`default-values/n-o-link/`](default-values/n-o-link/).

| what | rule | the default |
|---|---|---|
| zinc | a zinc ion with fewer than four Cys SG / His N ligands (within 2.8 Å) and at least two Zn–N(His) distances within 0.004 Å of 2.300 | 2.30 Å, sigma 0.03, applied to every such zinc; X-ray at 1.5 Å or better has a median of 2.034 Å |
| N–O pair | a protein N or O within 0.004 Å of 1.297 Å from a ligand O or N, with no bond record | 1.297 Å, sigma 0.01: the generic N–O bond value in cctbx, given to a lysine NZ and an ATP oxygen once they are within 2.20 Å. These atoms do not bond; a salt bridge sits near 2.8 Å |
| O–C pair | a protein O within 0.004 Å of 1.374 Å from a ligand C, with no bond record | 1.374 Å: the generic C–O bond value in the same cctbx table. Value match; the restraint has not been reproduced for this pair |
| chlorophyll-type Mg | the nearest outside nitrogen to each pigment magnesium (within 2.6 Å); the entry is listed when at least 3 and a quarter of them, or at least 2 and half of them, sit at 1.895 ± 0.010 or 2.302 ± 0.006 Å | 1.895 Å in cctbx builds 2020.10–2024.5, 2.300 Å from 2024.6, sigma 0.01; X-ray sits near 2.20 Å |

One row per zinc site, one row per entry and value for magnesium, one row per
N–O or O–C pair. A listed
distance is not by that fact wrong. It is the software's number, not the
map's.

Both filed weeks were re-run on 11 Nov 2026 with the element-pair geometry
rule and the full DEFAULT check; each brief says what changed.

From the 2026-11-11 release on, CONTACT measures to the correct symmetry image of
the ligand atom. Before that a contact made across a crystal symmetry
operator was listed with the distance to the unmoved copy.

## Running it

```
python3 release_audit.py 2026-11-04 [workdir]
```

Downloads that date's release from RCSB, scans it, writes the CSV and the
brief. Re-running is safe; the mmCIF cache is reused.

## What it is not

It flags candidates. It does not adjudicate them. A contact at 2.1 Å with no
record may be a missing bond, a disordered side chain, or a modelling error,
and the three are indistinguishable from coordinates alone.

The parity check is the exception. A CXXCH site either has its two thioethers
or it does not, and that is arithmetic rather than opinion.
