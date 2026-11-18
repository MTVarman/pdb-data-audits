#!/usr/bin/env python3
"""Recount the headline numbers of the pdb-data-audits dataset from its own tables.

Repository: github.com/MTVarman/pdb-data-audits  ·  Dataset: doi:10.5281/zenodo.23441847
Usage:  python3 verify_counts.py
Sits at the repository root, next to README.md. Standard library only.
Runs verify_strict.py as well when that file sits alongside.

A table that is absent is reported as a FAIL for its section and the remaining sections still run.

Counted elsewhere: the scan floor (16,827 sites / 2,568 entries) and the chemistry-consistent cut
(1,705 / 477) are rows of Table_S6_main_triage_sites.csv in the preprint supplement
(doi:10.5281/zenodo.23380304, S6_unrecorded_covalent/); the verify_counts.py shipped there recounts
them. The two-rulers scripts need a 506 MB archive index that is available on request.
"""
import csv, json, os, subprocess, sys
sys.stdout.reconfigure(line_buffering=True)

HERE = os.path.dirname(os.path.abspath(__file__))


class Missing(Exception):
    pass


def path(*parts):
    p = os.path.join(HERE, *parts)
    if not os.path.exists(p):
        raise Missing(p)
    return p


def rows(*parts):
    with open(path(*parts), newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh))


results = []
def check(label, got, expected):
    ok = got == expected
    results.append(ok)
    print(f'{"PASS" if ok else "FAIL"}  {label}: got {got}, expected {expected}')


def strict_table():
    return {r['site_id']: r for r in rows('unrecorded-covalent', 'strict-cut-491',
                                           'unrecorded_covalent_tight_cut_491_sites.csv')}


def remediation():
    strict = strict_table()
    ref = rows('unrecorded-covalent', 'remediation-history', 'refiled_check.csv')
    rem = {r['site_id'] for r in ref}
    check('refiled_check rows / distinct removed sites', (len(ref), len(rem)), (141, 135))
    check('removed sites that are in the strict table', sum(s in strict for s in rem), 135)
    check('entries with a removed strict site', len({strict[s]['entry'] for s in rem if s in strict}), 88)
    check('removed sites density-supported',
          sum(strict[s]['density_class'] == 'density_supported' for s in rem if s in strict), 98)
    dele = rows('unrecorded-covalent', 'remediation-history', 'deleted_rows_strict.csv')
    check('same 135 sites from the deleted-row table',
          {r['matched_site_id'] for r in dele if r['matched_site_id'].strip()} == rem, True)


def sars_cov_2():
    sc = rows('unrecorded-covalent', 'sars-cov-2', 'sars_cov2_unrecorded_covalent_appendix.csv')
    check('sites / entries', (len(sc), len({r['pdb_id'] for r in sc})), (73, 46))
    mp = [r for r in sc if r['residue'] == 'CYS 145']
    check('Mpro Cys145 sites density-supported',
          (sum(r['density_class'] == 'density_supported' for r in mp), len(mp)), (66, 66))


def benchmark_overlap():
    pb = rows('unrecorded-covalent', 'benchmark-overlap', 'posebusters_308.csv')
    with open(path('unrecorded-covalent', 'benchmark-overlap', 'posebusters_overlap_summary.json')) as fh:
        js = json.load(fh)
    check('PoseBusters pairs', len(pb), 308)
    check('PoseBusters exact pair hits among the 16,827 scan-floor sites', js['exact_pair_hits'], 0)
    rnp = rows('unrecorded-covalent', 'benchmark-overlap', 'strict_subset_in_runs_n_poses.csv')
    check("Runs N' Poses: strict-subset sites / entries", (len(rnp), len({r['entry'] for r in rnp})), (23, 7))


def lp_pdbbind():
    lc = rows('lp-pdbbind-label-conflicts', 'lp_pdbbind_label_conflicts.csv')
    twins = [r for r in lc if r['tier'] in ('no_substitution', 'few_substitutions')]
    check('label-conflict rows', len(lc), 37)
    check('identical-ligand twins labelled both ways (issue #8): entries', len({r['entry'] for r in twins}), 11)
    check('  of which train / test',
          (sum(r['split'] == 'train' for r in twins), sum(r['split'] == 'test' for r in twins)), (10, 1))
    lx = rows('unrecorded-covalent', 'lp-pdbbind-crosscheck', 'strict_vs_lp_27_entries.csv')
    ok = [r for r in lx if r['gemini_status'] in ('confirmed', 'plausible')]
    check('strict sites labelled covalent=False: sites / entries', (len(lx), len({r['entry'] for r in lx})), (38, 27))
    check('confirmed or plausible covalent: sites / entries', (len(ok), len({r['entry'] for r in ok})), (29, 20))


def reactive_residues():
    rr = rows('unrecorded-covalent', 'reactive-residues', 'removed_811_classified.csv')
    check('located sites (811 removed; 7 atoms not found)', len(rr), 804)
    check('density supported / not / unknown',
          tuple(sum(r['density'] == k for r in rr) for k in ('supported', 'not_supported', 'unknown')), (298, 122, 384))
    nc = rows('unrecorded-covalent', 'reactive-residues', 'non_cys_covalent_targets.csv')
    check('hand-picked density-supported non-Cys sites', len(nc), 23)
    check('KRAS G12D 8JHL among them', sum(r['entry'].upper() == '8JHL' for r in nc), 1)


def release_audit():
    for wk, exp in (('2026-11-04', (116, 45)), ('2026-11-11', (39, 20))):
        f = rows('release-audit', 'weeks', f'{wk}_findings.csv')
        check(f'{wk} findings / entries', (len(f), len({r["entry"] for r in f})), exp)


SECTIONS = [
    ('remediation-history: recorded, then removed by a revision', remediation),
    ('sars-cov-2', sars_cov_2),
    ('benchmark-overlap', benchmark_overlap),
    ('LP-PDBBind', lp_pdbbind),
    ('reactive-residues (archive-removed covalent sites)', reactive_residues),
    ('release-audit weekly findings', release_audit),
]

for title, fn in SECTIONS:
    print(f'-- {title} --')
    try:
        fn()
    except Missing as e:
        results.append(False)
        print(f'FAIL  missing {e}')
    print()

vs = os.path.join(HERE, 'verify_strict.py')
if os.path.exists(vs):
    print('-- verify_strict.py --')
    results.append(subprocess.call([sys.executable, vs]) == 0)
    print()

print(f'{sum(results)}/{len(results)} checks passed (verify_strict.py counted as one)')
sys.exit(0 if all(results) else 1)
