"""Figures for the note.  Usage: python3 -I figures.py data_dir out_dir
Figure 1 reads archive_residue_pairs.csv.gz (and boundary_exact_distances.csv for pairs on a boundary);
Figure 2 reads esters_4DL1.csv."""
import csv, sys, random
import gzip
def op(path):
    """Open a table that may be gzip-compressed."""
    return gzip.open(path, 'rt', newline='') if path.endswith('.gz') else open(path, newline='')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D, O = sys.argv[1], sys.argv[2]
BLUE, ORANGE = '#2a78d6', '#eb6834'          # recorded, not recorded (validated pair)
INK, INK2, GRID, BAND = '#0b0b0b', '#52514e', '#e3e2de', '#f0efec'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.edgecolor': GRID, 'axes.linewidth': 0.8,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.labelcolor': INK2, 'text.color': INK})
COV = {'C': 0.77, 'N': 0.70, 'O': 0.66, 'S': 1.04, 'P': 1.10}
BX = {}
for b in csv.DictReader(open(f'{D}/boundary_exact_distances.csv', newline='')):
    if b['source'] == 'archive':
        BX[(b['entry'], b['protein_chain'], b['protein_residue'], b['protein_seq'], b['protein_atom'], b['ligand_chain'], b['ligand_comp'], b['ligand_seq'], b['ligand_atom'])] = float(b['dist_exact'])
A = []
for r in csv.DictReader(op(f'{D}/archive_residue_pairs.csv.gz')):
    try: res = float(r['resolution']); dep = int(r['dep_year'])
    except Exception: continue
    if not r['method1'].startswith('X-RAY') or res > 2.5 or dep < 2014: continue
    a, b = sorted([r['prot_el'].capitalize(), r['lig_el'].capitalize()])
    d = BX.get((r['entry'], r['prot_chain'], r['prot_res'], r['prot_seq'], r['prot_atom'], r['lig_chain'], r['lig_comp'], r['lig_seq'], r['lig_atom']), float(r['dist']))
    A.append((r['prot_res'], f'{a}-{b}', d, r['recorded'] == '1' or r['disulf'] == '1'))
PANELS = [('Cys S to ligand C', ('CYS',), 'C-S', 2.00), ('Asp / Glu O to ligand C', ('ASP', 'GLU'), 'C-O', 1.55),
          ('Ser / Thr / Tyr O to ligand C', ('SER', 'THR', 'TYR'), 'C-O', 1.55), ('Lys N to ligand C', ('LYS',), 'C-N', 1.60),
          ('His N to ligand C', ('HIS',), 'C-N', 1.60), ('Cys S to ligand S', ('CYS',), 'S-S', 2.20)]
rng = random.Random(7)
fig, axes = plt.subplots(2, 3, figsize=(7.1, 3.9), dpi=220)
for ax, (title, residues, pair, upper) in zip(axes.flat, PANELS):
    a, b = pair.split('-'); pl = COV[a] + COV[b] + 0.2
    x0, x1 = upper - 0.13, upper + 0.23
    ax.axvspan(upper, pl, color=BAND, zorder=0, lw=0)
    pts = [(d, rec) for res, p, d, rec in A if res in residues and p == pair and x0 <= d <= x1]
    for want, col, y in ((True, BLUE, 1.0), (False, ORANGE, 0.0)):
        xs = [d for d, rec in pts if rec == want]
        ys = [y + rng.uniform(-0.30, 0.30) for _ in xs]
        ax.scatter(xs, ys, s=5.5 if len(pts) > 800 else 10, color=col, alpha=0.55 if len(pts) > 800 else 0.7, linewidths=0, zorder=2)
    ax.axvline(upper, color=INK, lw=1.1, zorder=3)
    ax.axvline(pl, color=INK2, lw=0.9, zorder=3)
    ax.text(upper - 0.004, 1.62, f'{upper:.2f}', ha='right', va='bottom', fontsize=7.5, color=INK, fontweight='bold')
    ax.text(pl + 0.004, 1.62, f'{pl:.2f}', ha='left', va='bottom', fontsize=7.5, color=INK2)
    ax.set_xlim(x0, x1); ax.set_ylim(-0.55, 1.55)
    ax.set_yticks([1.0, 0.0]); ax.set_yticklabels(['recorded', 'not\nrecorded'], fontsize=7.5)
    ax.tick_params(axis='y', length=0); ax.tick_params(axis='x', length=2.5, width=0.6, labelsize=7.5)
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    ax.set_title(title, fontsize=8.5, fontweight='bold', loc='left', pad=14, color=INK)
for ax in axes[1]: ax.set_xlabel('modelled distance (Å)', fontsize=8)
for ax in axes[:, 1:].flat: ax.set_yticklabels([])
fig.subplots_adjust(left=0.085, right=0.99, top=0.88, bottom=0.11, hspace=0.72, wspace=0.07)
fig.savefig(f'{O}/figure1.png', facecolor='white')
plt.close(fig)

# Figure 2: 4DL1, sixteen heme esters, one row each (table written by esters_4DL1.py)
pts = [(float(r['distance_A']), r['recorded'] == 'True', r['residue_chain'], r['close_contact_row']) for r in csv.DictReader(open(f'{D}/esters_4DL1.csv'))]
pts.sort()
n = len(pts)
fig, ax = plt.subplots(figsize=(7.1, 3.05), dpi=220)
ax.axvline(1.55, color=INK, lw=1.1, zorder=1)
for i, (d, rec, lab, row) in enumerate(pts):
    y = n - 1 - i
    ax.axhline(y, color=GRID, lw=0.5, zorder=0)
    ax.scatter([d], [y], s=30, color=BLUE if rec else ORANGE, edgecolors='white', linewidths=0.9, zorder=3)
ax.text(1.5503, n - 0.15, 'link-table limit 1.55 Å', ha='left', va='bottom', fontsize=7.5, fontweight='bold', color=INK)
lastrec = max(p for p in pts if p[1]); firstun = min(p for p in pts if not p[1])
ax.text(lastrec[0] - 0.0007, n - 1 - pts.index(lastrec), f'{lastrec[0]:.5f}', ha='right', va='center', fontsize=7.2, color=INK)
ax.text(firstun[0] + 0.0007, n - 1 - pts.index(firstun), f'{firstun[0]:.5f}', ha='left', va='center', fontsize=7.2, color=INK)
un = [n - 1 - i for i, p in enumerate(pts) if not p[1]]
xb = 1.5634
ax.plot([xb, xb + 0.0004, xb + 0.0004, xb], [max(un) + 0.3, max(un) + 0.3, min(un) - 0.3, min(un) - 0.3], color=INK2, lw=0.7)
ax.text(xb + 0.0009, (max(un) + min(un)) / 2, "rows 1 to 8\nof the entry's\nclose-contact\ntable", ha='left', va='center', fontsize=7.2, color=INK2)
ax.set_xlim(1.528, 1.5718); ax.set_ylim(-0.8, n + 0.6)
ax.set_xticks([1.53, 1.535, 1.54, 1.545, 1.55, 1.555, 1.56])
ax.set_yticks(range(n)); ax.set_yticklabels([pts[n - 1 - y][2] for y in range(n)], fontsize=6.8)
ax.tick_params(axis='y', length=0); ax.tick_params(axis='x', length=2.5, width=0.6, labelsize=7.5)
for sp in ('top', 'right', 'left'): ax.spines[sp].set_visible(False)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, markersize=5, label='covale record in _struct_conn'),
                   Line2D([], [], marker='o', ls='', color=ORANGE, markersize=5, label='no record')],
          loc='lower left', fontsize=7.2, frameon=False, handletextpad=0.3, borderaxespad=0.2)
ax.set_xlabel('modelled carboxylate O to heme methyl C distance (Å)', fontsize=8)
fig.subplots_adjust(left=0.105, right=0.985, top=0.95, bottom=0.15)
fig.savefig(f'{O}/figure2.png', facecolor='white')
print('figure 2 points', len(pts), 'recorded', sum(p[1] for p in pts), 'last recorded', lastrec[:3], 'first unrecorded', firstun[:3])
