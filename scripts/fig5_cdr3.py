#!/usr/bin/env python3
"""Figure 4 (v3): CDR3 — 4A length / 4B pooled AA heatmap (AVITI vs Nano)"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import LinearSegmentedColormap

from fig_utils import (get_cache, FIG_DIR, CHAINS, RABBITS, PLATFORMS, COLORS,
                       cdr3_lengths as _cdr3_lengths_from_summary, cdr3_cores)

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
os.makedirs(FIG_DIR3, exist_ok=True)

cache = get_cache()
DATA  = cache['data']

plt.rcParams.update({'font.family': 'Arial', 'axes.labelsize': 10,
                     'axes.titlesize': 12, 'xtick.labelsize': 8,
                     'ytick.labelsize': 8})

PLOT_ORDER = [(p,r) for r in RABBITS for p in PLATFORMS]
AAS        = list('ACDEFGHIKLMNPQRSTVWY')
POS_HALF   = 6
POSITIONS  = list(range(-POS_HALF, POS_HALF + 1))
C_AVITI    = '#1F6EA8'
C_NANO     = '#E07B39'

# ── CDR3 length ───────────────────────────────────────────────────────────────
def cdr3_lengths(platform, rabbit, chain):
    return _cdr3_lengths_from_summary(DATA[(platform, rabbit, chain)])

# ── CDR3 AA matrix (pooled across all rabbits for one platform) ───────────────
def cdr3_aa_matrix_pooled(platform, chain='IGH'):
    all_cores = []
    total_n   = 0
    for r in RABBITS:
        cores = cdr3_cores(DATA[(platform, r, chain)])
        total_n += len(cores)
        all_cores.extend(cores)
    if not all_cores:
        return None, 0, {aa: 0 for aa in AAS}
    counts = {p: {aa: 0 for aa in AAS} for p in POSITIONS}
    totpos = {p: 0 for p in POSITIONS}
    aa_total = {aa: 0 for aa in AAS}
    for seq in all_cores:
        L = len(seq)
        for i, aa in enumerate(seq):
            pos = i - (L - 1) // 2
            if pos in counts and aa in counts[pos]:
                counts[pos][aa] += 1
                totpos[pos]     += 1
                aa_total[aa]    += 1
    matrix = np.zeros((len(AAS), len(POSITIONS)))
    for pi, pos in enumerate(POSITIONS):
        if totpos[pos] > 0:
            for ai, aa in enumerate(AAS):
                matrix[ai, pi] = counts[pos][aa] / totpos[pos]
    return matrix, total_n, aa_total

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 18))
gs  = GridSpec(2, 1, figure=fig, hspace=0.55, height_ratios=[1.0, 1.4])

# ── 4A: CDR3 length (unchanged format, 4 lines per chain) ────────────────────
gs_a   = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], wspace=0.35)
LINE_MK = {('AVITI',1): ('-','o'), ('Nanopore',1): ('-','^'),
           ('AVITI',2): ('-','s'), ('Nanopore',2): ('-','D')}

XR_CDR3 = np.arange(3, 25)

for ci, chain in enumerate(CHAINS):
    ax = fig.add_subplot(gs_a[ci])
    for plat in PLATFORMS:
        col  = C_AVITI if plat == 'AVITI' else C_NANO
        curves = []
        for r in RABBITS:
            ls_arr = cdr3_lengths(plat, r, chain)
            if len(ls_arr) == 0:
                continue
            counts = pd.Series(ls_arr).value_counts().reindex(XR_CDR3, fill_value=0)
            curves.append((counts / counts.sum() * 100).values)
        if not curves:
            continue
        curves = np.array(curves)
        mean   = curves.mean(axis=0)
        std    = curves.std(axis=0)
        ax.plot(XR_CDR3, mean, color=col, lw=2.0, label=plat)
        ax.fill_between(XR_CDR3, mean - std, mean + std, color=col, alpha=0.18)
    ax.set_xlabel('CDR3 length (aa)', fontsize=10)
    ax.set_ylabel('Frequency (%)', fontsize=10)
    ax.set_title(chain, fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)
    if ci == 0:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI (mean ± SD)'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore (mean ± SD)')],
                  fontsize=8, frameon=False)

fig.text(0.01, gs[0].get_position(fig).y1 + 0.004,
         'A', fontsize=16, fontweight='bold', va='bottom')

# ── 4B: Pooled heatmap — AVITI | Nanopore ────────────────────────────────────
gs_b = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], wspace=0.35)

CMAP_MONO = LinearSegmentedColormap.from_list('mono_blue', ['white', '#053061'])

mats = [(cdr3_aa_matrix_pooled(p, 'IGH')) for p in PLATFORMS]

# 胺基酸列順序改依整體使用率（兩平台合併後的總出現次數）由高到低排序，
# 兩個 panel 共用同一個順序才能直接上下對照同一列是同一個胺基酸。
combined_totals = {aa: 0 for aa in AAS}
for mat, n, aa_total in mats:
    for aa in AAS:
        combined_totals[aa] += aa_total[aa]
AAS_SORTED = sorted(AAS, key=lambda aa: combined_totals[aa], reverse=True)
aa_reorder = [AAS.index(aa) for aa in AAS_SORTED]

for pi, (plat, (mat, n, aa_total)) in enumerate(zip(PLATFORMS, mats)):
    ax    = fig.add_subplot(gs_b[pi])
    col   = C_AVITI if plat == 'AVITI' else C_NANO
    title = f'{plat} (pooled, n={n:,})'
    if mat is None:
        ax.text(0.5, 0.5, 'No CDR3 data', transform=ax.transAxes,
                ha='center', va='center', fontsize=9, color='grey')
        ax.set_title(title, fontsize=10); continue
    mat_sorted = mat[aa_reorder, :]
    im = ax.imshow(mat_sorted, aspect='auto', cmap=CMAP_MONO,
                   vmin=0, vmax=0.7, interpolation='nearest')
    ax.set_yticks(range(len(AAS_SORTED)))
    ax.set_yticklabels(AAS_SORTED, fontsize=6)
    ax.set_xticks(range(len(POSITIONS)))
    ax.set_xticklabels([str(p) for p in POSITIONS], fontsize=7)
    ax.set_xlabel('CDR3 position (relative to center)', fontsize=8)
    ax.set_ylabel('Amino acid', fontsize=8)
    ax.set_title(title, fontsize=10, color=col, fontweight='bold')
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)

fig.text(0.01, gs[1].get_position(fig).y1 + 0.004,
         'B', fontsize=16, fontweight='bold', va='bottom')

out = os.path.join(FIG_DIR3, 'Figure5_CDR3_v4.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
