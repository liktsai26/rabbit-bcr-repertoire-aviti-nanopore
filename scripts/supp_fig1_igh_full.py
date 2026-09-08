#!/usr/bin/env python3
"""Supplementary Figure 1: full IGH gene usage (companion to Figure 3).
A) IGHV top-20 bar chart (Figure 3A used top-5 only).
B) IGHD top-15 bar chart (Figure 3B used top-5 only).
C) IGHJ, all genes (Figure 3C used top-5 only).
D) IGH V-gene scatter correlation using all genes (Figure 3E used top-5 only)."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import pearsonr

from fig_utils import (get_cache, gene_freq_pct, FIG_DIR,
                       RABBITS, PLATFORMS)

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
os.makedirs(FIG_DIR3, exist_ok=True)

cache = get_cache()
DATA  = cache['data']

plt.rcParams.update({'font.family': 'Arial', 'axes.labelsize': 10,
                     'axes.titlesize': 12, 'xtick.labelsize': 7,
                     'ytick.labelsize': 9})

C_AVITI = '#1F6EA8'
C_NANO  = '#E07B39'
W       = 0.35

def gene_freq(platform, rabbit, chain, col):
    return gene_freq_pct(DATA[(platform, rabbit, chain)], col)

vfreq = {(p,r): gene_freq(p, r, 'IGH', 'v_call') for p in PLATFORMS for r in RABBITS}
dfreq = {(p,r): gene_freq(p, r, 'IGH', 'd_call') for p in PLATFORMS for r in RABBITS}
jfreq = {(p,r): gene_freq(p, r, 'IGH', 'j_call') for p in PLATFORMS for r in RABBITS}

def top_n_by_aviti(fd, n):
    genes = set(fd[('AVITI',1)].index) | set(fd[('AVITI',2)].index)
    avg   = {g: (fd[('AVITI',1)].get(g,0) + fd[('AVITI',2)].get(g,0)) / 2 for g in genes}
    return sorted(avg, key=avg.get, reverse=True)[:n]

def sort_by_aviti(fd):
    genes = set().union(*[set(fd[(p,r)].index) for p in PLATFORMS for r in RABBITS])
    avg   = {g: (fd[('AVITI',1)].get(g,0) + fd[('AVITI',2)].get(g,0)) / 2 for g in genes}
    return sorted(avg, key=avg.get, reverse=True)

def bar2(ax, genes, get_val, legend=False):
    x = np.arange(len(genes))
    for pi, plat in enumerate(PLATFORMS):
        col  = C_AVITI if plat == 'AVITI' else C_NANO
        off  = -W/2 + pi * W
        vals_by_gene = [[get_val(plat, r, g) for r in RABBITS] for g in genes]
        means = [np.mean(v) for v in vals_by_gene]
        stds  = [np.std(v)  for v in vals_by_gene]
        ax.bar(x + off, means, W, color=col, alpha=0.85, label=plat,
               edgecolor='white', linewidth=0.3,
               yerr=stds, capsize=2.5,
               error_kw={'elinewidth': 0.8, 'ecolor': '#333333'})
        for r in RABBITS:
            ax.scatter(x + off, [get_val(plat, r, g) for g in genes],
                      color='black', s=18, zorder=4, edgecolors='white', linewidths=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Frequency (%)', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    if legend:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore')],
                  fontsize=8, frameon=False)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 17))
gs  = GridSpec(3, 1, figure=fig, hspace=0.60, height_ratios=[1.2, 1.0, 1.0])

# ── A: IGHV top 20 ────────────────────────────────────────────────────────────
ax_a    = fig.add_subplot(gs[0])
genes_v = top_n_by_aviti(vfreq, 20)
bar2(ax_a, genes_v, lambda p,r,g: vfreq[(p,r)].get(g,0), legend=True)
ax_a.set_title('IGHV (top 20)', fontsize=12)
fig.text(0.01, gs[0].get_position(fig).y1 + 0.004,
         'A', fontsize=16, fontweight='bold', va='bottom')

# ── B: IGHD top 15 | C: IGHJ all ──────────────────────────────────────────────
gs_bc = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], wspace=0.42)

ax_b    = fig.add_subplot(gs_bc[0])
genes_d = top_n_by_aviti(dfreq, 15)
bar2(ax_b, genes_d, lambda p,r,g: dfreq[(p,r)].get(g,0))
ax_b.set_title('IGHD (top 15)', fontsize=12)
pos_b = gs_bc[0].get_position(fig)
fig.text(pos_b.x0 - 0.01, pos_b.y1 + 0.004, 'B', fontsize=16, fontweight='bold', va='bottom')

ax_c    = fig.add_subplot(gs_bc[1])
genes_j = sort_by_aviti(jfreq)
bar2(ax_c, genes_j, lambda p,r,g: jfreq[(p,r)].get(g,0))
ax_c.set_title('IGHJ (all genes)', fontsize=12)
pos_c = gs_bc[1].get_position(fig)
fig.text(pos_c.x0 - 0.01, pos_c.y1 + 0.004, 'C', fontsize=16, fontweight='bold', va='bottom')

# ── D: Scatter — all V genes pooled across rabbits (centered, square) ────────
gs_d = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2], wspace=0.1)
ax_d  = fig.add_subplot(gs_d[0])
all_g = sorted(set().union(*[set(vfreq[(p,r)].index) for p in PLATFORMS for r in RABBITS]))
xall, yall = [], []
for rabbit in RABBITS:
    xall.extend([vfreq[('AVITI',    rabbit)].get(g, 0) for g in all_g])
    yall.extend([vfreq[('Nanopore', rabbit)].get(g, 0) for g in all_g])
ax_d.scatter(xall, yall, color='#555555', s=22, alpha=0.5, zorder=3)
valid = [(a,b) for a,b in zip(xall,yall) if a>0 or b>0]
if len(valid) >= 3:
    r_val, _ = pearsonr(*zip(*valid))
    ax_d.text(0.05, 0.90, f'r = {r_val:.3f}', transform=ax_d.transAxes,
              fontsize=9, color='#333333')
ax_max = max(xall + yall, default=0)
lim = ax_max * 1.1 or 1
ax_d.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.4)
ax_d.set_xlim(0, lim); ax_d.set_ylim(0, lim)
ax_d.set_xlabel('AVITI frequency (%)', fontsize=10)
ax_d.set_ylabel('Nanopore frequency (%)', fontsize=10)
ax_d.set_title('IGH (all V genes)', fontsize=12)
ax_d.set_aspect('equal', adjustable='box')
ax_d.spines[['top', 'right']].set_visible(False)
pos_d = gs_d[0].get_position(fig)
fig.text(pos_d.x0 - 0.01, pos_d.y1 + 0.004,
         'D', fontsize=16, fontweight='bold', va='bottom')

fig.suptitle(f'Supplementary Figure 1. Full IGH gene repertoire underlying Figure 3A/B/C/E (n={len(RABBITS)})',
             fontsize=13, fontweight='bold', y=1.01)

out = os.path.join(FIG_DIR3, 'SuppFigure1_IGH_FullGeneUsage.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
