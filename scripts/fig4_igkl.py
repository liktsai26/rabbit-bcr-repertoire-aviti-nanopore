#!/usr/bin/env python3
"""Figure 4 (v3): IgK only — top5 bar (A/B), top5 scatter (C), violin (D).
IGL moved entirely to Supplementary Figure 2 (see supp_fig2_igkl_full.py)."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import pearsonr

from fig_utils import (get_cache, gene_freq_pct, FIG_DIR,
                       RABBITS, PLATFORMS, COLORS, is_available)

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

vfreq = {(p,r): gene_freq(p, r, 'IGK', 'v_call') for p in PLATFORMS for r in RABBITS}
jfreq = {(p,r): gene_freq(p, r, 'IGK', 'j_call') for p in PLATFORMS for r in RABBITS}

def top_v_by_aviti(n=5):
    genes = set(vfreq[('AVITI',1)].index) | set(vfreq[('AVITI',2)].index)
    avg   = {g: (vfreq[('AVITI',1)].get(g,0)+vfreq[('AVITI',2)].get(g,0))/2
             for g in genes}
    return sorted(avg, key=avg.get, reverse=True)[:n]

def top_j_by_aviti(n=5):
    genes = set().union(*[set(jfreq[(p,r)].index) for p in PLATFORMS for r in RABBITS])
    avg   = {g: (jfreq[('AVITI',1)].get(g,0)+jfreq[('AVITI',2)].get(g,0))/2
             for g in genes}
    return sorted(avg, key=avg.get, reverse=True)[:n]

# V identity for panel D
vident = {}
for p in PLATFORMS:
    for r in RABBITS:
        vi = pd.Series(DATA[(p, r, 'IGK')]['v_identity'], dtype=float)
        if len(vi) and vi.max() <= 1.0:
            vi = vi * 100
        vident[(p, r)] = vi[(vi >= 80) & (vi <= 100)].values

# ── Plan B bar helper ─────────────────────────────────────────────────────────
def bar2(ax, genes, get_val, ylabel=True, legend=False):
    x = np.arange(len(genes))
    for pi, plat in enumerate(PLATFORMS):
        col  = C_AVITI if plat == 'AVITI' else C_NANO
        off  = -W/2 + pi * W
        rb   = [r for r in RABBITS if is_available(plat, r, 'IGK')]
        vals_by_gene = [[get_val(plat, r, g) for r in rb] for g in genes]
        means = [np.mean(v) for v in vals_by_gene]
        stds  = [np.std(v)  for v in vals_by_gene]
        ax.bar(x + off, means, W, color=col, alpha=0.85, label=plat,
               edgecolor='white', linewidth=0.3,
               yerr=stds, capsize=2.5,
               error_kw={'elinewidth': 0.8, 'ecolor': '#333333'})
        for r in rb:
            ax.scatter(x + off, [get_val(plat, r, g) for g in genes],
                      color='black', s=18, zorder=4, edgecolors='white', linewidths=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=7)
    if ylabel:
        ax.set_ylabel('Frequency (%)', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    if legend:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore')],
                  fontsize=8, frameon=False)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 13))
gs  = GridSpec(2, 1, figure=fig, hspace=0.35, height_ratios=[1.0, 1.0])

# ── A: IGKV top5 | B: IGKJ (all) ──────────────────────────────────────────────
gs_ab = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], wspace=0.38)

genes_v = top_v_by_aviti(5)
ax_a = fig.add_subplot(gs_ab[0])
bar2(ax_a, genes_v, lambda p,r,g: vfreq[(p,r)].get(g,0), legend=True)
ax_a.set_title('IGKV (top 5)', fontsize=12)
pos_a = gs_ab[0].get_position(fig)
fig.text(pos_a.x0 - 0.01, pos_a.y1 + 0.004, 'A', fontsize=16, fontweight='bold', va='bottom')

ax_b = fig.add_subplot(gs_ab[1])
bar2(ax_b, top_j_by_aviti(5), lambda p,r,g: jfreq[(p,r)].get(g,0))
ax_b.set_title('IGKJ (top 5)', fontsize=12)
pos_b = gs_ab[1].get_position(fig)
fig.text(pos_b.x0 - 0.01, pos_b.y1 + 0.004, 'B', fontsize=16, fontweight='bold', va='bottom')

# ── C: IGK scatter (top5 genes only) | D: IGK V identity violin ──────────────
gs_cd = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], wspace=0.38)

ax_c = fig.add_subplot(gs_cd[0])
rb   = [r for r in RABBITS if is_available('AVITI', r, 'IGK') and is_available('Nanopore', r, 'IGK')]
xall, yall = [], []
for rabbit in rb:
    xall.extend([vfreq[('AVITI',    rabbit)].get(g, 0) for g in genes_v])
    yall.extend([vfreq[('Nanopore', rabbit)].get(g, 0) for g in genes_v])
ax_c.scatter(xall, yall, color='#555555', s=22, alpha=0.5, zorder=3)
valid = [(a,b) for a,b in zip(xall,yall) if a>0 or b>0]
if len(valid) >= 3:
    r_val, _ = pearsonr(*zip(*valid))
    ax_c.text(0.05, 0.90, f'r = {r_val:.3f}', transform=ax_c.transAxes,
              fontsize=9, color='#333333')
ax_max = max(xall + yall, default=0)
lim = ax_max * 1.1 or 1
ax_c.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.4)
ax_c.set_xlim(0, lim); ax_c.set_ylim(0, lim)
ax_c.set_xlabel('AVITI frequency (%)', fontsize=10)
ax_c.set_ylabel('Nanopore frequency (%)', fontsize=10)
ax_c.set_title('IGK (top 5 V genes)', fontsize=12)
ax_c.set_aspect('equal', adjustable='box')
ax_c.spines[['top', 'right']].set_visible(False)
pos_c = gs_cd[0].get_position(fig)
fig.text(pos_c.x0 - 0.01, pos_c.y1 + 0.004, 'C', fontsize=16, fontweight='bold', va='bottom')

# ── D: IGK V identity violin (no per-rabbit dots) ────────────────────────────
ax_d = fig.add_subplot(gs_cd[1])
for pi, plat in enumerate(PLATFORMS):
    col    = C_AVITI if plat == 'AVITI' else C_NANO
    pooled = np.concatenate([vident[(plat, r)] for r in RABBITS])
    if len(pooled) < 5:
        continue
    parts = ax_d.violinplot([pooled], positions=[pi], showmedians=True,
                            showextrema=False, widths=0.7)
    for pc in parts['bodies']:
        pc.set_facecolor(col); pc.set_alpha(0.75); pc.set_edgecolor('none')
    parts['cmedians'].set_color('#222222'); parts['cmedians'].set_linewidth(2)

ax_d.set_xticks([0, 1])
ax_d.set_xticklabels(['AVITI', 'Nanopore'], fontsize=9)
ax_d.set_ylim(80, 100)
ax_d.set_title('IGK', fontsize=12)
ax_d.set_ylabel('V identity (%)', fontsize=10)
ax_d.spines[['top', 'right']].set_visible(False)
ax_d.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                     mpatches.Patch(color=C_NANO,  label='Nanopore')],
            fontsize=8, frameon=False)
pos_d = gs_cd[1].get_position(fig)
fig.text(pos_d.x0 - 0.01, pos_d.y1 + 0.004, 'D', fontsize=16, fontweight='bold', va='bottom')

out = os.path.join(FIG_DIR3, 'Figure4_IgKL_v4.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
