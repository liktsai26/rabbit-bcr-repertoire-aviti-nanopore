#!/usr/bin/env python3
"""Supplementary Figure 2: full IGK gene usage + all IGL results (companion to Figure 4).
A) IGKV top-20 bar chart (Figure 4A used top-5 only).
B) IGKJ, all genes (Figure 4B used top-5 only).
C) IGK V-gene scatter, all genes (Figure 4C used top-5 only).
D-G) IGL results — moved here entirely; Figure 4 (main) shows IGK only."""

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
                       RABBITS, PLATFORMS, is_available)

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
KL_CHAINS = ['IGK', 'IGL']

def gene_freq(platform, rabbit, chain, col):
    return gene_freq_pct(DATA[(platform, rabbit, chain)], col)

vfreq = {(p,r,ch): gene_freq(p, r, ch, 'v_call')
         for p in PLATFORMS for r in RABBITS for ch in KL_CHAINS}
jfreq = {(p,r,ch): gene_freq(p, r, ch, 'j_call')
         for p in PLATFORMS for r in RABBITS for ch in KL_CHAINS}

def top_v_by_aviti(chain, n=20):
    genes = set(vfreq[('AVITI',1,chain)].index) | set(vfreq[('AVITI',2,chain)].index)
    avg   = {g: (vfreq[('AVITI',1,chain)].get(g,0)+vfreq[('AVITI',2,chain)].get(g,0))/2
             for g in genes}
    return sorted(avg, key=avg.get, reverse=True)[:n]

def sort_j_by_aviti(chain):
    genes = set().union(*[set(jfreq[(p,r,chain)].index) for p in PLATFORMS for r in RABBITS])
    avg   = {g: (jfreq[('AVITI',1,chain)].get(g,0)+jfreq[('AVITI',2,chain)].get(g,0))/2
             for g in genes}
    return sorted(avg, key=avg.get, reverse=True)

# V identity for panel G (IGL only — IGK V identity is shown in Figure 4D)
vident = {}
for p in PLATFORMS:
    for r in RABBITS:
        vi = pd.Series(DATA[(p, r, 'IGL')]['v_identity'], dtype=float)
        if len(vi) and vi.max() <= 1.0:
            vi = vi * 100
        vident[(p, r)] = vi[(vi >= 80) & (vi <= 100)].values

def bar2(ax, genes, get_val, chain, ylabel=True, legend=False):
    x = np.arange(len(genes))
    for pi, plat in enumerate(PLATFORMS):
        col  = C_AVITI if plat == 'AVITI' else C_NANO
        off  = -W/2 + pi * W
        rb   = [r for r in RABBITS if is_available(plat, r, chain)]
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

def scatter_panel(ax, chain, genes=None, show_ylabel=True, title_suffix=''):
    rb = [r for r in RABBITS if is_available('AVITI', r, chain) and is_available('Nanopore', r, chain)]
    gl = genes if genes is not None else sorted(set().union(
        *[set(vfreq[(p,r,chain)].index) for p in PLATFORMS for r in rb]))
    xall, yall = [], []
    for rabbit in rb:
        xall.extend([vfreq[('AVITI',    rabbit, chain)].get(g, 0) for g in gl])
        yall.extend([vfreq[('Nanopore', rabbit, chain)].get(g, 0) for g in gl])
    ax.scatter(xall, yall, color='#555555', s=22, alpha=0.5, zorder=3)
    valid = [(a,b) for a,b in zip(xall,yall) if a>0 or b>0]
    ax_max = max(xall + yall, default=0)
    if len(valid) >= 3:
        r_val, _ = pearsonr(*zip(*valid))
        ax.text(0.05, 0.90, f'r = {r_val:.3f}', transform=ax.transAxes,
                fontsize=9, color='#333333')
    lim = ax_max * 1.1 or 1
    ax.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.4)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel('AVITI frequency (%)', fontsize=10)
    if show_ylabel:
        ax.set_ylabel('Nanopore frequency (%)', fontsize=10)
    ax.set_title(f'{chain}{title_suffix}', fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 24))
gs  = GridSpec(4, 1, figure=fig, hspace=0.45,
               height_ratios=[1.2, 1.0, 1.2, 1.0])

# ── Row 1: A(IGKV top20) | B(IGKJ all) ────────────────────────────────────────
gs_ab = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], wspace=0.38)

genes_kv = top_v_by_aviti('IGK', 20)
ax_a = fig.add_subplot(gs_ab[0])
bar2(ax_a, genes_kv, lambda p,r,g: vfreq[(p,r,'IGK')].get(g,0), 'IGK', legend=True)
ax_a.set_title('IGKV (top 20)', fontsize=12)
pos_a = gs_ab[0].get_position(fig)
fig.text(pos_a.x0 - 0.01, pos_a.y1 + 0.004, 'A', fontsize=16, fontweight='bold', va='bottom')

ax_b = fig.add_subplot(gs_ab[1])
bar2(ax_b, sort_j_by_aviti('IGK'), lambda p,r,g: jfreq[(p,r,'IGK')].get(g,0), 'IGK')
ax_b.set_title('IGKJ (all genes)', fontsize=12)
pos_b = gs_ab[1].get_position(fig)
fig.text(pos_b.x0 - 0.01, pos_b.y1 + 0.004, 'B', fontsize=16, fontweight='bold', va='bottom')

# ── Row 2: C(IGK scatter, all genes, centered/square) ────────────────────────
gs_c = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.1)
ax_c = fig.add_subplot(gs_c[0])
scatter_panel(ax_c, 'IGK', genes=None, title_suffix=' (all V genes)')
ax_c.set_aspect('equal', adjustable='box')
pos_c = gs_c[0].get_position(fig)
fig.text(pos_c.x0 - 0.01, pos_c.y1 + 0.004, 'C', fontsize=16, fontweight='bold', va='bottom')

# ── Row 3: D(IGLV top20) | E(IGLJ all) ────────────────────────────────────────
gs_de = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], wspace=0.38)

ax_d = fig.add_subplot(gs_de[0])
bar2(ax_d, top_v_by_aviti('IGL', 20), lambda p,r,g: vfreq[(p,r,'IGL')].get(g,0), 'IGL')
ax_d.set_title('IGLV (top 20)', fontsize=12)
pos_d = gs_de[0].get_position(fig)
fig.text(pos_d.x0 - 0.01, pos_d.y1 + 0.004, 'D', fontsize=16, fontweight='bold', va='bottom')

ax_e = fig.add_subplot(gs_de[1])
bar2(ax_e, sort_j_by_aviti('IGL'), lambda p,r,g: jfreq[(p,r,'IGL')].get(g,0), 'IGL')
ax_e.set_title('IGLJ (all genes)', fontsize=12)
pos_e = gs_de[1].get_position(fig)
fig.text(pos_e.x0 - 0.01, pos_e.y1 + 0.004, 'E', fontsize=16, fontweight='bold', va='bottom')

# ── Row 4: F(IGL scatter, all genes) | G(IGL V identity violin) ─────────────
gs_fg = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[3], wspace=0.38)

ax_f = fig.add_subplot(gs_fg[0])
scatter_panel(ax_f, 'IGL', genes=None, title_suffix=' (all V genes)')
ax_f.set_aspect('equal', adjustable='box')
pos_f = gs_fg[0].get_position(fig)
fig.text(pos_f.x0 - 0.01, pos_f.y1 + 0.004, 'F', fontsize=16, fontweight='bold', va='bottom')

ax_g = fig.add_subplot(gs_fg[1])
for pi, plat in enumerate(PLATFORMS):
    col    = C_AVITI if plat == 'AVITI' else C_NANO
    pooled = np.concatenate([vident[(plat, r)] for r in RABBITS])
    if len(pooled) < 5:
        continue
    parts = ax_g.violinplot([pooled], positions=[pi], showmedians=True,
                            showextrema=False, widths=0.7)
    for pc in parts['bodies']:
        pc.set_facecolor(col); pc.set_alpha(0.75); pc.set_edgecolor('none')
    parts['cmedians'].set_color('#222222'); parts['cmedians'].set_linewidth(2)

ax_g.set_xticks([0, 1])
ax_g.set_xticklabels(['AVITI', 'Nanopore'], fontsize=9)
ax_g.set_ylim(80, 100)
ax_g.set_title('IGL', fontsize=12)
ax_g.set_ylabel('V identity (%)', fontsize=10)
ax_g.spines[['top', 'right']].set_visible(False)
ax_g.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                     mpatches.Patch(color=C_NANO,  label='Nanopore')],
            fontsize=8, frameon=False)
pos_g = gs_fg[1].get_position(fig)
fig.text(pos_g.x0 - 0.01, pos_g.y1 + 0.004, 'G', fontsize=16, fontweight='bold', va='bottom')

fig.suptitle(f'Supplementary Figure 2. Full IGK gene repertoire underlying Figure 4A/B/C, '
             f'and complete IGL results (n={len(RABBITS)})',
             fontsize=13, fontweight='bold', y=1.005)

out = os.path.join(FIG_DIR3, 'SuppFigure2_IGKL_FullGeneUsage.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
