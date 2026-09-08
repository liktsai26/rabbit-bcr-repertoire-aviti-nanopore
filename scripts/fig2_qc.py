#!/usr/bin/env python3
"""Figure 2 (v3): QC — 2A sequencing depth table, 2B productive% bar+dots"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from fig_utils import (get_cache, productive_pct, FIG_DIR,
                       CHAINS, RABBITS, PLATFORMS, is_available)

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
os.makedirs(FIG_DIR3, exist_ok=True)

cache  = get_cache()
DATA   = cache['data']

plt.rcParams.update({'font.family': 'Arial', 'axes.labelsize': 10,
                     'axes.titlesize': 12, 'xtick.labelsize': 9,
                     'ytick.labelsize': 9})

C_AVITI = '#1F6EA8'
C_NANO  = '#E07B39'

def prod_pct(platform, rabbit, chain):
    if not is_available(platform, rabbit, chain):
        return np.nan
    return productive_pct(DATA[(platform, rabbit, chain)])

def reads_for(chain, platform):
    reads = []
    for r in RABBITS:
        if not is_available(platform, r, chain):
            continue
        reads.append(DATA[(platform, r, chain)]['total_count'])
    return reads

def fmt_reads_range(vals):
    # 大數字(AVITI, 十萬~百萬級)用 M 縮寫避免欄位擠不下；小數字(Nanopore, 千級)維持原樣
    if max(vals) >= 100_000:
        return f"{np.mean(vals)/1e6:.2f}M\n({min(vals)/1e6:.2f}–{max(vals)/1e6:.2f}M)"
    return f"{np.mean(vals):,.0f}\n({min(vals):,}–{max(vals):,})"

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 8.6))
gs  = GridSpec(2, 1, figure=fig, hspace=0.28, height_ratios=[1.15, 2.0])

# ── 2A: Sequencing depth table ────────────────────────────────────────────────
ax_tbl = fig.add_subplot(gs[0])
ax_tbl.axis('off')

rows = []
for ch in CHAINS:
    av_reads = reads_for(ch, 'AVITI')
    na_reads = reads_for(ch, 'Nanopore')
    depth_ratio = np.mean(av_reads) / np.mean(na_reads)
    rows.append([ch, fmt_reads_range(av_reads), fmt_reads_range(na_reads), f"{depth_ratio:,.0f}×"])

col_labels = ['Locus', 'AVITI reads\n(analyzed)', 'Nanopore reads\n(analyzed)', 'Depth\nratio']
tbl = ax_tbl.table(cellText=rows, colLabels=col_labels, loc='center', cellLoc='center',
                    bbox=[0.12, 0.10, 0.76, 0.90])
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.5)

for ci in range(len(col_labels)):
    cell = tbl[0, ci]
    cell.set_facecolor('#333333')
    cell.set_text_props(color='white', fontweight='bold')

col_bg = {1: '#DCE6F1', 2: '#FCE4D6'}
for ri in range(len(rows)):
    for ci in range(len(col_labels)):
        cell = tbl[ri + 1, ci]
        cell.set_facecolor(col_bg.get(ci, 'white'))
        if ci == 0:
            cell.set_text_props(fontweight='bold')

ax_tbl.text(0.5, 0.02, 'Values are mean (range) across rabbits; reads = post-QC/trim reads analyzed by IgBLAST '
                        '(raw pre-demux reads not shown).',
            transform=ax_tbl.transAxes, ha='center', fontsize=7.5, style='italic', color='#444444')
fig.text(0.01, gs[0].get_position(fig).y1 - 0.01,
         'A', fontsize=16, fontweight='bold', va='top')

# ── 2B: Productive % ──────────────────────────────────────────────────────────
gs_b = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.38)

for ci, chain in enumerate(CHAINS):
    ax = fig.add_subplot(gs_b[ci])
    for pi, plat in enumerate(PLATFORMS):
        col  = C_AVITI if plat == 'AVITI' else C_NANO
        vals = [prod_pct(plat, r, chain) for r in RABBITS]
        ax.bar(pi, np.nanmean(vals), 0.55, color=col, alpha=0.85, label=plat,
               edgecolor='white', linewidth=0.3,
               yerr=np.nanstd(vals), capsize=4,
               error_kw={'elinewidth': 1.1, 'ecolor': '#333333'})
        for v in vals:
            ax.scatter(pi, v, color='black', s=45, zorder=4,
                      edgecolors='white', linewidths=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['AVITI', 'Nanopore'], fontsize=9)
    ax.set_ylim(0, 115)
    ax.set_title(chain, fontsize=12)
    ax.set_ylabel('Productive reads (%)', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    if ci == 0:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore')],
                  fontsize=8, frameon=False)

fig.text(0.01, gs[1].get_position(fig).y1 + 0.004,
         'B', fontsize=16, fontweight='bold', va='bottom')

out = os.path.join(FIG_DIR3, 'Figure2_QC_v4.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
