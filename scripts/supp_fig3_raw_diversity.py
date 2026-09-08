#!/usr/bin/env python3
"""Supplementary Figure 3 (renumbered 2026-08-03, was Supplementary Figure 1):
Shannon H / Clonality (A) and clone-size distribution (B)
computed on RAW reads, with NO depth-matched subsampling — companion to Figure 6A/B,
which applies the same depth-matched subsampling used there. Shows what the naive
(uncorrected) comparison looks like, to make transparent how much of the apparent
AVITI/Nanopore difference in Figure 6A/B (before correction) was driven purely by
AVITI's much greater raw sequencing depth rather than true biological difference."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from fig_utils import (get_cache, CHAINS, RABBITS, PLATFORMS, clonotype_ids,
                       shannon_clonality, is_available)

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
os.makedirs(FIG_DIR3, exist_ok=True)

cache = get_cache()
DATA  = cache['data']

plt.rcParams.update({'font.family': 'Arial', 'axes.labelsize': 10,
                     'axes.titlesize': 12, 'xtick.labelsize': 8,
                     'ytick.labelsize': 9})

C_AVITI = '#1F6EA8'
C_NANO  = '#E07B39'
W       = 0.35

clono_raw = {}
for p in PLATFORMS:
    for r in RABBITS:
        for ch in CHAINS:
            clono_raw[(p, r, ch)] = clonotype_ids(DATA[(p, r, ch)])

# ── Panel A: raw Shannon H / Clonality (no subsampling) ────────────────────────
diversity_raw = {}
for p in PLATFORMS:
    for r in RABBITS:
        for ch in CHAINS:
            if not is_available(p, r, ch):
                diversity_raw[(p, r, ch)] = (np.nan, np.nan)
                continue
            diversity_raw[(p, r, ch)] = shannon_clonality(clono_raw[(p, r, ch)])

# ── Panel B: raw clone-size distribution (no subsampling) ──────────────────────
BINS_B1   = [1, 2, 11, 10**9]
LABELS_B1 = ['Singleton\n(=1)', 'Small\n(2–10)', 'Large\n(>10)']

def clone_3bin_raw(platform, rabbit, chain):
    if not is_available(platform, rabbit, chain):
        return [np.nan] * (len(BINS_B1) - 1)
    counts = pd.Series(clono_raw[(platform, rabbit, chain)]).value_counts()
    total  = len(counts)
    return [((counts >= lo) & (counts < hi)).sum() / total * 100 if total else 0
            for lo, hi in zip(BINS_B1[:-1], BINS_B1[1:])]

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14.5))
gs  = GridSpec(2, 1, figure=fig, hspace=0.55, height_ratios=[0.75, 0.8])

# Panel A
gs_a    = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], wspace=0.35)
METRICS = [('Shannon H', 0), ('Clonality', 1)]
x_pos   = np.arange(len(CHAINS))

for mi, (metric_name, metric_idx) in enumerate(METRICS):
    ax = fig.add_subplot(gs_a[mi])
    for pi, plat in enumerate(PLATFORMS):
        col  = C_AVITI if plat == 'AVITI' else C_NANO
        off  = -W/2 + pi * W
        vals_by_ch = [[diversity_raw[(plat, r, ch)][metric_idx] for r in RABBITS] for ch in CHAINS]
        means = [np.nanmean(v) for v in vals_by_ch]
        stds  = [np.nanstd(v)  for v in vals_by_ch]
        ax.bar(x_pos + off, means, W, color=col, alpha=0.85, label=plat,
               edgecolor='white', linewidth=0.3,
               yerr=stds, capsize=3,
               error_kw={'elinewidth': 0.9, 'ecolor': '#333333'})
        for r in RABBITS:
            vals = [diversity_raw[(plat, r, ch)][metric_idx] for ch in CHAINS]
            ax.scatter(x_pos + off, vals, color='black', s=22, zorder=4,
                      edgecolors='white', linewidths=0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(CHAINS, fontsize=10)
    ax.set_ylabel(metric_name, fontsize=10)
    ax.set_title(metric_name + '\n(raw reads, no depth correction)', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    if mi == 0:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore')],
                  fontsize=8, frameon=False)

fig.text(0.01, gs[0].get_position(fig).y1 + 0.006,
         'A', fontsize=16, fontweight='bold', va='bottom')

# Panel B
gs_b = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.35)

for ci, chain in enumerate(CHAINS):
    ax = fig.add_subplot(gs_b[ci])
    xb = np.arange(len(LABELS_B1))
    for pi, plat in enumerate(PLATFORMS):
        col        = C_AVITI if plat == 'AVITI' else C_NANO
        off        = -W/2 + pi * W
        all_r_vals = np.array([clone_3bin_raw(plat, r, chain) for r in RABBITS])
        means      = np.nanmean(all_r_vals, axis=0)
        stds       = np.nanstd(all_r_vals, axis=0)
        ax.bar(xb + off, means, W, color=col, alpha=0.85, label=plat,
               edgecolor='white', linewidth=0.3,
               yerr=stds, capsize=3,
               error_kw={'elinewidth': 0.9, 'ecolor': '#333333'})
        for r_vals in all_r_vals:
            ax.scatter(xb + off, r_vals, color='black', s=22, zorder=4,
                      edgecolors='white', linewidths=0.5)
    ax.set_xticks(xb)
    ax.set_xticklabels(LABELS_B1, fontsize=8)
    ax.set_ylabel('Clonotypes (%)', fontsize=10)
    ax.set_title(chain + '\n(raw reads, no depth correction)', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    if ci == 0:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore')],
                  fontsize=8, frameon=False)

fig.text(0.01, gs[1].get_position(fig).y1 + 0.006,
         'B', fontsize=16, fontweight='bold', va='bottom')

fig.suptitle('Supplementary Figure 3. Clonal diversity and clone-size distribution without\n'
             'depth correction — compare to Figure 6A/B (depth-matched subsampling)',
             fontsize=13, fontweight='bold', y=1.0)

out = os.path.join(FIG_DIR3, 'SuppFigure3_RawDiversity.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
