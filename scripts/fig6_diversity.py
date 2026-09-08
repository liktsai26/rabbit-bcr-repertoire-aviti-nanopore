#!/usr/bin/env python3
"""Figure 5 (v3): Diversity — Plan B bars (5A), clone size (5B), bar+dots (5C), rarefaction (5D)"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from fig_utils import (get_cache, FIG_DIR,
                       CHAINS, RABBITS, PLATFORMS, COLORS, LABELS,
                       clonotype_ids, shannon_clonality, is_available)

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

def subsample_shannon_clonality(clono_arr, target_n, n_iter=20, seed=42):
    rng = np.random.default_rng(seed)
    n = len(clono_arr)
    target_n = min(n, target_n)
    Hs, Cls = [], []
    for _ in range(n_iter):
        idx = rng.choice(n, size=target_n, replace=False)
        h, cl = shannon_clonality(clono_arr[idx])
        Hs.append(h); Cls.append(cl)
    return float(np.mean(Hs)), float(np.mean(Cls))

def subsample_clone_3bin(clono_arr, target_n, bins, n_iter=20, seed=42):
    rng = np.random.default_rng(seed)
    n = len(clono_arr)
    target_n = min(n, target_n)
    all_pcts = []
    for _ in range(n_iter):
        idx    = rng.choice(n, size=target_n, replace=False)
        counts = pd.Series(clono_arr[idx]).value_counts()
        total  = len(counts)
        pcts   = [((counts >= lo) & (counts < hi)).sum() / total * 100 if total else 0
                  for lo, hi in zip(bins[:-1], bins[1:])]
        all_pcts.append(pcts)
    return list(np.mean(all_pcts, axis=0))

# ── 5A: Diversity metrics (subsampled to equal depth) ─────────────────────────
print("Computing subsampled diversity (n_iter=20)...")
clono_raw = {}
for p in PLATFORMS:
    for r in RABBITS:
        for ch in CHAINS:
            clono_raw[(p, r, ch)] = clonotype_ids(DATA[(p, r, ch)])

diversity    = {}
target_depth = {}   # (rabbit, chain) -> shared subsampling depth, reused by 5B
for p in PLATFORMS:
    for r in RABBITS:
        for ch in CHAINS:
            # 兩平台配對比較才有意義；任一平台缺資料時整組標 NaN，
            # 避免另一平台的真實數據被強制拉去跟 0 配對
            if not (is_available('AVITI', r, ch) and is_available('Nanopore', r, ch)):
                diversity[(p, r, ch)] = (np.nan, np.nan)
                continue
            n_aviti  = len(clono_raw[('AVITI',    r, ch)])
            n_nano   = len(clono_raw[('Nanopore', r, ch)])
            target   = min(n_aviti, n_nano)
            target_depth[(r, ch)] = target
            diversity[(p, r, ch)] = subsample_shannon_clonality(
                clono_raw[(p, r, ch)], target)

# ── 5B: Clone size distribution ───────────────────────────────────────────────
BINS       = [1, 2, 3, 6, 11, 51, 101, 10**9]
BIN_LABELS = ['=1\n(singleton)', '=2', '3–5', '6–10', '11–50', '51–100', '>100']
GROUP_ORDER  = [('AVITI',1), ('Nanopore',1), ('AVITI',2), ('Nanopore',2)]
GROUP_LABELS = ['R1 AVITI', 'R1 Nano', 'R2 AVITI', 'R2 Nano']
GROUP_COLORS = [COLORS[k] for k in GROUP_ORDER]

def clone_size_dist(platform, rabbit, chain):
    df_p   = get_productive(DATA[(platform, rabbit, chain)])
    clono  = get_clonotypes(df_p)
    counts = pd.Series(clono).value_counts()
    total  = len(counts)
    result = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        n = int(((counts >= lo) & (counts < hi)).sum())
        result.append(n / total * 100 if total else 0)
    return result

# ── 5C: Top-100 clonotype overlap — IGH/IGK/IGL x 5 rabbits ───────────────────
# 注意：IGL 的 Nanopore repertoire 本身高度寡株化，部分兔子獨特 clonotype 總數
# 僅 ~116-356 個，Top-100 已涵蓋其中多數，並非只挑「最主力」的一小撮——這是
# IGL 生物學特性的真實反映，不是分析方法的問題，但解讀 IGL 這排時要注意。
TOP_N = 100

def topN(platform, rabbit, chain):
    clono = clonotype_ids(DATA[(platform, rabbit, chain)])
    return set(pd.Series(clono).value_counts().head(TOP_N).index)

overlap = {}
for ch in CHAINS:
    for r in RABBITS:
        if not (is_available('AVITI', r, ch) and is_available('Nanopore', r, ch)):
            overlap[(r, ch)] = None  # 缺資料，該格畫成「無資料」
            continue
        av   = topN('AVITI',    r, ch)
        nano = topN('Nanopore', r, ch)
        overlap[(r, ch)] = {'AVITI only': len(av - nano),
                             'Shared':     len(av & nano),
                             'Nano only':  len(nano - av)}

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 23))
gs  = GridSpec(3, 1, figure=fig, hspace=0.25, height_ratios=[0.75, 0.8, 1.05])

# ── 5A: Shannon H + Clonality — Plan B bars ──────────────────────────────────
gs_a    = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], wspace=0.35)
METRICS = [('Shannon H', 0), ('Clonality', 1)]
x_pos   = np.arange(len(CHAINS))

for mi, (metric_name, metric_idx) in enumerate(METRICS):
    ax = fig.add_subplot(gs_a[mi])
    for pi, plat in enumerate(PLATFORMS):
        col  = C_AVITI if plat == 'AVITI' else C_NANO
        off  = -W/2 + pi * W
        vals_by_ch = [[diversity[(plat, r, ch)][metric_idx] for r in RABBITS] for ch in CHAINS]
        means = [np.nanmean(v) for v in vals_by_ch]
        stds  = [np.nanstd(v)  for v in vals_by_ch]
        ax.bar(x_pos + off, means, W, color=col, alpha=0.85, label=plat,
               edgecolor='white', linewidth=0.3,
               yerr=stds, capsize=3,
               error_kw={'elinewidth': 0.9, 'ecolor': '#333333'})
        for r in RABBITS:
            vals = [diversity[(plat, r, ch)][metric_idx] for ch in CHAINS]
            ax.scatter(x_pos + off, vals, color='black', s=22, zorder=4,
                      edgecolors='white', linewidths=0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(CHAINS, fontsize=10)
    ax.set_ylabel(metric_name, fontsize=10)
    ax.set_title(metric_name + '\n(subsampled to equal depth)', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    if mi == 0:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore')],
                  fontsize=8, frameon=False)

fig.text(0.01, gs[0].get_position(fig).y1 + 0.004,
         'A', fontsize=16, fontweight='bold', va='bottom')

# ── 5B: Clone size — Plan B, 3 bins ──────────────────────────────────────────
BINS_B1   = [1, 2, 11, 10**9]
LABELS_B1 = ['Singleton\n(=1)', 'Small\n(2–10)', 'Large\n(>10)']

def clone_3bin(platform, rabbit, chain):
    # 比照 5A：任一平台缺資料時兩邊都標 NaN，且用兩平台共用的抽樣深度
    # （target_depth，等同 Nanopore 的深度）重新分類，避免深度較深的 AVITI
    # 用原始讀數算 Large 比例時，把單純的讀數深度優勢誤判成生物學上的克隆擴增。
    if not (is_available('AVITI', rabbit, chain) and is_available('Nanopore', rabbit, chain)):
        return [np.nan] * (len(BINS_B1) - 1)
    target = target_depth[(rabbit, chain)]
    return subsample_clone_3bin(clono_raw[(platform, rabbit, chain)], target, BINS_B1)

gs_b = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.35)

for ci, chain in enumerate(CHAINS):
    ax = fig.add_subplot(gs_b[ci])
    xb = np.arange(len(LABELS_B1))
    for pi, plat in enumerate(PLATFORMS):
        col       = C_AVITI if plat == 'AVITI' else C_NANO
        off       = -W/2 + pi * W
        all_r_vals = np.array([clone_3bin(plat, r, chain) for r in RABBITS])
        means     = np.nanmean(all_r_vals, axis=0)
        stds      = np.nanstd(all_r_vals, axis=0)
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
    ax.set_title(chain + '\n(subsampled to equal depth)', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    if ci == 0:
        ax.legend(handles=[mpatches.Patch(color=C_AVITI, label='AVITI'),
                           mpatches.Patch(color=C_NANO,  label='Nanopore')],
                  fontsize=8, frameon=False)

fig.text(0.01, gs[1].get_position(fig).y1 + 0.004,
         'B', fontsize=16, fontweight='bold', va='bottom')

# ── 5C: Venn diagrams — N 隻兔子(橫) x 3 chain(縱，IGH/IGK/IGL) ────────────────
gs_c = GridSpecFromSubplotSpec(len(CHAINS), len(RABBITS), subplot_spec=gs[2],
                                wspace=0.03, hspace=0.04)

def draw_venn(ax, av_only, shared, nano_only, show_platform_labels):
    ax.set_xlim(1.0, 8.9); ax.set_ylim(0.5, 5.9); ax.set_aspect('equal')
    ax.axis('off')
    c_av   = Circle((3.6, 3), 2.2, color=C_AVITI,  alpha=0.30, zorder=2)
    c_nano = Circle((6.4, 3), 2.2, color=C_NANO,    alpha=0.30, zorder=2)
    ax.add_patch(c_av); ax.add_patch(c_nano)
    ax.text(2.0, 3,   str(av_only),  ha='center', va='center',
            fontsize=12, fontweight='bold', zorder=3)
    ax.text(5.0, 3,   str(shared),   ha='center', va='center',
            fontsize=12, fontweight='bold', zorder=3)
    ax.text(8.0, 3,   str(nano_only),ha='center', va='center',
            fontsize=12, fontweight='bold', zorder=3)
    if show_platform_labels:
        ax.text(2.8, 5.6, 'AVITI',    ha='center', va='center',
                fontsize=8, color=C_AVITI,  fontweight='bold')
        ax.text(7.2, 5.6, 'Nanopore', ha='center', va='center',
                fontsize=8, color=C_NANO,   fontweight='bold')

def draw_no_data(ax):
    ax.set_xlim(1.0, 8.9); ax.set_ylim(0.5, 5.9)
    ax.axis('off')
    ax.text(5, 3, 'No data', ha='center', va='center',
            fontsize=9, color='grey', style='italic')

for ci, chain in enumerate(CHAINS):
    for ri, r in enumerate(RABBITS):
        ax_v = fig.add_subplot(gs_c[ci, ri])
        d = overlap[(r, chain)]
        if d is None:
            draw_no_data(ax_v)
        else:
            draw_venn(ax_v, d['AVITI only'], d['Shared'], d['Nano only'],
                      show_platform_labels=(ci == 0))
        if ci == 0:
            ax_v.set_title(f'Rabbit {r}', fontsize=10, fontweight='bold')
    # chain 標籤放在該列最左邊
    pos = gs_c[ci, 0].get_position(fig)
    fig.text(pos.x0 - 0.012, (pos.y0 + pos.y1) / 2, f'{chain}\nTop-{TOP_N}',
              ha='right', va='center', fontsize=10, fontweight='bold')

fig.text(0.01, gs[2].get_position(fig).y1 + 0.004,
         'C', fontsize=16, fontweight='bold', va='bottom')

out = os.path.join(FIG_DIR3, 'Figure6_Diversity_v4.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
