#!/usr/bin/env python3
"""Figure 7: Cross-animal Concordance — MA-plot (log2FC) / ICC.
Output: Figure7_Consistency_v4.png"""

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
                       CHAINS, RABBITS, PLATFORMS, COLORS,
                       shannon_clonality, is_available)

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
os.makedirs(FIG_DIR3, exist_ok=True)

cache = get_cache()
DATA  = cache['data']

plt.rcParams.update({'font.family': 'Arial', 'axes.labelsize': 10,
                     'axes.titlesize': 12, 'xtick.labelsize': 8,
                     'ytick.labelsize': 8})

CHAIN_MK  = {'IGH': 'o', 'IGK': 's', 'IGL': '^'}

def vgene_freq(platform, rabbit, chain):
    return gene_freq_pct(DATA[(platform, rabbit, chain)], 'v_call')

vfreq = {(p,r,c): vgene_freq(p,r,c)
         for p in PLATFORMS for r in RABBITS for c in CHAINS}

def rabbits_with_both(chain):
    """兩平台在該 chain 都有資料的兔子（排除如 Bov1_IGK 這種單邊缺口，避免假 0 污染配對統計）。"""
    return [r for r in RABBITS if is_available('AVITI', r, chain) and is_available('Nanopore', r, chain)]

def icc21(aviti_vals, nano_vals):
    ratings = np.column_stack([aviti_vals, nano_vals]).astype(float)
    n, k    = ratings.shape
    if n < 2: return np.nan
    grand = ratings.mean(); subj_m = ratings.mean(axis=1); rater_m = ratings.mean(axis=0)
    SSb = k*np.sum((subj_m-grand)**2); SSr = n*np.sum((rater_m-grand)**2)
    SSw = np.sum((ratings-subj_m[:,None])**2); SSe = SSw-SSr
    MSb = SSb/(n-1); MSe = SSe/((n-1)*(k-1)); MSr = SSr/(k-1)
    denom = MSb+(k-1)*MSe+k*(MSr-MSe)/n
    return float(np.clip((MSb-MSe)/denom, 0, 1)) if denom > 0 else np.nan

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 18))
gs  = GridSpec(2, 1, figure=fig, hspace=0.50,
               height_ratios=[1.0, 2.8])

# ── 6A: MA-plot（log2 ratio，取代原本的絕對差值 Bland-Altman）───────────────
# 原始 (Nano-AVITI) 版本有異方差問題：高頻基因的絕對差值天生比低頻基因大，
# 導致 LoA 被高頻基因拉寬、誤判成「高頻不準」。改用 log2 fold-change（MA-plot
# 慣例：A = 兩平台 log2 頻率的平均，M = log2 比值）可以消除這個尺度效應。
#
# 只有「兩平台都偵測到（頻率>0）」的基因才有 well-defined 的 fold-change；
# 只被一個平台偵測到的基因分母是 0、倍率是無限大，不能硬套 pseudocount 塞進
# 連續座標（那樣會做出一條假的對角線——一邊是常數 pseudocount、一邊隨真實值
# 變動，數學上必然排成直線，不是生物訊號）。這類基因改成用計數標註在圖上。

gs_a = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], wspace=0.38)

for ci, chain in enumerate(CHAINS):
    ax   = fig.add_subplot(gs_a[ci])
    rb    = rabbits_with_both(chain)
    all_g = sorted(set().union(*[set(vfreq[(p,r,chain)].index)
                                  for p in PLATFORMS for r in rb]))
    xall, yall = [], []
    n_only_av, n_only_nano = 0, 0
    for rabbit in rb:
        for g in all_g:
            av   = vfreq[('AVITI',    rabbit, chain)].get(g, 0)
            nano = vfreq[('Nanopore', rabbit, chain)].get(g, 0)
            if av > 0 and nano > 0:
                av_log, nano_log = np.log2(av), np.log2(nano)
                xall.append((av_log + nano_log) / 2)   # A：平均 log2 頻率
                yall.append(nano_log - av_log)          # M：log2(Nano/AVITI)
            elif av > 0:
                n_only_av += 1
            elif nano > 0:
                n_only_nano += 1

    ax.scatter(xall, yall, color='#555555', s=22, alpha=0.5, zorder=3)

    if yall:
        mean_m = np.mean(yall)
        sd_m   = np.std(yall)
        loa_hi, loa_lo = mean_m + 1.96*sd_m, mean_m - 1.96*sd_m
        ax.axhline(0, color='grey', lw=0.9, ls='--', alpha=0.6,
                   label='0 (no fold-change)')
        ax.axhline(mean_m, color='black', lw=1.1, ls='-',
                   label=f'Mean log2FC = {mean_m:+.2f} ({2**mean_m:.2f}×)')
        ax.axhline(loa_hi, color='black', lw=0.9, ls=':',
                   label=f'95% LoA = [{2**loa_lo:.2f}×, {2**loa_hi:.2f}×]')
        ax.axhline(loa_lo, color='black', lw=0.9, ls=':')
        ax.legend(fontsize=6.5, frameon=False, loc='upper right')

    ax.text(0.02, 0.02, f'AVITI only: {n_only_av} genes\nNanopore only: {n_only_nano} genes',
            transform=ax.transAxes, fontsize=6.5, color='grey', va='bottom', ha='left')

    ax.set_xlabel('A = mean log2(frequency %)\n[(log2 AVITI + log2 Nano)/2]', fontsize=8.5)
    if ci == 0:
        ax.set_ylabel('M = log2(Nano / AVITI)', fontsize=9)
    ax.set_title(chain, fontsize=12)
    ax.spines[['top','right']].set_visible(False)

fig.text(0.01, gs[0].get_position(fig).y1 + 0.004,
         'A', fontsize=16, fontweight='bold', va='bottom')

# ── 6B: ICC per V gene — 3 full-width rows ────────────────────────────────────
gs_d     = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1], hspace=1.0)
ICC_THRESH = 0.75
ICC_HI   = COLORS[('AVITI',    1)]
ICC_LO   = COLORS[('AVITI',    2)]

for ri, chain in enumerate(CHAINS):
    ax    = fig.add_subplot(gs_d[ri])
    rb    = rabbits_with_both(chain)
    all_g = sorted(set().union(*[set(vfreq[(p,r,chain)].index)
                                   for p in PLATFORMS for r in rb]))
    aviti_mean_freq = {g: np.mean([vfreq[('AVITI', r, chain)].get(g, 0) for r in rb])
                       for g in all_g}
    genes = sorted(all_g, key=lambda g: aviti_mean_freq[g], reverse=True)
    vals   = [0 if np.isnan(
                   icc21([vfreq[('AVITI',r,chain)].get(g,0) for r in rb],
                         [vfreq[('Nanopore',r,chain)].get(g,0) for r in rb]))
              else icc21([vfreq[('AVITI',r,chain)].get(g,0) for r in rb],
                         [vfreq[('Nanopore',r,chain)].get(g,0) for r in rb])
              for g in genes]
    bar_colors = [ICC_HI if v >= ICC_THRESH else ICC_LO for v in vals]
    x = np.arange(len(genes))
    ax.bar(x, vals, color=bar_colors, width=0.7, edgecolor='none')
    ax.axhline(ICC_THRESH, color='black', lw=1.2, ls='--', alpha=0.7)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=90, fontsize=6)
    ax.set_ylabel('ICC(2,1)', fontsize=9)
    ax.set_title(chain, fontsize=12)
    # x 用 axes-fraction（固定貼右緣）、y 用資料座標（跟著虛線），避免蓋到左側高頻基因的長條
    ax.text(0.995, ICC_THRESH + 0.025, 'ICC=0.75', transform=ax.get_yaxis_transform(),
            ha='right', fontsize=7, color='black',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=1.0))
    ax.spines[['top','right']].set_visible(False)
    if ri == 0:
        ax.legend(handles=[mpatches.Patch(color=ICC_HI, label=f'≥{ICC_THRESH}'),
                            mpatches.Patch(color=ICC_LO, label=f'<{ICC_THRESH}')],
                  fontsize=8, frameon=False, loc='lower right')

fig.text(0.01, gs[1].get_position(fig).y1 + 0.004,
         'B', fontsize=16, fontweight='bold', va='bottom')

out = os.path.join(FIG_DIR3, 'Figure7_Consistency_v4.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
