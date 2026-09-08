#!/usr/bin/env python3
"""Supplementary Figure 4 (renumbered 2026-08-03, was Supplementary Figure 2):
Cost efficiency of Nanopore vs AVITI.

AVITI pricing has a fixed run cost plus a per-rabbit variable cost (confirmed with
user 2026-07-11): $1,872 fixed + $115/rabbit (covers 3 amplicons/chains). Nanopore
is purely variable: $30/amplicon library, with failed libraries not charged.

This avoids the earlier "$/confirmed clone" framing, which was methodologically
asymmetric (Nanopore's cost was conditioned on independent confirmation by AVITI,
while AVITI's own Top-100 clones were unconditionally counted as "hits" against
itself). Cost STRUCTURE alone is a clean, symmetric comparison.

Panel A: absolute total cost vs. number of rabbits (both platforms), with the
         actual data point for THIS run's rabbit count marked (n=5 main run,
         n=3 for the N_3 subset — run this script once per BCR_RABBITS/
         BCR_FIG_SUBDIR context to get each version separately).
Panel B: cost ratio (AVITI/Nanopore) vs. number of rabbits — shows the advantage
         shrinking from the fixed-cost-dominated small-N regime toward the
         asymptotic marginal-cost ratio at large N.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from fig_utils import get_cache, RABBITS, is_available

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
os.makedirs(FIG_DIR3, exist_ok=True)

C_AVITI = '#1F6EA8'
C_NANO  = '#E07B39'

# ── Cost model (confirmed with user 2026-07-11) ───────────────────────────────
AVITI_FIXED       = 1872.0   # 固定跑一次定序的費用，不隨兔子數變動
AVITI_PER_RABBIT  = 115.0    # 每隻兔子(3 個 amplicon)的變動成本
NANO_PER_AMPLICON = 30.0     # 每個 amplicon library 的單價；失敗的 library 不計費
CHAINS_PER_RABBIT = 3

def aviti_cost(n_rabbits):
    return AVITI_FIXED + AVITI_PER_RABBIT * n_rabbits

def nano_cost_theoretical(n_rabbits):
    return NANO_PER_AMPLICON * CHAINS_PER_RABBIT * n_rabbits

# ── 實際發生的成本（有 1 個失敗 library：QYWFSM_Bov1_IGK，不計費）───────────────
cache = get_cache()
DATA  = cache['data']

def actual_nano_cost(rabbit_subset):
    n_success = sum(1 for r in rabbit_subset for ch in ['IGH', 'IGK', 'IGL']
                     if is_available('Nanopore', r, ch))
    return n_success * NANO_PER_AMPLICON

# 只標「這次執行所對應的兔子數」這一個點（n=5 主版跑一次、n=3 用 BCR_RABBITS
# 環境變數再跑一次，各自產生獨立檔案，不要兩個點擠在同一張圖上）。
N_THIS = len(RABBITS)
av_this = aviti_cost(N_THIS)
na_this = actual_nano_cost(RABBITS)
actual_points = [(N_THIS, av_this, na_this)]
print(f'N={N_THIS}: AVITI=${av_this:,.0f}  Nanopore(actual)=${na_this:,.0f}  ratio={av_this/na_this:.2f}x')

# ── Figure ────────────────────────────────────────────────────────────────────
# 用連續點(非整數)讓曲線平滑，尤其 Panel B 在小 N 時變化很陡，用整數點連線會看起來像
# 折線而不是平滑曲線。
N_range = np.linspace(1, 20, 400)
aviti_line = aviti_cost(N_range)
nano_line  = nano_cost_theoretical(N_range)

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.2))

# ── Panel A: absolute cost vs N ────────────────────────────────────────────────
ax_a.plot(N_range, aviti_line, color=C_AVITI, lw=2.2, label=r'AVITI (\$1,872 fixed + \$115/rabbit)')
ax_a.plot(N_range, nano_line,  color=C_NANO,  lw=2.2, label=r'Nanopore (\$30 × 3 amplicons/rabbit)')
n, av, na = actual_points[0]
ax_a.scatter([n], [av], color=C_AVITI, s=80, zorder=5, edgecolors='white', linewidths=1.3)
ax_a.scatter([n], [na], color=C_NANO,  s=80, zorder=5, edgecolors='white', linewidths=1.3)
# 標籤往右下方偏移，避開圖例（圖例固定在左上角）跟自己的線
ax_a.annotate(f'N={n} (actual): \\${av:,.0f}', (n, av), textcoords='offset points', xytext=(22, -6),
              fontsize=9, color=C_AVITI, fontweight='bold', va='top')
ax_a.annotate(f'N={n} (actual): \\${na:,.0f}', (n, na), textcoords='offset points', xytext=(22, -6),
              fontsize=9, color=C_NANO, fontweight='bold', va='top')
ax_a.set_xlabel('Number of rabbits sequenced', fontsize=10)
ax_a.set_ylabel('Total cost (USD)', fontsize=10)
ax_a.set_title(f'A. Total cost vs. study size (this study: n={N_THIS})', fontsize=11)
ax_a.spines[['top', 'right']].set_visible(False)
ax_a.legend(fontsize=8.5, frameon=False, loc='upper left')
ax_a.set_xlim(0, 24)

# ── Panel B: cost ratio vs N ───────────────────────────────────────────────────
ratio_line = aviti_line / nano_line
asymptote  = AVITI_PER_RABBIT / (NANO_PER_AMPLICON * CHAINS_PER_RABBIT)
ax_b.plot(N_range, ratio_line, color='#444444', lw=2.2)
ax_b.axhline(asymptote, color='#999999', lw=1.2, ls=':',
             label=f'Asymptote at large N: {asymptote:.2f}×\n(marginal cost ratio: \\${AVITI_PER_RABBIT:.0f} vs \\${CHAINS_PER_RABBIT*NANO_PER_AMPLICON:.0f} per rabbit)')
for n, av, na in actual_points:
    ratio_actual = av / na
    ax_b.scatter([n], [ratio_actual], color='#222222', s=70, zorder=5, edgecolors='white', linewidths=1.2)
    ax_b.annotate(f'N={n}: {ratio_actual:.1f}×', (n, ratio_actual), textcoords='offset points',
                  xytext=(8, 6), fontsize=9, fontweight='bold')
ax_b.set_xlabel('Number of rabbits sequenced', fontsize=10)
ax_b.set_ylabel('Cost ratio (AVITI / Nanopore)', fontsize=10)
ax_b.set_title(f'B. Nanopore\'s cost advantage shrinks as study size grows (this study: n={N_THIS})', fontsize=10.5)
ax_b.spines[['top', 'right']].set_visible(False)
ax_b.legend(fontsize=8, frameon=False)
ax_b.set_xlim(0, 24)

fig.subplots_adjust(top=0.78, bottom=0.14, wspace=0.35)
fig.suptitle(f"Supplementary Figure 4. Nanopore's cost advantage over AVITI is largest for small "
             f'cohorts and narrows as study size grows (n={N_THIS})',
             fontsize=13, fontweight='bold', y=0.99)
fig.text(0.5, 0.90,
         "AVITI's cost is dominated by a large fixed run cost, amortized over more rabbits as study size grows;\n"
         "Nanopore has no fixed cost, so its relative advantage is largest for small/pilot studies and narrows (but persists) at scale.",
         ha='center', fontsize=8.5, style='italic', color='#444444')

out = os.path.join(FIG_DIR3, 'SuppFigure4_CostEfficiency.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
