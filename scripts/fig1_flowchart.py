#!/usr/bin/env python3
"""Figure 1 (v3): Pipeline Flowchart — AVITI vs Nanopore"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import to_rgba

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
os.makedirs(FIG_DIR3, exist_ok=True)

C_AVITI  = '#1F6EA8'
C_NANO   = '#E07B39'
C_SHARED = '#444444'
ALPHA    = 0.12

TITLE_FS  = 11
PARAMS_FS = 7.5
LINESPACING = 1.4
PAD_X = 0.018   # 文字外圍留白（軸座標分數）
PAD_Y = 0.014
GAP   = 0.012   # 標題與參數文字間的間距

# ── Layout constants ──────────────────────────────────────────────────────────
X_AV  = 0.25    # AVITI column center (x)
X_NA  = 0.75    # Nanopore column center (x)
X_MID = 0.50    # shared column center (x)
ALW   = 2.0     # arrow linewidth

# Row Y centres (in [0,1] axes fraction, top→bottom)
Y_RAW    = 0.90
Y_DEMUX  = 0.735
Y_MERGE  = 0.575
Y_TRIM   = 0.415
Y_IGB    = 0.235
Y_AIRR   = 0.075

fig, ax = plt.subplots(1, 1, figsize=(12, 10))
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.axis('off')

_renderer = fig.canvas.get_renderer()
_inv      = ax.transAxes.inverted()

def _text_extent_axes(s, fontsize, fontweight='normal'):
    """量測一段文字（可能多行）在 axes-fraction 座標下的寬高。"""
    t = ax.text(0.5, 0.5, s, ha='center', va='center', fontsize=fontsize,
                fontweight=fontweight, linespacing=LINESPACING, alpha=0)
    bbox_px = t.get_window_extent(renderer=_renderer)
    (x0, y0) = _inv.transform((bbox_px.x0, bbox_px.y0))
    (x1, y1) = _inv.transform((bbox_px.x1, bbox_px.y1))
    t.remove()
    return abs(x1 - x0), abs(y1 - y0)

BOX_GEOM = {}  # (xc, yc) -> (bw, bh)，供 arrow() 查詢實際框高

def box(xc, yc, title, params, color, min_bw=0.0):
    tw, th = _text_extent_axes(title,  TITLE_FS,  fontweight='bold')
    pw, ph = _text_extent_axes(params, PARAMS_FS)

    content_w = max(tw, pw)
    content_h = th + GAP + ph
    bw = max(content_w / 2 + PAD_X, min_bw)
    bh = content_h / 2 + PAD_Y

    title_yc  = yc + content_h/2 - th/2
    params_yc = title_yc - th/2 - GAP - ph/2

    patch = FancyBboxPatch((xc - bw, yc - bh), 2*bw, 2*bh,
                           boxstyle='round,pad=0.006',
                           linewidth=2.0, edgecolor=color,
                           facecolor=to_rgba(color, ALPHA), zorder=2)
    ax.add_patch(patch)
    ax.text(xc, title_yc, title,
            ha='center', va='center', fontsize=TITLE_FS, fontweight='bold',
            color='#111111', zorder=3)
    ax.text(xc, params_yc, params,
            ha='center', va='center', fontsize=PARAMS_FS, color='#333333',
            zorder=3, linespacing=LINESPACING)

    BOX_GEOM[(xc, yc)] = (bw, bh)
    return bw, bh


def _bh_at(xc, yc):
    return BOX_GEOM[(xc, yc)][1]


ARROW_GAP = 0.014  # 箭頭與上下框之間各留的間距（軸座標分數），兩端對稱

def arrow(x, y_from, y_to, color='#333333'):
    bh_from = _bh_at(x, y_from)
    bh_to   = _bh_at(x, y_to)
    ax.annotate('', xy=(x, y_to + bh_to + ARROW_GAP),
                xytext=(x, y_from - bh_from - ARROW_GAP),
                arrowprops=dict(arrowstyle='->', color=color, lw=ALW,
                                mutation_scale=14))


def arrow_long(x, y_from, y_to, color='#333333'):
    """Arrow spanning two row levels (Nanopore skip Merge)."""
    arrow(x, y_from, y_to, color)


def converge_arrows(y_from, y_to):
    """Two arrows from Trim level converging into shared IgBLAST box."""
    bh_av  = _bh_at(X_AV, y_from)
    bh_na  = _bh_at(X_NA, y_from)
    bh_mid = _bh_at(X_MID, y_to)
    ax.annotate('', xy=(X_MID - 0.04, y_to + bh_mid + ARROW_GAP),
                xytext=(X_AV, y_from - bh_av - ARROW_GAP),
                arrowprops=dict(arrowstyle='->', color='#333333',
                                lw=ALW, mutation_scale=14))
    ax.annotate('', xy=(X_MID + 0.04, y_to + bh_mid + ARROW_GAP),
                xytext=(X_NA, y_from - bh_na - ARROW_GAP),
                arrowprops=dict(arrowstyle='->', color='#333333',
                                lw=ALW, mutation_scale=14))


# ── AVITI column: 先畫完所有框，箭頭最後統一畫（箭頭需要兩端框的實際大小）──────────
box(X_AV, Y_RAW,   'Raw FASTQ', 'Paired-end\n~12.7 M reads', C_AVITI)
box(X_AV, Y_DEMUX, 'Demux',
    'cutadapt\ne=0.1, overlap ≥15 bp\nIGH / IGK / IGL primers', C_AVITI)
box(X_AV, Y_MERGE, 'Merge',
    'FLASH2\noverlap 20–150 bp\nmismatch ≤0.1', C_AVITI)
box(X_AV, Y_TRIM,  'Trim',
    'cutadapt, e=0.1, Q≥20\nTSO + rev primers\n350–600 bp (IGH) / 320–600 bp (IGK/IGL)', C_AVITI)

arrow(X_AV, Y_RAW,   Y_DEMUX, C_AVITI)
arrow(X_AV, Y_DEMUX, Y_MERGE, C_AVITI)
arrow(X_AV, Y_MERGE, Y_TRIM,  C_AVITI)

# ── Nanopore column ────────────────────────────────────────────────────────────
box(X_NA, Y_RAW,   'Raw FASTQ', 'Long reads\n~5,400 reads', C_NANO)
box(X_NA, Y_DEMUX, 'Filter',
    'NanoFilt, Q≥10\n350–600 bp (IGH)\n320–600 bp (IGK/IGL)', C_NANO)
box(X_NA, Y_TRIM,  'Trim',
    'cutadapt, e=0.1\noverlap ≥15 bp\n--revcomp, --discard-untrimmed', C_NANO)

arrow(X_NA, Y_RAW, Y_DEMUX, C_NANO)
arrow_long(X_NA, Y_DEMUX, Y_TRIM, C_NANO)  # 跳過 Merge 那一列

# ── Shared steps ───────────────────────────────────────────────────────────────
box(X_MID, Y_IGB, 'IgBLAST',
    'igblastn, -outfmt 19\nrabbit germline DB\n-extend_align5end / 3end',
    C_SHARED)
box(X_MID, Y_AIRR, 'AIRR TSV  →  BCR Repertoire Analysis',
    'Productive reads | V/D/J assignment | junction_aa',
    C_SHARED)

converge_arrows(Y_TRIM, Y_IGB)
arrow(X_MID, Y_IGB, Y_AIRR, C_SHARED)

# ── Column headers ─────────────────────────────────────────────────────────────
ax.text(X_AV, 0.975, 'AVITI (Short-read)', ha='center', va='center',
        fontsize=13, fontweight='bold', color=C_AVITI)
ax.text(X_NA, 0.975, 'Nanopore (Long-read)', ha='center', va='center',
        fontsize=13, fontweight='bold', color=C_NANO)

fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
out = os.path.join(FIG_DIR3, 'Figure1_Flowchart_v4.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {out}')
