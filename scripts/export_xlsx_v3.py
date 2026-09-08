#!/usr/bin/env python3
"""Export underlying data for all v3 figures to Excel files."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from openpyxl import Workbook
from scipy.stats import gaussian_kde

from fig_utils import (get_cache, productive_pct, gene_freq_pct, clonotype_ids,
                       cdr3_lengths as fu_cdr3_lengths, cdr3_cores,
                       CHAINS, RABBITS, PLATFORMS, shannon_clonality, is_available)

ROOT     = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR3 = os.path.join(ROOT, "figures", os.environ.get("BCR_FIG_SUBDIR", ""))
cache    = get_cache()
DATA     = cache['data']
FASTAS   = cache['fasta']

AAS      = list('ACDEFGHIKLMNPQRSTVWY')
POS_HALF = 6
POSITIONS = list(range(-POS_HALF, POS_HALF + 1))

# ── Helpers ───────────────────────────────────────────────────────────────────
def gene_freq(platform, rabbit, chain, col):
    return gene_freq_pct(DATA[(platform, rabbit, chain)], col)

def prod_pct(platform, rabbit, chain):
    if not is_available(platform, rabbit, chain):
        return np.nan
    return productive_pct(DATA[(platform, rabbit, chain)])

def vident_arr(platform, rabbit, chain):
    vi = pd.Series(DATA[(platform, rabbit, chain)]['v_identity'], dtype=float)
    if len(vi) and vi.max() <= 1.0:
        vi = vi * 100
    return vi[(vi >= 80) & (vi <= 100)].values

def cdr3_lengths(platform, rabbit, chain):
    return fu_cdr3_lengths(DATA[(platform, rabbit, chain)])

# ── Depth-matched subsampling (must match fig6_diversity.py exactly, so the
#    Excel export reports the same numbers as Figure 6A/B) ────────────────────
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

clono_raw = {}
for p in PLATFORMS:
    for r in RABBITS:
        for ch in CHAINS:
            clono_raw[(p, r, ch)] = clonotype_ids(DATA[(p, r, ch)])

target_depth = {}   # (rabbit, chain) -> shared subsampling depth, needs both platforms
for r in RABBITS:
    for ch in CHAINS:
        if is_available('AVITI', r, ch) and is_available('Nanopore', r, ch):
            target_depth[(r, ch)] = min(len(clono_raw[('AVITI', r, ch)]),
                                         len(clono_raw[('Nanopore', r, ch)]))

def diversity_subsampled(platform, rabbit, chain):
    if (rabbit, chain) not in target_depth:
        return (np.nan, np.nan)
    target = target_depth[(rabbit, chain)]
    return subsample_shannon_clonality(clono_raw[(platform, rabbit, chain)], target)

def clone_3bin(platform, rabbit, chain):
    bins   = [1, 2, 11, 10**9]
    labels = ['Singleton(=1)', 'Small(2-10)', 'Large(>10)']
    if (rabbit, chain) not in target_depth:
        return {lbl: np.nan for lbl in labels}
    target = target_depth[(rabbit, chain)]
    vals = subsample_clone_3bin(clono_raw[(platform, rabbit, chain)], target, bins)
    return dict(zip(labels, vals))

TOP_N = 100

def topN_overlap(rabbit, chain):
    def topN(plat):
        clono = clonotype_ids(DATA[(plat, rabbit, chain)])
        return set(pd.Series(clono).value_counts().head(TOP_N).index)
    av = topN('AVITI'); nano = topN('Nanopore')
    return {'AVITI only': len(av-nano), 'Shared': len(av&nano), 'Nano only': len(nano-av)}

# ── Figure 1 ──────────────────────────────────────────────────────────────────
wb1 = Workbook(); wb1.active.title = '1A_data'
ws  = wb1.active
ws.append(['Chain','Platform','Rabbit','Total_reads','Productive_N','Productive_pct'])
for ch in CHAINS:
    for p in PLATFORMS:
        for r in RABBITS:
            s = DATA[(p,r,ch)]; tot = s['total_count']; pro = s['productive_count']
            ws.append([ch, p, r, tot, pro, round(pro/tot*100,2) if tot else 0])

ws2 = wb1.create_sheet('1B_data')
XR  = np.linspace(300, 700, 200)
ws2.append(['Chain','Platform','Seq_length_bp','Mean_density','Min_density','Max_density'])
for ch in CHAINS:
    for p in PLATFORMS:
        kdes = []
        for r in RABBITS:
            ls = np.array(FASTAS[(p,r,ch)]); ls = ls[(ls>=300)&(ls<=700)]
            if len(ls) >= 20:
                kdes.append(gaussian_kde(ls, bw_method='scott')(XR))
        if not kdes: continue
        kdes = np.array(kdes)
        for xi, x in enumerate(XR):
            ws2.append([ch, p, round(x,1),
                        round(kdes.mean(axis=0)[xi],6),
                        round(kdes.min(axis=0)[xi],6),
                        round(kdes.max(axis=0)[xi],6)])

wb1.save(os.path.join(FIG_DIR3, 'Figure2_QC_v4.xlsx'))
print('Saved Figure2_QC_v4.xlsx')

# ── Figure 2 ──────────────────────────────────────────────────────────────────
wb2 = Workbook(); wb2.active.title = '2A_IGHV'

def gene_sheet(wb, sheet_name, freq_dict, genes, chain=None):
    ws = wb.create_sheet(sheet_name) if sheet_name != wb.active.title else wb.active
    header = ['Gene'] + [f'AVITI_R{r}_pct' for r in RABBITS] + ['AVITI_mean'] \
                       + [f'Nano_R{r}_pct'  for r in RABBITS] + ['Nano_mean']
    ws.append(header)
    def rabbits_ok(plat):
        return [r for r in RABBITS if chain is None or is_available(plat, r, chain)]
    for g in genes:
        row = [g]
        for plat in ('AVITI', 'Nanopore'):
            rb   = rabbits_ok(plat)
            vals = []
            for r in RABBITS:
                fd = freq_dict[(plat,r)] if chain is None else freq_dict[(plat,r,chain)]
                vals.append(round(fd.get(g,0),4) if r in rb else '')
            row += vals
            mean_vals = [v for v,r in zip(vals, RABBITS) if r in rb]
            row.append(round(np.mean(mean_vals),4) if mean_vals else '')
        ws.append(row)

vf2 = {(p,r): gene_freq(p,r,'IGH','v_call') for p in PLATFORMS for r in RABBITS}
df2 = {(p,r): gene_freq(p,r,'IGH','d_call') for p in PLATFORMS for r in RABBITS}
jf2 = {(p,r): gene_freq(p,r,'IGH','j_call') for p in PLATFORMS for r in RABBITS}

def top_n(fd, n):
    genes = set(fd[('AVITI',1)].index)|set(fd[('AVITI',2)].index)
    avg   = {g:(fd[('AVITI',1)].get(g,0)+fd[('AVITI',2)].get(g,0))/2 for g in genes}
    return sorted(avg,key=avg.get,reverse=True)[:n]

def sort_desc(fd):
    genes = set().union(*[set(fd[(p,r)].index) for p in PLATFORMS for r in RABBITS])
    avg   = {g:(fd[('AVITI',1)].get(g,0)+fd[('AVITI',2)].get(g,0))/2 for g in genes}
    return sorted(avg,key=avg.get,reverse=True)

gene_sheet(wb2, '2A_IGHV',   vf2, top_n(vf2,20))
gene_sheet(wb2, '2B_IGHD',   df2, top_n(df2,15))
gene_sheet(wb2, '2C_IGHJ',   jf2, sort_desc(jf2))

ws2d = wb2.create_sheet('2D_Vident')
ws2d.append(['Platform','Rabbit','N_reads','Median_pct','Mean_pct'])
for p in PLATFORMS:
    pooled = np.concatenate([vident_arr(p,r,'IGH') for r in RABBITS])
    ws2d.append([p, 'Pooled', len(pooled), round(np.median(pooled),2), round(np.mean(pooled),2)])
    for r in RABBITS:
        vi = vident_arr(p,r,'IGH')
        ws2d.append([p, r, len(vi), round(np.median(vi),2) if len(vi) else '',
                     round(np.mean(vi),2) if len(vi) else ''])

ws2e = wb2.create_sheet('2E_scatter')
ws2e.append(['V_gene','AVITI_R1_pct','AVITI_R2_pct','Nano_R1_pct','Nano_R2_pct'])
all_g = sorted(set().union(*[set(vf2[(p,r)].index) for p in PLATFORMS for r in RABBITS]))
for g in all_g:
    ws2e.append([g, round(vf2[('AVITI',1)].get(g,0),4), round(vf2[('AVITI',2)].get(g,0),4),
                 round(vf2[('Nanopore',1)].get(g,0),4), round(vf2[('Nanopore',2)].get(g,0),4)])

wb2.save(os.path.join(FIG_DIR3, 'Figure3_IgH_v4.xlsx'))
print('Saved Figure3_IgH_v4.xlsx')

# ── Figure 3 ──────────────────────────────────────────────────────────────────
wb3 = Workbook(); wb3.active.title = '3A_IGKV'
KL  = ['IGK','IGL']
vf3 = {(p,r,ch): gene_freq(p,r,ch,'v_call') for p in PLATFORMS for r in RABBITS for ch in KL}
jf3 = {(p,r,ch): gene_freq(p,r,ch,'j_call') for p in PLATFORMS for r in RABBITS for ch in KL}

for sh, ch, col in [('3A_IGKV','IGK','v_call'), ('3B_IGKJ','IGK','j_call'),
                     ('3C_IGLV','IGL','v_call'), ('3D_IGLJ','IGL','j_call')]:
    fd    = vf3 if col=='v_call' else jf3
    genes_all = set(fd[('AVITI',1,ch)].index)|set(fd[('AVITI',2,ch)].index)
    avg   = {g:(fd[('AVITI',1,ch)].get(g,0)+fd[('AVITI',2,ch)].get(g,0))/2 for g in genes_all}
    if col=='v_call': genes = sorted(avg,key=avg.get,reverse=True)[:20]
    else:             genes = sorted(avg,key=avg.get,reverse=True)
    if sh == '3A_IGKV':
        wb3.active.title = sh
    gene_sheet(wb3, sh, fd, genes, chain=ch)

for sh, ch in [('3E_IGKscatter','IGK'),('3F_IGLscatter','IGL')]:
    ws = wb3.create_sheet(sh)
    header = ['V_gene'] + [f'AVITI_R{r}_pct' for r in RABBITS] + [f'Nano_R{r}_pct' for r in RABBITS]
    ws.append(header)
    rb_nano = [r for r in RABBITS if is_available('Nanopore', r, ch)]
    all_g3 = sorted(set().union(*[set(vf3[(p,r,ch)].index) for p in PLATFORMS for r in RABBITS]))
    for g in all_g3:
        row = [g] + [round(vf3[('AVITI',r,ch)].get(g,0),4) for r in RABBITS] \
                   + [round(vf3[('Nanopore',r,ch)].get(g,0),4) if r in rb_nano else '' for r in RABBITS]
        ws.append(row)

ws3g = wb3.create_sheet('3G_Vident')
ws3g.append(['Chain','Platform','Rabbit','N_reads','Median_pct','Mean_pct'])
for ch in KL:
    for p in PLATFORMS:
        pooled = np.concatenate([vident_arr(p,r,ch) for r in RABBITS])
        ws3g.append([ch,p,'Pooled',len(pooled),
                     round(np.median(pooled),2),round(np.mean(pooled),2)])
        for r in RABBITS:
            vi = vident_arr(p,r,ch)
            ws3g.append([ch,p,r,len(vi),
                         round(np.median(vi),2) if len(vi) else '',
                         round(np.mean(vi),2)   if len(vi) else ''])

wb3.save(os.path.join(FIG_DIR3, 'Figure4_IgKL_v4.xlsx'))
print('Saved Figure4_IgKL_v4.xlsx')

# ── Figure 4 ──────────────────────────────────────────────────────────────────
wb4 = Workbook(); wb4.active.title = '4A_CDR3len'
ws4a = wb4.active
XR_CDR3 = np.arange(3, 25)
ws4a.append(['Chain','Platform','CDR3_length_aa','Mean_freq_pct','Min_freq_pct','Max_freq_pct'])
for ch in CHAINS:
    for p in PLATFORMS:
        curves = []
        for r in RABBITS:
            ls = cdr3_lengths(p, r, ch)
            if len(ls) == 0: continue
            counts = pd.Series(ls).value_counts().reindex(XR_CDR3, fill_value=0)
            curves.append((counts / counts.sum() * 100).values)
        if not curves: continue
        curves = np.array(curves)
        for xi, length in enumerate(XR_CDR3):
            ws4a.append([ch, p, int(length),
                         round(curves.mean(axis=0)[xi],4),
                         round(curves.min(axis=0)[xi],4),
                         round(curves.max(axis=0)[xi],4)])

ws4b = wb4.create_sheet('4B_heatmap')
ws4b.append(['Platform','CDR3_position','Amino_acid','Frequency'])
for p in PLATFORMS:
    all_cores = []
    for r in RABBITS:
        all_cores.extend(cdr3_cores(DATA[(p, r, 'IGH')]))
    counts = {pos: {aa: 0 for aa in AAS} for pos in POSITIONS}
    totpos = {pos: 0 for pos in POSITIONS}
    for seq in all_cores:
        L = len(seq)
        for i, aa in enumerate(seq):
            pos = i - (L-1)//2
            if pos in counts and aa in counts[pos]:
                counts[pos][aa] += 1; totpos[pos] += 1
    for pos in POSITIONS:
        for aa in AAS:
            freq = counts[pos][aa]/totpos[pos] if totpos[pos]>0 else 0
            ws4b.append([p, pos, aa, round(freq,6)])

wb4.save(os.path.join(FIG_DIR3, 'Figure5_CDR3_v4.xlsx'))
print('Saved Figure5_CDR3_v4.xlsx')

# ── Figure 5 ──────────────────────────────────────────────────────────────────
wb5 = Workbook(); wb5.active.title = '5A_diversity'
ws5a = wb5.active
ws5a.append(['Chain','Platform','Rabbit','Shannon_H','Clonality',
             'Note: subsampled to equal AVITI/Nanopore depth (n_iter=20), matches Figure 6A'])
for ch in CHAINS:
    for p in PLATFORMS:
        for r in RABBITS:
            if (r, ch) not in target_depth:
                ws5a.append([ch, p, r, '', ''])
                continue
            H, cl = diversity_subsampled(p, r, ch)
            ws5a.append([ch, p, r, H, cl])

ws5b = wb5.create_sheet('5B_clonesize')
ws5b.append(['Chain','Platform','Rabbit','Singleton_pct','Small_pct','Large_pct',
             'Note: subsampled to equal AVITI/Nanopore depth (n_iter=20), matches Figure 6B'])
for ch in CHAINS:
    for p in PLATFORMS:
        for r in RABBITS:
            if (r, ch) not in target_depth:
                ws5b.append([ch, p, r, '', '', ''])
                continue
            d = clone_3bin(p, r, ch)
            ws5b.append([ch, p, r,
                         round(d['Singleton(=1)'],4),
                         round(d['Small(2-10)'],4),
                         round(d['Large(>10)'],4)])

ws5c = wb5.create_sheet('5C_overlap')
ws5c.append(['Chain','Rabbit','AVITI_only','Shared','Nano_only'])
for ch in CHAINS:
    for r in RABBITS:
        if not (is_available('AVITI', r, ch) and is_available('Nanopore', r, ch)):
            ws5c.append([ch, r, '', '', ''])
            continue
        d = topN_overlap(r, ch)
        ws5c.append([ch, r, d['AVITI only'], d['Shared'], d['Nano only']])

wb5.save(os.path.join(FIG_DIR3, 'Figure6_Diversity_v4.xlsx'))
print('Saved Figure6_Diversity_v4.xlsx')

# ── Figure 7 (Consistency — 從 n=5 cache 重新計算，不再複製舊 v2 檔案) ─────────
def rabbits_with_both(chain):
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

vf7 = {(p,r,ch): gene_freq(p,r,ch,'v_call') for p in PLATFORMS for r in RABBITS for ch in CHAINS}

wb7 = Workbook(); wb7.active.title = '7A_BlandAltman'
ws7a = wb7.active
ws7a.append(['Chain','Rabbit','V_gene','AVITI_pct','Nano_pct','Mean_pct','Diff_Nano_minus_AVITI'])
for ch in CHAINS:
    rb    = rabbits_with_both(ch)
    all_g = sorted(set().union(*[set(vf7[(p,r,ch)].index) for p in PLATFORMS for r in rb]))
    for r in rb:
        for g in all_g:
            av   = vf7[('AVITI',r,ch)].get(g,0)
            nano = vf7[('Nanopore',r,ch)].get(g,0)
            if av == 0 and nano == 0: continue
            ws7a.append([ch, r, g, round(av,4), round(nano,4),
                         round((av+nano)/2,4), round(nano-av,4)])

ws7b = wb7.create_sheet('7B_ICC')
ws7b.append(['Chain','V_gene','AVITI_mean_pct','ICC_2_1','N_rabbits_paired'])
for ch in CHAINS:
    rb    = rabbits_with_both(ch)
    all_g = sorted(set().union(*[set(vf7[(p,r,ch)].index) for p in PLATFORMS for r in rb]))
    aviti_mean_freq = {g: np.mean([vf7[('AVITI',r,ch)].get(g,0) for r in rb]) for g in all_g}
    genes = sorted(all_g, key=lambda g: aviti_mean_freq[g], reverse=True)
    for g in genes:
        icc = icc21([vf7[('AVITI',r,ch)].get(g,0) for r in rb],
                    [vf7[('Nanopore',r,ch)].get(g,0) for r in rb])
        ws7b.append([ch, g, round(aviti_mean_freq[g],4),
                     round(icc,4) if not np.isnan(icc) else '', len(rb)])

wb7.save(os.path.join(FIG_DIR3, 'Figure7_Consistency_v4.xlsx'))
print(f'Saved Figure7_Consistency_v4.xlsx (recomputed from n={len(RABBITS)} data)')

print('\nAll v3 xlsx files exported.')
