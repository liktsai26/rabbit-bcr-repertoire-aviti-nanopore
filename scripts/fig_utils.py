#!/usr/bin/env python3
"""Shared utilities, paths, and cache for BCR pipeline figures."""

import os, pickle, random
import numpy as np
import pandas as pd

ROOT       = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
FIG_DIR    = os.path.join(ROOT, "figures")
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# BCR_CACHE_NAME：覆蓋 cache 檔名（例如跑實驗版本時用 fig_cache_V4.pkl，避免覆蓋正式版 fig_cache.pkl）
CACHE_PATH = os.path.join(SCRIPT_DIR, os.environ.get("BCR_CACHE_NAME", "fig_cache.pkl"))
os.makedirs(FIG_DIR, exist_ok=True)

# BCR_NANO_RESULTS_SUBDIR：覆蓋 Nanopore AIRR/FASTA 的 results 子資料夾名稱（預設 "results"）。
# 只影響 Nanopore 路徑，AVITI 路徑不受影響——用於實驗性重跑（例如改 cutadapt 參數）時，
# 讓新結果寫到平行資料夾（如 "results_v4"），不覆蓋正式版資料。
NANO_RESULTS_SUBDIR = os.environ.get("BCR_NANO_RESULTS_SUBDIR", "results")

CHAINS    = ['IGH', 'IGK', 'IGL']
# 2026-07-16：n=3（rabbit 1-3）確認為最終正式版，預設改成 n=3；
# 可用環境變數 BCR_RABBITS（逗號分隔，如 "1,2,3,4,5"）覆蓋，供重新產生 n=5 存檔版用
_rabbits_override = os.environ.get('BCR_RABBITS')
RABBITS   = [int(x) for x in _rabbits_override.split(',')] if _rabbits_override else [1, 2, 3]
PLATFORMS = ['AVITI', 'Nanopore']

# rabbit -> Nanopore 資料夾與樣品前綴（AVITI 一律是 15651-LT-{rabbit}，命名規則一致不需對照表）
RABBIT_NANOPORE = {
    1: {'dir': 'MZPD8T_Nanopore', 'prefix': 'MZPD8T', 'label': 'WT'},
    2: {'dir': 'MZPD8T_Nanopore', 'prefix': 'MZPD8T', 'label': 'Immu'},
    3: {'dir': 'QYWFSM_Nanopore', 'prefix': 'QYWFSM', 'label': 'Bov1'},
    4: {'dir': 'QYWFSM_Nanopore', 'prefix': 'QYWFSM', 'label': 'Bov2'},
    5: {'dir': 'QYWFSM_Nanopore', 'prefix': 'QYWFSM', 'label': 'Bov3'},
}

# 已知缺失的資料組合（原始定序失敗，之後若補測序再從這裡移除）
# ('Nanopore', 3, 'IGK') 已於 2026-07-16 補測序（VVR5GH_1_Bov1-IGK.fastq，4,461 reads）並跑完 igblast，移除
KNOWN_MISSING = set()

COLS_NEED = [
    'productive', 'stop_codon', 'vj_in_frame',
    'v_call', 'd_call', 'j_call', 'v_identity', 'junction_aa',
]

# ── Colour palette ──────────────────────────────────────────────────────────────
# 每隻兔子一個 hue（categorical，CVD-safe 固定順序），AVITI=深色、Nanopore=淺色（~45% 向白混合）
COLORS = {
    ('AVITI',    1): '#1F6EA8', ('Nanopore', 1): '#7CB9E8',   # blue
    ('AVITI',    2): '#E07B39', ('Nanopore', 2): '#F5B68A',   # orange
    ('AVITI',    3): '#1BAF7A', ('Nanopore', 3): '#81D3B5',   # aqua
    ('AVITI',    4): '#4A3AA7', ('Nanopore', 4): '#9B92CE',   # violet
    ('AVITI',    5): '#E87BA4', ('Nanopore', 5): '#F2B6CC',   # magenta
}
LABELS = {
    ('AVITI',    r): f'R{r} AVITI'    for r in RABBITS
} | {
    ('Nanopore', r): f'R{r} Nanopore' for r in RABBITS
}
PLAT_COLORS = {'AVITI': '#1F6EA8', 'Nanopore': '#7CB9E8'}

# ── Path helpers ───────────────────────────────────────────────────────────────
def _nanopore_sample(rabbit, chain):
    info = RABBIT_NANOPORE[rabbit]
    d    = 'IGHV' if chain == 'IGH' else chain
    s    = f"{info['prefix']}_{info['label']}_{d}"
    return info['dir'], s, d

def airr_path(platform, rabbit, chain):
    if platform == 'AVITI':
        s = f"15651-LT-{rabbit}_AVITI_{chain}"
        return os.path.join(ROOT, f"15651-LT-{rabbit}_Aviti24",
                            "results", s, f"{s}_igblast.airr.tsv")
    nano_dir, s, _ = _nanopore_sample(rabbit, chain)
    return os.path.join(ROOT, nano_dir, NANO_RESULTS_SUBDIR, s, f"{s}_igblast.airr.tsv")

def fasta_path(platform, rabbit, chain):
    if platform == 'AVITI':
        s = f"15651-LT-{rabbit}_AVITI_{chain}"
        return os.path.join(ROOT, f"15651-LT-{rabbit}_Aviti24",
                            "results", s, f"{chain}_trimmed.fasta")
    nano_dir, s, d = _nanopore_sample(rabbit, chain)
    return os.path.join(ROOT, nano_dir, NANO_RESULTS_SUBDIR, s, f"{d}_trimmed.fasta")

def is_available(platform, rabbit, chain):
    """資料是否存在（排除已知缺失組合，並實際確認 AIRR 檔案存在）。"""
    if (platform, rabbit, chain) in KNOWN_MISSING:
        return False
    return os.path.exists(airr_path(platform, rabbit, chain))

# ── Gene call cleaning ─────────────────────────────────────────────────────────
def clean_gene_call(call_str):
    """'M93173|IGHV1S1*01|Oryctolagus,IGHV1S2*01' -> 'IGHV1S1*01'"""
    if pd.isna(call_str) or str(call_str).strip() == '':
        return None
    first = str(call_str).split(',')[0].strip()
    parts = first.split('|')
    return parts[1].strip() if len(parts) >= 2 else first

# ── Productive mask (3-tier fallback) ─────────────────────────────────────────
def productive_mask(df):
    def flag(col):
        return df[col].astype(str).str.strip() if col in df.columns \
               else pd.Series([''] * len(df), index=df.index)
    prod = flag('productive')
    if prod.isin(['T', 'F']).any():
        return prod == 'T'
    sc, vjif = flag('stop_codon'), flag('vj_in_frame')
    if vjif.isin(['T', 'F']).any():
        return (sc == 'F') & (vjif == 'T')
    return sc == 'F'

# ── FASTA length sampler ───────────────────────────────────────────────────────
def sample_fasta(path, n=50_000, seed=42):
    if not os.path.exists(path):
        print(f"  WARNING: FASTA not found: {path}")
        return []
    random.seed(seed)
    lengths, cur = [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur:
            lengths.append(cur)
    return random.sample(lengths, n) if len(lengths) > n else lengths

# ── Single file loader (chunked, compact summary — 8GB RAM 跑得動) ──────────────
# 不再把整份多欄位 DataFrame 存進 cache。改成邊讀邊算：基因 call 只留「已彙總的次數」
# （V/D/J 基因種類有限，彙總後很小），真正需要逐筆留著的只有 junction_aa（CDR3 序列，
# fig4 的胺基酸矩陣、fig5 的 clonotype 判定都需要精確內容）跟 v_identity（fig2/fig3 用，全 chain）。
# 這樣尖峰記憶體只等於「目前正在讀的這一個檔案」，而不是 30 個檔案的總和。
CHUNK_SIZE = 500_000

def _empty_summary():
    return {
        'total_count':      0,
        'productive_count': 0,
        'v_call_counts':    pd.Series(dtype=float),
        'd_call_counts':    pd.Series(dtype=float),
        'j_call_counts':    pd.Series(dtype=float),
        'v_identity':       np.array([], dtype=np.float32),
        'junction_aa':      np.array([], dtype=object),
    }

def _load_summary(platform, rabbit, chain):
    if not is_available(platform, rabbit, chain):
        print(f"  [{platform}] Rabbit {rabbit} {chain}  <- 缺資料，略過")
        return _empty_summary()
    path = airr_path(platform, rabbit, chain)
    print(f"  [{platform}] Rabbit {rabbit} {chain}  <- {os.path.basename(path)}（分批讀取中...）")
    avail = pd.read_csv(path, sep='\t', nrows=0).columns.tolist()
    use   = [c for c in COLS_NEED if c in avail]

    total, prod_total = 0, 0
    v_counts = pd.Series(dtype=float)
    d_counts = pd.Series(dtype=float)
    j_counts = pd.Series(dtype=float)
    v_ident_parts, junction_parts = [], []

    for chunk in pd.read_csv(path, sep='\t', usecols=use, low_memory=False,
                              dtype=str, chunksize=CHUNK_SIZE):
        total += len(chunk)
        for col in ['v_call', 'd_call', 'j_call']:
            if col in chunk.columns:
                chunk[col] = chunk[col].apply(clean_gene_call)
        if 'v_identity' in chunk.columns:
            chunk['v_identity'] = pd.to_numeric(chunk['v_identity'], errors='coerce')

        prod = chunk[productive_mask(chunk)]
        prod_total += len(prod)

        if 'v_call' in prod.columns:
            v_counts = v_counts.add(prod['v_call'].dropna().value_counts(), fill_value=0)
        if 'd_call' in prod.columns:
            d_counts = d_counts.add(prod['d_call'].dropna().value_counts(), fill_value=0)
        if 'j_call' in prod.columns:
            j_counts = j_counts.add(prod['j_call'].dropna().value_counts(), fill_value=0)
        if 'v_identity' in prod.columns:
            v_ident_parts.append(prod['v_identity'].dropna().values.astype(np.float32))
        if 'junction_aa' in prod.columns:
            junction_parts.append(prod['junction_aa'].values.astype(object))

    return {
        'total_count':      total,
        'productive_count': prod_total,
        'v_call_counts':    v_counts,
        'd_call_counts':    d_counts,
        'j_call_counts':    j_counts,
        'v_identity':       np.concatenate(v_ident_parts) if v_ident_parts else np.array([], dtype=np.float32),
        'junction_aa':      np.concatenate(junction_parts) if junction_parts else np.array([], dtype=object),
    }

# ── Cache build / load ─────────────────────────────────────────────────────────
def build_cache():
    print("=" * 60)
    print("Building compact summary cache from AIRR TSV files (chunked, low-memory)...")
    data, fasta = {}, {}
    missing = []
    for platform in PLATFORMS:
        for rabbit in RABBITS:
            for chain in CHAINS:
                key        = (platform, rabbit, chain)
                data[key]  = _load_summary(platform, rabbit, chain)
                if is_available(platform, rabbit, chain):
                    fasta[key] = sample_fasta(fasta_path(platform, rabbit, chain))
                else:
                    fasta[key] = []
                    missing.append(key)
    cache = {'data': data, 'fasta': fasta, 'missing': missing}
    with open(CACHE_PATH, 'wb') as f:
        pickle.dump(cache, f)
    print(f"Cache saved -> {CACHE_PATH}")
    if missing:
        print(f"缺資料組合（{len(missing)} 個）：{missing}\n")
    return cache

def get_cache():
    if os.path.exists(CACHE_PATH):
        print(f"Loading cache from {CACHE_PATH} ...")
        with open(CACHE_PATH, 'rb') as f:
            return pickle.load(f)
    return build_cache()

# ── Accessors on the compact summary dict ───────────────────────────────────────
def productive_pct(summary):
    tot = summary['total_count']
    return summary['productive_count'] / tot * 100 if tot else np.nan

def gene_freq_pct(summary, gene):
    """gene: 'v_call' / 'd_call' / 'j_call' -- % of productive reads per gene call."""
    counts = summary[f'{gene}_counts']
    total  = counts.sum()
    return counts / total * 100 if total else pd.Series(dtype=float)

def clonotype_ids(summary):
    """每條 productive read 的 clonotype ID：有 junction_aa 用 junction_aa，沒有的話用該筆在
    陣列中的位置產生唯一識別碼——效果等同舊版用 sequence_id 當 fallback（反正原本目的就是
    給每條沒有 CDR3 的 read 一個獨一無二的 ID，讓它自成一個 singleton clone）。"""
    junction = summary['junction_aa']
    if len(junction) == 0:
        return np.array([], dtype=object)
    ja = pd.Series(junction)
    has_cdr3 = ja.notna() & (ja.astype(str).str.strip() != '')
    result = ja.copy()
    if (~has_cdr3).any():
        result[~has_cdr3] = ['_NOJUNC_' + str(i) for i in np.flatnonzero((~has_cdr3).values)]
    return result.values

def cdr3_lengths(summary):
    """CDR3 AA 長度（扣頭尾各 1 個殘基，沿用原本定義）。"""
    junction = summary['junction_aa']
    if len(junction) == 0:
        return np.array([])
    ja = pd.Series(junction).dropna()
    ja = ja[ja.astype(str).str.strip() != '']
    return (ja.str.len() - 2).clip(lower=0).values

def cdr3_cores(summary):
    """CDR3 AA 序列去掉頭尾各 1 個殘基後的『核心』序列（fig4 胺基酸位置矩陣用）。"""
    junction = summary['junction_aa']
    if len(junction) == 0:
        return []
    ja = pd.Series(junction).dropna()
    ja = ja[ja.astype(str).str.strip() != '']
    cores = ja.apply(lambda x: x[1:-1] if len(x) > 2 else '')
    return cores[cores != ''].tolist()

def shannon_clonality(clono_arr):
    """Returns (Shannon H, Clonality)."""
    counts = pd.Series(clono_arr).value_counts().values.astype(float)
    total  = counts.sum()
    if total == 0:
        return 0.0, 0.0
    p = counts / total
    H = float(-np.sum(p * np.log2(p + 1e-300)))
    N = len(counts)
    clonality = float(1 - H / np.log2(N)) if N > 1 else 0.0
    return round(H, 4), round(max(0.0, clonality), 4)
