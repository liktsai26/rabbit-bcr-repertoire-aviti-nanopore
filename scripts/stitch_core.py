#!/usr/bin/env python3
"""
stitch_core.py
02g_stitch_qc_igh.py 與 04c_stitch_qc_igkl.py 共用的座標拼接+QC邏輯（原本兩份
各自複製一份，這裡合併成一份，避免IGH與IGK/IGL之後再出現差異）。

跟原版的差異（2026-10 code review 修正，見 tests/ 對應測試）：
1. R1、R2 各自的 v_germline 座標只對「各自的」主要V allele有意義。兩端被判給
   不同allele時（例如 IGHV1S40*01 296nt vs IGHV1S45*01 299nt，CDR1/CDR2長度不同），
   原版直接相減會讓切點錯開整數個密碼子，產生框內的重複/缺失（仍判為productive）。
   這裡先把兩條germline序列做全域比對，建立 alleleA 位置 -> alleleB 位置的對照，
   再換算切點。
2. 用 v_cigar 把 germline 位置換成 read 上的位置，read 本身相對germline的插入/缺失
   不再造成切點偏移。
3. AIRR 座標是 1-based 閉區間：R2 起點 = R1 終點 + 1 是「相鄰、無缺口」，原版
   gap 計算差1而把它排除；這裡修正，重疊長度也一併修正。
4. QC：在預測的對位上逐鹼基比對重疊區（不是在兩個視窗裡找最好的局部比對——
   局部比對對錯位不敏感，錯開3個鹼基的拼接也會拿到~100% identity），
   identity < QC_THRESHOLD 排除；比對過程出錯一律排除（原版是一律通過）。
5. 每個 chunk 的 FASTA 先寫暫存檔再改名；chunk_stats.csv 同一 chunk 重跑時取代
   舊列而不是重複追加；被排除的 read 與原因寫到 stitch_rejected/。
"""
import gzip
import os
import re
import sys

import pandas as pd
from Bio.Align import PairwiseAligner

QC_THRESHOLD = 0.90
# 重疊區小於這個長度時，可比對的鹼基太少，不做identity判斷（原版對極短重疊實際上
# 也幾乎都放行）；這類read在 chunk_stats.csv 的 short_overlap 欄位另外計數。
QC_MIN_OVERLAP = int(os.environ.get("BCR_STITCH_QC_MIN_OVERLAP", "10"))
# 允許的最小重疊（0 = 相鄰也接受，與原版「無缺口就保留」的意圖一致）
MIN_OVERLAP = int(os.environ.get("BCR_STITCH_MIN_OVERLAP", "0"))
# 兩條germline allele全域比對 identity 低於此值時，座標換算不可靠，排除
MIN_ALLELE_IDENTITY = 0.70

COORD_COLS = ['v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end']
R1_COLS = ['sequence_id', 'sequence', 'v_call', 'v_cigar'] + COORD_COLS
R2_COLS = ['sequence_id', 'sequence', 'v_call', 'v_cigar', 'productive', 'stop_codon', 'vj_in_frame'] + COORD_COLS
STATS_HEADER = ("chunk,r2_rows,r2_productive,paired_with_r1,nogap,qc_pass,"
                "short_overlap,shift_corrected,gap_rescued,different_allele,rej_gap,rej_qc,rej_unmappable,rej_error")

_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = ''.join(buf).upper()
    return seqs


def germline_to_query(cigar, v_seq_start, v_seq_end, v_germ_start, v_germ_end):
    """AIRR v_cigar -> {germline位置: query位置 or None(該germline鹼基在read中缺失)}，皆為1-based。
    CIGAR 與四個座標不一致時回傳 None。"""
    q = g = 0
    gmap = {}
    for n, op in _CIGAR_RE.findall(cigar or ''):
        n = int(n)
        if op == 'S':
            q += n
        elif op == 'N':
            g += n
        elif op in 'M=X':
            for _ in range(n):
                q += 1
                g += 1
                gmap[g] = q
        elif op == 'I':
            q += n
        elif op == 'D':
            for _ in range(n):
                g += 1
                gmap[g] = None
    aligned = [k for k, v in gmap.items() if v is not None]
    if not aligned:
        return None
    if (min(aligned) != v_germ_start or max(gmap) != v_germ_end
            or gmap[min(aligned)] != v_seq_start or gmap[max(aligned)] != v_seq_end):
        return None
    return gmap


class AlleleMapper:
    """alleleA 的 germline 位置 -> alleleB 的 germline 位置（全域比對，結果快取）。"""

    def __init__(self, germline_fasta):
        self.seqs = read_fasta(germline_fasta)
        self.cache = {}
        al = PairwiseAligner()
        al.mode = 'global'
        al.match_score, al.mismatch_score = 2, -1
        al.open_gap_score, al.extend_gap_score = -8, -1
        al.end_gap_score = -1
        self.aligner = al

    def mapping(self, a, b):
        """回傳 (dict a_pos->b_pos, identity)；a==b 時是恆等對照。找不到allele時回傳 (None, 0)。"""
        if a == b:
            if a not in self.seqs:
                return None, 0.0
            return None, 1.0  # None 代表恆等
        key = (a, b)
        if key in self.cache:
            return self.cache[key]
        if a not in self.seqs or b not in self.seqs:
            self.cache[key] = (None, 0.0)
            return self.cache[key]
        aln = self.aligner.align(self.seqs[a], self.seqs[b])[0]
        amap, same, cols = {}, 0, 0
        for (as_, ae), (bs, be) in zip(*aln.aligned):
            for k in range(ae - as_):
                amap[int(as_) + k + 1] = int(bs) + k + 1
                cols += 1
                same += self.seqs[a][as_ + k] == self.seqs[b][bs + k]
        ident = same / max(len(self.seqs[a]), len(self.seqs[b]))
        self.cache[key] = (amap, ident)
        return self.cache[key]


def _text(x):
    return x if isinstance(x, str) else ''


# 重疊區用哪一端的鹼基：'r2'（預設）= R2；'r1' = 原版作法（R1 到V終點，R2 接在後面）。
# R2 經過 Q20 3' trimming，R1 沒有；重疊區又剛好落在兩條read各自的3'尾端，所以用R2的鹼基
# 錯誤較少（模擬中重疊區裡的錯誤幾乎都來自R1尾端）。
OVERLAP_FROM = os.environ.get("BCR_STITCH_OVERLAP_FROM", "r2")
# 1 = 「germline座標上有缺口」但兩條read實際上仍有重疊（R1 的V比對常在read結束前就停了）時，
# 用整條R1與R2的無缺口比對驗證後救回。預設關閉（行為與原版相同：有缺口一律排除）。
RESCUE_GAP = os.environ.get("BCR_STITCH_RESCUE_GAP", "0") == "1"
SHIFT_SEARCH = 12
PHYSICAL_THRESHOLD = 0.85


def _ungapped_identity(s1, s2, q1, d):
    """R2 第 j 個鹼基對到 R1 第 j+d 個（1-based）時，R1[1..q1] 與 R2 重疊部分的 identity；
    重疊少於 QC_MIN_OVERLAP 時回傳 None。"""
    lo = max(1, d + 1)
    hi = min(q1, len(s2) + d)
    n = hi - lo + 1
    if n < QC_MIN_OVERLAP:
        return None
    same = sum(1 for i in range(lo, hi + 1) if s1[i - 1] == s2[i - d - 1])
    return same / n


def _best_register(s1, s2, q1, d0, threshold):
    """在 d0±SHIFT_SEARCH 內找 R1[1..q1] 與 R2 無缺口比對 identity 最高的對位。
    回傳 (best_d, best_identity, identity_at_d0)；重疊都太短時 best_d 為 None。"""
    cands = {}
    for d in range(d0 - SHIFT_SEARCH, d0 + SHIFT_SEARCH + 1):
        ident = _ungapped_identity(s1, s2, q1, d)
        if ident is not None:
            cands[d] = ident
    if not cands:
        return None, None, None
    best = max(cands, key=lambda d: (cands[d], -abs(d - d0)))
    if cands[best] < threshold:
        return None, cands[best], cands.get(d0)
    return best, cands[best], cands.get(d0)


def _join_register(s1, s2, d, r1_cut, r2_vstart):
    """固定對位 d（R1 位置 = R2 位置 + d）下的拼接。r1_cut = 依原版作法 R1 保留到的位置。"""
    if OVERLAP_FROM == 'r2':
        i = r2_vstart + d - 1          # R2 的V起點之前，R1 還要保留幾個鹼基
        if 0 <= i <= len(s1):
            return s1[:i] + s2[r2_vstart - 1:]
    return s1[:r1_cut] + s2[r1_cut - d:]


def stitch_pair(r1, r2, mapper):
    """r1/r2: dict 含 sequence, v_call, v_cigar, 四個座標(int)。
    回傳 (stitched_seq or None, reason, overlap, different_allele)；reason 為 'ok'（重疊區通過QC）、
    'ok_short'（相鄰或可比對鹼基 < QC_MIN_OVERLAP，未做QC）、'ok_shifted'（整段重疊比對修正了
    germline換算的對位）、'ok_rescued'（BCR_STITCH_RESCUE_GAP=1 時救回的缺口read）、
    'gap'、'qc'、'unmappable'。"""
    a = _text(r1['v_call']).split(',')[0]
    b = _text(r2['v_call']).split(',')[0]
    diff_allele = a != b
    m1 = germline_to_query(_text(r1['v_cigar']), r1['v_sequence_start'], r1['v_sequence_end'],
                           r1['v_germline_start'], r1['v_germline_end'])
    m2 = germline_to_query(_text(r2['v_cigar']), r2['v_sequence_start'], r2['v_sequence_end'],
                           r2['v_germline_start'], r2['v_germline_end'])
    if m1 is None or m2 is None:
        return None, 'unmappable', 0, diff_allele
    amap, ident = mapper.mapping(a, b)
    if ident < MIN_ALLELE_IDENTITY:
        return None, 'unmappable', 0, diff_allele

    def a2b(p):
        return p if amap is None else amap.get(p)

    a_end = r1['v_germline_end']
    # R2 起點換算到 alleleA 座標：找第一個能對到A的B位置
    b2a = {v: k for k, v in amap.items()} if amap is not None else None
    a_start2 = None
    for bp in range(r2['v_germline_start'], r2['v_germline_end'] + 1):
        ap = bp if b2a is None else b2a.get(bp)
        if ap is not None:
            a_start2 = ap
            break
    if a_start2 is None:
        return None, 'unmappable', 0, diff_allele
    overlap = a_end - a_start2 + 1
    s1, s2 = r1['sequence'], r2['sequence']
    v1_end, v2_start = r1['v_sequence_end'], r2['v_sequence_start']

    too_short = overlap < 0 or overlap < MIN_OVERLAP
    if too_short or overlap == 0:
        # 預測：R2 的V起點落在 R1 的第 v1_end + (a_start2 - a_end) 個位置
        d0 = v1_end + (a_start2 - a_end) - v2_start
        if too_short and not RESCUE_GAP:
            return None, 'gap', overlap, diff_allele
        # 相鄰（或開啟救回時的缺口）：用整條R1（含V比對終點之後的鹼基）與R2的實際重疊核對對位
        d, best_id, id0 = _best_register(s1, s2, len(s1), d0, PHYSICAL_THRESHOLD)
        if d is not None:
            if overlap < 0 and v2_start + d - 1 > v1_end:
                reason = 'ok_rescued'
            else:
                reason = 'ok' if d == d0 else 'ok_shifted'
            return s1[:v2_start + d - 1] + s2[v2_start - 1:], reason, overlap, diff_allele
        if too_short:
            return None, 'gap', overlap, diff_allele
        return s1[:v1_end] + s2[v2_start - 1:], 'ok_short', overlap, diff_allele

    # 在預測的對位上逐鹼基比對重疊區；記下「兩端都有鹼基」的最後一個A位置
    matches = compared = 0
    last = None
    for p in range(a_start2, a_end + 1):
        q1 = m1.get(p)
        bp = a2b(p)
        q2 = m2.get(bp) if bp is not None else None
        if q1 is None or q2 is None:
            continue
        compared += 1
        matches += s1[q1 - 1] == s2[q2 - 1]
        last = (q1, q2)
    if last is None:
        return None, 'unmappable', overlap, diff_allele
    q1, q2 = last
    d0 = q1 - q2
    if compared < QC_MIN_OVERLAP:
        # V比對區的重疊太短：改用整條R1（含V終點之後被soft-clip的鹼基）核對實際重疊
        d, best_id, id0 = _best_register(s1, s2, len(s1), d0, PHYSICAL_THRESHOLD)
        if d is not None and (d == d0 or id0 is None or best_id > id0 + 0.10):
            return _join_register(s1, s2, d, q1, v2_start), ('ok' if d == d0 else 'ok_shifted'), overlap, diff_allele
        return s1[:q1] + s2[q2:], 'ok_short', overlap, diff_allele
    # 兩條germline比對時，缺口若落在密碼子重複區（例如CDR2的AGTAGTAGT）位置本身就有歧義，
    # 換算出的對位可能錯開3的倍數：重疊區夠長時，在預測對位附近挑整段重疊區無缺口 identity 最高者。
    best_d, best_id, id0 = _best_register(s1, s2, q1, d0, QC_THRESHOLD)
    cigar_ok = matches / compared >= QC_THRESHOLD
    if best_d is not None and best_d != d0 and (id0 is None or best_id > id0 + 0.10 or not cigar_ok):
        return _join_register(s1, s2, best_d, q1, v2_start), 'ok_shifted', overlap, diff_allele
    if not cigar_ok:
        return None, 'qc', overlap, diff_allele
    if id0 is not None and id0 >= QC_THRESHOLD:
        # 整段重疊區在對位 d0 上都吻合：沿這條已驗證的對位拼接（OVERLAP_FROM=r2 時重疊區用R2）
        return _join_register(s1, s2, d0, q1, v2_start), 'ok', overlap, diff_allele
    # 重疊區內有read本身的插入/缺失（無缺口比對不成立、CIGAR逐位比對成立）：在最後一個共同位置切
    return s1[:q1] + s2[q2:], 'ok', overlap, diff_allele


def productive_mask(df):
    """逐列判斷：productive 有明確 T/F 就用它；空白時才用 stop_codon=F 且 vj_in_frame=T。
    （原版是看整個chunk裡「有沒有任何一列」有T/F來切換規則，同一列的結果會隨chunk內容改變。）"""
    def col(name):
        return df[name].astype(str).str.strip() if name in df.columns else pd.Series([''] * len(df), index=df.index)
    prod, sc, vjf = col('productive'), col('stop_codon'), col('vj_in_frame')
    explicit = prod.isin(['T', 'F'])
    return (explicit & (prod == 'T')) | (~explicit & (sc == 'F') & (vjf == 'T'))


def _to_int_rows(chunk):
    for c in COORD_COLS:
        chunk[c] = pd.to_numeric(chunk[c], errors='coerce')
    chunk = chunk.dropna(subset=COORD_COLS)
    for c in COORD_COLS:
        chunk[c] = chunk[c].astype(int)
    return chunk


def _replace_stats_row(stats_path, idx_str, line):
    rows = []
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            rows = [r for r in f.read().splitlines()[1:] if r and not r.startswith(idx_str + ',')]
    rows.append(line)
    rows.sort()
    tmp = stats_path + '.tmp'
    with open(tmp, 'w') as f:
        f.write(STATS_HEADER + '\n' + '\n'.join(rows) + '\n')
    os.replace(tmp, stats_path)


def run(label, r1_path, r2_path, out_dir, germline_v_fasta, next_hint):
    chunk_size = int(os.environ.get("BCR_CHUNK_SIZE", "500000"))
    chunk_dir = f"{out_dir}/chunks"
    stats_path = f"{out_dir}/chunk_stats.csv"
    rej_dir = f"{out_dir}/stitch_rejected"
    os.makedirs(chunk_dir, exist_ok=True)
    os.makedirs(rej_dir, exist_ok=True)
    for p in (r1_path, r2_path, germline_v_fasta):
        if not os.path.exists(p):
            sys.exit(f"[ERROR] 找不到 {p}")
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            if f.readline().strip() != STATS_HEADER:
                sys.exit(f"[ERROR] {stats_path} 是舊版欄位格式；請先移走舊的 chunk_stats.csv 與 chunks/*.stitch.done 再重跑")
    mapper = AlleleMapper(germline_v_fasta)
    complete_mark = f"{out_dir}/stitch.complete"
    if os.path.exists(complete_mark):
        os.remove(complete_mark)

    print(f"[{label}] Step 1: 建立R1查找表 <- {os.path.basename(r1_path)} ...")
    r1_lookup = {}
    n_r1_total = 0
    for chunk in pd.read_csv(r1_path, sep='\t', usecols=R1_COLS, dtype=str, chunksize=chunk_size,
                             keep_default_na=False, na_values=[''], engine='c'):
        n_r1_total += len(chunk)
        chunk = _to_int_rows(chunk)
        for row in chunk.itertuples(index=False):
            # v_call 只有幾十種，intern 後共用同一個字串物件，降低查找表的記憶體用量
            r1_lookup[row.sequence_id] = (row.sequence, sys.intern(_text(row.v_call)), row.v_cigar,
                                          int(row.v_sequence_start), int(row.v_sequence_end),
                                          int(row.v_germline_start), int(row.v_germline_end))
    print(f"  R1總筆數={n_r1_total:,}，有效V基因比對且建入查找表={len(r1_lookup):,}")

    print(f"[{label}] Step 2: 分塊處理R2-only <- {os.path.basename(r2_path)} ...")
    chunk_idx = 0
    n_r2_total = 0
    for chunk in pd.read_csv(r2_path, sep='\t', usecols=R2_COLS, dtype=str, chunksize=chunk_size,
                             keep_default_na=False, na_values=[''], engine='c'):
        idx_str = f"{chunk_idx:04d}"
        n_r2_total += len(chunk)
        done_mark = f"{chunk_dir}/chunk_{idx_str}.stitch.done"
        out_fasta = f"{chunk_dir}/chunk_{idx_str}.fasta"
        if os.path.exists(done_mark):
            print(f"  chunk {idx_str} 已完成，跳過")
            chunk_idx += 1
            continue

        n_r2_rows = len(chunk)
        prod = chunk[productive_mask(chunk)].copy()
        n_productive = len(prod)
        prod = _to_int_rows(prod)
        c = dict(paired=0, nogap=0, qc_pass=0, short=0, shifted=0, rescued=0, diff=0, gap=0, qc=0, unmappable=0, error=0)
        out_lines, rejected = [], []
        for row in prod.itertuples(index=False):
            r1v = r1_lookup.get(row.sequence_id)
            if r1v is None:
                continue
            c['paired'] += 1
            r1 = dict(zip(('sequence', 'v_call', 'v_cigar', 'v_sequence_start', 'v_sequence_end',
                           'v_germline_start', 'v_germline_end'), r1v))
            r2 = row._asdict()
            try:
                stitched, reason, overlap, diff = stitch_pair(r1, r2, mapper)
            except Exception as e:  # 出錯就排除（fail closed），並記錄原因
                stitched, reason, overlap, diff = None, 'error', 0, False
                rejected.append(f"{row.sequence_id}\terror:{type(e).__name__}\n")
            c['diff'] += diff
            if reason != 'gap' and reason != 'unmappable':
                c['nogap'] += 1
            if stitched is None:
                c[reason] += 1
                if reason != 'error':
                    rejected.append(f"{row.sequence_id}\t{reason}\t{overlap}\n")
                continue
            if reason == 'ok_short':
                c['short'] += 1
            if reason == 'ok_shifted':
                c['shifted'] += 1
            if reason == 'ok_rescued':
                c['rescued'] += 1
            c['qc_pass'] += 1
            out_lines.append(f">{row.sequence_id}\n{stitched}\n")

        tmp = out_fasta + '.tmp'
        with open(tmp, 'w') as f:
            f.writelines(out_lines)
        os.replace(tmp, out_fasta)
        with gzip.open(f"{rej_dir}/chunk_{idx_str}.tsv.gz", 'wt') as f:
            f.write("sequence_id\treason\toverlap\n")
            f.writelines(rejected)
        _replace_stats_row(stats_path, idx_str,
                           f"{idx_str},{n_r2_rows},{n_productive},{c['paired']},{c['nogap']},{c['qc_pass']},"
                           f"{c['short']},{c['shifted']},{c['rescued']},{c['diff']},{c['gap']},{c['qc']},{c['unmappable']},{c['error']}")
        with open(done_mark, 'w') as f:
            f.write('done')
        print(f"  chunk {idx_str}: R2 rows={n_r2_rows:,} productive={n_productive:,} 配對上R1={c['paired']:,} "
              f"無缺口={c['nogap']:,} QC通過={c['qc_pass']:,}（其中短重疊未QC={c['short']:,}；對位修正={c['shifted']:,}；缺口救回={c['rescued']:,}；"
              f"R1/R2不同allele={c['diff']:,}；排除 gap={c['gap']:,} qc={c['qc']:,} "
              f"無法對應={c['unmappable']:,} 錯誤={c['error']:,}）")
        chunk_idx += 1

    with open(complete_mark, 'w') as f:
        f.write(f"chunks={chunk_idx} r2_rows_total={n_r2_total}\n")
    print(f"[{label}] 拼接+QC完成，chunk FASTA已輸出至 {chunk_dir}")
    print(f"接下來執行：{next_hint}")
