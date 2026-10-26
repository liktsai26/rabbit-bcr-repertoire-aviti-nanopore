#!/usr/bin/env python3
"""
分輪次+深度上限的igblast執行器（2026-09-28，因應IGL測試遇到深層遞迴重跑問題設計）。

跟原本bash版igblast_verified()的差別：
- 原本：驗證失敗就對半切，遞迴下去，同一條路徑可以一路切到剩1條，每一層都要重新
  啟動一次igblastn（固定成本疊加，這是這次IGL測試變慢、process意外終止的主因）
- 現在：每一輪(round)裡，對半切最多到max_depth層就停止，把那一輪切到深度上限還是
  驗證不過的reads（用sequence_id比對，不是只看行數），收集起來當作下一輪的輸入，
  重新從正常大小的batch開始，而不是繼續往下切同一條路徑。
  這樣「頑固難處理」的reads會被集中到最後幾輪，數量越來越少，不會拖累整批。

用法：python igblast_round_retry.py <fasta> <out_airr.tsv> <germline_V> <germline_D> <germline_J>

2026-10 code review 修正：輸出列維持輸入順序（原本依 set 走訪順序，每次執行順序不同）；
欄位數跟表頭不符的殘缺列視為失敗重跑；igblastn 完全沒有輸出時直接報錯並附上 stderr
（原本會在最後 header 為 None 時丟 TypeError）；空輸入輸出空檔；不可處理ID清單的路徑
不再依賴 str.replace，檔尾補換行；igblastn 路徑與 stderr 紀錄可由 BCR_* 環境變數設定。
"""
import sys
import os
import subprocess
import tempfile
import shutil

FASTA, OUT, GV, GD, GJ = sys.argv[1:6]
MAX_DEPTH = 5
BATCH_SIZE = int(os.environ.get("BCR_BATCH_SIZE", "25000"))

ENV_PREFIX = os.environ.get("BCR_ENV_PREFIX", "/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env")
IGDATA = os.environ.get("BCR_IGDATA", f"{ENV_PREFIX}/share/igblast")
IGBLASTN = os.environ.get("BCR_IGBLASTN", f"{ENV_PREFIX}/bin/igblastn")
STDERR_LOG = os.environ.get("BCR_IGBLAST_STDERR_LOG")
AUX = f"{IGDATA}/optional_file/rabbit_gl.aux"
BLASTDB = f"{IGDATA}/database"

IGBLAST_ARGS = ["-organism", "rabbit", "-ig_seqtype", "Ig", "-outfmt", "19",
                "-num_threads", "2", "-extend_align5end", "-extend_align3end",
                "-auxiliary_data", AUX,
                # 資料庫名稱外加引號：路徑含空白時 BLAST 會把空白當成多個資料庫的分隔符號
                "-germline_db_V", f'"{GV}"', "-germline_db_D", f'"{GD}"', "-germline_db_J", f'"{GJ}"']


def read_fasta(path):
    seqs = {}
    cur_id, cur_seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(cur_seq)
    return seqs


def write_fasta(seqs, path):
    with open(path, "w") as f:
        for sid, seq in seqs.items():
            f.write(f">{sid}\n{seq}\n")


def run_igblastn(fasta_path, out_path):
    env = os.environ.copy()
    env["IGDATA"] = IGDATA
    env["BLASTDB"] = BLASTDB
    p = subprocess.run([IGBLASTN] + IGBLAST_ARGS + ["-query", fasta_path, "-out", out_path],
                       env=env, stderr=subprocess.PIPE, text=True)
    if p.stderr and STDERR_LOG:
        with open(STDERR_LOG, "a") as f:
            f.write(p.stderr)
    return p.returncode, p.stderr


def parse_airr(path):
    """回傳 (header_line, {sequence_id: full_row_line})"""
    rows = {}
    header = None
    if not os.path.exists(path):
        return None, rows
    with open(path) as f:
        header = f.readline().rstrip("\n") or None
        ncol = header.count("\t") if header else -1
        for line in f:
            if not line.endswith("\n") or line.count("\t") != ncol:
                continue  # 殘缺列（例如程序被中止時寫到一半）視為沒有輸出
            sid = line.split("\t", 1)[0]
            rows[sid] = line.rstrip("\n")
    return header, rows


def process_batch_capped(seqs, depth, tmpdir):
    """回傳 (成功的AIRR row list, 這批裡還沒成功的sequence_id set)。
    深度到上限或只剩1條就停止繼續切，把剩下驗證不過的丟回給呼叫者收集到下一輪。"""
    tmp_fasta = os.path.join(tmpdir, f"batch_{depth}_{len(seqs)}.fasta")
    tmp_out = os.path.join(tmpdir, f"batch_{depth}_{len(seqs)}.airr.tsv")
    write_fasta(seqs, tmp_fasta)
    rc, err = run_igblastn(tmp_fasta, tmp_out)
    header, rows = parse_airr(tmp_out)
    if header is None and depth == 0 and len(seqs) > 1 and rc != 0:
        sys.exit(f"[ERROR] igblastn 執行失敗且沒有任何輸出（exit {rc}），通常是設定問題：\n{err[-2000:]}")
    in_ids = set(seqs.keys())
    out_ids = set(rows.keys())
    missing = in_ids - out_ids

    os.remove(tmp_fasta)
    if os.path.exists(tmp_out):
        os.remove(tmp_out)

    if not missing:
        return header, [rows[i] for i in seqs], set()

    if depth >= MAX_DEPTH or len(seqs) <= 1:
        # 深度到上限（或只剩1條），不再往下切，成功的收下、失敗的丟回上層收集
        return header, [rows[i] for i in seqs if i in out_ids], missing

    # 對半切，繼續在這一輪內遞迴（但不超過max_depth）
    ids = list(seqs.keys())
    half = (len(ids) + 1) // 2
    part1 = {i: seqs[i] for i in ids[:half]}
    part2 = {i: seqs[i] for i in ids[half:]}
    h1, rows1, fail1 = process_batch_capped(part1, depth + 1, tmpdir)
    h2, rows2, fail2 = process_batch_capped(part2, depth + 1, tmpdir)
    header = header or h1 or h2
    return header, rows1 + rows2, fail1 | fail2


def main():
    all_seqs = read_fasta(FASTA)
    print(f"總序列數: {len(all_seqs):,}")
    if not all_seqs:
        open(OUT, "w").close()
        print("輸入沒有序列，輸出空檔")
        return
    pending = all_seqs
    round_num = 0
    header_final = None
    all_rows = []
    irrecoverable = set()

    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(FASTA))) as tmpdir:
        while pending:
            round_num += 1
            print(f"--- Round {round_num}: {len(pending):,} 條 ---")
            ids = list(pending.keys())
            leftover = {}
            n_batches = (len(ids) + BATCH_SIZE - 1) // BATCH_SIZE
            for bi in range(n_batches):
                batch_ids = ids[bi * BATCH_SIZE: (bi + 1) * BATCH_SIZE]
                batch = {i: pending[i] for i in batch_ids}
                header, rows, failed = process_batch_capped(batch, 0, tmpdir)
                if header:
                    header_final = header
                all_rows.extend(rows)
                for fid in failed:
                    leftover[fid] = pending[fid]
                print(f"  batch {bi}: {len(batch):,} 條，成功 {len(batch) - len(failed):,}，"
                      f"留到下一輪 {len(failed):,}")

            if len(leftover) == len(pending):
                print(f"Round {round_num} 完全沒有進展，剩下 {len(leftover):,} 條判定為不可處理，記錄後停止")
                irrecoverable = set(leftover.keys())
                break
            pending = leftover

    order = {sid: k for k, sid in enumerate(all_seqs)}
    all_rows.sort(key=lambda r: order.get(r.split("\t", 1)[0], len(order)))
    with open(OUT, "w") as f:
        if header_final is not None:
            f.write(header_final + "\n")
        for r in all_rows:
            f.write(r + "\n")

    print(f"完成：成功 {len(all_rows):,} / {len(all_seqs):,}，"
          f"不可處理 {len(irrecoverable):,} 條")
    if irrecoverable:
        base = OUT[:-len(".airr.tsv")] if OUT.endswith(".airr.tsv") else OUT
        irr_path = base + "_irrecoverable_ids.txt"
        with open(irr_path, "w") as f:
            f.write("".join(i + "\n" for i in sorted(irrecoverable)))
        print(f"不可處理的read ID清單: {irr_path}")


if __name__ == "__main__":
    main()
