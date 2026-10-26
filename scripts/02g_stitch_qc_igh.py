#!/usr/bin/env python3
"""
02g_stitch_qc_igh.py
接續02d(R2-only igblast)、02f(R1-only igblast)：用read ID把同一分子的R1、R2
配對起來，以igblast各自跑出來的座標為依據，把兩段序列接成一條完整V(D)J序列
（座標拼接法，Coordinate-guided Stitching），並驗證重疊區品質，只保留品質夠好的
拼接結果（2026-10 起不再用 Biopython 局部比對，見下方第5點）。

背景：R2-only單獨判定V基因，對IGHV1S40/1S45這類差異點集中在FR1/CDR1（R2讀
不到的區域）的近乎同源基因對，會系統性判錯（2026-09-24~25驗證，R2-only
判定比例跟Nanopore/Merge差好幾倍）。座標拼接法把R1的FR1資訊也接回來，
2026-09-25用LT-1/2/3各10萬條抽樣驗證過，結果跟Nanopore/Merge高度吻合。

流程：
1. 把該rabbit的R1-only AIRR TSV整個讀進記憶體，建成 sequence_id -> R1資訊 的
   查找表（R1資料量約400-500萬筆，8GB RAM機器上可負擔，比照paired_igh_cache.py
   的作法）
2. 分塊讀取R2-only AIRR TSV（500K rows/chunk），逐筆用sequence_id查R1
3. 排除：R1沒有資料、R2非productive、有缺口（R1終點跟R2起點之間有空隙）
4. 拼接：R1保留到v_sequence_end，R2從對應germline位置接續
   （sequence欄位在rev_comp=T時igblast已經轉成germline方向，不需要再手動
   revcomp——這是2026-09-25驗證時發現並修正過的重點，見討論記錄）
   2026-10 修正：R1/R2 被判給不同V allele時先換算兩條germline的座標，並用
   v_cigar 處理read內的插入/缺失（細節見 stitch_core.py）
5. QC：在預測對位上逐鹼基比對重疊區，identity<90%視為不可靠、排除
6. 輸出：拼接+QC通過的序列，寫成chunked FASTA供02h跑igblast

用法：python 02g_stitch_qc_igh.py <rabbit_number>
中斷後重跑：已完成的chunk會被跳過（chunk_XXX.stitch.done標記）。
"""

import os
import sys

from stitch_core import run

RABBIT = sys.argv[1] if len(sys.argv) > 1 else None
if not RABBIT:
    print("用法：python 02g_stitch_qc_igh.py <rabbit_number>")
    sys.exit(1)

BASE_DIR = os.environ.get("BCR_PIPELINE_ROOT", "/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline")
R1_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/IGH_R1only/15651-LT-{RABBIT}_AVITI_IGH_R1only.airr.tsv"
R2_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/IGH_R2only/15651-LT-{RABBIT}_AVITI_IGH_R2only.airr.tsv"
OUT_DIR = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/IGH_Stitched"
GERMLINE_V = f"{BASE_DIR}/reference/IMGT_reference_clean/IGHV.fasta"

run(f"Rabbit {RABBIT}", R1_PATH, R2_PATH, OUT_DIR, GERMLINE_V, f"bash 02h_stitch_igblast.sh {RABBIT}")
