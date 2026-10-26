#!/usr/bin/env python3
"""
04c_stitch_qc_igkl.py
接續04a(demux)、04b(R1-only/R2-only igblast)：座標拼接+QC，邏輯與02g共用
stitch_core.py（2026-10起不再各自複製一份），路徑改成IGK/IGL。

用法：python 04c_stitch_qc_igkl.py <rabbit_number> <IGK|IGL>
中斷後重跑：已完成的chunk會被跳過（chunk_XXX.stitch.done標記）。
"""

import os
import sys

from stitch_core import run

RABBIT = sys.argv[1] if len(sys.argv) > 1 else None
CHAIN = sys.argv[2] if len(sys.argv) > 2 else None
if not RABBIT or CHAIN not in ("IGK", "IGL"):
    print("用法：python 04c_stitch_qc_igkl.py <rabbit_number> <IGK|IGL>")
    sys.exit(1)

BASE_DIR = os.environ.get("BCR_PIPELINE_ROOT", "/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline")
R1_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/{CHAIN}_R1only/15651-LT-{RABBIT}_AVITI_{CHAIN}_R1only.airr.tsv"
R2_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/{CHAIN}_R2only/15651-LT-{RABBIT}_AVITI_{CHAIN}_R2only.airr.tsv"
OUT_DIR = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/{CHAIN}_Stitched"
GERMLINE_V = f"{BASE_DIR}/reference/IMGT_reference_clean/{CHAIN}V.fasta"

run(f"{CHAIN} Rabbit {RABBIT}", R1_PATH, R2_PATH, OUT_DIR, GERMLINE_V, f"bash 04d_stitch_igblast_igkl.sh {RABBIT} {CHAIN}")
