#!/bin/bash
# IGK全量座標拼接+QC pipeline正式執行入口（只跑IGK，IGL暫停）。
# 3隻兔子，每隻依序：demux(04a) -> R1-only igblast(04b) -> R2-only igblast(04b)
# -> 座標拼接+QC(04c) -> 第3次igblast(04d)。
# 中斷後重跑：各步驟自身都有chunk/子批次層級的.done標記，會自動跳過已完成部分。
set -e
SCRIPT_DIR="${BCR_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
RESULTS_DIR="${BASE_DIR}/results"
PYTHON="${BCR_PIPELINE_PYTHON:-${BCR_ENV_PREFIX:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env}/bin/python}"
CHAIN="IGK"

# 2026-10 code review：某隻兔子三份最終結果都有完成標記（.complete）就整隻跳過。原本只有
# 04b/04d 會跳過，04a/04c 沒有保護，IGL 清掉 chunks/ 後重跑會把已完成的兔子重新 demux、
# 重新拼接，chunk_stats.csv 也會重複追加。
rabbit_complete() {
    local base="${RESULTS_DIR}/15651-LT-$1_Aviti24" side
    for side in R1only R2only Stitched; do
        [ -f "${base}/${CHAIN}_${side}/15651-LT-$1_AVITI_${CHAIN}_${side}.airr.tsv.complete" ] || return 1
    done
    return 0
}

for RABBIT in 1 2 3; do
    if rabbit_complete "$RABBIT"; then
        echo "### Rabbit ${RABBIT} ${CHAIN} 三份最終結果都已完成（.complete），整隻跳過"
        continue
    fi
    echo ""
    echo "################################################################"
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step1 demux 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    echo "################################################################"
    bash "${SCRIPT_DIR}/04a_run_aviti_igkl_demux.sh" "$RABBIT" "$CHAIN"

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step2 R1-only igblast 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    bash "${SCRIPT_DIR}/04b_run_aviti_igkl_igblast.sh" "$RABBIT" "$CHAIN" r1only

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step2 R2-only igblast 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    bash "${SCRIPT_DIR}/04b_run_aviti_igkl_igblast.sh" "$RABBIT" "$CHAIN" r2only

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step3 座標拼接+QC 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    "$PYTHON" "${SCRIPT_DIR}/04c_stitch_qc_igkl.py" "$RABBIT" "$CHAIN"

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step4 拼接後第3次igblast 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    bash "${SCRIPT_DIR}/04d_stitch_igblast_igkl.sh" "$RABBIT" "$CHAIN"

    echo "### Rabbit ${RABBIT} ${CHAIN} 全部完成 $(date '+%Y-%m-%d %H:%M:%S')"
done

echo ""
echo "IGK_FULL_PIPELINE_ALL_DONE"
