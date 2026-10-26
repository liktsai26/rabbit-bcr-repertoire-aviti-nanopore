#!/bin/bash
# IGL全量座標拼接+QC pipeline正式執行入口。
# 3隻兔子，每隻依序：demux(04a) -> R1-only igblast(04b) -> R2-only igblast(04b)
# -> 座標拼接+QC(04c) -> 第3次igblast(04d)。
#
# 已知背景（見記憶bcr_pipeline_igl_r2only_retry_slowdown.md）：IGL R2-only/
# 拼接後第3次igblast會頻繁觸發igblastn 1.22.0本身的NULL pointer bug（已確認
# 是igblastn的問題，不是我們資料或pipeline的問題，查不到不犧牲準確度的修正
# 方式），round-based capped retry會自動處理、但速度慢很多，全量預估6-9天。
#
# 2026-10-02新增：跟IGK當初做法不同，這次**每隻兔子跑完Step4後，確認3份
# 最終合併檔存在且非空，就立即刪除該兔子的chunks/資料夾**，不等三隻兔子
# 都跑完才清——因為IGL全程耗時遠比IGK長，硬碟空間需要全程控制在「只疊加
# 當下這一隻」的量，不要累積到三隻份。
#
# 中斷後重跑：各步驟自身都有chunk/子批次層級的.done標記，會自動跳過已完成
# 部分；chunks 已被清掉的兔子由 rabbit_complete() 整隻跳過（2026-10 修正：原本只有
# 04b/04d 會跳過，04a/04c 會把已完成的兔子重跑一遍）。清理前要求三份最終檔都有
# .complete 完成標記；不可處理的read ID 與 chunk_stats.csv 都在 chunks/ 之外，不會被清掉。
# 設 BCR_KEEP_CHUNKS=1 可保留 chunks/ 不清理。
set -e
SCRIPT_DIR="${BCR_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
RESULTS_DIR="${BASE_DIR}/results"
PYTHON="${BCR_PIPELINE_PYTHON:-${BCR_ENV_PREFIX:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env}/bin/python}"
CHAIN="IGL"

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

cleanup_rabbit_chunks() {
    local rabbit="$1"
    local base="${RESULTS_DIR}/15651-LT-${rabbit}_Aviti24"
    local ok=1

    for side_dir in "${CHAIN}_R1only" "${CHAIN}_R2only" "${CHAIN}_Stitched"; do
        local final_airr
        final_airr="${base}/${side_dir}/15651-LT-${rabbit}_AVITI_${side_dir}.airr.tsv"
        # 2026-10：只看「存在且非空」無法分辨中斷時寫到一半的檔案，改成要求完成標記
        if [ ! -s "$final_airr" ] || [ ! -f "${final_airr}.complete" ]; then
            echo "  !! 警告：${side_dir} 的最終合併檔不存在、是空檔或沒有完成標記，不會清理 chunks/（${base}/${side_dir}）"
            ok=0
        fi
    done

    if [ "$ok" -eq 1 ]; then
        for side_dir in "${CHAIN}_R1only" "${CHAIN}_R2only" "${CHAIN}_Stitched"; do
            local chunks_dir="${base}/${side_dir}/chunks"
            if [ -d "$chunks_dir" ]; then
                local size
                size=$(du -sh "$chunks_dir" 2>/dev/null | cut -f1)
                rm -rf "$chunks_dir"
                echo "  已清理: ${side_dir}/chunks (${size})"
            fi
        done
    else
        echo "  !! Rabbit ${rabbit} 的chunks/保留不清，請人工確認後再手動處理"
    fi
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

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — 驗證最終檔案並清理chunks $(date '+%Y-%m-%d %H:%M:%S')"
    if [ "${BCR_KEEP_CHUNKS:-0}" = "1" ]; then
        echo "  BCR_KEEP_CHUNKS=1：保留 chunks/ 不清理"
    else
        cleanup_rabbit_chunks "$RABBIT"
    fi

    echo "### Rabbit ${RABBIT} ${CHAIN} 全部完成 $(date '+%Y-%m-%d %H:%M:%S')"
done

echo ""
echo "IGL_FULL_PIPELINE_ALL_DONE"
