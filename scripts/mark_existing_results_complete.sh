#!/bin/bash
# ============================================================
# mark_existing_results_complete.sh
# 一次性工具：給「舊版腳本已經跑完」的結果補上新版腳本使用的完成標記，
# 讓新版腳本把它們當成已完成而不重跑（新版會拒絕信任沒有標記的舊檔案）。
#
# 只在你確認這些舊結果完整時使用。標記內容會寫明「legacy, not verified」。
# 不會修改任何結果檔，只新增空的標記檔：
#   <stage>/<sample>.airr.tsv.complete、<R1only|R2only>/demux.complete、<Stitched>/stitch.complete
#
# 注意：若要用新版拼接重做 *_Stitched，請先把舊的 *_Stitched 資料夾移走再跑 02g/04c，
# 不要對它們使用本工具（或加 --skip-stitched）。
#
# 用法：bash mark_existing_results_complete.sh [--skip-stitched] [--dry-run]
#   結果目錄：${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}/results
# ============================================================
set -e
SKIP_STITCHED=0; DRY=0
for a in "$@"; do
    case "$a" in
        --skip-stitched) SKIP_STITCHED=1 ;;
        --dry-run) DRY=1 ;;
        *) echo "用法：bash mark_existing_results_complete.sh [--skip-stitched] [--dry-run]"; exit 1 ;;
    esac
done
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
RESULTS="${BASE_DIR}/results"
[ -d "$RESULTS" ] || { echo "找不到 ${RESULTS}"; exit 1; }
STAMP="legacy, not verified ($(date '+%Y-%m-%d'))"

mark() {
    if [ -e "$1" ]; then return 0; fi
    echo "  + $1"
    [ "$DRY" = "1" ] || echo "$STAMP" > "$1"
}

shopt -s nullglob
for final in "${RESULTS}"/15651-LT-*_Aviti24/*/15651-LT-*_AVITI_*.airr.tsv; do
    dir="$(dirname "$final")"
    case "$dir" in
        *_Stitched) [ "$SKIP_STITCHED" = "1" ] && continue; mark "${dir}/stitch.complete" ;;
        *_R1only|*_R2only) mark "${dir}/demux.complete" ;;
        *) continue ;;
    esac
    if [ ! -s "$final" ]; then
        echo "  ! 跳過空檔：$final"
        continue
    fi
    mark "${final}.complete"
done
[ "$DRY" = "1" ] && echo "（--dry-run：沒有寫入任何檔案）"
echo "完成。"
