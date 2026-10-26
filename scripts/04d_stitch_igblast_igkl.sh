#!/bin/bash
# ============================================================
# 04d_stitch_igblast_igkl.sh
# 接續04c_stitch_qc_igkl.py：對拼接+QC通過的chunk fasta跑第3次igblast，
# 取得完整V(D)J判讀。架構跟04b完全相同（chunk+子批次+checkpoint+
# round-based capped retry），只是輸入/輸出目錄換成{CHAIN}_Stitched。
#
# 用法：bash 04d_stitch_igblast_igkl.sh <rabbit> <IGK|IGL>
# 中斷後重跑：已完成的chunk/子批次會被跳過（.done標記）
# 2026-10 code review 修正：分塊/驗證/合併邏輯改由 igblast_chunk_lib.sh 共用；無法處理的
# read ID 現在保存在 <OUT_DIR>/irrecoverable/（原本寫在子批次資料夾，隨即被 rm -rf 刪掉）；
# igblast_round_retry.py 改用本腳本所在資料夾的那一份；路徑可用 BCR_* 環境變數覆寫。
# ============================================================
set -e
RABBIT="$1"; CHAIN="$2"
[ -z "$RABBIT" ] || [ -z "$CHAIN" ] && {
    echo "用法：bash 04d_stitch_igblast_igkl.sh <rabbit> <IGK|IGL>"; exit 1; }
case "$CHAIN" in IGK|IGL) ;; *) echo "CHAIN必須是IGK或IGL"; exit 1 ;; esac
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
ENV_PREFIX="${BCR_ENV_PREFIX:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env}"
IMGT="${BASE_DIR}/reference/IMGT_reference_clean"
PYTHON="${BCR_PIPELINE_PYTHON:-${ENV_PREFIX}/bin/python}"
ROUND_RETRY="${SCRIPT_DIR}/igblast_round_retry.py"
IGBLAST_BATCH="${BCR_IGBLAST_BATCH:-25000}"
IGBLAST_MODE=round

OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/${CHAIN}_Stitched"
CHUNK_DIR="${OUT_DIR}/chunks"
SAMPLE="15651-LT-${RABBIT}_AVITI_${CHAIN}_Stitched"
FINAL_AIRR="${OUT_DIR}/${SAMPLE}.airr.tsv"
STATS="${OUT_DIR}/${SAMPLE}_qc_stats.json"

GV="${IMGT}/${CHAIN}V.fasta"
GJ="${IMGT}/${CHAIN}J.fasta"
GD="${IMGT}/${CHAIN}J.fasta"   # IGK/IGL無D基因，借用J基因檔案佔位

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }
source "${SCRIPT_DIR}/igblast_chunk_lib.sh"

# 2026-09-30新增：最終合併檔已完成就直接跳過，不重新執行（避免chunk中間
# 檔案被清理後，重跑時「合併所有chunk」步驟報錯或覆蓋掉已完成的結果）。2026-10起以完成標記判斷。
UPSTREAM_MARKER="${OUT_DIR}/stitch.complete"
if final_guard; then
    info "最終合併檔已完成（${FINAL_AIRR}），視為此rabbit/chain已完成，略過重新執行。"
    info "若確定要強制重跑，請先手動搬移或刪除該檔案與 .complete 標記，並確認chunks/底下的中間檔案仍完整存在。"
    exit 0
fi
require_upstream "$UPSTREAM_MARKER" "04c_stitch_qc_igkl.py ${RABBIT} ${CHAIN}"

[ -f "$ROUND_RETRY" ] || error "找不到 ${ROUND_RETRY}"
[ -d "$CHUNK_DIR" ] || error "找不到 ${CHUNK_DIR}，請先執行 04c_stitch_qc_igkl.py ${RABBIT} ${CHAIN}"
FASTAS=("${CHUNK_DIR}"/chunk_*.fasta)
[ ${#FASTAS[@]} -gt 0 ] || error "找不到chunk fasta"
info "Rabbit ${RABBIT} ${CHAIN} Stitched：找到 ${#FASTAS[@]} 個 chunk fasta"

for FASTA in "${FASTAS[@]}"; do
    igblast_one_chunk "$FASTA"
done

step "合併所有 chunk"
merge_final_airr

"$PYTHON" - "$FINAL_AIRR" "$STATS" "$SAMPLE" "$CHAIN" "" <<'EOF'
import json, csv, sys
final_airr, stats, sample, chain, side = sys.argv[1:6]
n = 0; prod = 0
with open(final_airr) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        n += 1
        if row.get('productive') == 'T': prod += 1
out = {'sample': sample, 'platform': 'AVITI', 'chain': chain}

out.update({'method': 'coordinate-stitching + overlap QC + round-based retry igblast', 'reads_igblast': n, 'productive_count': prod})
with open(stats, 'w') as f:
    json.dump(out, f, indent=2)
print(f'productive={prod/n*100:.1f}%' if n else 'n=0')
EOF
mark_final_complete
info "✓ Rabbit ${RABBIT} ${CHAIN} Stitched igblast 全部完成：${FINAL_AIRR}"
