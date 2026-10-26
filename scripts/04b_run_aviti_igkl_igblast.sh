#!/bin/bash
# ============================================================
# 04b_run_aviti_igkl_igblast.sh
# 接續04a_run_aviti_igkl_demux.sh：對每個chunk fasta跑igblastn。
# 架構完全比照02d（chunk+子批次+checkpoint），但「驗證+重跑」機制改用
# round-based capped retry（igblast_round_retry.py，深度上限5層，超過
# 上限的reads收集起來，等這一輪全部batch跑完後，下一輪重新從頭批次
# 處理），取代02d原本的bash遞迴切半版本——10K測試證實IGL在遞迴式
# 重跑下會深切到10層以上、疑似觸發background process異常終止；
# round-based capped retry已在IGK/IGL各10K reads（r1only+r2only）
# 驗證可穩定收斂（IGK 0重試；IGL最多5輪收斂到個位數不可處理reads）。
#
# D-gene佔位符：IGK/IGL沒有D基因，但igblastn 1.22.0沒給-germline_db_D
# 會報錯，借用該chain自己的J基因fasta當佔位符（已驗證：d_call正確
# 維持空白，不產生假陽性）。
#
# 用法：bash 04b_run_aviti_igkl_igblast.sh <rabbit> <IGK|IGL> <r1only|r2only>
# 中斷後重跑：已完成的chunk/子批次會被跳過（.done標記）
# 2026-10 code review 修正：分塊/驗證/合併邏輯改由 igblast_chunk_lib.sh 共用；無法處理的
# read ID 現在保存在 <OUT_DIR>/irrecoverable/（原本寫在子批次資料夾，隨即被 rm -rf 刪掉）；
# igblast_round_retry.py 改用本腳本所在資料夾的那一份；路徑可用 BCR_* 環境變數覆寫。
# ============================================================
set -e
RABBIT="$1"; CHAIN="$2"; SIDE="$3"
[ -z "$RABBIT" ] || [ -z "$CHAIN" ] || [ -z "$SIDE" ] && {
    echo "用法：bash 04b_run_aviti_igkl_igblast.sh <rabbit> <IGK|IGL> <r1only|r2only>"; exit 1; }
case "$CHAIN" in IGK|IGL) ;; *) echo "CHAIN必須是IGK或IGL"; exit 1 ;; esac
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
ENV_PREFIX="${BCR_ENV_PREFIX:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env}"
IMGT="${BASE_DIR}/reference/IMGT_reference_clean"
PYTHON="${BCR_PIPELINE_PYTHON:-${ENV_PREFIX}/bin/python}"
ROUND_RETRY="${SCRIPT_DIR}/igblast_round_retry.py"
IGBLAST_BATCH="${BCR_IGBLAST_BATCH:-25000}"
IGBLAST_MODE=round

if [ "$SIDE" = "r1only" ]; then SIDE_TAG="R1only"; elif [ "$SIDE" = "r2only" ]; then SIDE_TAG="R2only"; else
    echo "SIDE必須是r1only或r2only"; exit 1; fi

OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/${CHAIN}_${SIDE_TAG}"
CHUNK_DIR="${OUT_DIR}/chunks"
SAMPLE="15651-LT-${RABBIT}_AVITI_${CHAIN}_${SIDE_TAG}"
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
UPSTREAM_MARKER="${OUT_DIR}/demux.complete"
if final_guard; then
    info "最終合併檔已完成（${FINAL_AIRR}），視為此rabbit/chain/side已完成，略過重新執行。"
    info "若確定要強制重跑，請先手動搬移或刪除該檔案與 .complete 標記，並確認chunks/底下的中間檔案仍完整存在。"
    exit 0
fi
require_upstream "$UPSTREAM_MARKER" "04a_run_aviti_igkl_demux.sh ${RABBIT} ${CHAIN}"

[ -f "$ROUND_RETRY" ] || error "找不到 ${ROUND_RETRY}"
[ -d "$CHUNK_DIR" ] || error "找不到 ${CHUNK_DIR}，請先執行 04a_run_aviti_igkl_demux.sh ${RABBIT} ${CHAIN}"
FASTAS=("${CHUNK_DIR}"/chunk_*.fasta)
[ ${#FASTAS[@]} -gt 0 ] || error "找不到chunk fasta"
info "Rabbit ${RABBIT} ${CHAIN} ${SIDE}：找到 ${#FASTAS[@]} 個 chunk fasta"

for FASTA in "${FASTAS[@]}"; do
    igblast_one_chunk "$FASTA"
done

step "合併所有 chunk"
merge_final_airr

"$PYTHON" - "$FINAL_AIRR" "$STATS" "$SAMPLE" "$CHAIN" "$SIDE" <<'EOF'
import json, csv, sys
final_airr, stats, sample, chain, side = sys.argv[1:6]
n = 0; prod = 0
with open(final_airr) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        n += 1
        if row.get('productive') == 'T': prod += 1
out = {'sample': sample, 'platform': 'AVITI', 'chain': chain}
out['side'] = side
out.update({'method': 'round-based capped retry igblast', 'reads_igblast': n, 'productive_count': prod})
with open(stats, 'w') as f:
    json.dump(out, f, indent=2)
print(f'productive={prod/n*100:.1f}%' if n else 'n=0')
EOF
mark_final_complete
info "✓ Rabbit ${RABBIT} ${CHAIN} ${SIDE} igblast 全部完成：${FINAL_AIRR}"
