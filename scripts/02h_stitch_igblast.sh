#!/bin/bash
# ============================================================
# 02h_stitch_igblast.sh
# 接續 02g_stitch_qc_igh.py：對每個拼接+QC通過的chunk fasta跑igblastn
# （完整V(D)J序列，判定跟merge/R2-only一樣完整的v_call/d_call/j_call/
# junction_aa/productive），逐塊執行、逐塊完成即算入最終檔案，最後合併
# 成一個AIRR TSV。
#
# 這是IGH的第三輪igblast（R1、R2各自先跑過一次，這裡是拼接後的序列
# 再跑一次），目的是用完整涵蓋FR1到J的序列做V基因判定，修正R2-only
# 看不到FR1/CDR1導致IGHV1S40/1S45這類近乎同源基因對判錯的問題
# （2026-09-24~25驗證，見討論記錄）。
#
# 用法：bash 02h_stitch_igblast.sh <rabbit_number>
# 長時間執行請用：
#   caffeinate -i bash 02h_stitch_igblast.sh 1 2>&1 | tee ~/Desktop/aviti_02h_stitch.log
#
# 中斷後重跑：已完成的 chunk 會被跳過（chunks/chunk_XXX.igblast.done 標記）。
# 2026-10 code review 修正：分塊/驗證/合併邏輯改由 igblast_chunk_lib.sh 共用（空 chunk、
# 中斷後重跑、合併不完整、無法處理的read紀錄等問題見該檔說明）；路徑可用 BCR_* 環境變數覆寫。
# ============================================================

RABBIT="$1"
[ -z "$RABBIT" ] && { echo "用法：bash 02h_stitch_igblast.sh <rabbit_number>"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
ENV_PREFIX="${BCR_ENV_PREFIX:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env}"
IMGT="${BASE_DIR}/reference/IMGT_reference_clean"
IGDATA="${BCR_IGDATA:-${ENV_PREFIX}/share/igblast}"
IGBLASTN="${BCR_IGBLASTN:-${ENV_PREFIX}/bin/igblastn}"
AUX="${IGDATA}/optional_file/rabbit_gl.aux"
BLASTDB="${IGDATA}/database"
THREADS=2
IGBLAST_BATCH="${BCR_IGBLAST_BATCH:-25000}"
PYTHON="${BCR_PIPELINE_PYTHON:-${ENV_PREFIX}/bin/python}"

OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/IGH_Stitched"
CHUNK_DIR="${OUT_DIR}/chunks"
SAMPLE="15651-LT-${RABBIT}_AVITI_IGH_Stitched"
FINAL_AIRR="${OUT_DIR}/${SAMPLE}.airr.tsv"
STATS="${OUT_DIR}/${SAMPLE}_qc_stats.json"

set -e
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }
source "${SCRIPT_DIR}/igblast_chunk_lib.sh"

# germline 資料庫名稱外加引號：路徑含空白時 BLAST 會把資料庫名稱裡的空白當成分隔符號
IGBLAST_ARGS=(
    -organism rabbit -ig_seqtype Ig
    -outfmt 19 -num_threads "$THREADS"
    -extend_align5end -extend_align3end
    -auxiliary_data "$AUX"
    -germline_db_V "\"${IMGT}/IGHV.fasta\""
    -germline_db_D "\"${IMGT}/IGHD.fasta\""
    -germline_db_J "\"${IMGT}/IGHJ.fasta\""
    -min_D_match 5
)

# 跟 02d/02f 完全相同的驗證＋切半重跑機制（igblast_chunk_lib.sh）。

# 2026-09-30新增：最終合併檔已完成就直接跳過，不重新執行（避免chunk中間檔案被清理後，
# 重跑時「合併所有chunk」步驟報錯或覆蓋掉已完成的結果）。2026-10起以完成標記判斷。
UPSTREAM_MARKER="${OUT_DIR}/stitch.complete"
if final_guard; then
    info "最終合併檔已完成（${FINAL_AIRR}），視為此rabbit已完成，略過重新執行。"
    info "若確定要強制重跑，請先手動搬移或刪除該檔案與 .complete 標記，並確認chunks/底下的中間檔案仍完整存在。"
    exit 0
fi
require_upstream "$UPSTREAM_MARKER" "02g_stitch_qc_igh.py ${RABBIT}"

step "環境確認"
[ -x "$IGBLASTN" ] || error "找不到 igblastn：$IGBLASTN"
[ -d "${IGDATA}/internal_data/rabbit" ] || error "rabbit internal_data 不存在"
[ -d "$CHUNK_DIR" ] || error "找不到 ${CHUNK_DIR}，請先執行 02g_stitch_qc_igh.py ${RABBIT}"
FASTAS=("${CHUNK_DIR}"/chunk_*.fasta)
[ ${#FASTAS[@]} -gt 0 ] || error "找不到任何 chunk fasta，請先執行 02g_stitch_qc_igh.py ${RABBIT}"
info "找到 ${#FASTAS[@]} 個 chunk fasta"

for FASTA in "${FASTAS[@]}"; do
    igblast_one_chunk "$FASTA"
done

step "合併所有 chunk 的 AIRR TSV"
merge_final_airr

"$PYTHON" - "$FINAL_AIRR" "$STATS" "$SAMPLE" <<'EOF'
import json, csv, sys
final_airr, stats, sample = sys.argv[1:4]
n = 0; prod = 0; has_j = 0
with open(final_airr) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        n += 1
        if row['productive'] == 'T': prod += 1
        if row['junction_aa']: has_j += 1
with open(stats, 'w') as f:
    json.dump({'sample': sample, 'platform': 'AVITI', 'chain': 'IGH',
               'method': 'Coordinate-guided Stitching (R1+R2 stitched by germline coordinates) + QC (overlap identity>=90%)',
               'reads_igblast': n, 'productive_count': prod, 'has_junction_aa': has_j}, f, indent=2)
print(f'productive={prod/n*100:.1f}%  has_junction_aa={has_j/n*100:.1f}%' if n else 'n=0')
EOF
mark_final_complete
info "✓ Rabbit ${RABBIT} IGH 座標拼接+QC 全部完成：${FINAL_AIRR}"
