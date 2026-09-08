#!/bin/bash
# ============================================================
# 02b_run_aviti_igblast.sh
# AVITI pipeline — Part 2: IgBLAST → Report
# 執行前請先完成 02a_run_aviti_qc.sh
# ============================================================

BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/lk-macbookair15/Desktop/Claude_code/BCR_pipeline}"
WORK_DIR="${BASE_DIR}/15651-LT-1_Aviti24"
IMGT="${BASE_DIR}/IMGT_reference_clean"
OUT_DIR="${WORK_DIR}/results"
PYTHON="${BCR_PIPELINE_PYTHON:-/Users/lk-macbookair15/miniconda3/envs/bcr_pipeline/bin/python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IGDATA="${BCR_IGDATA:-/Users/lk-macbookair15/miniconda3/envs/bcr_pipeline/share/igblast}"
IGBLASTN="${IGDATA}/bin/igblastn"
AUX="${IGDATA}/optional_file/rabbit_gl.aux"
BLASTDB="${IGDATA}/database"
THREADS=4
CHAIN_FILTER="${1:-}"   # 可選參數：IGH、IGK 或 IGL；不帶參數則跑全部

# ============================================================
set -e
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }

# ── 輸入驗證 ─────────────────────────────────────────────
if [ -n "$CHAIN_FILTER" ]; then
    [[ "$CHAIN_FILTER" =~ ^(IGH|IGK|IGL)$ ]] || \
        { echo -e "\033[0;31m[ERROR]\033[0m 無效的 chain：$CHAIN_FILTER（只接受 IGH / IGK / IGL）"; exit 1; }
fi

# ── 環境確認 ─────────────────────────────────────────────
step "環境確認"
[ -f "$IGBLASTN" ]                      || error "找不到 igblastn：$IGBLASTN"
[ -d "${IGDATA}/internal_data/rabbit" ] || error "rabbit internal_data 不存在"
[ -f "$AUX" ] && info "rabbit_gl.aux: ✓" || warn "rabbit_gl.aux 未找到"
info "IGDATA：$IGDATA"
info "BLASTDB：$BLASTDB"
info "igblastn：$($IGBLASTN -version 2>&1 | head -1)"

# ── 確認 Part 1 的 FASTA 都存在 ──────────────────────────
[ -n "$CHAIN_FILTER" ] && info "只處理 chain：$CHAIN_FILTER" || info "處理所有 chain：IGH IGK IGL"

for CHAIN in IGH IGK IGL; do
    [ -n "$CHAIN_FILTER" ] && [ "$CHAIN" != "$CHAIN_FILTER" ] && continue
    SAMPLE="15651-LT-1_AVITI_${CHAIN}"
    FASTA="${OUT_DIR}/${SAMPLE}/${CHAIN}_trimmed.fasta"
    [ -f "$FASTA" ] || error "找不到 FASTA：$FASTA\n請先執行 02a_run_aviti_qc.sh"
    SEQ_N=$(grep -c "^>" "$FASTA")
    info "${SAMPLE}：${SEQ_N} sequences ✓"
done

# ══════════════════════════════════════════════════════════
run_aviti_igblast() {
    local CHAIN="$1" LOCUS="$2"
    local SAMPLE="15651-LT-1_AVITI_${CHAIN}"
    local SDIR="${OUT_DIR}/${SAMPLE}"
    local FASTA="${SDIR}/${CHAIN}_trimmed.fasta"
    local AIRR="${SDIR}/${SAMPLE}_igblast.airr.tsv"
    local STATS="${SDIR}/${SAMPLE}_qc_stats.json"

    step "IgBLAST: ${SAMPLE}"

    local IGBLAST_ARGS=(
        -organism rabbit -ig_seqtype Ig
        -query "$FASTA" -out "$AIRR"
        -outfmt 19 -num_threads "$THREADS"
        -extend_align5end -extend_align3end
        -auxiliary_data "$AUX"
    )

    if [ "$LOCUS" = "IGH" ]; then
        IGBLAST_ARGS+=(
            -germline_db_V "${IMGT}/IGHV.fasta"
            -germline_db_D "${IMGT}/IGHD.fasta"
            -germline_db_J "${IMGT}/IGHJ.fasta"
            -min_D_match 5)
    elif [ "$LOCUS" = "IGK" ]; then
        IGBLAST_ARGS+=(
            -germline_db_V "${IMGT}/IGKV.fasta"
            -germline_db_J "${IMGT}/IGKJ.fasta")
    elif [ "$LOCUS" = "IGL" ]; then
        IGBLAST_ARGS+=(
            -germline_db_V "${IMGT}/IGLV.fasta"
            -germline_db_J "${IMGT}/IGLJ.fasta")
    fi

    IGDATA="$IGDATA" BLASTDB="$BLASTDB" "$IGBLASTN" "${IGBLAST_ARGS[@]}"

    local IGBLAST_N
    IGBLAST_N=$(tail -n +2 "$AIRR" | wc -l | tr -d ' ')
    info "  IgBLAST 輸出：${IGBLAST_N} sequences"

    # 更新 QC stats JSON（補入 igblast count）
    "$PYTHON" -c "
import json
with open('${STATS}') as f:
    data = json.load(f)
data['reads_igblast'] = ${IGBLAST_N}
with open('${STATS}','w') as f:
    json.dump(data, f, indent=2)
"
    info "報告生成中..."
    "$PYTHON" "${SCRIPT_DIR}/03_make_report.py" \
        --airr   "$AIRR" \
        --chain  "$LOCUS" \
        --sample "$SAMPLE" \
        --stats  "$STATS" \
        --output "${OUT_DIR}/${SAMPLE}.xlsx"
    info "✓ ${SAMPLE}.xlsx 完成"
}

[[ -z "$CHAIN_FILTER" || "$CHAIN_FILTER" == "IGH" ]] && run_aviti_igblast "IGH" "IGH"
[[ -z "$CHAIN_FILTER" || "$CHAIN_FILTER" == "IGK" ]] && run_aviti_igblast "IGK" "IGK"
[[ -z "$CHAIN_FILTER" || "$CHAIN_FILTER" == "IGL" ]] && run_aviti_igblast "IGL" "IGL"

echo ""
echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}  AVITI pipeline 完成！${NC}"
echo -e "${GREEN}  輸出：${OUT_DIR}${NC}"
echo -e "${GREEN}============================================${NC}"
ls -lh "${OUT_DIR}"/15651-LT-1_AVITI_*.xlsx 2>/dev/null || true
