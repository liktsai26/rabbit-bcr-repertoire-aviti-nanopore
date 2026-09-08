#!/bin/bash
IGDATA="${BCR_IGDATA:-/Users/lk-macbookair15/miniconda3/envs/bcr_pipeline/share/igblast}"
PYTHON="${BCR_PIPELINE_PYTHON:-/Users/lk-macbookair15/miniconda3/envs/bcr_pipeline/bin/python}"
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/lk-macbookair15/Desktop/Claude_code/BCR_pipeline}"
WORK_DIR="${BASE_DIR}/QYWFSM_Nanopore"
IMGT="${BASE_DIR}/IMGT_reference_clean"
OUT_DIR="${WORK_DIR}/results"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BLASTDB="${IGDATA}/database"
AUX="${IGDATA}/optional_file/rabbit_gl.aux"
IGBLASTN="${IGDATA}/bin/igblastn"
THREADS=4

# Bov1（原始 IgK 僅 29 條 reads 已跳過；VVR5GH_1 為重定序結果，4,461 reads，檔案位於專案根目錄）
BOV1_IGH_FASTQ="${WORK_DIR}/QYWFSM_1_Bov1-IgH.fastq"
BOV1_IGL_FASTQ="${WORK_DIR}/QYWFSM_3_Bov1-IgL.fastq"
BOV1_IGK_FASTQ="${BASE_DIR}/VVR5GH_1_Bov1-IGK.fastq"

# Bov2（原始 IgH 僅 9 條 reads 已跳過；Y6Y9DQ_1 為重定序結果，5,487 reads）
BOV2_IGH_FASTQ="${WORK_DIR}/Y6Y9DQ_1_rerun_Bov2-IgH_4.fastq"
BOV2_IGK_FASTQ="${WORK_DIR}/QYWFSM_5_Bov2-IgK.fastq"
BOV2_IGL_FASTQ="${WORK_DIR}/QYWFSM_6_Bov2-IgL.fastq"

# Bov3（三條 chain 齊全）
BOV3_IGH_FASTQ="${WORK_DIR}/QYWFSM_7_Bov3-IgH.fastq"
BOV3_IGK_FASTQ="${WORK_DIR}/QYWFSM_8_Bov3-IgK.fastq"
BOV3_IGL_FASTQ="${WORK_DIR}/QYWFSM_9_Bov3-IgL.fastq"

TSO="AAGCAGTGGTATCAACGCAGAGT"
IGH_REV1="CAGTGGGAAGACTGACGGAGCCTTAG"
IGH_REV2="CAGTGGGAAGACTGATGGAGCCTTAG"
IGK_REV1="TGGTGGGAAGAKGAGGACAGTAGG"
IGK_REV2="TGGTGGGAAGAKGAGGACACTAGG"
IGK_REV3="TGGTGGGAAGAKGAGGACAGAAGG"
IGL_REV1="CAAGGGGGCGACCACAGGCTGAC"
IGL_REV2="GTGAAGGAGTGACTACGGGTTGACC"
IGL_REV3="GAGGGGGTCACCGCGGGCTGAC"

LEN_MIN_IGH=350; LEN_MIN_IGK=320; LEN_MIN_IGL=320; LEN_MAX=600

set -e
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }
count_reads() { echo $(( $(wc -l < "$1") / 4 )); }
rc_iupac() {
    "$PYTHON" -c "
s = '$1'.upper()
comp = {'A':'T','T':'A','G':'C','C':'G','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D'}
print(''.join(comp.get(b,'N') for b in reversed(s)))
"
}

step "環境確認"
command -v NanoFilt &>/dev/null || error "找不到 NanoFilt"
command -v cutadapt &>/dev/null || error "找不到 cutadapt"
[ -f "$IGBLASTN" ]                      || error "找不到 igblastn：$IGBLASTN"
[ -d "${IGDATA}/internal_data/rabbit" ] || error "rabbit internal_data 不存在"
[ -f "$AUX" ] && info "rabbit_gl.aux: ✓" || warn "rabbit_gl.aux 未找到"
mkdir -p "$OUT_DIR"
info "IGDATA：$IGDATA"
info "igblastn：$($IGBLASTN -version 2>&1 | head -1)"

run_nano_chain() {
    local CHAIN="$1" LOCUS="$2" FASTQ="$3" LEN_MIN="$4" GROUP="$5"
    shift 5
    local REV_PRIMERS=("$@")
    local SAMPLE="QYWFSM_${GROUP}_${CHAIN}"
    local SDIR="${OUT_DIR}/${SAMPLE}"
    mkdir -p "$SDIR"

    step "處理 ${SAMPLE}"
    [ -f "$FASTQ" ] || error "找不到：$FASTQ"
    local RAW; RAW=$(count_reads "$FASTQ")
    info "  Raw reads：${RAW}"

    info "Step 1/4 — NanoFilt (-q 10, ${LEN_MIN}–${LEN_MAX} bp)"
    local QC="${SDIR}/${CHAIN}_qc.fastq"
    NanoFilt -q 10 -l "$LEN_MIN" --maxlength "$LEN_MAX" < "$FASTQ" > "$QC"
    local QC_N; QC_N=$(count_reads "$QC")
    info "  After QC：${QC_N} ($("$PYTHON" -c "print(f'{${QC_N}/${RAW}*100:.1f}%')"))"

    info "Step 2/4 — Adapter Trim (-e 0.1, --revcomp)"
    local TRIM="${SDIR}/${CHAIN}_trimmed.fastq"
    local CA_ARGS=(-g "$TSO")
    for REV in "${REV_PRIMERS[@]}"; do
        CA_ARGS+=(-a "$(rc_iupac "$REV")")
    done
    CA_ARGS+=(--revcomp -e 0.1 --overlap 15 --discard-untrimmed
              --minimum-length "$LEN_MIN" --maximum-length "$LEN_MAX" --quiet)
    cutadapt "${CA_ARGS[@]}" "$QC" -o "$TRIM"
    local TRIM_N; TRIM_N=$(count_reads "$TRIM")
    info "  After trim：${TRIM_N} ($("$PYTHON" -c "print(f'{${TRIM_N}/${QC_N}*100:.1f}%')"))"

    info "Step 3/4 — FASTQ → FASTA"
    local FASTA="${SDIR}/${CHAIN}_trimmed.fasta"
    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$TRIM" > "$FASTA"

    info "Step 4/4 — IgBLAST"
    local AIRR="${SDIR}/${SAMPLE}_igblast.airr.tsv"
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

    "$PYTHON" -c "
import json
with open('${SDIR}/${SAMPLE}_qc_stats.json','w') as f:
    json.dump({'sample':'${SAMPLE}','platform':'Nanopore','chain':'${LOCUS}',
               'reads_raw':${RAW},'reads_qc':${QC_N},
               'reads_trimmed':${TRIM_N},'reads_igblast':${IGBLAST_N}}, f, indent=2)
"
    info "報告生成中..."
    "$PYTHON" "${SCRIPT_DIR}/03_make_report.py" \
        --airr  "$AIRR" --chain "$LOCUS" --sample "$SAMPLE" \
        --stats "${SDIR}/${SAMPLE}_qc_stats.json" \
        --output "${OUT_DIR}/${SAMPLE}.xlsx"
    info "✓ ${SAMPLE}.xlsx 完成"
}

run_nano_chain "IGHV" "IGH" "$BOV1_IGH_FASTQ" "$LEN_MIN_IGH" "Bov1" "$IGH_REV1" "$IGH_REV2"
run_nano_chain "IGK"  "IGK" "$BOV1_IGK_FASTQ" "$LEN_MIN_IGK" "Bov1" "$IGK_REV1" "$IGK_REV2" "$IGK_REV3"
run_nano_chain "IGL"  "IGL" "$BOV1_IGL_FASTQ" "$LEN_MIN_IGL" "Bov1" "$IGL_REV1" "$IGL_REV2" "$IGL_REV3"
run_nano_chain "IGHV" "IGH" "$BOV2_IGH_FASTQ" "$LEN_MIN_IGH" "Bov2" "$IGH_REV1" "$IGH_REV2"
run_nano_chain "IGK"  "IGK" "$BOV2_IGK_FASTQ" "$LEN_MIN_IGK" "Bov2" "$IGK_REV1" "$IGK_REV2" "$IGK_REV3"
run_nano_chain "IGL"  "IGL" "$BOV2_IGL_FASTQ" "$LEN_MIN_IGL" "Bov2" "$IGL_REV1" "$IGL_REV2" "$IGL_REV3"
run_nano_chain "IGHV" "IGH" "$BOV3_IGH_FASTQ" "$LEN_MIN_IGH" "Bov3" "$IGH_REV1" "$IGH_REV2"
run_nano_chain "IGK"  "IGK" "$BOV3_IGK_FASTQ" "$LEN_MIN_IGK" "Bov3" "$IGK_REV1" "$IGK_REV2" "$IGK_REV3"
run_nano_chain "IGL"  "IGL" "$BOV3_IGL_FASTQ" "$LEN_MIN_IGL" "Bov3" "$IGL_REV1" "$IGL_REV2" "$IGL_REV3"

echo -e "\n${GREEN}============================================${NC}"
echo -e "${GREEN}  Nanopore pipeline 完成！輸出：${OUT_DIR}${NC}"
echo -e "${GREEN}============================================${NC}"
ls -lh "${OUT_DIR}"/QYWFSM_*.xlsx 2>/dev/null || true
