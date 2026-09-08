#!/bin/bash
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/lk-macbookair15/Desktop/Claude_code/BCR_pipeline}"
WORK_DIR="${BASE_DIR}/15651-LT-1_Aviti24"
R1="${WORK_DIR}/15651-LT-1_R1.fastq"
R2="${WORK_DIR}/15651-LT-1_R2.fastq"
OUT_DIR="${WORK_DIR}/results"
PYTHON="${BCR_PIPELINE_PYTHON:-/Users/lk-macbookair15/miniconda3/envs/bcr_pipeline/bin/python}"
THREADS=4
T7_TSO="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"
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
info()    { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error()   { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()    { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }
diskfree(){ df -h "$WORK_DIR" | awk 'NR==2{print $4}'; }
count_reads() { echo $(( $(wc -l < "$1") / 4 )); }
cleanup() { [ -f "$1" ] && rm -f "$1" && info "  🗑 刪除：$(basename "$1")"; }
rc_iupac() {
    "$PYTHON" -c "
s = '$1'.upper()
comp = {'A':'T','T':'A','G':'C','C':'G','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D'}
print(''.join(comp.get(b,'N') for b in reversed(s)))
"
}

step "環境確認"
command -v cutadapt &>/dev/null || error "找不到 cutadapt"
command -v flash2   &>/dev/null || error "找不到 flash2"
[ -f "$R1" ] || error "找不到 R1：$R1"
[ -f "$R2" ] || error "找不到 R2：$R2"
mkdir -p "$OUT_DIR"
DEMUX_DIR="${OUT_DIR}/demux"
mkdir -p "$DEMUX_DIR"
RAW_TOTAL=$(count_reads "$R1")
info "R1 raw reads：${RAW_TOTAL}"
info "空間：$(diskfree)"
info "開始：$(date '+%Y-%m-%d %H:%M:%S')"
START_TIME=$SECONDS

step "Step 1 / 4 — Demultiplex R2 (e=0.1)"
info "預估：45–60 分鐘"
info "→ IGHV R2..."
cutadapt -G "$IGH_REV1" -G "$IGH_REV2" -e 0.1 --overlap 15 \
    -o "${DEMUX_DIR}/IGH_s1_R1.fastq" -p "${DEMUX_DIR}/IGH_s1_R2.fastq" \
    --untrimmed-output "${DEMUX_DIR}/non_IGH_R1.fastq" \
    --untrimmed-paired-output "${DEMUX_DIR}/non_IGH_R2.fastq" \
    --quiet "$R1" "$R2"
IGH_S1=$(count_reads "${DEMUX_DIR}/IGH_s1_R1.fastq")
info "  IGH S1：${IGH_S1}"

info "→ IGK R2..."
cutadapt -G "$IGK_REV1" -G "$IGK_REV2" -G "$IGK_REV3" -e 0.1 --overlap 15 \
    -o "${DEMUX_DIR}/IGK_s1_R1.fastq" -p "${DEMUX_DIR}/IGK_s1_R2.fastq" \
    --untrimmed-output "${DEMUX_DIR}/non_IGK_R1.fastq" \
    --untrimmed-paired-output "${DEMUX_DIR}/non_IGK_R2.fastq" \
    --quiet "${DEMUX_DIR}/non_IGH_R1.fastq" "${DEMUX_DIR}/non_IGH_R2.fastq"
IGK_S1=$(count_reads "${DEMUX_DIR}/IGK_s1_R1.fastq")
info "  IGK S1：${IGK_S1}"
cleanup "${DEMUX_DIR}/non_IGH_R1.fastq"
cleanup "${DEMUX_DIR}/non_IGH_R2.fastq"

info "→ IGL R2..."
cutadapt -G "$IGL_REV1" -G "$IGL_REV2" -G "$IGL_REV3" -e 0.1 --overlap 15 \
    -o "${DEMUX_DIR}/IGL_s1_R1.fastq" -p "${DEMUX_DIR}/IGL_s1_R2.fastq" \
    --untrimmed-output "${DEMUX_DIR}/unassigned_R1.fastq" \
    --untrimmed-paired-output "${DEMUX_DIR}/unassigned_R2.fastq" \
    --quiet "${DEMUX_DIR}/non_IGK_R1.fastq" "${DEMUX_DIR}/non_IGK_R2.fastq"
IGL_S1=$(count_reads "${DEMUX_DIR}/IGL_s1_R1.fastq")
UN=$(count_reads "${DEMUX_DIR}/unassigned_R1.fastq")
info "  IGL S1：${IGL_S1}"
cleanup "${DEMUX_DIR}/non_IGK_R1.fastq"
cleanup "${DEMUX_DIR}/non_IGK_R2.fastq"
info "Unassigned：${UN}  空間：$(diskfree)"

step "Step 2 / 4 — Rescue Demultiplex R1 (e=0.1)"
info "預估：30–45 分鐘"
info "→ IGHV R1 rescue..."
cutadapt -g "$IGH_REV1" -g "$IGH_REV2" -e 0.1 --overlap 15 \
    -o "${DEMUX_DIR}/IGH_s2_trimR1.fastq" -p "${DEMUX_DIR}/IGH_s2_R2.fastq" \
    --untrimmed-output "${DEMUX_DIR}/non_IGH2_R1.fastq" \
    --untrimmed-paired-output "${DEMUX_DIR}/non_IGH2_R2.fastq" \
    --quiet "${DEMUX_DIR}/unassigned_R1.fastq" "${DEMUX_DIR}/unassigned_R2.fastq"
IGH_S2=$(count_reads "${DEMUX_DIR}/IGH_s2_trimR1.fastq")
info "  IGH rescued：${IGH_S2}"

info "→ IGK R1 rescue..."
cutadapt -g "$IGK_REV1" -g "$IGK_REV2" -g "$IGK_REV3" -e 0.1 --overlap 15 \
    -o "${DEMUX_DIR}/IGK_s2_trimR1.fastq" -p "${DEMUX_DIR}/IGK_s2_R2.fastq" \
    --untrimmed-output "${DEMUX_DIR}/non_IGK2_R1.fastq" \
    --untrimmed-paired-output "${DEMUX_DIR}/non_IGK2_R2.fastq" \
    --quiet "${DEMUX_DIR}/non_IGH2_R1.fastq" "${DEMUX_DIR}/non_IGH2_R2.fastq"
IGK_S2=$(count_reads "${DEMUX_DIR}/IGK_s2_trimR1.fastq")
info "  IGK rescued：${IGK_S2}"
cleanup "${DEMUX_DIR}/non_IGH2_R1.fastq"
cleanup "${DEMUX_DIR}/non_IGH2_R2.fastq"

info "→ IGL R1 rescue..."
cutadapt -g "$IGL_REV1" -g "$IGL_REV2" -g "$IGL_REV3" -e 0.1 --overlap 15 \
    -o "${DEMUX_DIR}/IGL_s2_trimR1.fastq" -p "${DEMUX_DIR}/IGL_s2_R2.fastq" \
    --untrimmed-output /dev/null --untrimmed-paired-output /dev/null \
    --quiet "${DEMUX_DIR}/non_IGK2_R1.fastq" "${DEMUX_DIR}/non_IGK2_R2.fastq"
IGL_S2=$(count_reads "${DEMUX_DIR}/IGL_s2_trimR1.fastq")
info "  IGL rescued：${IGL_S2}"
cleanup "${DEMUX_DIR}/non_IGK2_R1.fastq"
cleanup "${DEMUX_DIR}/non_IGK2_R2.fastq"
cleanup "${DEMUX_DIR}/unassigned_R1.fastq"
cleanup "${DEMUX_DIR}/unassigned_R2.fastq"

TOTAL_ASSIGNED=$(( IGH_S1 + IGK_S1 + IGL_S1 + IGH_S2 + IGK_S2 + IGL_S2 ))
info "Total assigned：${TOTAL_ASSIGNED} / ${RAW_TOTAL} ($("$PYTHON" -c "print(f'{${TOTAL_ASSIGNED}/${RAW_TOTAL}*100:.1f}%')"))"
info "空間：$(diskfree)"

run_aviti_chain() {
    local CHAIN="$1" LOCUS="$2" LEN_MIN="$3"
    shift 3
    local REV_PRIMERS=("$@")
    local SAMPLE="15651-LT-1_AVITI_${CHAIN}"
    local SDIR="${OUT_DIR}/${SAMPLE}"
    mkdir -p "$SDIR"

    step "Step 3+4: ${SAMPLE}  空間：$(diskfree)"
    local S1_N=$(count_reads "${DEMUX_DIR}/${CHAIN}_s1_R1.fastq")
    local S2_N=$(count_reads "${DEMUX_DIR}/${CHAIN}_s2_trimR1.fastq")
    info "  S1：${S1_N}  S2 rescued：${S2_N}"

    info "3a — Swap R1/R2 for rescued reads"
    cp "${DEMUX_DIR}/${CHAIN}_s2_R2.fastq"     "${SDIR}/${CHAIN}_rescued_newR1.fastq"
    cp "${DEMUX_DIR}/${CHAIN}_s2_trimR1.fastq" "${SDIR}/${CHAIN}_rescued_newR2.fastq"
    cleanup "${DEMUX_DIR}/${CHAIN}_s2_R2.fastq"
    cleanup "${DEMUX_DIR}/${CHAIN}_s2_trimR1.fastq"

    info "3b — Trim T7/TSO from R1"
    cutadapt -g "$T7_TSO" -e 0.1 --overlap 15 --quiet \
        -o "${SDIR}/${CHAIN}_s1_R1_trim.fastq" "${DEMUX_DIR}/${CHAIN}_s1_R1.fastq"
    cutadapt -g "$T7_TSO" -e 0.1 --overlap 15 --quiet \
        -o "${SDIR}/${CHAIN}_rescued_newR1_trim.fastq" "${SDIR}/${CHAIN}_rescued_newR1.fastq"
    cleanup "${DEMUX_DIR}/${CHAIN}_s1_R1.fastq"
    cleanup "${SDIR}/${CHAIN}_rescued_newR1.fastq"

    info "3c — Combine"
    cat "${SDIR}/${CHAIN}_s1_R1_trim.fastq" "${SDIR}/${CHAIN}_rescued_newR1_trim.fastq" > "${SDIR}/${CHAIN}_combined_R1.fastq"
    cat "${DEMUX_DIR}/${CHAIN}_s1_R2.fastq" "${SDIR}/${CHAIN}_rescued_newR2.fastq"      > "${SDIR}/${CHAIN}_combined_R2.fastq"
    COMBINED=$(count_reads "${SDIR}/${CHAIN}_combined_R1.fastq")
    info "  Combined：${COMBINED}"
    cleanup "${SDIR}/${CHAIN}_s1_R1_trim.fastq"
    cleanup "${SDIR}/${CHAIN}_rescued_newR1_trim.fastq"
    cleanup "${DEMUX_DIR}/${CHAIN}_s1_R2.fastq"
    cleanup "${SDIR}/${CHAIN}_rescued_newR2.fastq"
    info "  空間：$(diskfree)"

    info "4a — FLASH2 merge"
    flash2 --min-overlap 20 --max-overlap 150 --max-mismatch-density 0.1 \
        --threads "$THREADS" --output-directory "$SDIR" --output-prefix "${CHAIN}" \
        "${SDIR}/${CHAIN}_combined_R1.fastq" "${SDIR}/${CHAIN}_combined_R2.fastq" \
        > "${SDIR}/flash2.log" 2>&1
    local MERGED_FQ="${SDIR}/${CHAIN}.extendedFrags.fastq"
    [ -f "$MERGED_FQ" ] || error "FLASH2 失敗：${SDIR}/flash2.log"
    MERGED=$(count_reads "$MERGED_FQ")
    info "  Merged：${MERGED} ($("$PYTHON" -c "print(f'{${MERGED}/${COMBINED}*100:.1f}%')"))"
    cleanup "${SDIR}/${CHAIN}_combined_R1.fastq"
    cleanup "${SDIR}/${CHAIN}_combined_R2.fastq"
    cleanup "${SDIR}/${CHAIN}.notCombined_1.fastq"
    cleanup "${SDIR}/${CHAIN}.notCombined_2.fastq"
    cleanup "${SDIR}/${CHAIN}.hist"
    cleanup "${SDIR}/${CHAIN}.histogram"
    info "  空間：$(diskfree)"

    info "4b — Adapter Trim + QC"
    local TRIM="${SDIR}/${CHAIN}_trimmed.fastq"
    local CA_ARGS=(-g "$TSO")
    for REV in "${REV_PRIMERS[@]}"; do
        CA_ARGS+=(-a "$(rc_iupac "$REV")")
    done
    CA_ARGS+=(-e 0.1 --overlap 10 --quality-cutoff 20 --minimum-length "$LEN_MIN" --maximum-length "$LEN_MAX" --quiet)
    cutadapt "${CA_ARGS[@]}" "$MERGED_FQ" -o "$TRIM"
    TRIM_N=$(count_reads "$TRIM")
    info "  After trim：${TRIM_N}"
    cleanup "$MERGED_FQ"

    info "4c — FASTQ → FASTA"
    local FASTA="${SDIR}/${CHAIN}_trimmed.fasta"
    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$TRIM" > "$FASTA"
    cleanup "$TRIM"
    info "  FASTA：$(grep -c "^>" "$FASTA") sequences  空間：$(diskfree)"

    "$PYTHON" -c "
import json
with open('${SDIR}/${SAMPLE}_qc_stats.json','w') as f:
    json.dump({'sample':'${SAMPLE}','platform':'AVITI','chain':'${LOCUS}',
               'reads_raw':${RAW_TOTAL},'reads_s1':${S1_N},'reads_s2':${S2_N},
               'reads_combined':${COMBINED},'reads_merged':${MERGED},'reads_trimmed':${TRIM_N}}, f, indent=2)
"
    info "✓ ${SAMPLE} 完成"
}

run_aviti_chain "IGH" "IGH" "$LEN_MIN_IGH" "$IGH_REV1" "$IGH_REV2"
run_aviti_chain "IGK" "IGK" "$LEN_MIN_IGK" "$IGK_REV1" "$IGK_REV2" "$IGK_REV3"
run_aviti_chain "IGL" "IGL" "$LEN_MIN_IGL" "$IGL_REV1" "$IGL_REV2" "$IGL_REV3"

ELAPSED=$(( SECONDS - START_TIME ))
echo ""
echo -e "\033[0;32m============================================\033[0m"
echo -e "\033[0;32m  Part 1 完成！\033[0m"
echo -e "\033[0;32m============================================\033[0m"
info "結束：$(date '+%Y-%m-%d %H:%M:%S')  耗時：$(( ELAPSED/3600 ))h $(( (ELAPSED%3600)/60 ))m"
info "最終空間：$(diskfree)"
for CHAIN in IGH IGK IGL; do
    SAMPLE="15651-LT-1_AVITI_${CHAIN}"
    FASTA="${OUT_DIR}/${SAMPLE}/${CHAIN}_trimmed.fasta"
    [ -f "$FASTA" ] && info "${SAMPLE}：$(grep -c "^>" "$FASTA") sequences ✓" || warn "${SAMPLE}：FASTA 不存在！"
done
echo -e "\033[0;36m確認後執行：bash ~/Desktop/Tem/scripts/02b_run_aviti_igblast.sh\033[0m"
