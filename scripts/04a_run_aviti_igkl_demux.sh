#!/bin/bash
# ============================================================
# 04a_run_aviti_igkl_demux.sh
# 全量demux + trim：從raw fastq直接抽出IGK或IGL的R1角色(TSO-anchored)/
# R2角色(primer-anchored)reads，分別做對應trim，供後續座標拼接法用
# （R1-only igblast、R2-only igblast、拼接+QC）。
#
# 邏輯完全比照已驗證過的demux_test.sh（LT-1、1M raw pairs測試版）與
# 02a/02c/02e(IGH R1-only/R2-only)的demux慣例：cutadapt兩階段
# demux+rescue（primer在R2直接命中 / primer落在R1需要rescue+swap），
# R1角色只做T7/TSO trim，R2角色做Q20+150-300bp trim。
# 差異：這裡一次跑就同時產出R1角色與R2角色兩種chunk fasta（同一次
# demux分類結果分流輸出），不像02c/02e各自獨立重跑一次分類。
#
# 用法：bash 04a_run_aviti_igkl_demux.sh <rabbit_number> <IGK|IGL>
# 中斷後重跑：已完成的chunk會被跳過（chunks/chunk_XXX.done標記）
# 2026-10 code review 修正：路徑可用 BCR_* 環境變數覆寫；sed 讀到該 chunk 結尾就停（原本每個
# chunk 都把整個 FASTQ 從頭讀到尾）；記錄輸入檔指紋，輸入換了卻沿用舊的 .done 時停下來；
# chunk_stats.csv 同一 chunk 重跑時取代舊列；結尾統計在0筆時不再因 grep -c 回傳1而失敗。
# ============================================================
set -e
RABBIT="$1"
CHAIN="$2"
[ -z "$RABBIT" ] || [ -z "$CHAIN" ] && { echo "用法：bash 04a_run_aviti_igkl_demux.sh <rabbit_number> <IGK|IGL>"; exit 1; }

BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
ENV_PREFIX="${BCR_ENV_PREFIX:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env}"
R1="${BASE_DIR}/samples/15651-LT-${RABBIT}_R1.fastq"
R2="${BASE_DIR}/samples/15651-LT-${RABBIT}_R2.fastq"
OUT_DIR1="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/${CHAIN}_R1only"
OUT_DIR2="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/${CHAIN}_R2only"
CHUNK_DIR1="${OUT_DIR1}/chunks"
CHUNK_DIR2="${OUT_DIR2}/chunks"
TMP_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/${CHAIN}_demux_tmp"
PYTHON="${BCR_PIPELINE_PYTHON:-${ENV_PREFIX}/bin/python}"
export PATH="${ENV_PREFIX}/bin:${PATH}"

CHUNK_RAW_PAIRS="${BCR_CHUNK_RAW_PAIRS:-1500000}"   # 跟02c/02e相同慣例

T7_TSO="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"
TSO="AAGCAGTGGTATCAACGCAGAGT"

if [ "$CHAIN" = "IGK" ]; then
    REV1="TGGTGGGAAGAKGAGGACAGTAGG"
    REV2="TGGTGGGAAGAKGAGGACACTAGG"
    REV3="TGGTGGGAAGAKGAGGACAGAAGG"
elif [ "$CHAIN" = "IGL" ]; then
    REV1="CAAGGGGGCGACCACAGGCTGAC"
    REV2="GTGAAGGAGTGACTACGGGTTGACC"
    REV3="GAGGGGGTCACCGCGGGCTGAC"
else
    echo "CHAIN必須是IGK或IGL"; exit 1
fi

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }
count_reads() { echo $(( $(wc -l < "$1") / 4 )); }

# 輸入檔指紋（檔名、大小、第一筆與最後一筆read名稱、chunk大小）。已經有 .done 標記但指紋不同，代表輸入換過，
# 舊的 chunk 不能沿用。
input_fingerprint() {
    local f
    for f in "$R1" "$R2"; do
        printf '%s\t%s\t%s\t%s\n' "$(basename "$f")" "$(wc -c < "$f" | tr -d ' ')" "$(head -1 "$f")" "$(tail -n 4 "$f" | head -1)"
    done
    echo "CHUNK_RAW_PAIRS=${CHUNK_RAW_PAIRS}"
}
check_fingerprint() {
    local fp_file="$1" now
    now=$(input_fingerprint)
    if [ -f "$fp_file" ]; then
        [ "$(cat "$fp_file")" = "$now" ] || error "輸入檔或chunk設定跟上次不同（見 ${fp_file}）。若要重跑，請先移走 $(dirname "$fp_file") 底下舊的 chunks/ 與 chunk_stats.csv。"
    else
        echo "$now" > "$fp_file"
    fi
}
# 取代 chunk_stats.csv 中同一個 chunk 的舊列（中斷後重跑不會重複計數）
replace_stats_row() {
    local stats="$1" idx="$2" line="$3"
    { head -1 "$stats"; tail -n +2 "$stats" | grep -v "^${idx}," || true; echo "$line"; } > "${stats}.tmp"
    mv "${stats}.tmp" "$stats"
}
rc_iupac() {
    "$PYTHON" -c "
s = '$1'.upper()
comp = {'A':'T','T':'A','G':'C','C':'G','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D'}
print(''.join(comp.get(b,'N') for b in reversed(s)))
"
}
TSO_RC=$(rc_iupac "$TSO")

step "環境確認：Rabbit ${RABBIT} ${CHAIN}"
[ -f "$R1" ] || error "找不到 R1：$R1"
[ -f "$R2" ] || error "找不到 R2：$R2"
command -v cutadapt &>/dev/null || error "找不到 cutadapt"
mkdir -p "$CHUNK_DIR1" "$CHUNK_DIR2" "$TMP_DIR"
[ -f "${OUT_DIR1}/chunk_stats.csv" ] || echo "chunk,raw_pairs,demux,r1only_trimmed,r2only_trimmed" > "${OUT_DIR1}/chunk_stats.csv"
check_fingerprint "${OUT_DIR1}/demux.input_fingerprint"
# 已經完整跑完（有 demux.complete、輸入指紋相同）就整個跳過：IGL 會清掉 chunks/，若照 .done 判斷會整份重做，
# 而重做又會讓下游以為上游更新了而跟著重跑。要強制重做請刪除 demux.complete。
if [ -f "${OUT_DIR1}"/demux.complete ] && [ -f "${OUT_DIR2}"/demux.complete ]; then
    info "demux 先前已完整完成（demux.complete），略過；若要重做請先刪除 demux.complete"
    exit 0
fi
rm -f "${OUT_DIR1}"/demux.complete   # 跑完全部 chunk 後才重新寫入；下游 igblast 步驟會檢查它
rm -f "${OUT_DIR2}"/demux.complete   # 跑完全部 chunk 後才重新寫入；下游 igblast 步驟會檢查它

TOTAL_PAIRS=$(( $(wc -l < "$R1") / 4 ))
NUM_CHUNKS=$(( (TOTAL_PAIRS + CHUNK_RAW_PAIRS - 1) / CHUNK_RAW_PAIRS ))
info "raw total ${TOTAL_PAIRS} pairs，切成 ${NUM_CHUNKS} 塊（每塊上限 ${CHUNK_RAW_PAIRS} pairs）"

for (( i=0; i<NUM_CHUNKS; i++ )); do
    IDX=$(printf "%03d" "$i")
    DONE_MARK="${CHUNK_DIR1}/chunk_${IDX}.done"
    R1_FASTA="${CHUNK_DIR1}/chunk_${IDX}.fasta"
    R2_FASTA="${CHUNK_DIR2}/chunk_${IDX}.fasta"

    if [ -f "$DONE_MARK" ]; then
        info "chunk ${IDX} 已完成，跳過"
        continue
    fi

    step "Chunk ${IDX} / $(printf '%03d' $((NUM_CHUNKS-1)))"
    START_LINE=$(( i * CHUNK_RAW_PAIRS * 4 + 1 ))
    END_LINE=$(( (i+1) * CHUNK_RAW_PAIRS * 4 ))
    sed -n "${START_LINE},${END_LINE}p;${END_LINE}q" "$R1" > "${TMP_DIR}/raw_R1.fastq"
    sed -n "${START_LINE},${END_LINE}p;${END_LINE}q" "$R2" > "${TMP_DIR}/raw_R2.fastq"
    N_RAW=$(count_reads "${TMP_DIR}/raw_R1.fastq")
    if [ "$N_RAW" -eq 0 ]; then
        info "chunk ${IDX} 沒有資料（已到檔案結尾），結束"
        rm -f "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"
        break
    fi
    info "  raw: ${N_RAW} pairs"

    cutadapt -G "$REV1" -G "$REV2" -G "$REV3" -e 0.1 --overlap 15 \
        -o "${TMP_DIR}/s1_R1.fastq" -p "${TMP_DIR}/s1_R2.fastq" \
        --untrimmed-output "${TMP_DIR}/non_R1.fastq" \
        --untrimmed-paired-output "${TMP_DIR}/non_R2.fastq" \
        --quiet "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"
    rm -f "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"

    cutadapt -g "$REV1" -g "$REV2" -g "$REV3" -e 0.1 --overlap 15 \
        -o "${TMP_DIR}/s2_trimR1.fastq" -p "${TMP_DIR}/s2_R2.fastq" \
        --untrimmed-output /dev/null --untrimmed-paired-output /dev/null \
        --quiet "${TMP_DIR}/non_R1.fastq" "${TMP_DIR}/non_R2.fastq"
    rm -f "${TMP_DIR}/non_R1.fastq" "${TMP_DIR}/non_R2.fastq"

    # R2角色（primer-anchored side）= s1_R2 + s2_trimR1(rescued)
    cat "${TMP_DIR}/s1_R2.fastq" "${TMP_DIR}/s2_trimR1.fastq" > "${TMP_DIR}/combined_R2.fastq"
    # R1角色（TSO-anchored side）= s1_R1 + s2_R2(rescued)
    cat "${TMP_DIR}/s1_R1.fastq" "${TMP_DIR}/s2_R2.fastq" > "${TMP_DIR}/combined_R1.fastq"
    rm -f "${TMP_DIR}/s1_R1.fastq" "${TMP_DIR}/s1_R2.fastq" "${TMP_DIR}/s2_R2.fastq" "${TMP_DIR}/s2_trimR1.fastq"
    N_DEMUX=$(count_reads "${TMP_DIR}/combined_R2.fastq")
    info "  ${CHAIN} demux: ${N_DEMUX} reads"

    cutadapt -g "$T7_TSO" -e 0.1 --overlap 15 --minimum-length 1 --quiet \
        "${TMP_DIR}/combined_R1.fastq" -o "${TMP_DIR}/r1_trimmed.fastq"
    rm -f "${TMP_DIR}/combined_R1.fastq"
    N_R1=$(count_reads "${TMP_DIR}/r1_trimmed.fastq")

    cutadapt -a "$TSO_RC" -e 0.1 --overlap 10 --quality-cutoff 20 \
        --minimum-length 150 --maximum-length 300 --quiet \
        "${TMP_DIR}/combined_R2.fastq" -o "${TMP_DIR}/r2_trimmed.fastq"
    rm -f "${TMP_DIR}/combined_R2.fastq"
    N_R2=$(count_reads "${TMP_DIR}/r2_trimmed.fastq")

    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "${TMP_DIR}/r1_trimmed.fastq" > "$R1_FASTA"
    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "${TMP_DIR}/r2_trimmed.fastq" > "$R2_FASTA"
    rm -f "${TMP_DIR}/r1_trimmed.fastq" "${TMP_DIR}/r2_trimmed.fastq"

    replace_stats_row "${OUT_DIR1}/chunk_stats.csv" "$IDX" "${IDX},${N_RAW},${N_DEMUX},${N_R1},${N_R2}"
    touch "$DONE_MARK"
    touch "${CHUNK_DIR2}/chunk_${IDX}.done"
    info "  ✓ chunk ${IDX} 完成：raw=${N_RAW} demux=${N_DEMUX} R1only=${N_R1} R2only=${N_R2}"
done

date "+%Y-%m-%d %H:%M:%S" > "${OUT_DIR1}"/demux.complete
date "+%Y-%m-%d %H:%M:%S" > "${OUT_DIR2}"/demux.complete
rmdir "$TMP_DIR" 2>/dev/null || true
step "Demux + trim 完成"
N_R1_TOTAL=$(cat "${CHUNK_DIR1}"/chunk_*.fasta 2>/dev/null | grep -c "^>" || true)
N_R2_TOTAL=$(cat "${CHUNK_DIR2}"/chunk_*.fasta 2>/dev/null | grep -c "^>" || true)
info "Rabbit ${RABBIT} ${CHAIN}：R1-only總計 ${N_R1_TOTAL} 條，R2-only總計 ${N_R2_TOTAL} 條"
