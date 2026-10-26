#!/bin/bash
# ============================================================
# 02c_run_aviti_r2only_igh.sh
# AVITI pipeline — IGH CDR3 專用替代路徑：跳過 FLASH2 overlap merge，
# 100% 直接用 R2（primer-anchored side）判讀 CDR3。
#
# 背景：現行 02a/02b 只用 FLASH2 merge 成功的 read pairs（約 25-27%）跑
# igblast，merge 成功與否跟 insert 總長（進而跟 CDR3 長度）相關，導致
# 長 CDR3 分子被系統性排除、官方 AVITI CDR3 長度曲線偏短。驗證過程見
# discussion_notes/aviti_r2only_preanalysis_20260915/。
#
# 跟 02a_run_aviti_qc.sh 的差異（其餘demux/primer/igblast邏輯不變）：
#   - 完全不呼叫 FLASH2，100% reads 都送進下游（不是只有merge成功的~25-27%）
#   - R1 只用來做 rescue+swap 判斷方向，trim後不使用（R1 讀不到 CDR3，
#     見前述驗證：R1-only junction_aa 判讀率僅 0.2%）
#   - QC/trim 門檻改成 R2 專用版本：quality-cutoff 20、
#     minimum-length 150、maximum-length 300（已用實測數據驗證此門檻
#     不會切掉CDR3——連卡在150bp下限的read，CDR3判讀完整率仍有98.4%）、
#     新增 -a 去除 revcomp(TSO) read-through（短insert讀穿保護，對應
#     02a對merged序列做的兩端adapter清理邏輯）
#   - 只處理 IGH（IGK/IGL 沿用現行 02a/02b 官方 merged-based 結果）
#
# 用法：bash 02c_run_aviti_r2only_igh.sh <rabbit_number>
#   例：bash 02c_run_aviti_r2only_igh.sh 1
# 長時間執行請用（跟 archive/scripts_archive/run_02b.sh 相同慣例，
# 這批腳本本身不內建 caffeinate，由呼叫時外加）：
#   caffeinate -i bash 02c_run_aviti_r2only_igh.sh 1 2>&1 | tee ~/Desktop/aviti_02c_r1.log
#
# 中斷後重跑：已完成的 chunk 會被跳過（chunks/chunk_XXX.done 標記），
# 只會繼續處理尚未完成的部分，不會整個重跑。
# 2026-10 code review 修正：路徑可用 BCR_* 環境變數覆寫；sed 讀到該 chunk 結尾就停（原本每個
# chunk 都把整個 FASTQ 從頭讀到尾）；記錄輸入檔指紋，輸入換了卻沿用舊的 .done 時停下來；
# chunk_stats.csv 同一 chunk 重跑時取代舊列；結尾統計在0筆時不再因 grep -c 回傳1而失敗。
# ============================================================

RABBIT="$1"
[ -z "$RABBIT" ] && { echo "用法：bash 02c_run_aviti_r2only_igh.sh <rabbit_number>"; exit 1; }

BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
ENV_PREFIX="${BCR_ENV_PREFIX:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env}"
R1="${BASE_DIR}/samples/15651-LT-${RABBIT}_R1.fastq"
R2="${BASE_DIR}/samples/15651-LT-${RABBIT}_R2.fastq"
OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/IGH_R2only"
CHUNK_DIR="${OUT_DIR}/chunks"
TMP_DIR="${OUT_DIR}/tmp"
PYTHON="${BCR_PIPELINE_PYTHON:-${ENV_PREFIX}/bin/python}"

CHUNK_RAW_PAIRS="${BCR_CHUNK_RAW_PAIRS:-1500000}"   # 每塊處理的原始 raw read pairs 數量（demux前，~13.5M/1.5M≈9塊）

T7_TSO="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"
TSO="AAGCAGTGGTATCAACGCAGAGT"
IGH_REV1="CAGTGGGAAGACTGACGGAGCCTTAG"
IGH_REV2="CAGTGGGAAGACTGATGGAGCCTTAG"

set -e
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
comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
print(''.join(comp.get(b,'N') for b in reversed(s)))
"
}
TSO_RC=$(rc_iupac "$TSO")

step "環境確認"
[ -f "$R1" ] || error "找不到 R1：$R1"
[ -f "$R2" ] || error "找不到 R2：$R2"
command -v cutadapt &>/dev/null || error "找不到 cutadapt（記得先 conda activate bcr_pipeline）"
mkdir -p "$CHUNK_DIR" "$TMP_DIR"
[ -f "${OUT_DIR}/chunk_stats.csv" ] || echo "chunk,raw_pairs,igh_demux,qc_trimmed" > "${OUT_DIR}/chunk_stats.csv"
check_fingerprint "${OUT_DIR}/demux.input_fingerprint"
# 已經完整跑完（有 demux.complete、輸入指紋相同）就整個跳過：IGL 會清掉 chunks/，若照 .done 判斷會整份重做，
# 而重做又會讓下游以為上游更新了而跟著重跑。要強制重做請刪除 demux.complete。
if [ -f "${OUT_DIR}"/demux.complete ]; then
    info "demux 先前已完整完成（demux.complete），略過；若要重做請先刪除 demux.complete"
    exit 0
fi
rm -f "${OUT_DIR}"/demux.complete   # 跑完全部 chunk 後才重新寫入；下游 igblast 步驟會檢查它

TOTAL_LINES=$(wc -l < "$R1")
TOTAL_PAIRS=$(( TOTAL_LINES / 4 ))
NUM_CHUNKS=$(( (TOTAL_PAIRS + CHUNK_RAW_PAIRS - 1) / CHUNK_RAW_PAIRS ))
info "Rabbit ${RABBIT}：raw total ${TOTAL_PAIRS} pairs，切成 ${NUM_CHUNKS} 塊（每塊上限 ${CHUNK_RAW_PAIRS} pairs）"

for (( i=0; i<NUM_CHUNKS; i++ )); do
    IDX=$(printf "%03d" "$i")
    DONE_MARK="${CHUNK_DIR}/chunk_${IDX}.done"
    OUT_FASTA="${CHUNK_DIR}/chunk_${IDX}.fasta"

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

    # Step1：demux R2（direct，chain-specific reverse primer 直接命中在raw R2上）
    cutadapt -G "$IGH_REV1" -G "$IGH_REV2" -e 0.1 --overlap 15 \
        -o "${TMP_DIR}/s1_R1.fastq" -p "${TMP_DIR}/s1_R2.fastq" \
        --untrimmed-output "${TMP_DIR}/non_R1.fastq" \
        --untrimmed-paired-output "${TMP_DIR}/non_R2.fastq" \
        --quiet "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"
    rm -f "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"

    # Step2：rescue（primer其實落在raw R1上，swap後才是正確的R2端）
    cutadapt -g "$IGH_REV1" -g "$IGH_REV2" -e 0.1 --overlap 15 \
        -o "${TMP_DIR}/s2_trimR1.fastq" -p "${TMP_DIR}/s2_R2.fastq" \
        --untrimmed-output /dev/null --untrimmed-paired-output /dev/null \
        --quiet "${TMP_DIR}/non_R1.fastq" "${TMP_DIR}/non_R2.fastq"
    rm -f "${TMP_DIR}/non_R1.fastq" "${TMP_DIR}/non_R2.fastq"

    # Step3：彙整兩個來源的「正確方向R2」（s1_R2已經primer-trim好；
    # s2_trimR1是rescue時被trim掉primer的R1，swap後當作R2用）
    # R1（s1_R1 / s2_R2）在R2-only流程中不使用，直接丟棄，不做TSO trim。
    cat "${TMP_DIR}/s1_R2.fastq" "${TMP_DIR}/s2_trimR1.fastq" > "${TMP_DIR}/combined_R2.fastq"
    rm -f "${TMP_DIR}/s1_R1.fastq" "${TMP_DIR}/s1_R2.fastq" "${TMP_DIR}/s2_R2.fastq" "${TMP_DIR}/s2_trimR1.fastq"
    N_IGH=$(count_reads "${TMP_DIR}/combined_R2.fastq")
    info "  IGH demux: ${N_IGH} reads"

    # R2-only QC/trim（跳過FLASH2，這是跟02a唯一的實質差異）
    cutadapt -a "$TSO_RC" -e 0.1 --overlap 10 --quality-cutoff 20 \
        --minimum-length 150 --maximum-length 300 --quiet \
        "${TMP_DIR}/combined_R2.fastq" -o "${TMP_DIR}/trimmed.fastq"
    rm -f "${TMP_DIR}/combined_R2.fastq"
    N_TRIM=$(count_reads "${TMP_DIR}/trimmed.fastq")
    info "  QC/trim 後: ${N_TRIM} reads"

    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "${TMP_DIR}/trimmed.fastq" > "$OUT_FASTA"
    rm -f "${TMP_DIR}/trimmed.fastq"

    replace_stats_row "${OUT_DIR}/chunk_stats.csv" "$IDX" "${IDX},${N_RAW},${N_IGH},${N_TRIM}"
    touch "$DONE_MARK"
    info "  ✓ chunk ${IDX} 完成：${N_RAW} raw → ${N_IGH} IGH demux → ${N_TRIM} QC/trim 後"
done

date "+%Y-%m-%d %H:%M:%S" > "${OUT_DIR}"/demux.complete
rmdir "$TMP_DIR" 2>/dev/null || true
step "Demux + R2-only trim 完成"
TOTAL_FASTA=$(cat "${CHUNK_DIR}"/chunk_*.fasta 2>/dev/null | grep -c "^>" || true)
N_CHUNKS_DONE=$(ls "${CHUNK_DIR}"/chunk_*.fasta 2>/dev/null | wc -l | tr -d ' ')
info "Rabbit ${RABBIT} IGH R2-only 總計：${TOTAL_FASTA} 條序列，${N_CHUNKS_DONE} 個 chunk fasta"
info "接下來執行：bash 02d_run_aviti_r2only_igblast.sh ${RABBIT}"
