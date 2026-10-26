# ============================================================
# igblast_chunk_lib.sh
# 02d/02f/02h/04b/04d 與 Nanopore 腳本共用的 igblast 分塊執行、驗證、合併邏輯
# （原本每支腳本各複製一份，2026-10 code review 時發現各份已經不一致，合併成一份）。
# 由呼叫端 source；需要呼叫端先設定：
#   IGBLASTN IGDATA BLASTDB IGBLAST_ARGS(陣列)      — bash 切半重跑模式（IGBLAST_MODE=halving）
#   PYTHON ROUND_RETRY GV GD GJ                     — round-based retry 模式（IGBLAST_MODE=round）
#   CHUNK_DIR OUT_DIR FINAL_AIRR IGBLAST_BATCH，以及 info/warn/error/step 四個訊息函式
#
# 跟原版的差異：
#   - grep -c 在0筆時回傳1，原版在 set -e 下會讓空的 chunk 直接中止腳本（02d/02f/04b），
#     02h/04d 的 "|| echo 0" 則會得到 "0\n0" 使判斷失效；這裡一律用 count_fasta()。
#   - igblastn 非0結束（例如被系統砍掉）也走切半重跑，不再讓 set -e 直接結束整支腳本；
#     stderr 存到 ${OUT_DIR}/igblastn_stderr.log 而不是丟掉。
#   - 合併時表頭取第一個非空檔案，空檔/只有表頭的子批次不會吃掉下一個檔案的第一筆資料。
#   - 子批次切分先寫完再留 .split.done，中斷後重跑不會沿用切到一半的子批次而漏序列。
#   - 無法處理的 read ID 存到 ${OUT_DIR}/irrecoverable/（不在 chunks/ 底下，不會被清掉），
#     最後彙整成 <sample>_irrecoverable_ids.txt。
#   - 最終合併先寫 .tmp、核對「輸出筆數 + 不可處理筆數 = 輸入筆數」後才改名，並寫
#     <final>.complete 完成標記；「最終檔存在就跳過」改成檢查完成標記，避免把中斷時
#     寫到一半的檔案當成已完成。
# ============================================================
shopt -s nullglob

count_fasta() {
    local n
    n=$(grep -c '^>' "$1" 2>/dev/null || true)
    echo "${n:-0}"
}

airr_rows() {
    local n=0
    if [ -s "$1" ]; then n=$(( $(wc -l < "$1") - 1 )); fi
    [ "$n" -lt 0 ] && n=0
    echo "$n"
}

append_airr() {
    # append_airr <src> <dst>：dst 還是空的就連表頭一起寫入，否則只接資料列
    [ -s "$1" ] || return 0
    if [ -s "$2" ]; then tail -n +2 "$1" >> "$2"; else cat "$1" > "$2"; fi
}

split_fasta() {
    # split_fasta <input.fasta> <out_prefix> <batch_size>  (每筆序列剛好2行：header+seq)
    awk -v batch="$3" -v prefix="$2" '
    /^>/ {
        if (n % batch == 0) {
            if (fi > 0) close(outfile)
            outfile = sprintf("%s_%04d.fasta", prefix, fi)
            fi++
        }
        n++
    }
    { print > outfile }
    ' "$1"
}

# igblastn 1.22.0 已知問題：單次呼叫累積處理夠多序列後會對特定序列丟出
# CEXCEPTION(NULL pointer)，該序列之後排隊的部分被靜默丟棄且exit code仍是0。
# 核對輸入輸出行數，對不上就切半重跑；切到剩1條仍失敗就記錄到 IRR 檔。
igblast_verified() {
    local FASTA="$1" OUT="$2" IRR="$3"
    local N_IN N_OUT TMP_OUT="${OUT}.attempt"
    N_IN=$(count_fasta "$FASTA")
    rm -f "$TMP_OUT"
    IGDATA="$IGDATA" BLASTDB="$BLASTDB" "$IGBLASTN" "${IGBLAST_ARGS[@]}" \
        -query "$FASTA" -out "$TMP_OUT" 2>>"${IGBLAST_STDERR_LOG:-${OUT_DIR}/igblastn_stderr.log}" || true
    [ -f "$TMP_OUT" ] || : > "$TMP_OUT"
    N_OUT=$(airr_rows "$TMP_OUT")

    if [ "$N_OUT" -eq "$N_IN" ]; then
        mv "$TMP_OUT" "$OUT"
        return 0
    fi

    if [ "$N_IN" -le 1 ]; then
        warn "  無法處理的序列（igblastn持續失敗），記錄為缺失並跳過：$(grep '^>' "$FASTA")"
        grep '^>' "$FASTA" | cut -c2- | awk '{print $1}' >> "$IRR"
        head -1 "$TMP_OUT" > "$OUT"
        rm -f "$TMP_OUT"
        return 0
    fi

    warn "  批次驗證失敗（輸入${N_IN}條，輸出${N_OUT}條），切半重跑：$(basename "$FASTA")"
    local HALF=$(( (N_IN + 1) / 2 ))
    local SPLIT_DIR HALF_FASTA
    SPLIT_DIR=$(mktemp -d "$(dirname "$FASTA")/retry.XXXXXX")
    split_fasta "$FASTA" "${SPLIT_DIR}/part" "$HALF"
    local MERGED="${OUT}.merged"
    : > "$MERGED"
    for HALF_FASTA in "${SPLIT_DIR}"/part_*.fasta; do
        igblast_verified "$HALF_FASTA" "${HALF_FASTA}.airr.tsv" "$IRR"
        append_airr "${HALF_FASTA}.airr.tsv" "$MERGED"
    done
    rm -rf "$SPLIT_DIR"
    mv "$MERGED" "$OUT"
    rm -f "$TMP_OUT"
}

# 最終檔已完成（有完成標記）就回傳0（呼叫端據此跳過）。
# 最終檔存在但沒有完成標記：chunks/ 還在就從 chunks 重新合併；chunks/ 已被清掉則停下來，
# 請使用者確認（可能是上次合併中斷留下的半截檔案）。
final_guard() {
    if [ -s "$FINAL_AIRR" ] && [ -f "${FINAL_AIRR}.complete" ]; then
        if [ -n "${UPSTREAM_MARKER:-}" ] && [ "$UPSTREAM_MARKER" -nt "${FINAL_AIRR}.complete" ]; then
            warn "上游步驟在最終檔完成之後又重跑過（${UPSTREAM_MARKER} 較新），重新檢查各 chunk 並重新合併"
            rm -f "${FINAL_AIRR}.complete"
            return 1
        fi
        return 0
    fi
    if [ -e "$FINAL_AIRR" ] && [ ! -d "$CHUNK_DIR" ]; then
        error "${FINAL_AIRR} 存在但沒有完成標記（.complete），且 chunks/ 已不在，無法確認是否完整（可能是舊版腳本產生，或上次合併被中斷）。確認無誤後可執行：touch \"${FINAL_AIRR}.complete\"；否則請移走該檔並從上游重跑。"
    fi
    return 1
}

# 上游步驟（demux 或拼接）沒有完整跑完就停下來：原版只用 glob 抓現有的 chunk FASTA，上游中斷
# （或 `| tee` 把上游的錯誤碼吃掉）時會只處理一部分 chunk，之後又被當成已完成。
require_upstream() {
    [ -f "$1" ] && return 0
    if [ "${BCR_SKIP_UPSTREAM_CHECK:-0}" = "1" ]; then
        warn "略過上游完成檢查（BCR_SKIP_UPSTREAM_CHECK=1）：$1"
        return 0
    fi
    error "上游步驟尚未完整跑完（找不到 $1）。請先跑完：$2；若確定上游結果完整（例如由舊版腳本產生），可設 BCR_SKIP_UPSTREAM_CHECK=1 再執行。"
}

igblast_one_chunk() {
    local FASTA="$1"
    local IDX AIRR DONE_MARK SUB_DIR N_SEQ SUB_FASTA SUB_ID SUB_AIRR SUB_DONE SUB_N SUB_OUT_N IRR N_IRR
    IDX=$(basename "$FASTA" .fasta | sed 's/chunk_//')
    AIRR="${CHUNK_DIR}/chunk_${IDX}.airr.tsv"
    DONE_MARK="${CHUNK_DIR}/chunk_${IDX}.igblast.done"
    SUB_DIR="${CHUNK_DIR}/chunk_${IDX}_sub"

    N_SEQ=$(count_fasta "$FASTA")
    if [ -f "$DONE_MARK" ]; then
        # .done 裡記錄當時的輸入筆數；上游重跑使筆數改變時重做這個 chunk（舊版留下的空白 .done 照舊信任）
        if [ ! -s "$DONE_MARK" ] || [ "$(cat "$DONE_MARK")" = "$N_SEQ" ]; then
            rm -rf "$SUB_DIR"
            info "chunk ${IDX} igblast 已完成，跳過"
            return 0
        fi
        warn "chunk ${IDX} 的輸入筆數（${N_SEQ}）與上次完成時（$(cat "$DONE_MARK")）不同，重做"
        rm -f "$DONE_MARK" "$AIRR"
        rm -rf "$SUB_DIR"
    fi
    if [ "$N_SEQ" -eq 0 ]; then
        info "chunk ${IDX} 沒有任何序列，標記完成並跳過"
        rm -f "$AIRR"
        echo 0 > "$DONE_MARK"
        return 0
    fi

    step "IgBLAST chunk ${IDX}（${N_SEQ} 條序列，切成 ${IGBLAST_BATCH}/批）"
    if [ ! -f "${SUB_DIR}/.split.done" ]; then
        rm -rf "$SUB_DIR"
        mkdir -p "$SUB_DIR"
        split_fasta "$FASTA" "${SUB_DIR}/sub" "$IGBLAST_BATCH"
        touch "${SUB_DIR}/.split.done"
    fi
    local SUB_FASTAS=("${SUB_DIR}"/sub_*.fasta)
    info "  分成 ${#SUB_FASTAS[@]} 個子批次"
    mkdir -p "${OUT_DIR}/irrecoverable"

    for SUB_FASTA in "${SUB_FASTAS[@]}"; do
        SUB_ID=$(basename "$SUB_FASTA" .fasta)
        SUB_AIRR="${SUB_DIR}/${SUB_ID}.airr.tsv"
        SUB_DONE="${SUB_DIR}/${SUB_ID}.done"
        IRR="${OUT_DIR}/irrecoverable/chunk_${IDX}_${SUB_ID}.txt"
        if [ -f "$SUB_DONE" ] && [ -f "$SUB_AIRR" ]; then
            info "    ${SUB_ID} 已完成，跳過"
            continue
        fi
        SUB_N=$(count_fasta "$SUB_FASTA")
        rm -f "$IRR"
        if [ "${IGBLAST_MODE:-halving}" = "round" ]; then
            info "    跑 ${SUB_ID}（${SUB_N} 條，round-based retry）..."
            BCR_IGBLAST_STDERR_LOG="${OUT_DIR}/igblastn_stderr.log" \
                "$PYTHON" "$ROUND_RETRY" "$SUB_FASTA" "$SUB_AIRR" "$GV" "$GD" "$GJ"
            if [ -f "${SUB_DIR}/${SUB_ID}_irrecoverable_ids.txt" ]; then
                cp "${SUB_DIR}/${SUB_ID}_irrecoverable_ids.txt" "$IRR"
            fi
        else
            info "    跑 ${SUB_ID}（${SUB_N} 條）..."
            igblast_verified "$SUB_FASTA" "$SUB_AIRR" "$IRR"
        fi
        SUB_OUT_N=$(airr_rows "$SUB_AIRR")
        N_IRR=0
        [ -f "$IRR" ] && N_IRR=$(wc -l < "$IRR" | tr -d ' ')
        if [ $(( SUB_OUT_N + N_IRR )) -ne "$SUB_N" ]; then
            error "    ${SUB_ID} 筆數對不上：輸入${SUB_N}條 ≠ 輸出${SUB_OUT_N}條 + 不可處理${N_IRR}條"
        fi
        [ "$N_IRR" -eq 0 ] || warn "    ${SUB_ID} 有 ${N_IRR} 條無法處理，ID 記錄於 ${IRR}"
        touch "$SUB_DONE"
    done

    # 合併這個 chunk 底下所有子批次的 AIRR TSV（先寫暫存檔，完成後才改名）
    : > "${AIRR}.tmp"
    for SUB_AIRR in "${SUB_DIR}"/sub_*.airr.tsv; do
        append_airr "$SUB_AIRR" "${AIRR}.tmp"
    done
    mv "${AIRR}.tmp" "$AIRR"
    echo "$N_SEQ" > "$DONE_MARK"   # 先標記完成再刪子批次：刪到一半被中斷時，重跑不會用殘缺的子批次覆蓋完整結果
    rm -rf "$SUB_DIR"
    info "  ✓ chunk ${IDX} 完成：$(airr_rows "$AIRR") rows"
}

# 合併所有 chunk → FINAL_AIRR（先寫 .tmp，核對筆數後才改名）。完成標記由呼叫端在統計
# 檔寫完後呼叫 mark_final_complete 寫入。
merge_final_airr() {
    local TMP="${FINAL_AIRR}.tmp" A FASTA N_IN=0 N_OUT N_IRR
    local IRR_ALL="${FINAL_AIRR%.airr.tsv}_irrecoverable_ids.txt"
    : > "$TMP"
    for A in "${CHUNK_DIR}"/chunk_*.airr.tsv; do
        append_airr "$A" "$TMP"
    done
    if [ ! -s "$TMP" ]; then
        rm -f "$TMP"
        error "沒有任何 chunk 產生 AIRR TSV 結果（所有 chunk 都是空的？），請檢查上游輸出"
    fi
    for FASTA in "${CHUNK_DIR}"/chunk_*.fasta; do
        N_IN=$(( N_IN + $(count_fasta "$FASTA") ))
    done
    # 逐檔 cat：nullglob 下沒有任何 .txt 時，`cat 目錄/*.txt` 會變成不帶參數的 cat 而去讀 stdin（卡住）
    : > "${IRR_ALL}.tmp"
    for A in "${OUT_DIR}"/irrecoverable/*.txt; do
        cat "$A" >> "${IRR_ALL}.tmp"
    done
    sort -u "${IRR_ALL}.tmp" > "$IRR_ALL"
    rm -f "${IRR_ALL}.tmp"
    N_IRR=$(wc -l < "$IRR_ALL" | tr -d ' ')
    N_OUT=$(airr_rows "$TMP")
    if [ $(( N_OUT + N_IRR )) -ne "$N_IN" ]; then
        error "合併後筆數對不上：輸入${N_IN}條 ≠ 輸出${N_OUT}條 + 不可處理${N_IRR}條（未改名，暫存檔保留在 ${TMP}）"
    fi
    mv "$TMP" "$FINAL_AIRR"
    MERGE_SUMMARY="input=${N_IN} rows=${N_OUT} irrecoverable=${N_IRR}"
    info "合併完成：${FINAL_AIRR}（${N_OUT} rows；不可處理 ${N_IRR} 條）"
}

mark_final_complete() {
    echo "${MERGE_SUMMARY:-}" > "${FINAL_AIRR}.complete"
}
