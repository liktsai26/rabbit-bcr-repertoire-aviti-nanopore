#!/bin/bash
# Fetch public test inputs and build the IgBLAST germline databases used by the tests.
#   - IMGT/V-QUEST reference directory, Oryctolagus cuniculus IG (V/D/J, gapped)
#   - IMGT/GENE-DB rabbit constant genes (IGHG, IGKC, IGLC) to place primers realistically
#   - OAS unpaired data unit Banerjee_2017/SRR5060324_Heavy_Bulk (real rabbit IgH, 71,806 seqs)
# Usage: bash tests/fetch_test_data.sh [WORK_DIR]   (default: ./work)
# Needs: curl, an env with makeblastdb + edit_imgt_file.pl (see tests/README.md).
set -euo pipefail
WORK="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/work}"
ENV_BIN="${BCR_TEST_ENV_BIN:-${WORK}/env/bin}"
DL="${WORK}/downloads"
REF="${WORK}/reference/IMGT_reference_clean"
mkdir -p "${DL}/imgt" "${DL}/constant" "${DL}/oas" "$REF"

for g in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ; do
    [ -s "${DL}/imgt/${g}.fasta" ] || curl -fsSL --retry 3 -o "${DL}/imgt/${g}.fasta" \
        "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Oryctolagus_cuniculus/IG/${g}.fasta"
done
for g in IGHG IGKC IGLC; do
    [ -s "${DL}/constant/q_14.1_${g}.html" ] || curl -fsSL --retry 3 -o "${DL}/constant/q_14.1_${g}.html" \
        "https://www.imgt.org/genedb/GENElect?query=14.1+${g}&species=Oryctolagus+cuniculus"
done
OAS_URL="https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Banerjee_2017/csv/SRR5060324_Heavy_Bulk.csv.gz"
[ -s "${DL}/oas/SRR5060324_Heavy_Bulk.csv.gz" ] || curl -fsSL --retry 3 -o "${DL}/oas/SRR5060324_Heavy_Bulk.csv.gz" "$OAS_URL"

# Clean IMGT FASTA the standard IgBLAST way (allele-name headers, IMGT gaps removed) and build DBs.
for g in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ; do
    "${ENV_BIN}/perl" "${ENV_BIN}/edit_imgt_file.pl" "${DL}/imgt/${g}.fasta" > "${REF}/${g}.fasta" 2>/dev/null
    "${ENV_BIN}/makeblastdb" -parse_seqids -dbtype nucl -in "${REF}/${g}.fasta" -out "${REF}/${g}.fasta" > /dev/null
done
( cd "$DL" && find . -type f \( -name '*.fasta' -o -name '*.gz' -o -name '*.html' \) -print0 | sort -z | xargs -0 shasum -a 256 ) > "${DL}/SHA256SUMS"
echo "test data ready under ${WORK}"
