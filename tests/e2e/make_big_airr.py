"""Replicate a sandbox's R1-only/R2-only IGH AIRR rows (needed columns only) to N rows with unique
IDs, to measure memory/time of the stitching step at a realistic scale.
Usage: python make_big_airr.py SRC_RESULTS_DIR DST_BASE N"""
import csv, os, sys
csv.field_size_limit(sys.maxsize)
src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
cols = ["sequence_id", "sequence", "v_call", "v_cigar", "productive", "stop_codon", "vj_in_frame",
        "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end"]
for side in ("R1only", "R2only"):
    rows = list(csv.DictReader(open(f"{src}/15651-LT-1_Aviti24/IGH_{side}/15651-LT-1_AVITI_IGH_{side}.airr.tsv"), delimiter="\t"))
    od = f"{dst}/results/15651-LT-9_Aviti24/IGH_{side}"; os.makedirs(od, exist_ok=True)
    with open(f"{od}/15651-LT-9_AVITI_IGH_{side}.airr.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(cols)
        for k in range(n):
            r = rows[k % len(rows)]
            w.writerow([f"{r['sequence_id']}#{k // len(rows)}"] + [r[c] for c in cols[1:]])
