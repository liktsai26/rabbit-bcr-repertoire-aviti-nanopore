"""Measure, on EXISTING R1-only / R2-only AIRR files, how often the original 02g/04c join differs
from the fixed join (stitch_core.py). No IgBLAST run and no truth needed, so it can be pointed at
the production results to see the real size of the junction-shift problem.

Usage:
  python compare_stitch_methods.py R1ONLY.airr.tsv R2ONLY.airr.tsv GERMLINE_V.fasta [--max-rows N]
Prints: pairs compared, how many joins differ, the length differences (multiples of 3 indicate the
allele-coordinate shift), and the allele pairs responsible.
"""
import argparse
import collections
import csv
import os
import sys

csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "scripts"))
import stitch_core as sc  # noqa: E402

COLS = ("v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end")


def original_join(r1, r2):
    """02g_stitch_qc_igh.py lines 132-144 (upstream revision 26339d3), QC omitted."""
    gap = r2["v_germline_start"] - r1["v_germline_end"]
    if gap > 0:
        return None
    r2_cut = int(r2["v_sequence_start"] + (r1["v_germline_end"] + 1 - r2["v_germline_start"]))
    return r1["sequence"][:r1["v_sequence_end"]] + r2["sequence"][max(r2_cut - 1, 0):]


def load(path, need_productive=False, max_rows=None):
    out = {}
    with open(path, newline="") as f:
        for i, r in enumerate(csv.DictReader(f, delimiter="\t")):
            if max_rows and i >= max_rows:
                break
            if need_productive and r.get("productive") != "T":
                continue
            try:
                for c in COLS:
                    r[c] = int(r[c])
            except (ValueError, TypeError):
                continue
            out[r["sequence_id"]] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("r1")
    ap.add_argument("r2")
    ap.add_argument("germline_v")
    ap.add_argument("--max-rows", type=int, default=None)
    a = ap.parse_args()
    r1 = load(a.r1, max_rows=a.max_rows)
    r2 = load(a.r2, need_productive=True, max_rows=a.max_rows)
    mapper = sc.AlleleMapper(a.germline_v)
    n = same = differ = 0
    deltas = collections.Counter()
    pairs = collections.Counter()
    for sid, b in r2.items():
        x = r1.get(sid)
        if x is None:
            continue
        old = original_join(x, b)
        new, reason, _, _ = sc.stitch_pair(x, b, mapper)
        if old is None or new is None:
            continue
        n += 1
        if old == new:
            same += 1
            continue
        differ += 1
        d = len(old) - len(new)
        deltas[d] += 1
        if d:
            pairs[(x["v_call"].split(",")[0], b["v_call"].split(",")[0])] += 1
    print(f"pairs stitched by both methods: {n}")
    if not n:
        return
    shifted = sum(v for k, v in deltas.items() if k)
    print(f"identical joins: {same} ({100 * same / n:.1f}%)")
    print(f"length differs (original minus fixed): {shifted} ({100 * shifted / n:.1f}%)  "
          f"{ {k: v for k, v in deltas.most_common(10) if k} }")
    print(f"same length, different bases (overlap now taken from R2): {deltas.get(0, 0)}")
    print("allele pairs (R1 call, R2 call) behind the length differences:")
    for (p, q), c in pairs.most_common(10):
        print(f"  {p:14s} {q:14s} {c}")


if __name__ == "__main__":
    main()
