"""Compare pipeline outputs in a sandbox against the simulation truth.

Usage: python evaluate.py --truth truth_molecules.tsv --read-truth read_truth.tsv
                          --results SANDBOX/BCR/results --out metrics.json [--label NAME]
Per chain (summed over rabbits) it reports stage yields, stitched-sequence correctness
against the true molecule ('ACATGGG' + cDNA up to the primer), V-gene / junction_aa accuracy
of R2-only vs stitched calls, and how often wrong joins are still called productive.
"""
import argparse
import collections
import glob
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sim"))
from common import gene, read_airr, read_tsv  # noqa: E402

from Bio.Align import PairwiseAligner  # noqa: E402

_al = PairwiseAligner()
_al.mode = "global"
_al.match_score, _al.mismatch_score = 1, -1
_al.open_gap_score, _al.extend_gap_score = -3, -1
_al.end_gap_score = 0


def mol_of(seq_id):
    parts = seq_id.split(":")
    return parts[5] if len(parts) > 5 else seq_id


def classify(stitched, expected):
    if stitched == expected:
        return "exact", 0, 0
    if len(stitched) == len(expected):
        return "substitutions_only", sum(a != b for a, b in zip(stitched, expected)), 0
    # semi-global: stitched may start a few bases late / end early only via trimming, so count internal gaps
    aln = _al.align(stitched, expected)[0]
    ins = dele = 0
    coords = aln.aligned
    (s_blocks, e_blocks) = coords
    for k in range(1, len(s_blocks)):
        ds = s_blocks[k][0] - s_blocks[k - 1][1]
        de = e_blocks[k][0] - e_blocks[k - 1][1]
        ins += ds
        dele += de
    if ins == 0 and dele == 0:
        return "truncated_or_extended", 0, len(stitched) - len(expected)
    return "internal_indel", ins + dele, len(stitched) - len(expected)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--read-truth", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    truth = {r["mol_id"]: r for r in read_tsv(a.truth)}
    reads = {r["read_id"]: r for r in read_tsv(a.read_truth)}
    out = {"label": a.label}
    for chain, letter in (("IGH", "H"), ("IGK", "K"), ("IGL", "L")):
        m = collections.Counter()
        art_examples = []
        n_truth = sum(1 for t in truth.values() if t["chain"] == chain)
        m["truth_molecules"] = n_truth
        for rabbit_dir in sorted(glob.glob(os.path.join(a.results, "15651-LT-*_Aviti24"))):
            rabbit = os.path.basename(rabbit_dir).split("-")[2].split("_")[0]
            pre = os.path.join(rabbit_dir, f"{chain}_{{}}", f"15651-LT-{rabbit}_AVITI_{chain}_{{}}.airr.tsv")
            r2p, r1p, stp = pre.format("R2only", "R2only"), pre.format("R1only", "R1only"), pre.format("Stitched", "Stitched")
            r2 = {x["sequence_id"]: x for x in read_airr(r2p)} if os.path.exists(r2p) else {}
            r1 = {x["sequence_id"]: x for x in read_airr(r1p)} if os.path.exists(r1p) else {}
            st = read_airr(stp) if os.path.exists(stp) else []
            m["r2only_rows"] += len(r2)
            m["r1only_rows"] += len(r1)
            m["stitched_rows"] += len(st)
            for sid, x in r2.items():
                mid = mol_of(sid)
                if mid[0] != letter:
                    m["r2only_wrong_chain_or_junk"] += 1
                    continue
                t = truth[mid]
                m["r2only_productive"] += x["productive"] == "T"
                m["r2only_v_gene_correct"] += gene(x["v_call"]) == gene(t["truth_v_call"])
                m["r2only_junction_aa_correct"] += x["junction_aa"] == t["truth_junction_aa"]
            for x in st:
                sid = x["sequence_id"]
                mid = mol_of(sid)
                if mid[0] != letter:
                    m["stitched_wrong_chain_or_junk"] += 1
                    continue
                t = truth[mid]
                rt = reads.get(sid, {})
                expected = "ACATGGG" + t["cdna"]
                cls, n, dlen = classify(x["sequence"].upper(), expected)
                read_indels = int(rt.get("tso_side_indels", 0) or 0) + int(rt.get("primer_side_indels", 0) or 0)
                m["stitched_" + cls] += 1
                vg_ok = gene(x["v_call"]) == gene(t["truth_v_call"])
                m["stitched_v_gene_correct"] += vg_ok
                m["stitched_v_allele_correct"] += x["v_call"].split(",")[0] == t["truth_v_call"].split(",")[0]
                m["stitched_junction_aa_correct"] += x["junction_aa"] == t["truth_junction_aa"]
                m["stitched_productive"] += x["productive"] == "T"
                r1x, r2x = r1.get(sid), r2.get(sid)
                differ = bool(r1x and r2x and r1x["v_call"].split(",")[0] != r2x["v_call"].split(",")[0])
                m["stitched_r1r2_primary_v_differ"] += differ
                if cls == "internal_indel" and read_indels == 0:
                    m["stitch_artifact_indel_no_read_indel"] += 1
                    m["stitch_artifact_and_productive_T"] += x["productive"] == "T"
                    m["stitch_artifact_and_r1r2_v_differ"] += differ
                    m[f"stitch_artifact_len_delta_{dlen:+d}"] += 1
                    if len(art_examples) < 5:
                        art_examples.append({"read": sid, "len_delta": dlen, "productive": x["productive"],
                                             "r1_v": r1x["v_call"] if r1x else "", "r2_v": r2x["v_call"] if r2x else "",
                                             "truth_v": t["truth_v_call"], "stitched_v": x["v_call"]})
                if cls == "internal_indel" and read_indels > 0:
                    m["stitched_indel_from_read_indel"] += 1
        d = dict(m)
        good = d.get("stitched_rows", 0) - d.get("stitched_wrong_chain_or_junk", 0)
        if good:
            for k in ("stitched_exact", "stitched_substitutions_only", "stitch_artifact_indel_no_read_indel",
                      "stitched_v_gene_correct", "stitched_junction_aa_correct", "stitched_productive",
                      "stitched_r1r2_primary_v_differ"):
                d[k + "_pct"] = round(100 * d.get(k, 0) / good, 2)
        r2good = d.get("r2only_rows", 0) - d.get("r2only_wrong_chain_or_junk", 0)
        if r2good:
            d["r2only_v_gene_correct_pct"] = round(100 * d.get("r2only_v_gene_correct", 0) / r2good, 2)
            d["r2only_junction_aa_correct_pct"] = round(100 * d.get("r2only_junction_aa_correct", 0) / r2good, 2)
        if n_truth:
            d["stitched_yield_pct_of_truth"] = round(100 * good / n_truth, 2)
        d["artifact_examples"] = art_examples
        out[chain] = d
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    for chain in ("IGH", "IGK", "IGL"):
        d = out[chain]
        print(chain, {k: v for k, v in d.items() if k.endswith("_pct") or k in (
            "truth_molecules", "stitched_rows", "stitch_artifact_indel_no_read_indel", "stitch_artifact_and_productive_T")})


if __name__ == "__main__":
    main()
