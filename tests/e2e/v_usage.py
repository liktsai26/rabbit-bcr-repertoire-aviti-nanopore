"""V-gene usage of stitched output vs the truth molecules (per chain): total variation distance
and per-gene observed/true frequency ratios for the most common genes.
Usage: python v_usage.py TRUTH_TSV LABEL=RESULTS_DIR [LABEL=RESULTS_DIR ...]"""
import collections, glob, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sim"))
from common import gene, read_airr, read_tsv  # noqa: E402
truth = read_tsv(sys.argv[1])
runs = [a.split("=", 1) for a in sys.argv[2:]]
for chain in ("IGH", "IGK", "IGL"):
    t = collections.Counter(gene(r["truth_v_call"]) for r in truth if r["chain"] == chain)
    tn = sum(t.values())
    top = [g for g, _ in t.most_common(5)]
    print(f"== {chain}  truth n={tn}  " + "  ".join(f"{g}={t[g]/tn:.3f}" for g in top))
    for label, res in runs:
        c = collections.Counter(); n = 0
        for p in glob.glob(os.path.join(res, "*", f"{chain}_Stitched", f"*_{chain}_Stitched.airr.tsv")):
            for r in read_airr(p):
                if r["sequence_id"].split(":")[5][0] != chain[2] or r["productive"] != "T":
                    continue
                c[gene(r["v_call"])] += 1; n += 1
        if not n:
            continue
        tvd = 0.5 * sum(abs(c[g] / n - t[g] / tn) for g in set(c) | set(t))
        print(f"  {label:10s} n={n:5d} TVD={tvd:.3f}  " + "  ".join(f"{g}:{(c[g]/n)/(t[g]/tn):.2f}x" for g in top))
