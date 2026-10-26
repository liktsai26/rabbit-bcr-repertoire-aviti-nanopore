"""Stitched yield by true junction_aa length (does the gap filter drop long CDR3s?).
Usage: python cdr3_bias.py TRUTH_TSV CHAIN LABEL=RESULTS_DIR ..."""
import collections, glob, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sim"))
from common import read_airr, read_tsv  # noqa: E402
truth = {r["mol_id"]: r for r in read_tsv(sys.argv[1])}
chain = sys.argv[2]
bins = [(0, 13), (14, 15), (16, 17), (18, 19), (20, 99)]
def b(n):
    return next(f"{lo}-{hi}" for lo, hi in bins if lo <= n <= hi)
tot = collections.Counter(b(len(t["truth_junction_aa"])) for t in truth.values() if t["chain"] == chain)
print("bin        " + "  ".join(f"{lo}-{hi}:{tot[f'{lo}-{hi}']}" for lo, hi in bins))
for arg in sys.argv[3:]:
    label, res = arg.split("=", 1)
    seen = set()
    for p in glob.glob(os.path.join(res, "*", f"{chain}_Stitched", f"*_{chain}_Stitched.airr.tsv")):
        seen |= {r["sequence_id"].split(":")[5] for r in read_airr(p)}
    got = collections.Counter(b(len(truth[m]["truth_junction_aa"])) for m in seen if m in truth and truth[m]["chain"] == chain)
    print(f"{label:10s} " + "  ".join(f"{lo}-{hi}:{got[f'{lo}-{hi}']/max(1,tot[f'{lo}-{hi}']):.2f}" for lo, hi in bins))
