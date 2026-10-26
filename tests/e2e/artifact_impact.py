"""Quantify downstream impact of stitch-junction artifacts: v_identity shift, CDR2/FR3 changes,
alignment indel counts, relative to IgBLAST on the true molecule."""
import collections, glob, os, statistics, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sim"))
from common import read_airr, read_tsv, run_igblast  # noqa: E402
truth_tsv, results, workdir = sys.argv[1:4]
truth = {r["mol_id"]: r for r in read_tsv(truth_tsv)}
os.makedirs(workdir, exist_ok=True)
for chain in ("IGH", "IGK", "IGL"):
    rows = []
    for p in glob.glob(os.path.join(results, "*", f"{chain}_Stitched", f"*_{chain}_Stitched.airr.tsv")):
        rows += read_airr(p)
    rows = [r for r in rows if r["sequence_id"].split(":")[5][0] == chain[2]]
    # truth annotation with full columns
    fa = os.path.join(workdir, f"{chain}_truth.fasta")
    with open(fa, "w") as f:
        for m in {r["sequence_id"].split(":")[5] for r in rows}:
            f.write(f">{m}\nACATGGG{truth[m]['cdna']}\n")
    out = os.path.join(workdir, f"{chain}_truth.airr.tsv")
    run_igblast(fa, out, chain, threads=6)
    t = {r["sequence_id"]: r for r in read_airr(out)}
    groups = collections.defaultdict(list)
    for r in rows:
        tm = t[r["sequence_id"].split(":")[5]]
        exp_len = len(tm["sequence"])
        cls = "len_changed" if len(r["sequence"]) != exp_len else "len_same"
        groups[cls].append((r, tm))
    print(f"== {chain}")
    for cls, lst in groups.items():
        dv = [float(r["v_identity"]) - float(tm["v_identity"]) for r, tm in lst if r["v_identity"] and tm["v_identity"]]
        cdr2 = sum(r["cdr2_aa"] != tm["cdr2_aa"] for r, tm in lst)
        fwr3 = sum(r["fwr3_aa"] != tm["fwr3_aa"] for r, tm in lst)
        aln_gap = sum(("D" in r["v_cigar"] or "I" in r["v_cigar"]) and not ("D" in tm["v_cigar"] or "I" in tm["v_cigar"]) for r, tm in lst)
        seq_aa = sum(r["sequence_alignment_aa"] != tm["sequence_alignment_aa"] for r, tm in lst)
        print(f"  {cls}: n={len(lst)}  mean v_identity shift={statistics.mean(dv):+.2f} pts (median {statistics.median(dv):+.2f})"
              f"  cdr2_aa differs={cdr2}  fwr3_aa differs={fwr3}  new V-alignment indel={aln_gap}  full aa seq differs={seq_aa}")
