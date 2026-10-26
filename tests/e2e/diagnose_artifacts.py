"""List stitched reads whose length differs from the truth without any read indel, with the
reason stitch_core assigns, to see which code path still produces junction artifacts."""
import collections, sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(here), "sim")); sys.path.insert(0, here)
truth_tsv, read_truth, results, scripts = sys.argv[1:5]
sys.path.insert(0, scripts)
from common import read_airr, read_tsv  # noqa
from evaluate import classify  # noqa
import stitch_core as sc  # noqa
truth = {r["mol_id"]: r for r in read_tsv(truth_tsv)}
reads = {r["read_id"]: r for r in read_tsv(read_truth)}
conv = lambda r: {k: (int(r[k]) if k in sc.COORD_COLS else r[k]) for k in r}  # noqa
for chain in ("IGH", "IGK", "IGL"):
    mapper = sc.AlleleMapper(os.path.join(os.environ.get("BCR_REF", ""), f"{chain}V.fasta"))
    cnt = collections.Counter(); allr = collections.Counter()
    for rab in "123":
        base = f"{results}/15651-LT-{rab}_Aviti24"
        r1 = {x["sequence_id"]: x for x in read_airr(f"{base}/{chain}_R1only/15651-LT-{rab}_AVITI_{chain}_R1only.airr.tsv")}
        r2 = {x["sequence_id"]: x for x in read_airr(f"{base}/{chain}_R2only/15651-LT-{rab}_AVITI_{chain}_R2only.airr.tsv")}
        for x in read_airr(f"{base}/{chain}_Stitched/15651-LT-{rab}_AVITI_{chain}_Stitched.airr.tsv"):
            sid = x["sequence_id"]; mid = sid.split(":")[5]
            if mid[0] != chain[2]: continue
            s, reason, ov, diff = sc.stitch_pair(conv(r1[sid]), conv(r2[sid]), mapper)
            allr[(reason, diff)] += 1
            cls, n, d = classify(x["sequence"].upper(), "ACATGGG" + truth[mid]["cdna"])
            rt = reads[sid]
            if cls == "internal_indel" and int(rt["tso_side_indels"]) + int(rt["primer_side_indels"]) == 0:
                cnt[(reason, diff, d)] += 1
    print(chain, "artifacts by (reason, diff_allele, len_delta):", dict(cnt.most_common()))
    print(chain, "all accepted by (reason, diff_allele):", dict(allr.most_common()))
