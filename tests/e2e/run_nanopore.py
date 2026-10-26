"""Run 02_run_nanopore.sh from SCRIPTS on simulated Nanopore reads and score it against truth.

Usage: python run_nanopore.py --scripts DIR --sim SIM_DIR --sandbox NEW_DIR
Scores per sample: reads surviving each step, V gene / junction_aa / productive agreement with
IgBLAST on the true molecule, and how many trimmed reads still carry TSO or primer sequence.
"""
import argparse
import collections
import json
import os
import shutil
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(here), "sim"))
from common import ENV_PREFIX, PRIMERS, REF_DIR, TSO, gene, rc, read_airr, read_fasta, read_tsv  # noqa: E402


def has_approx(seq, motif, max_mm=3):
    m = len(motif)
    for i in range(0, len(seq) - m + 1):
        mm = 0
        for a, b in zip(seq[i:i + m], motif):
            if a != b and b not in "KMRYSWN":
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--sim", required=True)
    ap.add_argument("--sandbox", required=True)
    ap.add_argument("--truth", required=True)
    a = ap.parse_args()
    if os.path.exists(a.sandbox):
        sys.exit("sandbox exists")
    base = os.path.join(a.sandbox, "BCR")
    os.makedirs(os.path.join(base, "reference"))
    os.symlink(REF_DIR, os.path.join(base, "reference", "IMGT_reference_clean"))
    os.symlink(os.path.abspath(os.path.join(a.sim, "samples")), os.path.join(base, "samples"))
    shutil.copytree(a.scripts, os.path.join(base, "script"))
    env = dict(os.environ, BCR_PIPELINE_ROOT=base, BCR_PIPELINE_PYTHON=os.path.join(ENV_PREFIX, "bin", "python"),
               BCR_IGDATA=os.path.join(ENV_PREFIX, "share", "igblast"),
               PATH=os.path.join(ENV_PREFIX, "bin") + os.pathsep + os.environ["PATH"])
    with open(os.path.join(a.sandbox, "run.log"), "w") as lf:
        p = subprocess.run(["/bin/bash", os.path.join(base, "script", "02_run_nanopore.sh")], stdout=lf,
                           stderr=subprocess.STDOUT, env=env)
    truth = {r["mol_id"]: r for r in read_tsv(a.truth)}
    out = {"exit": p.returncode}
    res = os.path.join(base, "results", "MZPD8T_Nanopore", "results")
    for d in sorted(os.listdir(res)) if os.path.isdir(res) else []:
        sd = os.path.join(res, d)
        chain = "IGH" if "IGHV" in d else d.split("_")[-1]
        m = collections.Counter()
        fa = os.path.join(sd, f"{d.split('_')[-1]}_trimmed.fasta")
        seqs = read_fasta(fa) if os.path.exists(fa) else {}
        m["trimmed_reads"] = len(seqs)
        for sid, s in seqs.items():
            head, tail = s[:60], s[-60:]
            if has_approx(head, TSO) or has_approx(head, rc(TSO)) or has_approx(tail, TSO) or has_approx(tail, rc(TSO)):
                m["residual_TSO"] += 1
            if any(has_approx(tail, rc(pr)) or has_approx(head, pr) for pr in PRIMERS[chain]):
                m["residual_primer"] += 1
        airr = os.path.join(sd, f"{d}_igblast.airr.tsv")
        rows = read_airr(airr) if os.path.exists(airr) else []
        m["airr_rows"] = len(rows)
        for r in rows:
            t = truth[r["sequence_id"].split(":")[1]]
            m["v_gene_correct"] += gene(r["v_call"]) == gene(t["truth_v_call"])
            m["junction_aa_correct"] += r["junction_aa"] == t["truth_junction_aa"]
            m["productive"] += r["productive"] == "T"
        out[d] = dict(m)
    with open(os.path.join(a.sandbox, "metrics.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
