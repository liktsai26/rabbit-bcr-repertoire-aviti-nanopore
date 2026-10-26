"""Simulate Nanopore-like FASTQ for the MZPD8T samples (02_run_nanopore.sh layout).

read = [0-15 nt flank] + TSO + 'ACATGGG' + cDNA + RC(primer) + [0-15 nt flank], random strand,
ONT-like errors (substitution / insertion / deletion, default 1.5% / 1% / 1.5%). 5% of reads lose
the TSO end and 5% lose the primer end (incomplete molecules).

Usage: python simulate_nanopore.py TRUTH_TSV OUT_DIR [--per-sample 300] [--seed S]
Writes OUT_DIR/samples/MZPD8T_Nanopore/MZPD8T_{1..6}_*.fastq and OUT_DIR/nano_truth.tsv
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import PRIMERS, TSO, rc, read_tsv, resolve_iupac, write_tsv  # noqa: E402

SAMPLES = [("1_WT-IGHV", "IGH"), ("2_WT-IGK", "IGK"), ("3_WT-IGL", "IGL"),
           ("4_Immu-IGHV", "IGH"), ("5_Immu-IGK", "IGK"), ("6_Immu-IGL", "IGL")]


def ont_errors(seq, rng, sub, ins, dele):
    out, quals = [], []
    for c in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + sub:
            out.append(rng.choice([b for b in "ACGT" if b != c]))
            quals.append(rng.randint(4, 12))
        else:
            out.append(c)
            quals.append(max(5, min(40, int(rng.gauss(22, 5)))))
        if rng.random() < ins:
            out.append(rng.choice("ACGT"))
            quals.append(rng.randint(3, 10))
    return "".join(out), "".join(chr(q + 33) for q in quals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("truth")
    ap.add_argument("out")
    ap.add_argument("--per-sample", type=int, default=300)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--sub", type=float, default=0.015)
    ap.add_argument("--ins", type=float, default=0.01)
    ap.add_argument("--del", dest="dele", type=float, default=0.015)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    mols = read_tsv(a.truth)
    by_chain = {c: [m for m in mols if m["chain"] == c] for c in ("IGH", "IGK", "IGL")}
    sdir = os.path.join(a.out, "samples", "MZPD8T_Nanopore")
    os.makedirs(sdir, exist_ok=True)
    rows = []
    for sname, chain in SAMPLES:
        pick = rng.sample(by_chain[chain], min(a.per_sample, len(by_chain[chain])))
        with open(os.path.join(sdir, f"MZPD8T_{sname}.fastq"), "w") as f:
            for k, m in enumerate(pick):
                primer = resolve_iupac(PRIMERS[chain][int(m["primer_idx"])], rng)
                left = TSO + "ACATGGG" if rng.random() > 0.05 else ""
                right = rc(primer) if rng.random() > 0.05 else ""
                mol = left + m["cdna"] + right
                flank = lambda: "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 15)))  # noqa: E731
                read = flank() + mol + flank()
                strand = "+" if rng.random() < 0.5 else "-"
                if strand == "-":
                    read = rc(read)
                seq, qual = ont_errors(read, rng, a.sub, a.ins, a.dele)
                rid = f"{sname}:{m['mol_id']}:{k}"
                f.write(f"@{rid} runid=sim\n{seq}\n+\n{qual}\n")
                rows.append({"read_id": rid, "sample": sname, "chain": chain, "mol_id": m["mol_id"], "strand": strand,
                             "has_tso": int(bool(left)), "has_primer": int(bool(right))})
    write_tsv(os.path.join(a.out, "nano_truth.tsv"), rows, list(rows[0].keys()))
    print(f"{len(rows)} nanopore reads written to {sdir}")


if __name__ == "__main__":
    main()
