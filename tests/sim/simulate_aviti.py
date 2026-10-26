"""Simulate AVITI-like 2x300 paired FASTQ for rabbits 1..3 from truth_molecules.tsv.

Library model (matches the demux logic in 02c/02e/04a):
  insert (sense) = T7_TSO + 'ACATGGG' + cDNA(5'UTR+leader+V(D)J+C spacer) + RC(reverse primer)
  R1 = insert[:L], R2 = RC(insert)[:L]; ~30% of pairs are swapped (primer on R1) to exercise the
  rescue/swap branch. Errors: position-dependent substitutions drawn from the simulated Phred
  qualities, plus rare 1-nt indels (rate configurable). All chains are mixed in one FASTQ pair, with
  4% random junk pairs and 1% primer dimers.

Usage: python simulate_aviti.py TRUTH_TSV OUT_DIR [--indel-rate 2e-5] [--seed S] [--rabbits 3]
Writes OUT_DIR/samples/15651-LT-{r}_R1.fastq / _R2.fastq and OUT_DIR/read_truth.tsv
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import PRIMERS, T7_TSO, rc, read_tsv, resolve_iupac, write_tsv  # noqa: E402

L = 300
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def quals(rng):
    qend = rng.uniform(8, 36)
    out = []
    for i in range(L):
        q = 38 - (38 - qend) * (i / (L - 1)) ** 3 + rng.gauss(0, 1.5)
        out.append(max(2, min(41, int(round(q)))))
    return out


def corrupt(seq, qs, indel_rate, rng):
    s = list(seq)
    n_sub = n_ind = 0
    for i, q in enumerate(qs[:len(s)]):
        if rng.random() < 10 ** (-q / 10):
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
            n_sub += 1
    if indel_rate > 0:
        out = []
        for c in s:
            r = rng.random()
            if r < indel_rate / 2:
                n_ind += 1          # deletion
                continue
            out.append(c)
            if r > 1 - indel_rate / 2:
                out.append(rng.choice("ACGT"))
                n_ind += 1          # insertion
        s = out
    return "".join(s), n_sub, n_ind


def make_read(template, rng, indel_rate):
    # take a little more than L so indels do not shorten the read
    raw = (template + ADAPTER + "".join(rng.choice("ACGT") for _ in range(L)))[:L + 20]
    qs = quals(rng)
    seq, n_sub, n_ind = corrupt(raw, qs + [qs[-1]] * 20, indel_rate, rng)
    seq = seq[:L]
    return seq, "".join(chr(q + 33) for q in qs[:len(seq)]), n_sub, n_ind


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("truth")
    ap.add_argument("out")
    ap.add_argument("--indel-rate", type=float, default=2e-5)
    ap.add_argument("--swap", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--rabbits", type=int, default=3)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    mols = read_tsv(a.truth)
    os.makedirs(os.path.join(a.out, "samples"), exist_ok=True)
    per_rabbit = {r: [] for r in range(1, a.rabbits + 1)}
    for i, m in enumerate(mols):
        per_rabbit[i % a.rabbits + 1].append(m)
    truth_rows = []
    for rabbit, lst in per_rabbit.items():
        pairs = []
        for m in lst:
            primer = resolve_iupac(PRIMERS[m["chain"]][int(m["primer_idx"])], rng)
            insert = T7_TSO + "ACATGGG" + m["cdna"] + rc(primer)
            pairs.append((m["mol_id"], insert, "mol"))
        n_junk = max(1, len(lst) // 25)
        for j in range(n_junk):
            pairs.append((f"J{rabbit}{j:04d}", "".join(rng.choice("ACGT") for _ in range(450)), "junk"))
        for j in range(max(1, len(lst) // 100)):
            p = resolve_iupac(rng.choice(sum(PRIMERS.values(), [])), rng)
            pairs.append((f"D{rabbit}{j:04d}", T7_TSO + "ACATGGG" + rc(p), "dimer"))
        rng.shuffle(pairs)
        f1 = open(os.path.join(a.out, "samples", f"15651-LT-{rabbit}_R1.fastq"), "w")
        f2 = open(os.path.join(a.out, "samples", f"15651-LT-{rabbit}_R2.fastq"), "w")
        for k, (mid, insert, kind) in enumerate(pairs):
            name = f"AV:1:FCSIM:1:{1101 + k // 5000}:{mid}:{k}"
            r1, q1, s1, i1 = make_read(insert, rng, a.indel_rate)
            r2, q2, s2, i2 = make_read(rc(insert), rng, a.indel_rate)
            swapped = rng.random() < a.swap
            if swapped:
                r1, q1, r2, q2 = r2, q2, r1, q1
            f1.write(f"@{name} 1:N:0:ACGTACGT\n{r1}\n+\n{q1}\n")
            f2.write(f"@{name} 2:N:0:ACGTACGT\n{r2}\n+\n{q2}\n")
            truth_rows.append({"read_id": name, "rabbit": rabbit, "mol_id": mid, "kind": kind, "swapped": int(swapped),
                               "tso_side_subs": s2 if swapped else s1, "tso_side_indels": i2 if swapped else i1,
                               "primer_side_subs": s1 if swapped else s2, "primer_side_indels": i1 if swapped else i2})
        f1.close()
        f2.close()
        print(f"rabbit {rabbit}: {len(pairs)} pairs", flush=True)
    write_tsv(os.path.join(a.out, "read_truth.tsv"), truth_rows, list(truth_rows[0].keys()))


if __name__ == "__main__":
    main()
