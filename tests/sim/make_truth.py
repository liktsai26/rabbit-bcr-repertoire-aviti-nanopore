"""Build the truth molecule set for the end-to-end simulation.

IGH : real rabbit IgH V(D)J regions from OAS (Banerjee 2017, SRR5060324), re-annotated with
      IgBLAST 1.22.0 against the same IMGT reference the pipeline uses; only rows whose V
      alignment starts at germline position 1 (full FR1) are used.
IGK/IGL : synthetic rearrangements of IMGT F alleles with SHM-like substitutions,
      gene-conversion-like tracts copied from other V alleles (aligned by IMGT position),
      occasional codon indels in CDR1/CDR2, V/J trimming and short N additions.
Every molecule is the sense cDNA  5'UTR + leader + V(D)J + constant-region spacer
(no TSO, no primer). Truth calls = IgBLAST on that cDNA. Only productive molecules are kept.

Usage: python make_truth.py OUT_DIR [--igh N] [--igk N] [--igl N] [--seed S]
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (C_SPACER, LEADERS, PRIMER_CHOICE, WORK, read_airr, run_igblast,  # noqa: E402
                    write_tsv)

IMGT_DIR = os.path.join(WORK, "downloads", "imgt")
OAS = os.path.join(WORK, "downloads", "oas", "SRR5060324_Heavy_Bulk.csv.gz")
CDR_IMGT = [(27, 38), (56, 65)]


def read_imgt_gapped(path):
    recs, h, s = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if h:
                recs[h[0]] = (h[1], "".join(s).upper())
            f = line[1:].split("|")
            h, s = (f[1], f[3]), []
        elif line:
            s.append(line)
    if h:
        recs[h[0]] = (h[1], "".join(s).upper())
    return recs


def utr(rng):
    n = rng.randint(20, 90)
    return "".join(rng.choice("ACGT") for _ in range(n))


def build_igh(n, rng, tmp):
    import pandas as pd
    df = pd.read_csv(OAS, skiprows=1, usecols=["sequence", "v_call"], low_memory=False)
    df = df[~df["sequence"].str.contains("N")]
    df["g"] = df["v_call"].str.split(",").str[0].str.split("*").str[0]
    # keep every rare-gene sequence (up to 60 per gene) plus a random sample of the rest,
    # so the near-identical IGHV1S40/1S45 majority and the rarer genes are both represented
    parts = [grp.sample(min(len(grp), 60), random_state=rng.randint(0, 10**6)) for _, grp in df.groupby("g")]
    rare = pd.concat(parts)
    rest = df.drop(rare.index).sample(min(len(df) - len(rare), int(n * 1.6)), random_state=rng.randint(0, 10**6))
    pool = pd.concat([rare, rest])
    fa = os.path.join(tmp, "igh_pool.fasta")
    with open(fa, "w") as f:
        for i, s in enumerate(pool["sequence"]):
            f.write(f">p{i}\n{s}\n")
    out = os.path.join(tmp, "igh_pool.airr.tsv")
    run_igblast(fa, out, "IGH", threads=6)
    mols = []
    for r in read_airr(out):
        if r["productive"] != "T" or r["v_germline_start"] != "1" or not r["j_call"] or not r["j_sequence_end"]:
            continue
        seq = r["sequence"]  # germline orientation
        vdj = seq[int(r["v_sequence_start"]) - 1:int(r["j_sequence_end"])]
        mols.append(vdj)
    rng.shuffle(mols)
    return mols[:int(n * 1.15)]


def mutate_light(vgap, donors, rng):
    """vgap: IMGT-gapped V; returns ungapped mutated V (gene conversion + SHM + rare codon indel)."""
    seq = list(vgap)
    if donors and rng.random() < 0.6:  # gene-conversion-like tract
        d = rng.choice(donors)
        start = rng.randint(1, 100)
        length = rng.randint(8, 40)
        for p in range(start, min(start + length, len(seq), len(d))):
            if seq[p] != "." and d[p] != ".":
                seq[p] = d[p]
    for p in range(len(seq)):
        if seq[p] == ".":
            continue
        imgt = p + 1
        in_cdr = any(a <= imgt <= b for a, b in CDR_IMGT)
        if rng.random() < (0.08 if in_cdr else 0.03):
            seq[p] = rng.choice([b for b in "ACGT" if b != seq[p]])
    v = "".join(c for c in seq if c != ".")
    if rng.random() < 0.08:  # codon indel in CDR1/CDR2 region (approx ungapped 75-200)
        pos = rng.randrange(75, min(200, len(v) - 3), 3)
        if rng.random() < 0.5:
            v = v[:pos] + v[pos + 3:]
        else:
            v = v[:pos] + "".join(rng.choice("ACGT") for _ in range(3)) + v[pos:]
    return v


def build_light(chain, n, rng, tmp):
    vrec = read_imgt_gapped(os.path.join(IMGT_DIR, f"{chain}V.fasta"))
    jrec = read_imgt_gapped(os.path.join(IMGT_DIR, f"{chain}J.fasta"))
    vf = [k for k, (func, _) in vrec.items() if func == "F"]
    jf = [k for k, (func, _) in jrec.items() if func in ("F", "ORF")]
    allv = [s for _, s in vrec.values()]
    cands = []
    for _ in range(int(n * 3)):
        va = rng.choice(vf)
        v = mutate_light(vrec[va][1], [d for d in allv if d != vrec[va][1]], rng)
        v = v[:len(v) - rng.randint(0, 4)]
        j = jrec[rng.choice(jf)][1].replace(".", "")
        j = j[rng.randint(0, 4):]
        nadd = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 3)))
        cands.append(v + nadd + j)
    fa = os.path.join(tmp, f"{chain}_pool.fasta")
    with open(fa, "w") as f:
        for i, s in enumerate(cands):
            f.write(f">p{i}\n{LEADERS[chain]}{s}\n")
    out = os.path.join(tmp, f"{chain}_pool.airr.tsv")
    run_igblast(fa, out, chain, threads=6)
    ok = [cands[int(r["sequence_id"][1:])] for r in read_airr(out) if r["productive"] == "T"]
    rng.shuffle(ok)
    return ok[:int(n * 1.15)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--igh", type=int, default=4500)
    ap.add_argument("--igk", type=int, default=3000)
    ap.add_argument("--igl", type=int, default=1800)
    ap.add_argument("--seed", type=int, default=20261022)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(a.out, exist_ok=True)
    tmp = os.path.join(a.out, "tmp")
    os.makedirs(tmp, exist_ok=True)
    vdjs = {"IGH": build_igh(a.igh, rng, tmp), "IGK": build_light("IGK", a.igk, rng, tmp),
            "IGL": build_light("IGL", a.igl, rng, tmp)}
    # assemble cDNA and annotate truth with IgBLAST on the whole cDNA
    rows = []
    for chain, lst in vdjs.items():
        fa = os.path.join(tmp, f"{chain}_cdna.fasta")
        recs = []
        for i, vdj in enumerate(lst):
            mol_id = f"{chain[2]}{i:05d}"
            u, lead, sp = utr(rng), LEADERS[chain], rng.choice(C_SPACER[chain])
            cdna = u + lead + vdj + sp
            recs.append({"mol_id": mol_id, "chain": chain, "cdna": cdna, "vdj_start": len(u) + len(lead),
                         "vdj_end": len(u) + len(lead) + len(vdj), "primer_idx": rng.choice(PRIMER_CHOICE[chain])})
        with open(fa, "w") as f:
            for r in recs:
                f.write(f">{r['mol_id']}\n{r['cdna']}\n")
        out = os.path.join(tmp, f"{chain}_cdna.airr.tsv")
        run_igblast(fa, out, chain, threads=6)
        ann = {r["sequence_id"]: r for r in read_airr(out)}
        kept = 0
        target = {"IGH": a.igh, "IGK": a.igk, "IGL": a.igl}[chain]
        for r in recs:
            t = ann.get(r["mol_id"])
            if not t or t["productive"] != "T" or t["rev_comp"] != "F" or kept >= target:
                continue
            kept += 1
            for k in ("v_call", "d_call", "j_call", "junction", "junction_aa", "productive", "v_germline_start",
                      "v_germline_end", "v_sequence_start", "v_sequence_end", "v_cigar", "v_identity"):
                r["truth_" + k] = t[k]
            rows.append(r)
        print(f"{chain}: {kept} truth molecules", flush=True)
    fields = ["mol_id", "chain", "cdna", "vdj_start", "vdj_end", "primer_idx"] + [
        "truth_" + k for k in ("v_call", "d_call", "j_call", "junction", "junction_aa", "productive",
                               "v_germline_start", "v_germline_end", "v_sequence_start", "v_sequence_end",
                               "v_cigar", "v_identity")]
    write_tsv(os.path.join(a.out, "truth_molecules.tsv"), rows, fields)


if __name__ == "__main__":
    main()
