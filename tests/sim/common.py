"""Shared helpers for the review test harness (simulation + evaluation).

Everything here is test-only code: it never touches the pipeline scripts under
scripts/ except by running them as subprocesses.
"""
import csv
import os
import subprocess
import sys

csv.field_size_limit(sys.maxsize)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
WORK = os.environ.get("BCR_TEST_WORK", os.path.join(REPO, "work"))
ENV_PREFIX = os.environ.get("BCR_TEST_ENV", os.path.join(WORK, "env"))
REF_DIR = os.path.join(WORK, "reference", "IMGT_reference_clean")
IGDATA = os.path.join(ENV_PREFIX, "share", "igblast")

# Library constants copied from the pipeline scripts (02c/04a/02_run_nanopore.sh).
T7_TSO = "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"
TSO = "AAGCAGTGGTATCAACGCAGAGT"
PRIMERS = {
    "IGH": ["CAGTGGGAAGACTGACGGAGCCTTAG", "CAGTGGGAAGACTGATGGAGCCTTAG"],
    "IGK": ["TGGTGGGAAGAKGAGGACAGTAGG", "TGGTGGGAAGAKGAGGACACTAGG", "TGGTGGGAAGAKGAGGACAGAAGG"],
    "IGL": ["CAAGGGGGCGACCACAGGCTGAC", "GTGAAGGAGTGACTACGGGTTGACC", "GAGGGGGTCACCGCGGGCTGAC"],
}
# Constant-region bases between the J/C boundary and the primer site (IMGT GENE-DB,
# split-codon N resolved to G): IGHG CH1[0:7], IGKC1/IGKC2 [0:15], IGLC [0:1].
C_SPACER = {
    "IGH": ["GGGCAAC"],
    "IGK": ["GGTGATCCAGTTGCA"],
    "IGL": ["G"],
}
# Primer index usable for each spacer/isotype in the simulation.
PRIMER_CHOICE = {"IGH": [0, 1], "IGK": [0, 2], "IGL": [0, 1, 2]}

LEADERS = {
    # rabbit VH leader METGLRWLLLVAVLKGVQC (seen in OAS Banerjee_2017 reads)
    "IGH": "ATGGAGACTGGGCTGCGCTGGCTTCTCCTGGTCGCTGTGCTCAAAGGTGTCCAGTGT",
    # rabbit kappa leader MDTRAPTQLLGLLLLWLPGARC
    "IGK": "ATGGACACGAGGGCCCCCACTCAGCTGCTGGGGCTCCTGCTGCTCTGGCTCCCAGGTGCCAGATGT",
    # lambda-like leader MAWTPLLLGLLSHCTGSLS
    "IGL": "ATGGCCTGGACCCCTCTCCTCCTCGGCCTCCTCTCTCACTGCACAGGTTCTCTCTCC",
}

_COMP = str.maketrans("ACGTNKMRYSWBVDHacgtnkmryswbvdh", "TGCANMKYRSWVBHDtgcanmkyrswvbhd")


def rc(seq):
    return seq.translate(_COMP)[::-1]


def resolve_iupac(seq, rng):
    choices = {"K": "GT", "M": "AC", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "N": "ACGT"}
    return "".join(rng.choice(choices[c]) if c in choices else c for c in seq)


def igblast_args(chain, threads=4):
    gl = {"IGH": ("IGHV", "IGHD", "IGHJ"), "IGK": ("IGKV", "IGKJ", "IGKJ"), "IGL": ("IGLV", "IGLJ", "IGLJ")}[chain]
    args = ["-organism", "rabbit", "-ig_seqtype", "Ig", "-outfmt", "19", "-num_threads", str(threads),
            "-extend_align5end", "-extend_align3end",
            "-auxiliary_data", os.path.join(IGDATA, "optional_file", "rabbit_gl.aux"),
            "-germline_db_V", os.path.join(REF_DIR, gl[0] + ".fasta"),
            "-germline_db_D", os.path.join(REF_DIR, gl[1] + ".fasta"),
            "-germline_db_J", os.path.join(REF_DIR, gl[2] + ".fasta")]
    if chain == "IGH":
        args += ["-min_D_match", "5"]
    return args


def run_igblast(fasta, out, chain, threads=4):
    env = dict(os.environ, IGDATA=IGDATA, BLASTDB=os.path.join(IGDATA, "database"))
    cmd = [os.path.join(ENV_PREFIX, "bin", "igblastn")] + igblast_args(chain, threads) + ["-query", fasta, "-out", out]
    subprocess.run(cmd, env=env, check=True)


def read_airr(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def gene(call):
    """'IGHV1S40*01,IGHV1S45*01' -> 'IGHV1S40' (first listed gene)."""
    if not call:
        return ""
    return call.split(",")[0].split("*")[0]


def allele(call):
    return call.split(",")[0] if call else ""
