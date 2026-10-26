"""Regression tests: one test per confirmed bug in the upstream scripts.

Run against this repo's scripts (default) or the unmodified upstream scripts with
BCR_TEST_SCRIPTS=/path/to/upstream/scripts to see each bug reproduce.
All shell scripts run under /bin/bash (3.2 on macOS, as on the authors' machine).
"""
import glob
import os
import subprocess
import sys


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from conftest import (ENV_PREFIX, REPO, SCRIPTS, SIM_SAMPLES, copy_fastq_subset,  # noqa: E402
                      needs_sim, needs_tools)

sys.path.insert(0, os.path.join(REPO, "tests", "sim"))
from common import PRIMERS, TSO, rc, read_airr, read_tsv  # noqa: E402

TRUTH = os.path.join(REPO, "work", "sim", "truth", "truth_molecules.tsv")
pytestmark = [needs_tools, needs_sim]


def truth_seqs(chain, n):
    rows = [r for r in read_tsv(TRUTH) if r["chain"] == chain][:n]
    return [(f"{chain}_{i}", r["cdna"]) for i, r in enumerate(rows)]


def write_fasta(path, recs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for name, s in recs:
            f.write(f">{name}\n{s}\n")


def mark_upstream_done(chunk_dir):
    """The fixed IgBLAST stages require the upstream stage's completion marker."""
    os.makedirs(chunk_dir, exist_ok=True)
    out_dir = os.path.dirname(chunk_dir)
    name = "stitch.complete" if out_dir.endswith("_Stitched") else "demux.complete"
    open(os.path.join(out_dir, name), "w").write("test\n")


def rows(path):
    return len(read_airr(path)) if os.path.exists(path) else -1


def files_mentioning(root, text):
    hits = []
    for p in glob.glob(os.path.join(root, "**", "*"), recursive=True):
        if os.path.isfile(p) and not p.endswith((".fasta", ".airr.tsv")):
            try:
                if text in open(p, errors="ignore").read():
                    hits.append(p)
            except OSError:
                pass
    return hits


# ---------------------------------------------------------------- IgBLAST stages

def test_02d_empty_chunk_does_not_abort(sandbox):
    """K2: an empty chunk FASTA made 02d exit silently (grep -c returns 1 under set -e)."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), [])
    write_fasta(os.path.join(d, "chunk_001.fasta"), truth_seqs("IGH", 3))
    p = sandbox.run("02d_run_aviti_r2only_igblast.sh", 1)
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv")) == 3


def test_02h_empty_chunk_does_not_abort(sandbox):
    """K1: 02h's `grep -c || echo 0` yields "0\\n0", the skip branch never runs, head -1 "" fails."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_Stitched", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_0000.fasta"), truth_seqs("IGH", 3))
    write_fasta(os.path.join(d, "chunk_0001.fasta"), [])
    p = sandbox.run("02h_stitch_igblast.sh", 1)
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(sandbox.path("15651-LT-1_Aviti24", "IGH_Stitched", "15651-LT-1_AVITI_IGH_Stitched.airr.tsv")) == 3


def test_04b_empty_chunk_does_not_abort(sandbox):
    d = sandbox.path("15651-LT-1_Aviti24", "IGK_R1only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGK", 3))
    write_fasta(os.path.join(d, "chunk_001.fasta"), [])
    p = sandbox.run("04b_run_aviti_igkl_igblast.sh", 1, "IGK", "r1only")
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(sandbox.path("15651-LT-1_Aviti24", "IGK_R1only", "15651-LT-1_AVITI_IGK_R1only.airr.tsv")) == 3


def test_partial_final_file_is_not_treated_as_complete(sandbox):
    """K3: an interrupted final merge leaves a non-empty partial FINAL_AIRR; `[ -s ]` accepted it."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    seqs = truth_seqs("IGH", 6)
    write_fasta(os.path.join(d, "chunk_000.fasta"), seqs[:3])
    write_fasta(os.path.join(d, "chunk_001.fasta"), seqs[3:])
    final = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv")
    assert sandbox.run("02d_run_aviti_r2only_igblast.sh", 1).returncode == 0
    assert rows(final) == 6
    lines = open(final).read().splitlines(True)
    with open(final, "w") as f:          # what a kill during the merge loop leaves behind
        f.writelines(lines[:2])
    for marker in glob.glob(final + ".complete"):
        os.remove(marker)
    p = sandbox.run("02d_run_aviti_r2only_igblast.sh", 1)
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(final) == 6


def test_split_resume_does_not_drop_sequences(sandbox):
    """K16: a sub-batch directory left by an interrupted split was reused as if complete."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    seqs = truth_seqs("IGH", 6)
    write_fasta(os.path.join(d, "chunk_000.fasta"), seqs)
    write_fasta(os.path.join(d, "chunk_000_sub", "sub_0000.fasta"), seqs[:2])  # split killed half-way
    p = sandbox.run("02d_run_aviti_r2only_igblast.sh", 1)
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv")) == 6


def test_nonzero_igblastn_exit_is_retried(fake_sandbox, tmp_path):
    """K15: igblastn exiting non-zero (e.g. killed under memory pressure) aborted 02d under set -e."""
    d = fake_sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGH", 4))
    p = fake_sandbox.run("02d_run_aviti_r2only_igblast.sh", 1,
                         extra_env={"FAKE_FAIL_ONCE": str(tmp_path / "failed_once")})
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(fake_sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv")) == 4


def test_irrecoverable_ids_survive_in_04b(fake_sandbox):
    """K4: igblast_round_retry.py wrote the unprocessable read IDs into SUB_DIR, which 04b deletes."""
    d = fake_sandbox.path("15651-LT-1_Aviti24", "IGK_R1only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGK", 4))
    p = fake_sandbox.run("04b_run_aviti_igkl_igblast.sh", 1, "IGK", "r1only", extra_env={"FAKE_DROP_IDS": "IGK_2"})
    assert p.returncode == 0, p.stdout[-2000:]
    out_dir = fake_sandbox.path("15651-LT-1_Aviti24", "IGK_R1only")
    assert rows(os.path.join(out_dir, "15651-LT-1_AVITI_IGK_R1only.airr.tsv")) == 3
    assert files_mentioning(out_dir, "IGK_2"), "irrecoverable read ID was not kept anywhere"


def test_irrecoverable_ids_recorded_by_bash_retry(fake_sandbox):
    """Bash halving retry (02d/02f/02h): the dropped read was only printed to the console."""
    d = fake_sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGH", 4))
    p = fake_sandbox.run("02d_run_aviti_r2only_igblast.sh", 1, extra_env={"FAKE_DROP_IDS": "IGH_1"})
    assert p.returncode == 0, p.stdout[-2000:]
    out_dir = fake_sandbox.path("15651-LT-1_Aviti24", "IGH_R2only")
    assert rows(os.path.join(out_dir, "15651-LT-1_AVITI_IGH_R2only.airr.tsv")) == 3
    assert files_mentioning(out_dir, "IGH_1")


# ---------------------------------------------------------------- retry engine

def _patched_engine(tmp_path):
    """The upstream engine hardcodes /Users/newuser paths; point it at the test env."""
    src = open(os.path.join(SCRIPTS, "igblast_round_retry.py")).read()
    src = src.replace("/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env", ENV_PREFIX)
    p = tmp_path / "igblast_round_retry.py"
    p.write_text(src)
    return str(p)


def test_round_retry_output_order_is_deterministic(tmp_path):
    """K14: rows were emitted in set-iteration order, which changes with PYTHONHASHSEED."""
    fasta = str(tmp_path / "in.fasta")
    write_fasta(fasta, truth_seqs("IGK", 30))
    eng = _patched_engine(tmp_path)
    outs = []
    for seed in (1, 2):
        out = str(tmp_path / f"o{seed}.airr.tsv")
        ref = os.path.join(REPO, "work", "reference", "IMGT_reference_clean")
        env = dict(os.environ, PYTHONHASHSEED=str(seed), BCR_ENV_PREFIX=ENV_PREFIX)
        subprocess.run([os.path.join(ENV_PREFIX, "bin", "python"), eng, fasta, out, f"{ref}/IGKV.fasta",
                        f"{ref}/IGKJ.fasta", f"{ref}/IGKJ.fasta"], check=True, capture_output=True, env=env)
        outs.append([r["sequence_id"] for r in read_airr(out)])
    assert outs[0] == outs[1] == [f"IGK_{i}" for i in range(30)]


def test_round_retry_empty_input(tmp_path):
    """K14: an empty input crashed with TypeError (header None)."""
    fasta = str(tmp_path / "empty.fasta")
    open(fasta, "w").close()
    eng = _patched_engine(tmp_path)
    ref = os.path.join(REPO, "work", "reference", "IMGT_reference_clean")
    p = subprocess.run([os.path.join(ENV_PREFIX, "bin", "python"), eng, fasta, str(tmp_path / "o.airr.tsv"),
                        f"{ref}/IGKV.fasta", f"{ref}/IGKJ.fasta", f"{ref}/IGKJ.fasta"],
                       capture_output=True, text=True, env=dict(os.environ, BCR_ENV_PREFIX=ENV_PREFIX))
    assert p.returncode == 0, p.stderr[-1500:]


# ---------------------------------------------------------------- demux

def _demux_inputs(sandbox, n=60):
    copy_fastq_subset(SIM_SAMPLES, sandbox.samples, 1, n)


def test_demux_rerun_does_not_duplicate_stats(sandbox):
    """chunk_stats.csv was appended again for a chunk that is re-processed after an interruption."""
    _demux_inputs(sandbox, 2500)   # 3 chunks with the shrunken chunk size
    assert sandbox.run("02c_run_aviti_r2only_igh.sh", 1).returncode == 0
    os.remove(sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks", "chunk_001.done"))
    os.remove(sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "demux.complete"))   # an interrupted run has none
    assert sandbox.run("02c_run_aviti_r2only_igh.sh", 1).returncode == 0
    with open(sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunk_stats.csv")) as f:
        idx = [line.split(",")[0] for line in f.read().splitlines()[1:]]
    assert sorted(idx) == ["000", "001", "002"], idx


def test_demux_detects_changed_input(sandbox):
    """GPT-reported: stale .done markers are trusted after the FASTQ is replaced."""
    _demux_inputs(sandbox, 60)
    assert sandbox.run("02c_run_aviti_r2only_igh.sh", 1).returncode == 0
    for side in ("R1", "R2"):   # replace with a different set of reads
        src = os.path.join(SIM_SAMPLES, f"15651-LT-2_{side}.fastq")
        with open(src) as f, open(os.path.join(sandbox.samples, f"15651-LT-1_{side}.fastq"), "w") as g:
            g.writelines(f.readlines()[:240])
    p = sandbox.run("02c_run_aviti_r2only_igh.sh", 1)
    assert p.returncode != 0 and "跟上次不同" in p.stdout


def test_demux_with_no_matching_reads_exits_zero(sandbox):
    """K18: zero demuxed reads made the final `grep -c` return 1 under set -e (exit 1 after success)."""
    for side in ("R1", "R2"):
        with open(os.path.join(sandbox.samples, f"15651-LT-1_{side}.fastq"), "w") as g:
            for i in range(20):
                g.write(f"@r{i} {side[1]}:N:0:A\n{'ACGT' * 60}\n+\n{'I' * 240}\n")
    p = sandbox.run("02c_run_aviti_r2only_igh.sh", 1)
    assert p.returncode == 0, p.stdout[-1500:]


def test_sed_stops_reading_after_chunk(sandbox):
    """K17: each chunk re-read the whole FASTQ; the fixed sed quits at the chunk end."""
    txt = open(os.path.join(SCRIPTS, "02c_run_aviti_r2only_igh.sh"), encoding="utf-8").read()
    assert '${END_LINE}q' in txt


# ---------------------------------------------------------------- IGL entry rerun

def test_run_igl_rerun_skips_finished_rabbits(sandbox):
    """K5: after the IGL wrapper deleted chunks/, a rerun re-demuxed and re-stitched finished rabbits."""
    for r in (1, 2, 3):
        copy_fastq_subset(SIM_SAMPLES, sandbox.samples, r, 120)
    p = sandbox.run("run_igl_full_pipeline.sh")
    assert p.returncode == 0, p.stdout[-2000:]
    stats = sandbox.path("15651-LT-1_Aviti24", "IGL_R1only", "chunk_stats.csv")
    before = open(stats).read()
    p = sandbox.run("run_igl_full_pipeline.sh")
    assert p.returncode == 0, p.stdout[-2000:]
    assert "raw:" not in p.stdout, "demux was redone for a finished rabbit"
    assert open(stats).read() == before


# ---------------------------------------------------------------- stitching (scientific)

def _stitch_artifact_rate(sandbox, chain="IGH", n_pairs=400):
    copy_fastq_subset(SIM_SAMPLES, sandbox.samples, 1, n_pairs)
    for s in ("02c_run_aviti_r2only_igh.sh", "02d_run_aviti_r2only_igblast.sh",
              "02e_run_aviti_r1only_igh.sh", "02f_run_aviti_r1only_igblast.sh", "02g_stitch_qc_igh.py"):
        p = sandbox.run(s, 1)
        assert p.returncode == 0, (s, p.stdout[-2000:])
    truth = {r["mol_id"]: r for r in read_tsv(TRUTH)}
    reads = {r["read_id"]: r for r in read_tsv(os.path.join(os.path.dirname(SIM_SAMPLES), "read_truth.tsv"))}
    n = bad = 0
    for fa in glob.glob(sandbox.path("15651-LT-1_Aviti24", "IGH_Stitched", "chunks", "chunk_*.fasta")):
        name = None
        for line in open(fa):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                continue
            mid = name.split(":")[5]
            rt = reads[name]
            if mid[0] != "H" or int(rt["tso_side_indels"]) + int(rt["primer_side_indels"]):
                continue
            n += 1
            bad += len(line) != len("ACATGGG" + truth[mid]["cdna"])
    return n, bad


def test_stitched_length_matches_true_molecule(sandbox):
    """K8: R1/R2 assigned to different V alleles (e.g. IGHV1S40*01 296 nt vs IGHV1S45*01 299 nt)
    produced in-frame codon duplications/deletions at the stitch point (~60% of IgH reads)."""
    n, bad = _stitch_artifact_rate(sandbox)
    assert n > 150
    assert bad / n < 0.02, f"{bad}/{n} stitched reads have a length different from the true molecule"


# ---------------------------------------------------------------- Nanopore

def test_nanopore_trims_both_adapters(sandbox):
    """K10: with -g TSO and -a RC(primer), cutadapt removed only one adapter per read."""
    cdna = [r for r in read_tsv(TRUTH) if r["chain"] == "IGH"][0]["cdna"]
    sdir = os.path.join(sandbox.base, "samples", "MZPD8T_Nanopore")
    os.makedirs(sdir, exist_ok=True)
    read = "GGTAC" + TSO + "ACATGGG" + cdna + rc(PRIMERS["IGH"][1]) + "TTGCA"
    names = ["1_WT-IGHV", "2_WT-IGK", "3_WT-IGL", "4_Immu-IGHV", "5_Immu-IGK", "6_Immu-IGL"]
    for nm in names:
        with open(os.path.join(sdir, f"MZPD8T_{nm}.fastq"), "w") as f:
            if nm.endswith("IGHV"):
                for k, s in enumerate((read, rc(read))):
                    f.write(f"@n{k}\n{s}\n+\n{'5' * len(s)}\n")
    env = {"BCR_PIPELINE_ROOT": sandbox.base}
    p = sandbox.run("02_run_nanopore.sh", extra_env=env)
    trimmed = os.path.join(sandbox.base, "results", "MZPD8T_Nanopore", "results", "MZPD8T_WT_IGHV", "IGHV_trimmed.fasta")
    assert os.path.exists(trimmed), p.stdout[-2000:]
    seqs = [line.strip() for line in open(trimmed) if not line.startswith(">")]
    assert seqs and all(s in ("ACATGGG" + cdna, cdna) for s in seqs), seqs


# ---------------------------------------------------------------- issues reported by the second review round

def test_igblast_stage_refuses_unfinished_upstream(sandbox):
    """F1: 02d globbed whatever chunk FASTAs existed, so an interrupted 02c (or one whose error
    was hidden by `| tee`) produced a partial result that was later treated as complete."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGH", 3))   # no demux.complete
    p = sandbox.run("02d_run_aviti_r2only_igblast.sh", 1)
    assert p.returncode != 0
    assert not os.path.exists(sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv"))


def test_bash_retry_survives_zero_byte_leaf(fake_sandbox):
    """A1/F2: when the isolated failing read gives a 0-byte igblastn output (real 1.22.0 behaviour),
    the halving merge lost the AIRR header or dropped an extra good read."""
    d = fake_sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGH", 5))
    p = fake_sandbox.run("02d_run_aviti_r2only_igblast.sh", 1, extra_env={"FAKE_DROP_IDS": "IGH_0", "FAKE_ZERO_BYTE": "1"})
    final = fake_sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv")
    assert p.returncode == 0, p.stdout[-2000:]
    assert open(final).readline().startswith("sequence_id\t")
    assert sorted(r["sequence_id"] for r in read_airr(final)) == ["IGH_1", "IGH_2", "IGH_3", "IGH_4"]


def test_kill_during_subbatch_cleanup_does_not_truncate(sandbox):
    """F3: a kill while `rm -rf SUB_DIR` was running left sub .done files without their AIRR; the
    rerun skipped them and overwrote the complete chunk AIRR with a truncated one."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    seqs = truth_seqs("IGH", 320)                       # 2 sub-batches with the shrunken batch size (300)
    write_fasta(os.path.join(d, "chunk_000.fasta"), seqs)
    out_dir = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only")
    final = os.path.join(out_dir, "15651-LT-1_AVITI_IGH_R2only.airr.tsv")
    assert sandbox.run("02d_run_aviti_r2only_igblast.sh", 1).returncode == 0
    lines = open(os.path.join(d, "chunk_000.airr.tsv")).read().splitlines(True)
    # state left by the original scripts when killed half-way through `rm -rf chunk_000_sub`
    sub = os.path.join(d, "chunk_000_sub")
    write_fasta(os.path.join(sub, "sub_0000.fasta"), seqs[:300])
    write_fasta(os.path.join(sub, "sub_0001.fasta"), seqs[300:])
    for name in ("sub_0000.done", "sub_0001.done"):
        open(os.path.join(sub, name), "w").close()
    open(os.path.join(sub, "sub_0001.airr.tsv"), "w").writelines([lines[0]] + lines[301:])
    for p in (os.path.join(d, "chunk_000.igblast.done"), final, final + ".complete"):
        if os.path.exists(p):
            os.remove(p)
    p = sandbox.run("02d_run_aviti_r2only_igblast.sh", 1)
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(final) == 320


def test_04b_rejects_heavy_chain(sandbox):
    """F4: 04b/04d accepted CHAIN=IGH, used IGHJ as the D database and wrote to 02d's file names."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGH", 3))
    p = sandbox.run("04b_run_aviti_igkl_igblast.sh", 1, "IGH", "r2only")
    assert p.returncode != 0
    assert not os.path.exists(sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv"))


def test_space_in_pipeline_root(tmp_path):
    """GI7: BLAST splits an unquoted -germline_db_* value at spaces (e.g. 'Google Drive/...')."""
    from conftest import Sandbox
    sb = Sandbox(tmp_path, root_name="with space")
    d = sb.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGH", 3))
    p = sb.run("02d_run_aviti_r2only_igblast.sh", 1)
    assert p.returncode == 0, p.stdout[-2000:]
    assert rows(sb.path("15651-LT-1_Aviti24", "IGH_R2only", "15651-LT-1_AVITI_IGH_R2only.airr.tsv")) == 3


def _nanopore_sample(sandbox, n_reads, offset=0):
    """WT-IGHV gets n_reads reads starting at truth row `offset`; every other sample gets 3 reads of
    its own chain so that no sample is empty (an empty sample aborts the original script, A2)."""
    sdir = os.path.join(sandbox.base, "samples", "MZPD8T_Nanopore")
    os.makedirs(sdir, exist_ok=True)
    truth = read_tsv(TRUTH)
    for nm in ["1_WT-IGHV", "2_WT-IGK", "3_WT-IGL", "4_Immu-IGHV", "5_Immu-IGK", "6_Immu-IGL"]:
        chain = "IGH" if nm.endswith("IGHV") else nm[-3:]
        rows_ = [r for r in truth if r["chain"] == chain]
        rows_ = rows_[offset:offset + n_reads] if nm == "1_WT-IGHV" else rows_[:3]
        with open(os.path.join(sdir, f"MZPD8T_{nm}.fastq"), "w") as f:
            for k, r in enumerate(rows_):
                primer = PRIMERS[chain][int(r["primer_idx"])].replace("K", "G")
                s = TSO + "ACATGGG" + r["cdna"] + rc(primer)
                f.write(f"@n{offset + k}\n{s}\n+\n{'5' * len(s)}\n")


def test_nanopore_rerun_detects_swapped_input(sandbox):
    """A7: the per-sample 'already done' check ignored the input file, so a re-sequenced FASTQ
    placed under the same sample name was silently skipped."""
    env = {"BCR_PIPELINE_ROOT": sandbox.base}
    _nanopore_sample(sandbox, 3)
    assert sandbox.run("02_run_nanopore.sh", extra_env=env).returncode == 0
    _nanopore_sample(sandbox, 6, offset=10)
    sandbox.run("02_run_nanopore.sh", extra_env=env)
    airr = os.path.join(sandbox.base, "results", "MZPD8T_Nanopore", "results", "MZPD8T_WT_IGHV",
                        "MZPD8T_WT_IGHV_igblast.airr.tsv")
    assert rows(airr) == 6


def test_stage_does_not_read_stdin(sandbox):
    """Found while re-testing the fix: with nullglob on and no irrecoverable reads,
    `cat dir/*.txt` became a bare `cat` and blocked on stdin (a terminal) at the final merge."""
    d = sandbox.path("15651-LT-1_Aviti24", "IGH_R2only", "chunks")
    mark_upstream_done(d)
    write_fasta(os.path.join(d, "chunk_000.fasta"), truth_seqs("IGH", 3))
    p = subprocess.Popen(["/bin/bash", os.path.join(sandbox.sdir, "02d_run_aviti_r2only_igblast.sh"), "1"],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=sandbox.env)
    try:
        p.wait(timeout=180)
    except subprocess.TimeoutExpired:
        p.kill()
        raise AssertionError("stage blocked waiting on stdin")
    finally:
        p.stdin.close()
    assert p.returncode == 0
