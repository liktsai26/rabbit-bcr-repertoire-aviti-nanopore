"""Unit tests for scripts/stitch_core.py on hand-built AIRR-like rows (no IgBLAST needed).

Germline allele A is 120 nt; allele B is A with 3 nt deleted at position 31-33, the same kind of
length difference as IGHV1S45*01 (299 nt) vs IGHV1S40*01 (296 nt). The molecule follows A.
R1 is aligned to A and R2 (which cannot see position 31-33) to B.
"""
import importlib
import os
import random
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from conftest import SCRIPTS  # noqa: E402

sys.path.insert(0, SCRIPTS)
stitch_core = pytest.importorskip("stitch_core")

rng = random.Random(5)
A = "".join(rng.choice("ACGT") for _ in range(120))
B = A[:30] + A[33:]
MOL = A


@pytest.fixture
def mapper(tmp_path):
    fa = tmp_path / "V.fasta"
    fa.write_text(f">ALLELE_A*01\n{A}\n>ALLELE_B*01\n{B}\n")
    return stitch_core.AlleleMapper(str(fa))


def row(seq, call, cigar, vs, ve, gs, ge):
    return dict(sequence=seq, v_call=call, v_cigar=cigar, v_sequence_start=vs, v_sequence_end=ve,
                v_germline_start=gs, v_germline_end=ge)


def r1_row(end=80):
    return row(MOL[:end], "ALLELE_A*01", f"{end}M", 1, end, 1, end)


def r2_row(mol_start):
    """R2 covering molecule positions mol_start..120 (1-based), aligned to allele B."""
    seq = MOL[mol_start - 1:]
    b_start = mol_start - 3          # every position after 33 is 3 smaller on B
    return row(seq, "ALLELE_B*01", f"{b_start - 1}N{len(seq)}M", 1, len(seq), b_start, b_start + len(seq) - 1)


def test_germline_to_query_handles_clips_skips_and_indels():
    m = stitch_core.germline_to_query("5S2N10M3D4M2I6M4S", 6, 27, 3, 25)
    assert m[3] == 6 and m[12] == 15 and m[13] is None and m[15] is None
    assert m[16] == 16 and m[19] == 19 and m[20] == 22 and m[25] == 27
    assert stitch_core.germline_to_query("10M", 1, 11, 1, 10) is None   # coordinates disagree


def test_different_alleles_with_length_difference_stitch_exactly(mapper):
    """K8: the original subtracts coordinates of two different alleles -> 3-nt deletion here."""
    s, reason, overlap, diff = stitch_core.stitch_pair(r1_row(), r2_row(60), mapper)
    assert diff and reason.startswith("ok")
    assert s == MOL


def test_adjacent_reads_are_kept(mapper):
    """K6: R2 starting right after R1's last base is 'no gap' in 1-based closed coordinates."""
    s, reason, overlap, _ = stitch_core.stitch_pair(r1_row(), r2_row(81), mapper)
    assert overlap == 0 and s == MOL


def test_true_gap_is_rejected(mapper):
    s, reason, overlap, _ = stitch_core.stitch_pair(r1_row(), r2_row(90), mapper)
    assert s is None and reason == "gap" and overlap < 0


def test_unrelated_overlap_fails_qc(mapper):
    """K9: the original local-alignment QC passes overlaps that do not match at the joined register."""
    r2 = r2_row(50)
    junk = "".join(rng.choice("ACGT") for _ in range(31))
    r2["sequence"] = junk + r2["sequence"][31:]
    s, reason, overlap, _ = stitch_core.stitch_pair(r1_row(), r2, mapper)
    assert s is None and reason == "qc"


def test_short_positive_overlap_below_min_is_not_duplicated(mapper, monkeypatch):
    """Found while re-testing: with MIN_OVERLAP>0 and gap rescue on, a too-short positive overlap
    must not fall through to the 'adjacent' join (which would duplicate the overlapping bases)."""
    monkeypatch.setattr(stitch_core, "MIN_OVERLAP", 20)
    monkeypatch.setattr(stitch_core, "RESCUE_GAP", True)
    r1 = r1_row()
    r1["sequence"] = MOL[:80] + "".join(rng.choice("ACGT") for _ in range(20))   # unverifiable tail
    s, reason, overlap, _ = stitch_core.stitch_pair(r1, r2_row(76), mapper)
    assert 0 < overlap < 20
    assert s is None or s == MOL


@pytest.mark.parametrize("prod,stop,vj,expected", [
    ("T", "F", "T", True), ("F", "F", "T", False), ("", "F", "T", True), ("", "T", "T", False),
])
def test_productive_mask_is_per_row(prod, stop, vj, expected):
    """GPT-reported: the fallback rule used to switch per chunk depending on other rows."""
    import pandas as pd
    df = pd.DataFrame({"productive": [prod, "T"], "stop_codon": [stop, "F"], "vj_in_frame": [vj, "T"]})
    assert bool(stitch_core.productive_mask(df).iloc[0]) is expected


def test_module_reloads_with_env_overrides(monkeypatch):
    monkeypatch.setenv("BCR_STITCH_OVERLAP_FROM", "r1")
    mod = importlib.reload(stitch_core)
    assert mod.OVERLAP_FROM == "r1"
    monkeypatch.delenv("BCR_STITCH_OVERLAP_FROM")
    importlib.reload(stitch_core)
