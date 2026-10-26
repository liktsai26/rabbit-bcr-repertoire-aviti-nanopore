# Rabbit BCR Repertoire — Analysis Scripts (AVITI vs Nanopore, n=3)

This folder contains the analysis code, final figures, and environment
documentation for the rabbit IgG BCR repertoire comparison (AVITI short-read
vs Nanopore long-read sequencing, rabbits 1–3), prepared for manuscript
submission and collaborator review.

**Scope**: this folder contains code and small summary/figure files only.
It does **not** contain raw FASTQ files or IgBLAST AIRR TSV output (each can
be several GB per sample/locus) — those stay in the private working
directory and are excluded from version control. This code is for
**reviewing analysis logic and reproducing figures from intermediate
data**, not for a one-command rerun from raw reads without that
underlying data.

## Environment setup

All analysis was run in a single conda environment to keep versions
consistent across every step. See `docs/` for the full version record:

- `Software_and_Computational_Environment.docx` — Methods-ready paragraph + Table S1 (English)
  (at the repository root in this revision, not under `docs/`)
- `environment.yml` — machine-readable export (`conda env export --no-builds`), recreate with:
  ```
  conda env create -f environment.yml
  ```
  > Review note (2026-10): this export does not solve as-is on macOS (osx-arm64:
  > `sqlite=3.51.2` conflicts with `libsqlite>=3.52.0` required by `igblast=1.22.0`, and several
  > `defaults`-channel builds are rejected under strict channel priority; osx-64: no
  > `libffi=3.4.8`). It also lacks `chopper`, which the Nanopore scripts call. A minimal,
  > solvable environment with the key pinned versions is in `tests/environment.test.yml`.

Key tool versions used by the scripts in this folder (verified 2026-09-08,
unchanged since first recorded 2026-08-03):
Python 3.10.20 · cutadapt 5.2 · NanoFilt 2.8.0 (note: the Nanopore scripts actually call
`chopper`, version not recorded) · igblastn 1.22.0 (BLAST+ 2.16.0) ·
pandas 2.3.3 · numpy 2.2.6 · scipy 1.15.3 · matplotlib 3.10.9 · biopython 1.87 ·
openpyxl 3.1.5. (`environment.yml` and the `docs/` version records describe
the full analysis environment, including python-docx, which is used only by
the internal report-compilation script not included here.)

### Path configuration

`scripts/fig_utils.py` (and the figure/table scripts that each also define
their own `ROOT`) resolve the pipeline's root directory as:

```python
ROOT = os.environ.get("BCR_PIPELINE_ROOT", os.path.expanduser("~/Desktop/Claude_code/BCR_pipeline"))
```

To run these scripts on your own machine, set `BCR_PIPELINE_ROOT` to point
at your local copy of the full `BCR_pipeline` working directory (the one
containing `15651-LT-*_Aviti24/`, `MZPD8T_Nanopore/`, `QYWFSM_Nanopore/`,
`IMGT_reference_clean/`, and `figures/`), e.g.:

```bash
export BCR_PIPELINE_ROOT=/path/to/your/BCR_pipeline
python3 fig2_qc.py
```

The two Nanopore upstream scripts (`02_run_nanopore.sh`,
`02_run_nanopore_qywfsm.sh`) resolve three more machine-specific paths the
same way — override them if your setup differs:

| Variable | Default | What it points to |
|---|---|---|
| `BCR_PIPELINE_ROOT` | `~/Desktop/Claude_code/BCR_pipeline` | pipeline root (same variable as above) |
| `BCR_PIPELINE_PYTHON` | `.../envs/bcr_pipeline/bin/python` | the conda env's Python |
| `BCR_IGDATA` | `.../envs/bcr_pipeline/share/igblast` | IgBLAST data dir (database/, optional_file/) |

```bash
export BCR_PIPELINE_ROOT=/path/to/your/BCR_pipeline
export BCR_PIPELINE_PYTHON=/path/to/your/conda/envs/bcr_pipeline/bin/python
export BCR_IGDATA=/path/to/your/conda/envs/bcr_pipeline/share/igblast
bash 02_run_nanopore.sh
```

The AVITI upstream scripts (IGH/IGK/IGL pipelines, `02c`–`04d`,
`igblast_round_retry.py`, `run_igk_full_pipeline.sh`, `run_igl_full_pipeline.sh`)
accept the same kind of overrides since the 2026-10 review; the defaults are
still the original absolute paths, so behaviour on the original machine is unchanged:

| Variable | Default | What it points to |
|---|---|---|
| `BCR_PIPELINE_ROOT` | `/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline` | contains `samples/`, `reference/IMGT_reference_clean/`, `results/` |
| `BCR_ENV_PREFIX` | `.../miniforge3/envs/igblast-env` | conda env prefix (`bin/python`, `bin/igblastn`, `share/igblast`) |
| `BCR_PIPELINE_PYTHON`, `BCR_IGDATA`, `BCR_IGBLASTN` | derived from `BCR_ENV_PREFIX` | individual overrides |
| `BCR_KEEP_CHUNKS` | `0` | `1` keeps IGL `chunks/` instead of deleting them after each rabbit |
| `BCR_CHUNK_RAW_PAIRS`, `BCR_IGBLAST_BATCH`, `BCR_CHUNK_SIZE`, `BCR_BATCH_SIZE` | 1500000 / 25000 / 500000 / 25000 | chunk sizes (only changed by the tests) |

`04b`/`04d` now call the `igblast_round_retry.py` that sits next to them (the
original path `${BASE_DIR}/script/…` only existed on the original machine).

### Known limitation: gene order can vary between runs on tied values

`sort_j_by_aviti()` / `top_v_by_aviti()`-style helpers in `supp_fig2_igkl_full.py`
(and the same `set()`-based ranking pattern in `fig3_igh.py`, `fig4_igkl.py`,
`supp_fig1_igh_full.py`, `fig7_consistency.py`, `export_xlsx_v3.py`, `fig_utils.py`)
rank genes using a Python `set` before sorting by average frequency. When two
genes tie exactly, Python's string-hash randomization means their left/right
order on the x-axis can differ between separate script runs — this already
happened for `IGKJ1-1*03` vs `IGKJ1-5*01` in `SuppFigure2_IGKL_FullGeneUsage.png`. This
only swaps the position of tied categories; it does not change any underlying
value, count, or statistic. Not fixed in this release — noted here so it
isn't mistaken for a data error if you regenerate the figure and a gene pair
lands in a different order than the shipped PNG.

## 2026-10 code review changes (summary)

See `tests/README.md` for how each change was tested (public IMGT reference, real rabbit IgH
sequences from OAS, simulated AVITI/Nanopore reads with known truth, bash 3.2).

* **Stitching (`02g`/`04c` → shared `stitch_core.py`)**: germline coordinates of R1 and R2 are
  only comparable when both reads were assigned the same V allele. When they differ (most IgH
  pairs, e.g. IGHV1S40\*01 296 nt vs IGHV1S45\*01 299 nt) the cut point shifted by whole
  codons, producing in-frame duplications/deletions that stay `productive=T`. The fix maps
  coordinates between the two germline alleles, uses `v_cigar` for read indels, fixes the
  1-based adjacency off-by-one, and checks the overlap at the predicted register instead of with
  a shift-insensitive local alignment (exceptions now reject instead of pass).
* **IgBLAST stages (`02d`/`02f`/`02h`/`04b`/`04d` → shared `igblast_chunk_lib.sh`)**: empty
  chunks no longer abort; a non-zero igblastn exit is retried; interrupted sub-batch splits are
  redone; unprocessable read IDs are kept in `<OUT_DIR>/irrecoverable/` and
  `<sample>_irrecoverable_ids.txt`; final merges are written to a temp file, checked
  (`rows + irrecoverable = input`) and marked with `<final>.complete`, which replaces the
  "non-empty file means done" check.
* **Demux (`02c`/`02e`/`04a`)**: `sed` stops at the chunk end, an input fingerprint stops reuse
  of stale `.done` markers, `chunk_stats.csv` rows are replaced rather than appended on rerun,
  empty R1 inserts are dropped, the final count no longer fails on zero reads.
* **Entry points**: finished rabbits (all three `.complete` markers) are skipped, so a rerun after
  the IGL clean-up no longer re-demuxes and re-stitches them; clean-up requires the markers.
* **Nanopore**: `cutadapt -n 2` so both the TSO and the primer are removed when both are present
  (the default removes only one adapter per read).

### Using the fixed scripts on results produced by the old scripts

Script names, arguments, output file names and `qc_stats.json` keys are unchanged. The new
scripts, however, only trust results that carry a completion marker, so on an existing results
folder do this once:

```bash
# 1. keep the old stitched results for comparison
for d in "$BCR_PIPELINE_ROOT"/results/15651-LT-*_Aviti24/IG?_Stitched; do mv "$d" "${d}_old"; done
# 2. mark the existing demux / R1-only / R2-only results as complete (only if you trust them)
bash scripts/mark_existing_results_complete.sh --skip-stitched
# 3. re-run only stitching + the third IgBLAST pass
python scripts/02g_stitch_qc_igh.py 1 && bash scripts/02h_stitch_igblast.sh 1     # each rabbit
bash scripts/run_igk_full_pipeline.sh && bash scripts/run_igl_full_pipeline.sh     # skips demux/R1/R2
```

Without step 2 the new stages stop with a message instead of trusting unmarked old files
(`BCR_SKIP_UPSTREAM_CHECK=1` also bypasses the upstream check). To measure how much the stitching
fix changes the existing data before re-running anything, use
`tests/e2e/compare_stitch_methods.py R1only.airr.tsv R2only.airr.tsv IGHV.fasta`.

## Scripts

### Upstream: demuxing, QC, IgBLAST

AVITI (short-read) analysis uses a **coordinate-guided stitching** method:
R1 (TSO-anchored side) and R2 (primer-anchored side) are each IgBLAST'd
independently, then stitched into one full-length V(D)J sequence using each
read's own IgBLAST-reported coordinates (not sequence overlap), and the
stitched sequence is IgBLAST'd a third time for the final call. This
replaced an earlier FLASH2-merge-based method after R2-only calls were found
to systematically misassign near-identical V genes whose distinguishing
positions fall in FR1/CDR1 (outside R2's coverage) — not included in this
release.

**IGH** (no single entry-point script; run per rabbit, in order):

| Script | Purpose |
|---|---|
| `02c_run_aviti_r2only_igh.sh` | Demux + QC raw R1/R2 FASTQ → IGH "R2-role" reads (primer-anchored, covers CDR3) |
| `02d_run_aviti_r2only_igblast.sh` | IgBLAST the R2-only chunks (main CDR3/D/J source) |
| `02e_run_aviti_r1only_igh.sh` | Demux the same read pairs' "R1-role" reads (TSO-anchored, covers FR1–CDR1/FR2), IDs matched to `02c`'s output |
| `02f_run_aviti_r1only_igblast.sh` | IgBLAST the R1-only chunks (stop-codon check + coverage coordinates; not a CDR3/D/J source) |
| `02g_stitch_qc_igh.py` | Coordinate-guided stitching of R1+R2 IgBLAST output into full-length V(D)J; Biopython local-alignment QC on the overlap region |
| `02h_stitch_igblast.sh` | Third-pass IgBLAST on the stitched sequences — final production result |

```bash
for RABBIT in 1 2 3; do
    bash 02c_run_aviti_r2only_igh.sh ${RABBIT}
    bash 02d_run_aviti_r2only_igblast.sh ${RABBIT}
    bash 02e_run_aviti_r1only_igh.sh ${RABBIT}
    bash 02f_run_aviti_r1only_igblast.sh ${RABBIT}
    python 02g_stitch_qc_igh.py ${RABBIT}
    bash 02h_stitch_igblast.sh ${RABBIT}
done
```

**IGK / IGL** (same coordinate-stitching method, parameterized by chain;
each has a single full-pipeline entry point):

| Script | Purpose |
|---|---|
| `04a_run_aviti_igkl_demux.sh` | Demux raw FASTQ into both R1-role and R2-role chunks in one pass |
| `04b_run_aviti_igkl_igblast.sh` | IgBLAST the R1-only / R2-only chunks, using round-based capped retry (see below) |
| `04c_stitch_qc_igkl.py` | Coordinate-guided stitching + QC (same logic as `02g_stitch_qc_igh.py`, parameterized for IGK/IGL) |
| `04d_stitch_igblast_igkl.sh` | Third-pass IgBLAST on the stitched sequences — final production result |
| `igblast_round_retry.py` | Shared retry engine (called by `04b`/`04d`): works around a reproducible `igblastn` 1.22.0 `NULL pointer` crash tied to small `-word_size` on a subset of IGLV alleles, by retrying failed batches in capped-depth rounds rather than recursing on a single branch to exhaustion (which was found to risk killing the background process on long IGL runs) |
| `run_igk_full_pipeline.sh` | IGK entry point — runs `04a`→`04b`(r1only)→`04b`(r2only)→`04c`→`04d` for rabbits 1–3 in sequence |
| `run_igl_full_pipeline.sh` | IGL entry point — same as above, plus deletes each rabbit's intermediate chunk files once its final output is verified, to bound disk usage during the multi-day IGL run |

```bash
bash run_igk_full_pipeline.sh
bash run_igl_full_pipeline.sh
```

**Nanopore + per-sample report**:

| Script | Purpose |
|---|---|
| `02_run_nanopore.sh` | Nanopore QC/trim + IgBLAST pipeline (rabbits 1–2, MZPD8T run) — long reads, no stitching needed |
| `02_run_nanopore_qywfsm.sh` | Nanopore QC/trim + IgBLAST pipeline (rabbits 3–5, QYWFSM run) |
| `03_make_report.py` | Per-sample/per-locus Excel QC + repertoire report |

An earlier FLASH2-merge-based AVITI method (`02a_run_aviti_qc.sh` /
`02b_run_aviti_igblast.sh`) is **not included** in this release — it was
superseded by the coordinate-stitching pipeline above for all three loci
(IGH in 2026-09-28, IGK/IGL in 2026-10-15) and its intermediate output has
since been deleted from the working directory.

### Shared library

> Not included in this revision: `03_make_report.py`, `fig_utils.py` and all figure/table
> scripts listed below were removed from the repository on 2026-10-19 (commit c8f81f9), and
> `figures/` was removed as well. The tables are kept for reference.

| Script | Purpose |
|---|---|
| `fig_utils.py` | Shared paths, constants, and the data cache (`fig_cache.pkl`) used by every figure/table script below |

### Main figures

| Script | Output |
|---|---|
| `fig1_flowchart.py` | `Figure1_Flowchart_v6.png` |
| `fig2_qc.py` | `Figure2_QC_v4.png` |
| `fig3_igh.py` | `Figure3_IgH_v4.png` |
| `fig4_igkl.py` | `Figure4_IgKL_v4.png` |
| `fig5_cdr3.py` | `Figure5_CDR3_v4.png` |
| `fig6_diversity.py` | `Figure6_Diversity_v4.png` |
| `fig7_consistency.py` | `Figure7_Consistency_v4.png` |

### Supplementary figures

| Script | Output |
|---|---|
| `supp_fig1_igh_full.py` | `SuppFigure1_IGH_FullGeneUsage.png` |
| `supp_fig2_igkl_full.py` | `SuppFigure2_IGKL_FullGeneUsage.png` |
| `supp_fig3_raw_diversity.py` | `SuppFigure3_RawDiversity.png` |
| `supp_fig4_cost_efficiency.py` | `SuppFigure4_CostEfficiency.png` |

(A fifth supplementary figure, top-50 clonotype overlap, is tracked
separately and intentionally not included in this batch.)

### Tables / data export

| Script | Output |
|---|---|
| `export_xlsx_v3.py` | Underlying numeric data behind each main/supplementary figure (Excel) |

Two scripts from the private working directory are intentionally not
included here:
- A Word report-compilation script used only for internal manuscript
  drafting — has no bearing on the analysis logic.
- `supp_table1.py` (intended Supplementary Table 1 — per-sample
  sequencing/QC statistics) — currently non-functional against the present
  `fig_utils.py` cache schema and missing rabbit 3 from its sample list; it
  has never produced output under the current pipeline, so it is left out
  of this release rather than shared broken.

## Figures included in this folder

`figures/` contains the final n=3 (rabbits 1–3) versions matching the table
above: `Figure1_Flowchart_v6.png`, `Figure2_QC_v4.png` through
`Figure7_Consistency_v4.png`, and `SuppFigure1`–`SuppFigure4`.
