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

All analysis was run in a single conda environment (`bcr_pipeline`) to keep
versions consistent across every step. See `docs/` for the full version
record:

- `docs/Software_and_Computational_Environment.docx` — Methods-ready paragraph + Table S1 (English)
- `environment.yml` — machine-readable export (`conda env export --no-builds`), recreate with:
  ```
  conda env create -f environment.yml
  ```

Key tool versions used by the scripts in this folder (verified 2026-09-08,
unchanged since first recorded 2026-08-03):
Python 3.10.20 · cutadapt 5.2 · NanoFilt 2.8.0 · igblastn 1.22.0 (BLAST+ 2.16.0) ·
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

The four upstream shell scripts (`02_run_nanopore.sh`, `02_run_nanopore_qywfsm.sh`,
`02a_run_aviti_qc.sh`, `02b_run_aviti_igblast.sh`) resolve three more
machine-specific paths the same way — override them if your setup differs:

| Variable | Default | What it points to |
|---|---|---|
| `BCR_PIPELINE_ROOT` | `~/Desktop/Claude_code/BCR_pipeline` | pipeline root (same variable as above) |
| `BCR_PIPELINE_PYTHON` | `.../envs/bcr_pipeline/bin/python` | the `bcr_pipeline` conda env's Python |
| `BCR_IGDATA` | `.../envs/bcr_pipeline/share/igblast` | IgBLAST data dir (database/, optional_file/) |

```bash
export BCR_PIPELINE_ROOT=/path/to/your/BCR_pipeline
export BCR_PIPELINE_PYTHON=/path/to/your/conda/envs/bcr_pipeline/bin/python
export BCR_IGDATA=/path/to/your/conda/envs/bcr_pipeline/share/igblast
bash 02a_run_aviti_qc.sh
```

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

## Scripts

### Upstream: QC, trimming, IgBLAST (shell + per-sample report)

| Script | Purpose |
|---|---|
| `02_run_nanopore.sh` | Nanopore QC/trim + IgBLAST pipeline (rabbits 1–2, MZPD8T run) |
| `02_run_nanopore_qywfsm.sh` | Nanopore QC/trim + IgBLAST pipeline (rabbits 3–5, QYWFSM run) |
| `02a_run_aviti_qc.sh` | AVITI short-read QC/adapter trimming (cutadapt, `-e 0.1`) |
| `02b_run_aviti_igblast.sh` | AVITI IgBLAST V(D)J assignment |
| `03_make_report.py` | Per-sample/per-locus Excel QC + repertoire report |

### Shared library

| Script | Purpose |
|---|---|
| `fig_utils.py` | Shared paths, constants, and the data cache (`fig_cache.pkl`, not included here) used by every figure/table script below |

### Main figures

| Script | Output |
|---|---|
| `fig1_flowchart.py` | `Figure1_Flowchart_v4.png` |
| `fig2_qc.py` | `Figure2_QC_v4.png` |
| `fig3_igh.py` | `Figure3_IgH_v4.png` |
| `fig4_igkl.py` | `Figure4_IgKL_v4.png` |
| `fig5_cdr3.py` | `Figure5_CDR3_v4.png` |
| `fig6_diversity.py` | `Figure6_Diversity_v4.png` |
| `fig7_consistency.py` | `Figure7_Consistency_v4.png` |

Script numbering now matches the figure numbering directly (renamed
2026-09-08 from the working-copy names `fig_flowchart_v3.py`,
`fig1_qc_v3.py`, ..., `fig6_consistency_V07162026.py`, which were off by
one against the manuscript figure numbers).

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
above: `Figure1_Flowchart_v4.png` through `Figure7_Consistency_v4.png`, and
`SuppFigure1`–`SuppFigure4`.
