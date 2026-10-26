# Tests for the upstream pipeline

Everything here runs locally with public data only; no study FASTQ/AIRR files are needed.
All shell entry points are executed with `/bin/bash` (3.2 on macOS, the shell the scripts were
written for).

## 1. Environment and public reference data

```bash
conda env create -p work/env -f tests/environment.test.yml
bash tests/fetch_test_data.sh            # IMGT rabbit IG reference + constant genes + OAS rabbit IgH
```

`fetch_test_data.sh` downloads into `work/downloads/` (checksums in `work/downloads/SHA256SUMS`)
and builds IgBLAST databases in `work/reference/IMGT_reference_clean/` with IgBLAST's own
`edit_imgt_file.pl` + `makeblastdb -parse_seqids`.

## 2. Simulated data with known truth

```bash
work/env/bin/python tests/sim/make_truth.py work/sim/truth
work/env/bin/python tests/sim/simulate_aviti.py work/sim/truth/truth_molecules.tsv work/sim/aviti_main
work/env/bin/python tests/sim/simulate_aviti.py work/sim/truth/truth_molecules.tsv work/sim/aviti_indel --indel-rate 5e-4 --seed 8
work/env/bin/python tests/sim/simulate_nanopore.py work/sim/truth/truth_molecules.tsv work/sim/nano_main
```

* IgH molecules are **real rabbit IgH V(D)J sequences** (OAS, Banerjee et al. 2017, SRR5060324),
  so near-identical genes such as IGHV1S40/IGHV1S45, somatic mutation, gene conversion and
  real CDR3s are represented. IgK/IgL molecules are synthetic rearrangements of IMGT alleles
  with SHM, gene-conversion-like tracts and occasional codon indels.
* Reads follow the library layout the demux code expects (T7-TSO on R1, chain reverse primer
  at the start of R2, primers located on the real IMGT constant genes), with position-dependent
  substitution errors, rare indels, 30% swapped pairs, junk pairs and primer dimers.
* "Truth" calls are IgBLAST 1.22.0 on the clean molecule with the same parameters.

## 3. Regression tests (one per confirmed bug)

```bash
work/env/bin/python -m pytest tests/unit -q
# the same tests against the unmodified upstream scripts (all are expected to fail):
BCR_TEST_SCRIPTS=/path/to/upstream/scripts work/env/bin/python -m pytest tests/unit -q
```

`tests/conftest.py` builds a sandbox per test (scripts copied, hardcoded upstream paths
rewritten, chunk sizes shrunk so multi-chunk/resume paths run on small inputs) and provides a
fake `igblastn` wrapper that reproduces the igblastn 1.22.0 silent truncation and non-zero exits.

## 4. End-to-end accuracy and stress runs

```bash
work/env/bin/python tests/e2e/run_pipeline.py --scripts scripts --samples work/sim/aviti_main/samples \
    --sandbox work/runs/fix_main --small-chunks
work/env/bin/python tests/e2e/evaluate.py --truth work/sim/truth/truth_molecules.tsv \
    --read-truth work/sim/aviti_main/read_truth.tsv --results work/runs/fix_main/BCR/results --out work/runs/fix_main/metrics.json
work/env/bin/python tests/e2e/v_usage.py work/sim/truth/truth_molecules.tsv fixed=work/runs/fix_main/BCR/results
work/env/bin/python tests/e2e/artifact_impact.py work/sim/truth/truth_molecules.tsv work/runs/fix_main/BCR/results work/runs/fix_main/impact
work/env/bin/python tests/e2e/run_nanopore.py --scripts scripts --sim work/sim/nano_main --sandbox work/runs/nano_fix \
    --truth work/sim/truth/truth_molecules.tsv
work/env/bin/python tests/e2e/chaos_resume.py --scripts scripts --samples work/sim/aviti_main/samples \
    --sandbox work/runs/chaos_fix_igk --entry igk --kills 8
```

`compare_stitch_methods.py R1only.airr.tsv R2only.airr.tsv IGHV.fasta` needs no truth or IgBLAST
run: it compares the original and the fixed join on existing R1-only/R2-only AIRR files (e.g. the
production results) and reports how many joins change length and which allele pairs cause it.

`rerun_stitch.py` re-runs only stitching + the third IgBLAST pass on an existing sandbox, and
`diagnose_artifacts.py` lists which stitching code path produced any remaining junction artifacts.
