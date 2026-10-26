"""Re-run only the stitching + third IgBLAST pass (02g/02h, 04c/04d) of SCRIPTS on the
R1-only/R2-only results of an existing sandbox, into a NEW sandbox (nothing is deleted).

Usage: python rerun_stitch.py --scripts DIR --from SANDBOX --sandbox NEW_SANDBOX --samples DIR
"""
import argparse
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as rp  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--scripts", required=True)
ap.add_argument("--from", dest="src", required=True)
ap.add_argument("--sandbox", required=True)
ap.add_argument("--samples", required=True)
a = ap.parse_args()
if os.path.exists(a.sandbox):
    sys.exit(f"{a.sandbox} exists; choose a new sandbox name")
base, sdir = rp.prepare(a.scripts, a.samples, a.sandbox, True)
for d in glob.glob(os.path.join(a.src, "BCR", "results", "*", "*_R[12]only")):
    dst = os.path.join(base, "results", os.path.basename(os.path.dirname(d)), os.path.basename(d))
    shutil.copytree(d, dst, ignore=shutil.ignore_patterns("chunks"))
env = rp.env_for(base, True)
rows = []
py = os.path.join(rp.ENV_PREFIX, "bin", "python")
for r in "123":
    rp.step(a.sandbox, f"IGH_r{r}_02g", [py, os.path.join(sdir, "02g_stitch_qc_igh.py"), r], env, rows)
    rp.step(a.sandbox, f"IGH_r{r}_02h", ["/bin/bash", os.path.join(sdir, "02h_stitch_igblast.sh"), r], env, rows)
    for ch in ("IGK", "IGL"):
        rp.step(a.sandbox, f"{ch}_r{r}_04c", [py, os.path.join(sdir, "04c_stitch_qc_igkl.py"), r, ch], env, rows)
        rp.step(a.sandbox, f"{ch}_r{r}_04d", ["/bin/bash", os.path.join(sdir, "04d_stitch_igblast_igkl.sh"), r, ch], env, rows)
