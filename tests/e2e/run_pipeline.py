"""Run the AVITI pipeline scripts end-to-end in an isolated sandbox.

The original scripts hardcode /Users/newuser/... paths, so they are copied into
SANDBOX/script/ with those prefixes rewritten (the fixed scripts read BCR_* env vars
instead; both are handled). --small-chunks shrinks chunk constants so that multi-chunk,
multi-sub-batch and resume paths are exercised with a few thousand reads.

Usage:
  python run_pipeline.py --scripts DIR --samples DIR --sandbox DIR [--small-chunks]
                         [--chains IGH,IGK,IGL] [--rabbits 1,2,3]
Every step's exit code, wall time and log go to SANDBOX/run_log.tsv and SANDBOX/logs/.
Uses /bin/bash (3.2 on macOS, as on the authors' machine).
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sim"))
from common import ENV_PREFIX, REF_DIR  # noqa: E402

ORIG_BASE = "/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline"
ORIG_ENV = "/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env"
SMALL = {"CHUNK_RAW_PAIRS": "1000", "IGBLAST_BATCH": "300", "CHUNK_SIZE": "700", "BATCH_SIZE": "300"}


def prepare(scripts, samples, sandbox, small):
    base = os.path.join(sandbox, "BCR")
    sdir = os.path.join(base, "script")
    os.makedirs(sdir, exist_ok=True)
    os.makedirs(os.path.join(base, "reference"), exist_ok=True)
    os.makedirs(os.path.join(base, "results"), exist_ok=True)
    os.makedirs(os.path.join(sandbox, "logs"), exist_ok=True)
    link = os.path.join(base, "reference", "IMGT_reference_clean")
    if not os.path.exists(link):
        os.symlink(REF_DIR, link)
    sl = os.path.join(base, "samples")
    if not os.path.exists(sl):
        os.symlink(os.path.abspath(samples), sl)
    for fn in os.listdir(scripts):
        src = os.path.join(scripts, fn)
        if not os.path.isfile(src):
            continue
        txt = open(src, encoding="utf-8").read()
        txt = txt.replace(ORIG_BASE, base).replace(ORIG_ENV, ENV_PREFIX)
        if small:
            txt = re.sub(r"^(CHUNK_RAW_PAIRS)=1500000", r"\1=" + SMALL["CHUNK_RAW_PAIRS"], txt, flags=re.M)
            txt = re.sub(r"^(IGBLAST_BATCH)=25000", r"\1=" + SMALL["IGBLAST_BATCH"], txt, flags=re.M)
            txt = re.sub(r"^(CHUNK_SIZE) = 500_000", r"\1 = " + SMALL["CHUNK_SIZE"], txt, flags=re.M)
            txt = re.sub(r"^(BATCH_SIZE) = 25000", r"\1 = " + SMALL["BATCH_SIZE"], txt, flags=re.M)
        with open(os.path.join(sdir, fn), "w", encoding="utf-8") as f:
            f.write(txt)
        shutil.copymode(src, os.path.join(sdir, fn))
    return base, sdir


def env_for(base, small):
    env = dict(os.environ)
    env["PATH"] = os.path.join(ENV_PREFIX, "bin") + os.pathsep + env.get("PATH", "")
    env.update({"BCR_PIPELINE_ROOT": base, "BCR_PIPELINE_PYTHON": os.path.join(ENV_PREFIX, "bin", "python"),
                "BCR_IGDATA": os.path.join(ENV_PREFIX, "share", "igblast"), "BCR_ENV_PREFIX": ENV_PREFIX,
                "BCR_SCRIPT_DIR": os.path.join(base, "script")})
    if small:
        env.update({"BCR_" + k: v for k, v in SMALL.items()})
    return env


def step(sandbox, name, cmd, env, log_rows):
    t0 = time.time()
    logp = os.path.join(sandbox, "logs", name + ".log")
    with open(logp, "w") as lf:
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)
    dt = time.time() - t0
    log_rows.append((name, p.returncode, f"{dt:.1f}"))
    print(f"{name}\texit={p.returncode}\t{dt:.1f}s", flush=True)
    return p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--sandbox", required=True)
    ap.add_argument("--small-chunks", action="store_true")
    ap.add_argument("--chains", default="IGH,IGK,IGL")
    ap.add_argument("--rabbits", default="1,2,3")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    base, sdir = prepare(a.scripts, a.samples, a.sandbox, a.small_chunks)
    env = env_for(base, a.small_chunks)
    py = os.path.join(ENV_PREFIX, "bin", "python")
    rows = []
    chains = a.chains.split(",")
    tag = a.tag
    if "IGH" in chains:
        for r in a.rabbits.split(","):
            for s in ["02c_run_aviti_r2only_igh.sh", "02d_run_aviti_r2only_igblast.sh",
                      "02e_run_aviti_r1only_igh.sh", "02f_run_aviti_r1only_igblast.sh"]:
                step(a.sandbox, f"{tag}IGH_r{r}_{s[:3]}", ["/bin/bash", os.path.join(sdir, s), r], env, rows)
            step(a.sandbox, f"{tag}IGH_r{r}_02g", [py, os.path.join(sdir, "02g_stitch_qc_igh.py"), r], env, rows)
            step(a.sandbox, f"{tag}IGH_r{r}_02h", ["/bin/bash", os.path.join(sdir, "02h_stitch_igblast.sh"), r], env, rows)
    if "IGK" in chains:
        step(a.sandbox, f"{tag}IGK_full", ["/bin/bash", os.path.join(sdir, "run_igk_full_pipeline.sh")], env, rows)
    if "IGL" in chains:
        step(a.sandbox, f"{tag}IGL_full", ["/bin/bash", os.path.join(sdir, "run_igl_full_pipeline.sh")], env, rows)
    with open(os.path.join(a.sandbox, "run_log.tsv"), "a") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
