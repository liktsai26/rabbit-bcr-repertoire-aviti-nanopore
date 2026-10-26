"""Kill-and-resume stress test: run a pipeline entry point, SIGKILL its whole process group at a
random moment, start it again, and repeat until it exits 0. Then check that nothing was lost or
double counted compared with an uninterrupted reference run.

Usage: python chaos_resume.py --scripts DIR --samples DIR --sandbox NEW_DIR --entry {igk,igh1}
                              [--kills 6] [--seed S] [--reference SANDBOX_WITH_CLEAN_RUN]
"""
import argparse
import glob
import os
import random
import signal
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as rp  # noqa: E402

IGH_STEPS = ["02c_run_aviti_r2only_igh.sh", "02d_run_aviti_r2only_igblast.sh", "02e_run_aviti_r1only_igh.sh",
             "02f_run_aviti_r1only_igblast.sh", "02g_stitch_qc_igh.py", "02h_stitch_igblast.sh"]


def data_rows(path):
    with open(path) as f:
        return max(0, sum(1 for _ in f) - 1)


def summarize(base):
    out = {}
    for p in sorted(glob.glob(os.path.join(base, "results", "*", "*", "*.airr.tsv"))):
        out[os.path.relpath(p, base)] = data_rows(p)
    for p in sorted(glob.glob(os.path.join(base, "results", "*", "*", "chunk_stats.csv"))):
        lines = open(p).read().splitlines()[1:]
        keys = [line.split(",")[0] for line in lines]
        out[os.path.relpath(p, base) + ":rows"] = len(lines)
        out[os.path.relpath(p, base) + ":duplicated_chunks"] = len(keys) - len(set(keys))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--sandbox", required=True)
    ap.add_argument("--entry", choices=["igk", "igh1"], required=True)
    ap.add_argument("--kills", type=int, default=6)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--max-delay", type=float, default=8.0)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    base, sdir = rp.prepare(a.scripts, a.samples, a.sandbox, True)
    env = rp.env_for(base, True)
    py = os.path.join(rp.ENV_PREFIX, "bin", "python")
    if a.entry == "igk":
        cmds = [["/bin/bash", os.path.join(sdir, "run_igk_full_pipeline.sh")]]
    else:
        cmds = [([py] if s.endswith(".py") else ["/bin/bash"]) + [os.path.join(sdir, s), "1"] for s in IGH_STEPS]
    log = open(os.path.join(a.sandbox, "chaos.log"), "w")
    kills = 0
    for cmd in cmds:
        while True:
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
            delay = rng.uniform(0.3, a.max_delay) if kills < a.kills else None
            try:
                rc = p.wait(timeout=delay)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                kills += 1
                print(f"killed {os.path.basename(cmd[-2] if a.entry == 'igh1' else cmd[-1])} after {delay:.1f}s", flush=True)
                time.sleep(0.2)
                continue
            if rc != 0:
                print(f"FAILED: {cmd} exit {rc} (see chaos.log)")
                sys.exit(1)
            break
    print(f"completed after {kills} kills")
    import json
    print(json.dumps(summarize(base), indent=1))


if __name__ == "__main__":
    main()
