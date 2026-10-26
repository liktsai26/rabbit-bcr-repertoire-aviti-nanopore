"""Shared fixtures for the regression tests.

BCR_TEST_SCRIPTS selects which scripts directory is tested (default: this repo's scripts/),
so the same tests can be pointed at the unmodified upstream scripts to demonstrate each bug.
Tests that need IgBLAST/cutadapt/the IMGT reference are skipped when the test env is missing
(see tests/README.md and tests/fetch_test_data.sh).
"""
import os
import shutil
import stat
import subprocess
import sys

import pytest

TESTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(TESTS, "sim"))
sys.path.insert(0, os.path.join(TESTS, "e2e"))
import run_pipeline as rp  # noqa: E402
from common import ENV_PREFIX, REF_DIR  # noqa: E402

REPO = os.path.dirname(TESTS)
SCRIPTS = os.environ.get("BCR_TEST_SCRIPTS", os.path.join(REPO, "scripts"))
HAVE_TOOLS = os.path.exists(os.path.join(ENV_PREFIX, "bin", "igblastn")) and os.path.exists(
    os.path.join(REF_DIR, "IGHV.fasta.nin"))
needs_tools = pytest.mark.skipif(not HAVE_TOOLS, reason="IgBLAST test env / reference not installed")

FAKE_IGBLASTN = r'''#!/bin/bash
# Test double around the real igblastn: FAKE_DROP_IDS (comma list) simulates the
# igblastn 1.22.0 NULL-pointer behaviour (output stops before that read, exit code 0,
# every time — so the read can never be processed);
# FAKE_FAIL_ONCE=<marker file> makes the first call exit 1 without output;
# FAKE_ZERO_BYTE=1 makes a call whose first read is a dropped read write a 0-byte file.
REAL="__REAL__"
args=("$@"); out=""; q=""
for ((i=0; i<${#args[@]}; i++)); do
  [ "${args[$i]}" = "-out" ] && out="${args[$((i+1))]}"
  [ "${args[$i]}" = "-query" ] && q="${args[$((i+1))]}"
done
if [ -n "$FAKE_FAIL_ONCE" ] && [ ! -f "$FAKE_FAIL_ONCE" ]; then
  touch "$FAKE_FAIL_ONCE"; echo "fake crash" >&2; exit 1
fi
"$REAL" "$@" || exit $?
first=$(grep -m1 '^>' "$q" | cut -c2- | awk '{print $1}')
if [ -n "$FAKE_ZERO_BYTE" ] && [ -n "$FAKE_DROP_IDS" ] && echo ",$FAKE_DROP_IDS," | grep -q ",$first,"; then
  : > "$out"; exit 0     # real igblastn 1.22.0: CEXCEPTION in batch 0 -> 0-byte output, exit 0
fi
if [ -n "$FAKE_DROP_IDS" ] && [ -n "$out" ]; then
  awk -v ids="$FAKE_DROP_IDS" 'BEGIN{n=split(ids,a,",");for(i=1;i<=n;i++)d[a[i]]=1}
       NR==1{print;next} {split($0,f,"\t"); if(f[1] in d){exit} print}' "$out" > "$out.cut" && mv "$out.cut" "$out"
fi
exit 0
'''


def make_fake_env(tmp_path):
    """A fake env prefix: real python/cutadapt/... via symlinks, igblastn replaced by the test double."""
    fake = tmp_path / "fakeenv"
    (fake / "bin").mkdir(parents=True)
    (fake / "share").mkdir()
    os.symlink(os.path.join(ENV_PREFIX, "share", "igblast"), fake / "share" / "igblast")
    for name in os.listdir(os.path.join(ENV_PREFIX, "bin")):
        if name != "igblastn":
            os.symlink(os.path.join(ENV_PREFIX, "bin", name), fake / "bin" / name)
    real = os.path.join(ENV_PREFIX, "share", "igblast", "bin", "igblastn")
    p = fake / "bin" / "igblastn"
    p.write_text(FAKE_IGBLASTN.replace("__REAL__", real))
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return str(fake)


class Sandbox:
    def __init__(self, tmp_path, env_prefix=None, small=True, root_name="sb"):
        self.root = str(tmp_path / root_name)
        self.samples = str(tmp_path / "samples")
        os.makedirs(self.samples, exist_ok=True)
        self.env_prefix = env_prefix or ENV_PREFIX
        saved = rp.ENV_PREFIX
        rp.ENV_PREFIX = self.env_prefix
        try:
            self.base, self.sdir = rp.prepare(SCRIPTS, self.samples, self.root, small)
            self.env = rp.env_for(self.base, small)
        finally:
            rp.ENV_PREFIX = saved
        self.env["BCR_ENV_PREFIX"] = self.env_prefix
        self.env["BCR_PIPELINE_PYTHON"] = os.path.join(self.env_prefix, "bin", "python")
        self.env["BCR_IGDATA"] = os.path.join(self.env_prefix, "share", "igblast")
        self.env["PATH"] = os.path.join(self.env_prefix, "bin") + os.pathsep + self.env["PATH"]
        self.results = os.path.join(self.base, "results")

    def run(self, script, *args, extra_env=None):
        env = dict(self.env, **(extra_env or {}))
        cmd = (["/bin/bash"] if script.endswith(".sh") else [os.path.join(self.env_prefix, "bin", "python")])
        p = subprocess.run(cmd + [os.path.join(self.sdir, script)] + [str(a) for a in args],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
        return p

    def path(self, *parts):
        return os.path.join(self.results, *parts)


@pytest.fixture
def sandbox(tmp_path):
    return Sandbox(tmp_path)


@pytest.fixture
def fake_sandbox(tmp_path):
    return Sandbox(tmp_path, env_prefix=make_fake_env(tmp_path))


def copy_fastq_subset(src_dir, dst_dir, rabbit, n_pairs):
    for side in ("R1", "R2"):
        src = os.path.join(src_dir, f"15651-LT-{rabbit}_{side}.fastq")
        with open(src) as f, open(os.path.join(dst_dir, f"15651-LT-{rabbit}_{side}.fastq"), "w") as g:
            for i, line in enumerate(f):
                if i >= 4 * n_pairs:
                    break
                g.write(line)


SIM_SAMPLES = os.path.join(REPO, "work", "sim", "aviti_main", "samples")
needs_sim = pytest.mark.skipif(not os.path.exists(os.path.join(SIM_SAMPLES, "15651-LT-1_R1.fastq")),
                               reason="simulated reads not generated (tests/sim/*.py)")


def rmtree(p):
    shutil.rmtree(p, ignore_errors=True)
