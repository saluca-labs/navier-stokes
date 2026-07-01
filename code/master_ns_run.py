#!/usr/bin/env python
r"""master_ns_run.py — ONE comprehensive, resumable run of the whole NS geometric-depletion campaign,
capturing everything into a single downloadable report: MASTER_RESULTS.md (in the persistent dir).

Drives the existing, validated tools SERIALLY (JHTDB = one job at a time) — it does NOT re-implement
the physics. Per dataset:

  STAGE 0  recon            — mounted datasets, package versions.
  STAGE 1  METHOD GATE      — fd_vs_spectral.py across the gradient-order ladder
                              (fd4noint -> fd6noint -> fd8noint -> m2q8): does the extreme-tail
                              depletion D/mean order-CONVERGE, and to what percentile depth?
  STAGE 2  PRODUCTION TAILS — sciserver_cutout.py at high statistics: mean + top 5/1/0.1/0.01/0.001%.

The report puts isotropic8192 (Re_lambda~1300) next to isotropic32768 (Re_lambda~2500) so the
cross-Re trend is readable in one place.

RESUMABLE: master_progress.json marks finished (dataset, step); re-running skips them and appends.
Every sub-step's full stdout is captured verbatim. Walk away; download MASTER_RESULTS.md when done.

RUN (in the SciServer container, in the dir holding fd_vs_spectral.py + sciserver_cutout.py):
  export JHTDB_TOKEN=YOUR_REGISTERED_JHTDB_TOKEN
  nohup python master_ns_run.py > master_run.log 2>&1 &
  tail -f master_run.log

Config via env (all optional):
  DATASETS=isotropic8192,isotropic32768   LADDER=fd4noint,fd6noint,fd8noint,m2q8
  LADDER_NBOX=300   PROD_B=128   PROD_NBOX=120   STEP_TIMEOUT=21600
"""
import os, sys, json, time, subprocess, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
PERSIST = Path("/home/idies/workspace/Storage/saluca/persistent")
OUTDIR = PERSIST if PERSIST.is_dir() else HERE
REPORT = OUTDIR / "MASTER_RESULTS.md"
PROGRESS = OUTDIR / "master_progress.json"
PY = sys.executable

TOKEN = os.environ.get("JHTDB_TOKEN", "")
DATASETS = [d for d in os.environ.get("DATASETS", "isotropic8192,isotropic32768").split(",") if d.strip()]
LADDER = [m for m in os.environ.get("LADDER", "fd4noint,fd6noint,fd8noint,m2q8").split(",") if m.strip()]
LADDER_NBOX = os.environ.get("LADDER_NBOX", "300")
PROD_B = os.environ.get("PROD_B", "128")
PROD_NBOX = os.environ.get("PROD_NBOX", "120")
STEP_TIMEOUT = int(os.environ.get("STEP_TIMEOUT", "21600"))   # 6 h per step ceiling


def log(msg):
    print(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def append(md):
    with open(REPORT, "a", encoding="utf-8") as f:
        f.write(md + "\n")


def progress():
    try:
        return json.loads(PROGRESS.read_text())
    except Exception:
        return {}


def mark(key):
    p = progress(); p[key] = "done"; PROGRESS.write_text(json.dumps(p, indent=2))


def step(key, title, cmd, env_extra):
    """Run one sub-step; capture stdout+stderr; append to the report; mark progress (resumable)."""
    if progress().get(key) == "done":
        log(f"SKIP (already done): {key}")
        return
    pretty = " ".join(f"{k}={v}" for k, v in env_extra.items() if k != "JHTDB_TOKEN") + " " + " ".join(cmd)
    log(f"RUN  {key}  ::  {pretty}")
    append(f"\n### {title}\n\n```text\n$ {pretty}")
    env = {**os.environ, **env_extra}
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=str(HERE), env=env, capture_output=True, text=True, timeout=STEP_TIMEOUT)
        out = (r.stdout or "").rstrip()
        if r.stderr and r.stderr.strip():
            out += "\n[stderr]\n" + r.stderr.strip()[-4000:]
        if r.returncode != 0:
            out += f"\n[exit code {r.returncode}]"
    except subprocess.TimeoutExpired:
        out = f"STEP TIMED OUT after {STEP_TIMEOUT}s"
    except Exception as e:
        out = f"STEP FAILED: {type(e).__name__}: {e}"
    dt = time.time() - t0
    append(out + f"\n```\n_completed in {dt:.0f}s_\n")
    mark(key)
    log(f"DONE {key}  ({dt:.0f}s)")


def main():
    if not TOKEN:
        sys.exit("Set JHTDB_TOKEN first:  export JHTDB_TOKEN=YOUR_REGISTERED_JHTDB_TOKEN")
    if not REPORT.exists():
        append("# Navier–Stokes geometric-depletion — MASTER RUN\n")
        append(f"_Started {datetime.datetime.now():%Y-%m-%d %H:%M} · token …{TOKEN[-8:]} · "
               f"datasets {DATASETS} · ladder {LADDER} · LADDER_NBOX={LADDER_NBOX} · "
               f"PROD B={PROD_B} NBOX={PROD_NBOX}_\n")
        append("**What this is.** D = ⟨ω·Sω⟩/⟨|ω|²|S|⟩ conditioned on |ω| percentile. Stage 1 is the "
               "gradient-ORDER convergence gate (is the extreme tail trustworthy, and to what depth). "
               "Stage 2 is the high-statistics tail per dataset. Read them together: Stage 1 says which "
               "operator/depth to trust; Stage 2 gives the numbers; 8192-vs-32768 at the converged depth "
               "is the Re trend.\n")

    step("recon", "Stage 0 — environment recon",
         ["bash", "-lc",
          "echo 'mounted datasets:'; ls /home/idies/workspace/turbulence-ceph/sciserver-turbulence 2>/dev/null | grep -v checkpoint; "
          "python -c \"import numpy,giverny; print('numpy',numpy.__version__,'| giverny',getattr(giverny,'__version__','?'))\" 2>/dev/null"],
         {})

    for ds in DATASETS:
        append(f"\n---\n## {ds}\n")
        # STAGE 1 — method gate: the order ladder
        for m in LADDER:
            step(f"{ds}:gate:{m}", f"Stage 1 · {ds} · method gate — GRAD_METHOD={m}",
                 [PY, "fd_vs_spectral.py", "--nbox", LADDER_NBOX],
                 {"DATASET": ds, "GRAD_METHOD": m, "JHTDB_TOKEN": TOKEN})
        # STAGE 2 — production tails (high statistics, 4th-order cutout)
        step(f"{ds}:prod", f"Stage 2 · {ds} · production tails — cutout B={PROD_B}, NBOX={PROD_NBOX}",
             [PY, "sciserver_cutout.py"],
             {"DATASET": ds, "B": PROD_B, "NBOX": PROD_NBOX, "JHTDB_TOKEN": TOKEN})

    append(f"\n---\n## DONE — {datetime.datetime.now():%Y-%m-%d %H:%M}\n")
    append("**How to read this report:**\n"
           "1. **Method gate (Stage 1):** per dataset, compare the `JHTDB-<method>` tail ratios across "
           "fd4 → fd6 → fd8 → m2q8. If they PLATEAU (fd6≈fd8≈m2q8), the tail is order-converged — that "
           "plateau is the physical D/mean, claimable to the depth where it stops drifting. If they keep "
           "DRIFTING at top0.01%/0.001%, the DNS doesn't resolve that depth → back off to the converged "
           "percentile. (The `FD-4th` row is `sciserver_cutout.py`'s operator — if it matches the converged "
           "reference, the high-statistics Stage-2 numbers are trustworthy as-is; if not, only its shallower "
           "tails are.)\n"
           "2. **Production (Stage 2):** the high-N tails. Trust them to the depth Stage 1 says is converged.\n"
           "3. **The NS trend:** isotropic8192 (Re_λ~1300) vs isotropic32768 (Re_λ~2500), tail D/mean at the "
           "converged depth. Does de-depletion (tail/mean > 1) STRENGTHEN with Re? That is the result.\n")
    log(f"ALL DONE -> {REPORT}")


if __name__ == "__main__":
    main()
