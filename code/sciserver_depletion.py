#!/usr/bin/env python
r"""sciserver_depletion.py  —  RUN INSIDE THE SciServer Compute container (_salucaTurbulence).

The geometric-depletion extreme-tail measurement, now with the rate-limit GONE: inside SciServer the JHTDB
DNS is mounted read-only at /home/idies/workspace/turbulence-ceph and `giverny` reads it locally (no web
service). With 100 GB RAM we can pull MILLIONS of velocity-gradient points and resolve the top 0.01% / 0.001%
extreme tail of  <w.Sw>/<|w|^2|S|>  with real statistics — the quantity our local TG arc (Re<=3200) and the
60k web-API pool (Re_lambda~433, no de-depletion) reached but couldn't push past.

THE PRIZE: if the high-Re isotropic sets (isotropic4096 ~Re_lambda 610, isotropic8192 ~Re_lambda 1300) are
mounted, this lands us squarely in the Buaria-Pumir extreme-DE-depletion band — the open frontier.

HOW TO RUN (in the container):
  Option A (terminal):  put this file in ~/workspace/Storage/saluca/persistent/ and:
        python sciserver_depletion.py
  Option B (notebook):  paste this whole file into one Python-3 cell and run it.

Tunables via env:  NPTS (default 2_000_000), CHUNK (200_000), DATASETS (comma list; default auto+high-Re).
Durable: raw (w2,P,den) + summary saved to PERSIST dir after every dataset; safe to re-run / resume per-set.
"""
import os, sys, time, traceback
import numpy as np

# ---------- config ----------
NPTS  = int(os.environ.get("NPTS", "2000000"))
CHUNK = int(os.environ.get("CHUNK", "4096"))   # testing token caps getData at 4096 pts/call (client-side policy)
PERSIST = "/home/idies/workspace/Storage/saluca/persistent"
if not os.path.isdir(PERSIST):
    PERSIST = os.getcwd()
OUT = os.path.join(PERSIST, "depletion_out")
os.makedirs(OUT, exist_ok=True)
SCRATCH = "/home/idies/workspace/Temporary/saluca/scratch"
WORK = SCRATCH if os.path.isdir(SCRATCH) else os.path.join(OUT, "_giv")
os.makedirs(WORK, exist_ok=True)
TOKEN = os.environ.get("JHTDB_TOKEN", "edu.jhu.pha.turbulence.testing-201406")  # local mode ignores quota

# ---------- recon: what data is mounted, which package is here ----------
CEPH = "/home/idies/workspace/turbulence-ceph"
print("=" * 78)
print("RECON")
print("  persist out :", OUT)
print("  scratch work:", WORK)
print("  ceph mount  :", CEPH, "->", "PRESENT" if os.path.isdir(CEPH) else "MISSING")
mounted = []
if os.path.isdir(CEPH):
    try:
        mounted = sorted(os.listdir(CEPH))
        print("  mounted entries:")
        for d in mounted:
            print("     ", d)
            sub = os.path.join(CEPH, d)
            if os.path.isdir(sub) and not d.startswith("."):
                try:
                    kids = sorted(os.listdir(sub))
                    for k in kids[:80]:
                        print("         /", k)
                    if len(kids) > 80:
                        print(f"         ... (+{len(kids)-80} more)")
                except Exception:
                    pass
    except Exception as e:
        print("  (could not list ceph mount:", e, ")")

pkg = None
for name in ("giverny", "givernylocal"):
    try:
        __import__(name)
        pkg = name
        break
    except Exception:
        continue
print("  giverny package:", pkg)
if pkg is None:
    print("\n  No giverny package importable. Try:  pip install giverny   (or givernylocal)")
    sys.exit(1)

if pkg == "giverny":
    from giverny.turbulence_dataset import turb_dataset
    from giverny.turbulence_toolkit import getData
else:
    from givernylocal.turbulence_dataset import turb_dataset
    from givernylocal.turbulence_toolkit import getData
print("=" * 78, flush=True)

L = 2 * np.pi

# Per-dataset VALID TIMES. JHTDB datasets differ: isotropic1024coarse uses CONTINUOUS float time (0..10.056);
# the big high-Re sets use INTEGER snapshot indices (isotropic8192 -> [1..6], isotropic32768 -> [1]). Passing
# the wrong type/value errors out. Sample across all available times (independent realizations) for statistics.
TIME_MAP = {
    "isotropic1024coarse": [round(float(t), 1) for t in np.linspace(0.0, 10.0, 11)],  # Re_lambda~433 (web)
    "isotropic1024fine":   [round(float(t), 1) for t in np.linspace(0.0, 10.0, 11)],
    "isotropic4096":       [1, 2, 3, 4, 5],          # if present
    "isotropic8192":       [1, 2, 3, 4, 5, 6],        # Re_lambda~1300  PRIZE (6 snapshots)
    "isotropic32768":      [1],                       # largest DNS, single snapshot
}

def times_for(title):
    return TIME_MAP.get(title, [1])   # default integer-1 for unknown high-Re sets

env_ds = os.environ.get("DATASETS", "").strip()
if env_ds:
    DATASETS = [(d.strip(), times_for(d.strip())) for d in env_ds.split(",") if d.strip()]
else:
    DATASETS = [(t, times_for(t)) for t in ("isotropic8192", "isotropic32768", "isotropic1024coarse")]


def depletion(w2, P, den):
    w2, P, den = map(np.asarray, (w2, P, den))
    out = {"mean": float(P.sum() / (den.sum() + 1e-30)), "n": int(len(w2))}
    order = np.argsort(w2)
    for q in (0.05, 0.01, 0.001, 0.0001):
        k = max(1, int(q * len(w2)))
        idx = order[-k:]
        out[q] = (float(P[idx].sum() / (den[idx].sum() + 1e-30)), int(k))
    return out


def grad_to_depletion(arr, n):
    A = np.asarray(arr).reshape(n, 3, 3)               # A[:,i,j] = d u_i / d x_j
    S = 0.5 * (A + A.transpose(0, 2, 1))
    w = np.stack([A[:, 2, 1] - A[:, 1, 2],
                  A[:, 0, 2] - A[:, 2, 0],
                  A[:, 1, 0] - A[:, 0, 1]], -1)
    w2 = (w * w).sum(-1)
    P  = np.einsum("ki,kij,kj->k", w, S, w)
    Smag = np.sqrt(2 * (S * S).sum((-1, -2)))
    return w2, P, w2 * Smag


def run_dataset(title, times):
    print(f"\n########## {title}  (target {NPTS} pts, chunk {CHUNK}) ##########", flush=True)
    raw = os.path.join(OUT, f"pool_{title}.npz")
    if os.path.exists(raw):
        z = np.load(raw)
        w2, P, den, b = list(z["w2"]), list(z["P"]), list(z["den"]), int(z["b"])
        print(f"  resumed {len(w2)} pts ({b} batches)", flush=True)
    else:
        w2, P, den, b = [], [], [], 0
    try:
        ds = turb_dataset(dataset_title=title, output_path=WORK, auth_token=TOKEN)
    except Exception as e:
        print(f"  SKIP {title}: turb_dataset failed: {type(e).__name__}: {str(e)[:160]}")
        return None
    while len(w2) < NPTS:
        t = times[b % len(times)]
        rng = np.random.RandomState(7000 + b)
        pts = rng.rand(CHUNK, 3) * L
        t0 = time.time()
        try:
            res = getData(ds, "velocity", t, "none", "fd4lag4", "gradient", pts)  # t type preserved (int vs float)
            arr = np.asarray(getattr(res, "result", res))
            if arr.dtype == object or arr.ndim == 1:
                arr = np.asarray([np.asarray(r) for r in res]).reshape(CHUNK, -1)
            a, p, c = grad_to_depletion(arr, CHUNK)
        except Exception as e:
            print(f"  batch {b} FAILED: {type(e).__name__}: {str(e)[:140]}")
            traceback.print_exc()
            break
        w2 += list(a); P += list(p); den += list(c); b += 1
        np.savez(raw, w2=np.asarray(w2), P=np.asarray(P), den=np.asarray(den), b=b)
        d = depletion(w2, P, den)
        print(f"  batch {b:3d} t={t:4.1f} ({time.time()-t0:5.1f}s) N={d['n']:8d} "
              f"mean={d['mean']:.4f} top1%={d[0.01][0]:.4f} top0.1%={d[0.001][0]:.4f} "
              f"top0.01%={d[0.0001][0]:.4f}", flush=True)
    d = depletion(w2, P, den) if w2 else None
    if d:
        print(f"\n  == {title} FINAL ==  N={d['n']}")
        print(f"     MEAN      = {d['mean']:.4f}")
        for q in (0.05, 0.01, 0.001, 0.0001):
            v, k = d[q]
            print(f"     top {q*100:7.4f}% = {v:.4f}  ({k} pts)")
    return d


def main():
    summary = {}
    for title, times in DATASETS:
        try:
            d = run_dataset(title, times)
            if d:
                summary[title] = d
        except Exception as e:
            print(f"  {title} crashed: {e}")
            traceback.print_exc()
    # write a compact human summary
    sp = os.path.join(OUT, "SUMMARY.txt")
    with open(sp, "w") as f:
        f.write("geometric depletion <w.Sw>/<|w|^2|S|>  — SciServer local (no rate limit)\n")
        f.write(f"NPTS target {NPTS}, chunk {CHUNK}, package {pkg}\n\n")
        f.write(f"{'dataset':22s} {'N':>9s} {'mean':>7s} {'top1%':>7s} {'top.1%':>7s} {'top.01%':>8s} {'top.001%':>9s}\n")
        for t, d in summary.items():
            f.write(f"{t:22s} {d['n']:9d} {d['mean']:7.4f} {d[0.01][0]:7.4f} {d[0.001][0]:7.4f} "
                    f"{d[0.0001][0]:8.4f} {'-':>9s}\n")
        f.write("\nDISCRIMINATOR: tail < mean = regularizing (bulk-like or more depleted); "
                "tail >> mean = extreme DE-depletion (singularity-relevant). Watch the trend with Re_lambda:\n"
                "isotropic1024coarse~433  <  isotropic4096~610  <  isotropic8192~1300.\n")
    print("\nwrote", sp)
    print("\n================ CROSS-DATASET SUMMARY ================")
    for t, d in summary.items():
        print(f"  {t:22s} N={d['n']:8d} mean={d['mean']:.4f} top1%={d[0.01][0]:.4f} "
              f"top0.1%={d[0.001][0]:.4f} top0.01%={d[0.0001][0]:.4f}")


main()
