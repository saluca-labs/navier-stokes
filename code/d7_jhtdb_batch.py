r"""d7_jhtdb_batch.py — DEFINITIVE-as-possible extreme-tail geometric-depletion measurement on real
high-Re DNS (JHTDB isotropic1024coarse, Re_lambda~433), by BATCHING the rate-limited public testing token.

The open NS question our D1-D5 Taylor-Green arc reached but couldn't pass: at our Re (<=3200) the well-
resolved EXTREME tail (top 0.01% |omega|) of the geometric depletion <w.Sw>/<|w|^2|S|> strengthens toward
the mean (toward regularity). Buaria-Pumir report extreme DE-depletion only at far higher Re. JHTDB's
isotropic1024coarse sits at Re_lambda~433 (>> TG Re=3200) — so its extreme tail is the discriminating
measurement. ns_jhtdb.py proved the pipeline at 2000 pts (mean 0.126, top1% 0.113 — no de-depletion). A
real top-0.1% tail needs ~1e5+ points; the testing token caps each call small, so we ACCUMULATE:

  - many small point-batches (CHUNK each), across SLICES independent times (fresh turbulence realizations),
  - persist (w2, P, den) to disk after every batch  -> resumable / survives the session,
  - recompute mean + top {5,1,0.1,0.01}% tails on the running pool each time.

Durable: re-running resumes from the saved pool. Stop anytime; the pool on disk is the result.

  python d7_jhtdb_batch.py [TARGET_POINTS=60000] [CHUNK=2000]
"""
import sys, os, time, numpy as np

TARGET = int(sys.argv[1]) if len(sys.argv) > 1 else 60000
CHUNK  = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
TOKEN  = "edu.jhu.pha.turbulence.testing-201406"
DS     = "isotropic1024coarse"
L      = 2 * np.pi
TIMES  = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]   # independent realizations
POOL   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_d7_pool.npz")

from givernylocal.turbulence_dataset import turb_dataset
from givernylocal.turbulence_toolkit import getData

ds = turb_dataset(dataset_title=DS, output_path="./_giv", auth_token=TOKEN)

# ---- resume from disk if present ---------------------------------------------
if os.path.exists(POOL):
    z = np.load(POOL)
    w2, P, den, nb = list(z["w2"]), list(z["P"]), list(z["den"]), int(z["nb"])
    print(f"resumed pool: {len(w2)} pts, {nb} batches done", flush=True)
else:
    w2, P, den, nb = [], [], [], 0


def depletion(w2a, Pa, dena):
    w2a, Pa, dena = map(np.asarray, (w2a, Pa, dena))
    out = {"mean": Pa.sum() / (dena.sum() + 1e-30), "n": len(w2a)}
    order = np.argsort(w2a)
    for q in (0.05, 0.01, 0.001, 0.0001):
        k = max(1, int(q * len(w2a)))
        idx = order[-k:]
        out[q] = (Pa[idx].sum() / (dena[idx].sum() + 1e-30), k)
    return out


def query(npts, t, seed):
    pts = (np.random.RandomState(seed).rand(npts, 3) * L)
    res = getData(ds, "velocity", float(t), "none", "fd4lag4", "gradient", pts)
    arr = np.asarray(getattr(res, "result", res))
    if arr.dtype == object or arr.ndim == 1:
        arr = np.asarray([np.asarray(r) for r in res]).reshape(npts, -1)
    A = arr.reshape(npts, 3, 3)                      # A[:,i,j] = d u_i / d x_j
    S = 0.5 * (A + A.transpose(0, 2, 1))
    w = np.stack([A[:, 2, 1] - A[:, 1, 2],
                  A[:, 0, 2] - A[:, 2, 0],
                  A[:, 1, 0] - A[:, 0, 1]], -1)
    w2b = (w * w).sum(-1)
    Pb  = np.einsum("ki,kij,kj->k", w, S, w)
    Smag = np.sqrt(2 * (S * S).sum((-1, -2)))
    return w2b, Pb, w2b * Smag


print(f"target {TARGET} pts, chunk {CHUNK}, {len(TIMES)} time-slices, Re_lambda~433\n", flush=True)
while len(w2) < TARGET:
    t = TIMES[nb % len(TIMES)]
    seed = 1000 + nb                       # distinct points each batch
    t0 = time.time()
    try:
        a, b, c = query(CHUNK, t, seed)
    except Exception as e:
        print(f"  batch {nb} FAILED ({type(e).__name__}: {str(e)[:80]}); backing off 20s", flush=True)
        time.sleep(20); continue
    w2 += list(a); P += list(b); den += list(c); nb += 1
    np.savez(POOL, w2=np.asarray(w2), P=np.asarray(P), den=np.asarray(den), nb=nb)
    d = depletion(w2, P, den)
    msg = f"  batch {nb:3d} t={t:4.1f} ({time.time()-t0:4.1f}s)  N={d['n']:6d}  mean={d['mean']:.3f}"
    for q in (0.05, 0.01, 0.001):
        msg += f"  top{q*100:g}%={d[q][0]:.3f}"
    print(msg, flush=True)

# ---- final report ------------------------------------------------------------
d = depletion(w2, P, den)
print("\n================ d7 FINAL (JHTDB isotropic1024coarse, Re_lambda~433) ================")
print(f"  pooled points : {d['n']}  across {len(TIMES)} independent time-slices")
print(f"  MEAN depletion <w.Sw>/<|w|^2|S|> = {d['mean']:.4f}")
for q in (0.05, 0.01, 0.001, 0.0001):
    val, k = d[q]
    print(f"  top {q*100:6.3f}% |w| tail        = {val:.4f}   ({k} pts)")
print("\n  TG arc (well-resolved): extremes -> mean toward regularity through Re=3200.")
print("  Discriminator: is the high-Re extreme tail MORE depleted (regularizing) or LESS (de-depletion)?")
