r"""block_uncertainty.py — honest error bars on the production tail ratios that respect within-cutout
spatial correlation. The pool_cutout_*.npz arrays are concatenated box-by-box, so CONTIGUOUS blocks are
groups of whole cutouts. We split each pool into K contiguous blocks, compute the conditional depletion
ratio D(q)/D_mean independently in each block, and report mean +/- t_{K-1}*SEM across blocks.

NAMING + MULTIPLIER CORRECTION (v0.6, 2026-07-01): this is a BATCH-MEANS (block-means) partition
estimator, not a resample (nothing is resampled); earlier revisions called it a "block resample" and used
the normal 1.96 multiplier, which is ~7% too narrow at K=20 — the correct multiplier is Student t_{K-1}
(2.093 at K=20). It still does NOT pretend the millions of correlated points inside a box are independent;
a point bootstrap would give falsely-tight intervals (the paper warns against exactly this). Note the
block-mean-of-ratios is a slightly different estimand from the pooled ratio-of-sums (each block's top-q
threshold is shallower than the pool's); the paper reports both.

  python block_uncertainty.py <pool.npz> <label> [K]
"""
import sys, numpy as np

# t_{K-1} 97.5% quantiles, df 1..40 (fallback if scipy absent)
_T = [12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.36, 2.31, 2.26, 2.23, 2.20, 2.18, 2.16, 2.14, 2.13,
      2.12, 2.11, 2.10, 2.09, 2.093, 2.08, 2.07, 2.07, 2.06, 2.06, 2.06, 2.05, 2.05, 2.05, 2.04,
      2.04, 2.04, 2.03, 2.03, 2.03, 2.03, 2.03, 2.02, 2.02, 2.02]


def tmult(k):
    try:
        from scipy import stats
        return float(stats.t.ppf(0.975, k - 1))
    except Exception:
        return _T[min(k - 2, len(_T) - 1)] if k > 1 else float("nan")

path = sys.argv[1]
label = sys.argv[2] if len(sys.argv) > 2 else path
K = int(sys.argv[3]) if len(sys.argv) > 3 else 20

z = np.load(path)
w2, P, den = z["w2"], z["P"], z["den"]
N = len(w2)
# split into K contiguous blocks (= groups of whole cutouts, correlation preserved)
edges = np.linspace(0, N, K + 1).astype(np.int64)


def ratio(w2b, Pb, denb, q):
    mean = Pb.sum() / (denb.sum() + 1e-30)
    k = max(1, int(q * len(w2b)))
    o = np.argpartition(w2b, -k)[-k:]
    tail = Pb[o].sum() / (denb[o].sum() + 1e-30)
    return tail / mean, mean


qs = [0.01, 0.001]
means = []
rat = {q: [] for q in qs}
for i in range(K):
    a, b = edges[i], edges[i + 1]
    w2b, Pb, denb = w2[a:b], P[a:b], den[a:b]
    _, m = ratio(w2b, Pb, denb, 0.01)
    means.append(m)
    for q in qs:
        r, _ = ratio(w2b, Pb, denb, q)
        rat[q].append(r)

print(f"\n==== {label} ====")
print(f"  N={N}  K={K} contiguous blocks (~{N//K:,} pts each)")
mm = np.array(means)
print(f"  D_mean        = {mm.mean():.4f} +/- {mm.std(ddof=1):.4f} (block sd)   SEM={mm.std(ddof=1)/np.sqrt(K):.4f}")
tm = tmult(K)
for q in qs:
    r = np.array(rat[q])
    sem = r.std(ddof=1) / np.sqrt(K)
    lo, hi = r.mean() - tm * sem, r.mean() + tm * sem
    tag = "top%.3g%%" % (q * 100)
    flag = "  <-- straddles 1.0" if lo < 1.0 < hi else ("  (resolved <1)" if hi < 1.0 else "  (resolved >1)")
    print(f"  {tag:9} D/mean = {r.mean():.4f} +/- {r.std(ddof=1):.4f} (block sd)  "
          f"95% t_{{{K-1}}}-CI[{lo:.4f},{hi:.4f}]{flag}")
