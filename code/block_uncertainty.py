r"""block_uncertainty.py — honest error bars on the production tail ratios that respect within-cutout
spatial correlation. The pool_cutout_*.npz arrays are concatenated box-by-box, so CONTIGUOUS blocks are
groups of whole 128^3 cutouts. We split each pool into K contiguous blocks, compute the conditional
depletion ratio D(q)/D_mean independently in each block, and report mean +/- std and SEM across blocks.

This is a block (not point) resample: it does NOT pretend the millions of correlated points inside a box
are independent. Point bootstrap would give falsely-tight intervals (the paper warns against exactly this).

  python block_uncertainty.py <pool.npz> <label> [K]
"""
import sys, numpy as np

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
for q in qs:
    r = np.array(rat[q])
    sem = r.std(ddof=1) / np.sqrt(K)
    lo, hi = r.mean() - 1.96 * sem, r.mean() + 1.96 * sem
    tag = "top%.3g%%" % (q * 100)
    flag = "  <-- straddles 1.0" if lo < 1.0 < hi else ("  (resolved <1)" if hi < 1.0 else "  (resolved >1)")
    print(f"  {tag:9} D/mean = {r.mean():.4f} +/- {r.std(ddof=1):.4f} (block sd)  "
          f"95%CI[{lo:.4f},{hi:.4f}]{flag}")
