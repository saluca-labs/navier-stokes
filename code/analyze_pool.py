r"""analyze_pool.py — rich diagnostics on a depletion pool_<dataset>.npz produced by sciserver_depletion.py.

Beyond the SUMMARY quantiles: the full tail/mean RATIO curve (the discriminator), the negative-production
fraction, and production concentration. Run on each dataset's pool and compare across Re_lambda.

  python analyze_pool.py path\to\pool_isotropic8192.npz [label]
"""
import sys, numpy as np

path = sys.argv[1]
label = sys.argv[2] if len(sys.argv) > 2 else path
z = np.load(path)
w2, P, den = z["w2"], z["P"], z["den"]
b = int(z["b"]) if "b" in z else -1
N = len(w2)
mean = P.sum() / (den.sum() + 1e-30)

print(f"\n==== {label} ====")
print(f"  N={N}  batches={b}")
print(f"  MEAN depletion <w.Sw>/<|w|^2|S|> = {mean:.4f}")
print(f"  {'quantile':>10s} {'depletion':>10s} {'ratio/mean':>11s} {'npts':>7s}")
order = np.argsort(w2)
for q in (0.10, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001):
    k = max(1, int(q * N))
    idx = order[-k:]
    val = P[idx].sum() / (den[idx].sum() + 1e-30)
    flag = "  <- NOISE" if k < 30 else ""
    print(f"  {q*100:9.3f}% {val:10.4f} {val/mean:11.3f} {k:7d}{flag}")
print(f"  fraction NEGATIVE production (w.Sw<0) = {(P < 0).mean():.3f}")
k1 = max(1, int(0.01 * N))
print(f"  top 1% enstrophy carries {P[order[-k1:]].sum()/P.sum()*100:.1f}% of total production")
print(f"  |w|^2: mean={w2.mean():.1f}  max={w2.max():.1f}  (intermittency: max/mean={w2.max()/w2.mean():.0f})")
print("\n  READ: ratio/mean falling into the tail = regularizing (extremes MORE depleted).")
print("        ratio/mean rising above 1 in the tail = extreme DE-depletion (singularity-relevant).")
