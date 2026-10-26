r"""plot_crossRe.py — Fig 2: conditional geometric-depletion ratio D(q)/D_mean vs vorticity percentile for the
two high-Reynolds JHTDB datasets, with block-resampled 95% bands. Claimed depth = top-0.1%; deeper points are
drawn faded (resolution-limited, §3.2). Saves paper/depletion_crossRe.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Directory holding the two pool_cutout_*.npz files (set NS_POOL_DIR; defaults to the current directory).
DL = os.environ.get("NS_POOL_DIR", ".")
SETS = [
    (os.path.join(DL, "pool_cutout_isotropic8192(2).npz"),  r"isotropic8192  ($R_\lambda\approx1300$, $k_{max}\eta\approx1.4$)", "#1f77b4"),
    (os.path.join(DL, "pool_cutout_isotropic32768(2).npz"), r"isotropic32768 ($R_\lambda=2500$, $k_{max}\eta\approx2.8$)",     "#d62728"),
]
QS = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001]
CLAIM = 0.001  # top-0.1% : claimed depth; deeper = faded
K = 20


def block_ratio(w2, P, den, q, K):
    N = len(w2); edges = np.linspace(0, N, K + 1).astype(np.int64); vals = []
    for i in range(K):
        a, b = edges[i], edges[i + 1]
        Pb, denb, w2b = P[a:b], den[a:b], w2[a:b]
        mean = Pb.sum() / (denb.sum() + 1e-30)
        k = max(1, int(q * len(w2b)))
        o = np.argpartition(w2b, -k)[-k:]
        vals.append((Pb[o].sum() / (denb[o].sum() + 1e-30)) / mean)
    v = np.array(vals); sem = v.std(ddof=1) / np.sqrt(K)
    return v.mean(), 1.96 * sem


fig, ax = plt.subplots(figsize=(7.6, 5.0))
for path, label, c in SETS:
    z = np.load(path); w2, P, den = z["w2"], z["P"], z["den"]
    m = np.array([block_ratio(w2, P, den, q, K) for q in QS])
    mean, half = m[:, 0], m[:, 1]
    x = np.array(QS) * 100
    solid = np.array(QS) >= CLAIM
    ax.fill_between(x, mean - half, mean + half, color=c, alpha=0.18, lw=0)
    ax.plot(x[solid], mean[solid], "-o", color=c, lw=2, ms=5, label=label)
    ax.plot(x[~solid], mean[~solid], ":o", color=c, lw=1.3, ms=4, mfc="white", alpha=0.7)
    del w2, P, den, z
ax.axhline(1.0, ls="--", c="0.4", lw=1)
ax.text(0.06, 1.005, "de-depletion threshold", color="0.4", fontsize=8, va="bottom")
ax.axvspan(QS[-2] * 100, QS[0] * 100 * 1.6, alpha=0)  # padding
ax.axvline(CLAIM * 100, color="0.7", ls=":", lw=1)
ax.text(CLAIM * 100 * 1.25, 0.715, "claimed depth\n(top-0.1%)", fontsize=7.5, color="0.45", va="bottom", ha="right")
ax.text(0.00018 * 100, 0.80, "resolution-\nlimited (§3.2)", fontsize=7.5, color="0.55", va="center", ha="center")
ax.set_xscale("log"); ax.invert_xaxis()
ax.set_xlabel(r"vorticity percentile: top-$q\%$ of $|\omega|^2$  (deeper tail $\rightarrow$)")
ax.set_ylabel(r"$D(q)/D_{\mathrm{mean}}$  (ratio $<1$ = depleted / regularizing)")
ax.set_title(r"Extreme-tail geometric depletion weakens with Reynolds number")
ax.set_ylim(0.70, 1.04); ax.legend(loc="lower left", fontsize=8.5, framealpha=0.9)
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper", "depletion_crossRe.png")
fig.savefig(out, dpi=150)
print("saved", out)
