r"""plot_depletion.py — the money figure: geometric-depletion tail/mean ratio vs quantile, per dataset/Re_lambda.

Each curve = how the depletion <w.Sw>/<|w|^2|S|> in the most-intense-vorticity tail compares to the bulk mean.
Below the dashed line at 1.0 => extremes MORE depleted than bulk (regularizing). Crossing ABOVE 1.0 in the
tail => extreme DE-depletion (Buaria-Pumir; singularity-relevant). The scientific claim is the TREND with
Re_lambda: 433 (regularizing) -> 1300 / largest-DNS (does it invert?).

  python plot_depletion.py POOL1.npz "label1@Re433" POOL2.npz "label2@Re1300" ...
"""
import sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QS = np.array([0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001])
args = sys.argv[1:]
pairs = [(args[i], args[i + 1]) for i in range(0, len(args) - 1, 2)]
if not pairs:
    print("usage: python plot_depletion.py POOL.npz LABEL [POOL2.npz LABEL2 ...]"); sys.exit(1)

fig, ax = plt.subplots(figsize=(8.5, 6))
for path, label in pairs:
    z = np.load(path)
    w2, P, den = z["w2"], z["P"], z["den"]
    N = len(w2); mean = P.sum() / (den.sum() + 1e-30)
    order = np.argsort(w2)
    ratio, xs = [], []
    for q in QS:
        k = max(1, int(q * N))
        if k < 20:            # don't plot bins too sparse to trust
            continue
        v = P[order[-k:]].sum() / (den[order[-k:]].sum() + 1e-30)
        ratio.append(v / mean); xs.append(q * 100)
    ax.plot(xs, ratio, "o-", lw=2, ms=6, label=f"{label}  (mean={mean:.3f}, N={N//1000}k)")

ax.axhline(1.0, ls="--", color="0.4", lw=1.3)
ax.text(ax.get_xlim()[0], 1.005, "de-depletion threshold (extremes LESS depleted than bulk)",
        fontsize=8, color="0.35", va="bottom")
ax.set_xscale("log")
ax.invert_xaxis()
ax.set_xlabel("vorticity percentile  (top X% of |omega|, into the extreme tail -->)")
ax.set_ylabel(r"depletion ratio  $\langle\omega\cdot S\omega\rangle_{tail}/\langle\omega\cdot S\omega\rangle_{mean}$")
ax.set_title("Geometric depletion in the extreme vorticity tail vs Reynolds number\n"
             "(below 1 = regularizing;  crossing above 1 = extreme de-depletion)")
ax.legend(loc="best", fontsize=9)
ax.grid(alpha=0.3)
out = "depletion_ratio_vs_Re.png"
fig.tight_layout(); fig.savefig(out, dpi=140)
print("wrote", out)
