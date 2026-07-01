r"""ns_jhtdb.py — run OUR geometric-depletion diagnostics on real high-Re DNS (Johns Hopkins Turbulence DB).

The frontier route: instead of simulating high Re ourselves (resolution wall), analyze turbulence JHU
already simulated. isotropic1024coarse is forced isotropic turbulence at Re_lambda ~ 433 — well past our
Taylor-Green Re=3200. Query the velocity GRADIENT at K random points, form S and omega, and measure the
geometric depletion  <omega.S.omega>/<|omega|^2|S|>  (mean + extreme tail), exactly as in D2-D5.

NOTE: the public 'testing' token is rate-limited to small queries; a robust extreme-tail (top 0.01% needs
~1e6 points) requires a free registered token + batched download. This proves the pipeline.

  python ns_jhtdb.py [K]
"""
import sys, numpy as np

K = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
TOKEN = "edu.jhu.pha.turbulence.testing-201406"
DS = "isotropic1024coarse"

from givernylocal.turbulence_dataset import turb_dataset
from givernylocal.turbulence_toolkit import getData

ds = turb_dataset(dataset_title=DS, output_path="./_giv", auth_token=TOKEN)
L = 2 * np.pi
pts = (np.random.RandomState(0).rand(K, 3) * L)               # float64
print(f"querying {K} velocity-gradient points from {DS} (Re_lambda~433)...", flush=True)
res = getData(ds, "velocity", 0.0, "none", "fd4lag4", "gradient", pts)   # velocity gradient via operator
arr = np.asarray(getattr(res, "result", res))                 # (K,9) duxdx,duxdy,duxdz,duydx,...
if arr.dtype == object or arr.ndim == 1:
    arr = np.asarray([np.asarray(r) for r in res]).reshape(K, -1)
A = arr.reshape(K, 3, 3)                                       # A[:,i,j] = d u_i / d x_j

S = 0.5 * (A + A.transpose(0, 2, 1))          # strain
w = np.stack([A[:, 2, 1] - A[:, 1, 2], A[:, 0, 2] - A[:, 2, 0], A[:, 1, 0] - A[:, 0, 1]], -1)  # vorticity
w2 = (w * w).sum(-1)
P = np.einsum("ki,kij,kj->k", w, S, w)        # omega . S . omega (enstrophy production, per point)
Smag = np.sqrt(2 * (S * S).sum((-1, -2)))
den = w2 * Smag

mean_depl = P.mean() / (den.mean() + 1e-30)
print(f"\n  geometric depletion on real Re_lambda~433 DNS ({K} points):")
print(f"    MEAN  <w.Sw>/<|w|^2|S|> = {mean_depl:.3f}")
for q in [0.05, 0.01]:
    k = max(1, int(q * K)); idx = np.argsort(w2)[-k:]
    print(f"    top {q*100:.0f}% |w| tail        = {P[idx].sum()/(den[idx].sum()+1e-30):.3f}  ({k} pts)")
print(f"\n  vs our Taylor-Green DNS: mean ~0.10 across Re=400-3200.")
print("  (small-sample tail is noisy; robust extreme tail needs a registered token + ~1e6 pts.)")
