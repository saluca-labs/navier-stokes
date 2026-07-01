r"""fd_vs_spectral.py — FD-vs-spectral cross-check of the geometric-depletion tail (run in SciServer).

Paired comparison on identical grid points: for each box we compute the depletion (w2, P=omega.S.omega,
den=|omega|^2|S|) two ways and accumulate both pools, point-matched:
  (a) FD path     — 4th-order central differences on a velocity cutout (what sciserver_cutout.py does).
  (b) SPECTRAL    — JHTDB spectral velocity-gradient (getData) at the SAME interior grid points.

Because it is paired, it detects an FD-vs-tail bias with far less data than an absolute extreme-tail
measurement (the testing-token 4096/call cap is fine: a 16^3 box has 12^3=1728 interior pts < 4096).

The giverny point-gradient API is the one unknown; this script PROBES it first (prints the getData
signature + tries one call) and tells you exactly what to fix if the guess is wrong.

  python fd_vs_spectral.py --probe            # just confirm the spectral-gradient API works
  python fd_vs_spectral.py --nbox 40          # full paired comparison over 40 boxes (~70k pts each path)
"""
from __future__ import annotations
import os, sys, json, argparse, inspect
import numpy as np
from giverny.turbulence_dataset import turb_dataset
from giverny.turbulence_toolkit import getCutout
try:
    from giverny.turbulence_toolkit import getData
except Exception:
    getData = None

DS = os.environ.get("DATASET", "isotropic8192")
N = {"isotropic8192": 8192, "isotropic32768": 32768, "isotropic1024coarse": 1024}.get(DS, 8192)
B = 16
dx = 2 * np.pi / N
TOKEN = os.environ.get("JHTDB_TOKEN", "edu.jhu.pha.turbulence.testing-201406")
# isotropic8192: use only the Re_lambda~1300 snapshots (giverny times 1-5; time 6 = snap 5 = Re~610)
TIMES = {"isotropic8192": [1, 2, 3, 4, 5]}.get(DS, [1])
ds = turb_dataset(dataset_title=DS, output_path="/tmp/_giv", auth_token=TOKEN)


def d4(f, ax):
    return (-np.roll(f, -2, ax) + 8 * np.roll(f, -1, ax) - 8 * np.roll(f, 1, ax) + np.roll(f, 2, ax)) / (12 * dx)


def cutout_vel(t, c):
    xyzt = np.array([[c[0], c[0] + B - 1], [c[1], c[1] + B - 1], [c[2], c[2] + B - 1], [t, t]], dtype=np.int64)
    res = getCutout(ds, "velocity", xyzt, np.array([1, 1, 1, 1], dtype=np.int64))
    da = res[list(res.data_vars)[0]]
    return np.squeeze(np.asarray(da.values, np.float64)).transpose(2, 1, 0, 3)   # (x,y,z,comp)


def depl_from_A(A):
    """A[...,i,j] = du_i/dx_j -> (w2, P=omega.S.omega, den=|omega|^2|S|)."""
    S = 0.5 * (A + np.swapaxes(A, -1, -2))
    w = np.stack([A[..., 2, 1] - A[..., 1, 2], A[..., 0, 2] - A[..., 2, 0], A[..., 1, 0] - A[..., 0, 1]], -1)
    w2 = (w * w).sum(-1)
    P = np.einsum("...i,...ij,...j->...", w, S, w)
    Sm = np.sqrt(2 * (S * S).sum((-1, -2)))
    return w2, P, w2 * Sm


def fd_box(t, c):
    fld = cutout_vel(t, c)                         # (B,B,B,3) in (x,y,z,comp)
    A = np.empty(fld.shape[:3] + (3, 3))
    for i in range(3):
        for j in range(3):
            A[..., i, j] = d4(fld[..., i], j)
    A = A[2:-2, 2:-2, 2:-2]                         # interior (B-4)^3
    w2, P, den = depl_from_A(A)
    g = np.arange(2, B - 2)                         # interior offsets; grid index (1-based) = c[ax] + g
    gx = (c[0] + g - 1) * dx; gy = (c[1] + g - 1) * dx; gz = (c[2] + g - 1) * dx
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")   # axes match A's (x,y,z) order
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], -1).astype(np.float64)
    return w2.ravel(), P.ravel(), den.ravel(), pts


def _grad_to_A(arr, npts):
    """Reshape a getData gradient return (npts x 9) into A[...,i,j]=du_i/dx_j.
    JHTDB order: [duxdx,duxdy,duxdz, duydx,duydy,duydz, duzdx,duzdy,duzdz]."""
    a = np.asarray(arr, np.float64)
    a = a.reshape(npts, 9) if a.size == npts * 9 else a.reshape(npts, -1)[:, :9]
    return a.reshape(npts, 3, 3)


# confirmed giverny convention: getData(ds,'velocity',t,'none',<method>,'gradient',pts).
# valid isotropic8192 gradient methods: m1q4, m2q8, fd4noint, fd6noint, fd8noint, fd4lag4.
# at grid points the 'noint' methods are exact FD (no interpolation); fd8noint (8th order) is the
# discriminating high-order reference; fd4noint (4th order) should ~match our own on-grid FD (sanity).
METHOD = os.environ.get("GRAD_METHOD", "fd8noint")


def jhtdb_A_at(pts, t, method=None):
    """JHTDB velocity gradient at points (npts,3) at time t via spatial method `method` -> A (npts,3,3)."""
    method = method or METHOD
    res = getData(ds, "velocity", t, "none", method, "gradient", pts)
    arr = res
    if hasattr(arr, "data_vars"):
        arr = arr[list(arr.data_vars)[0]]
    if hasattr(arr, "values"):
        arr = arr.values
    return _grad_to_A(arr, pts.shape[0])


def report(tag, w2, P, den):
    w2 = np.asarray(w2); P = np.asarray(P); den = np.asarray(den)
    mean = P.sum() / (den.sum() + 1e-30)
    o = np.argsort(w2)
    line = f"  {tag:9} N={len(w2)} Dmean={mean:.4f}"
    for q in (0.05, 0.01, 0.001, 0.0001):
        k = max(1, int(q * len(w2)))
        line += f"  top{q*100:g}%={P[o[-k:]].sum()/(den[o[-k:]].sum()+1e-30)/mean:.3f}"
    print(line, flush=True)
    return mean


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--nbox", type=int, default=40)
    a = ap.parse_args()
    if getData is None:
        print("getData not importable from giverny.turbulence_toolkit — check the package."); sys.exit(2)
    print(f"dataset={DS} N={N}  getData signature: {inspect.signature(getData)}", flush=True)

    print(f"reference JHTDB method = {METHOD} (ours = on-grid 4th-order FD)", flush=True)
    if a.probe:
        c = np.array([100, 200, 300]); t = TIMES[0]
        _, _, _, pts = fd_box(t, c)
        try:
            A = jhtdb_A_at(pts[:3], t)
            print(f"PROBE OK ({METHOD})\n  A[0]=\n{A[0]}")
        except Exception as e:
            print(f"PROBE FAILED ({METHOD}): {type(e).__name__}: {str(e)[:200]}\n"
                  f"  getData signature: {inspect.signature(getData)}")
        return

    rng = np.random.RandomState(2024)
    W2f = []; Pf = []; Df = []; W2s = []; Ps = []; Ds = []; corr = []
    for n in range(a.nbox):
        t = TIMES[n % len(TIMES)]
        c = rng.randint(1, N - B, size=3)
        try:
            w2f, pf, df, pts = fd_box(t, c)
            A = jhtdb_A_at(pts, t)
            w2s, ps, dss = depl_from_A(A)
        except Exception as e:
            print(f"box {n} FAILED: {type(e).__name__}: {str(e)[:160]}"); break
        W2f += list(w2f); Pf += list(pf); Df += list(df)
        W2s += list(w2s); Ps += list(ps); Ds += list(dss)
        r = np.corrcoef(w2f, w2s)[0, 1] if len(w2f) > 2 else float("nan")
        corr.append(r)
        if (n + 1) % 10 == 0:
            print(f"  box {n+1}/{a.nbox}  per-box corr(w2_FD,w2_spec)={r:.4f}", flush=True)
    if not W2f:
        return
    print(f"\n=== differentiation-order convergence ({DS}): our 4th-order FD vs JHTDB {METHOD} ===")
    print(f"  point-match sanity: mean per-box corr(w2_FD, w2_ref) = {np.nanmean(corr):.4f} (want ~1.0)")
    mf = report("FD-4th", W2f, Pf, Df)
    ms = report(f"JHTDB-{METHOD}", W2s, Ps, Ds)
    print(f"  Dmean ratio FD/ref = {mf/ms:.3f}")
    print("  -> if the top-q% ratios MATCH, the extreme-tail gradients are order-converged (on-grid 4th-order")
    print("     FD is adequate; the paper's gradient caveat can be removed). If they DIFFER, the tail is")
    print("     under-resolved by 4th-order FD and the higher-order reference should be used.")


if __name__ == "__main__":
    main()
