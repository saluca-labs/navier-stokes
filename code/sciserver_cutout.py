#!/usr/bin/env python
r"""sciserver_cutout.py — FAST geometric-depletion on high-Re JHTDB via CONTIGUOUS CUTOUTS (run in SciServer).

Random-point getData on the 8192^3 / 32768^3 cubes is pathological (random I/O ~3.5 min / 4096 pts). The
efficient access for a locally-mounted DNS is a CONTIGUOUS sub-volume: one sequential read of a B^3 box gives
~B^3 grid points at once. We then compute the velocity-gradient tensor A_ij = du_i/dx_j ON the grid with
4th-order central differences (interior only), form S=sym(A), omega=curl, and the geometric depletion
<w.Sw>/<|w|^2|S|> — mean + extreme tails — exactly the D2-D5 / sciserver_depletion diagnostic, but at millions
of points per box instead of 4096 per slow call.

  python sciserver_cutout.py            # defaults: isotropic8192, B=128, 8 boxes
  DATASET=isotropic32768 B=256 NBOX=4 python sciserver_cutout.py
  DATASET=isotropic8192 B=64 NBOX=1 python sciserver_cutout.py     # <-- quick probe first (validate + time)

Durable: appends each box's (w2,P,den) to pool_cutout_<dataset>.npz in the persistent dir after every box.
"""
import os, sys, time, numpy as np

DATASET = os.environ.get("DATASET", "isotropic8192")
B       = int(os.environ.get("B", "128"))      # cube edge in grid cells
NBOX    = int(os.environ.get("NBOX", "8"))     # how many boxes to sample
TOKEN   = os.environ.get("JHTDB_TOKEN", "edu.jhu.pha.turbulence.testing-201406")
PERSIST = "/home/idies/workspace/Storage/saluca/persistent"
if not os.path.isdir(PERSIST):
    PERSIST = os.getcwd()
OUT = os.path.join(PERSIST, "depletion_out"); os.makedirs(OUT, exist_ok=True)
SCRATCH = "/home/idies/workspace/Temporary/saluca/scratch"
WORK = SCRATCH if os.path.isdir(SCRATCH) else os.path.join(OUT, "_giv"); os.makedirs(WORK, exist_ok=True)

NGRID = {"isotropic1024coarse": 1024, "isotropic1024fine": 1024, "isotropic4096": 4096,
         "isotropic8192": 8192, "isotropic32768": 32768}.get(DATASET, 8192)
# isotropic8192: getCutout time t -> snapshot (t-1) via JHTDB timeIndexShift=-1; snapshots 0-4 (times 1-5)
# are Re_lambda~1200-1300, snapshot 5 (time 6) is Re_lambda~610 -> EXCLUDE time 6 for a clean high-Re sample.
TIMES = {"isotropic8192": [1, 2, 3, 4, 5], "isotropic32768": [1], "isotropic4096": [1]}.get(
    DATASET, [round(float(t), 1) for t in np.linspace(0.0, 10.0, 11)])
LBOX = 2 * np.pi
dx = LBOX / NGRID

# package
pkg = None
for nm in ("giverny", "givernylocal"):
    try:
        __import__(nm); pkg = nm; break
    except Exception:
        continue
if pkg == "giverny":
    from giverny.turbulence_dataset import turb_dataset
    from giverny.turbulence_toolkit import getCutout
else:
    from givernylocal.turbulence_dataset import turb_dataset
    from givernylocal.turbulence_toolkit import getCutout
print(f"package={pkg}  dataset={DATASET}  N={NGRID}  dx={dx:.3e}  box=B^3={B}^3  nbox={NBOX}", flush=True)

ds = turb_dataset(dataset_title=DATASET, output_path=WORK, auth_token=TOKEN)


def d4(f, axis):
    """4th-order central derivative along axis (np.roll wraps; trim 2 cells/side after)."""
    return (-np.roll(f, -2, axis) + 8 * np.roll(f, -1, axis)
            - 8 * np.roll(f, 1, axis) + np.roll(f, 2, axis)) / (12 * dx)


def cutout_field(t, corner):
    """fetch B^3 velocity cutout at 1-based grid corner; return (B,B,B,3) float array, spatial axes
    aligned to velocity components. giverny: getCutout(cube,var,xyzt_axes_ranges,xyzt_strides).
    JHTDB cutouts return [(t,)z,y,x, comp] with comp in (x,y,z) order -> reverse the 3 spatial axes."""
    c = corner
    xyzt = np.array([[int(c[0]), int(c[0] + B - 1)],
                     [int(c[1]), int(c[1] + B - 1)],
                     [int(c[2]), int(c[2] + B - 1)],
                     [int(t), int(t)]], dtype=np.int64)    # giverny checks .dtype -> needs np int arrays
    strides = np.array([1, 1, 1, 1], dtype=np.int64)
    res = getCutout(ds, "velocity", xyzt, strides)        # returns xarray.Dataset dims (zcoor,ycoor,xcoor,values=3)
    dv = list(res.data_vars)                              # velocity lives in a single data variable
    da = res[dv[0]] if dv else res.to_array()
    arr = np.squeeze(np.asarray(da.values, dtype=np.float64))   # (z,y,x,comp), comp in (x,y,z) order
    print(f"    [cutout vars={dv} raw {da.shape} -> squeezed {arr.shape}]", flush=True)
    cax = [i for i, s in enumerate(arr.shape) if s == 3]
    if not cax:
        raise ValueError(f"no length-3 component axis in cutout shape {arr.shape}")
    arr = np.moveaxis(arr, cax[-1], -1)                   # (s0,s1,s2,3)
    if arr.ndim != 4:
        raise ValueError(f"unexpected cutout ndim {arr.ndim}, shape {arr.shape}")
    arr = arr.transpose(2, 1, 0, 3)                       # reverse spatial -> (x,y,z,comp)
    return arr


def depletion_of_box(field):
    """A_ij = du_i/dx_j via 4th-order FD; S, omega, depletion; trim 2-cell halo. returns w2,P,den (flat)."""
    A = np.empty(field.shape[:3] + (3, 3))
    for i in range(3):
        A[..., i, 0] = d4(field[..., i], 0)
        A[..., i, 1] = d4(field[..., i], 1)
        A[..., i, 2] = d4(field[..., i], 2)
    s = (slice(2, -2),) * 3
    A = A[s]                                   # valid interior (B-4)^3
    S = 0.5 * (A + np.swapaxes(A, -1, -2))
    w = np.stack([A[..., 2, 1] - A[..., 1, 2],
                  A[..., 0, 2] - A[..., 2, 0],
                  A[..., 1, 0] - A[..., 0, 1]], -1)
    w2 = (w * w).sum(-1)
    P = np.einsum("...i,...ij,...j->...", w, S, w)
    Smag = np.sqrt(2 * (S * S).sum((-1, -2)))
    return w2.ravel(), P.ravel(), (w2 * Smag).ravel()


def report(w2, P, den, tag):
    w2, P, den = map(np.asarray, (w2, P, den))
    mean = P.sum() / (den.sum() + 1e-30)
    order = np.argsort(w2)
    print(f"\n==== {tag}  N={len(w2)} ====")
    print(f"  MEAN = {mean:.4f}")
    for q in (0.05, 0.01, 0.001, 0.0001, 1e-5):
        k = max(1, int(q * len(w2))); v = P[order[-k:]].sum() / (den[order[-k:]].sum() + 1e-30)
        fl = "  <-NOISE" if k < 30 else ""
        print(f"  top {q*100:8.4f}% = {v:.4f}  ratio/mean={v/mean:5.3f}  ({k}){fl}")
    print(f"  neg-production frac = {(P < 0).mean():.3f}")
    return mean


pool = os.path.join(OUT, f"pool_cutout_{DATASET}.npz")
if os.path.exists(pool):
    z = np.load(pool); W2, PP, DN, nb = list(z["w2"]), list(z["P"]), list(z["den"]), int(z["nb"])
    print(f"resumed {len(W2)} pts, {nb} boxes")
else:
    W2, PP, DN, nb = [], [], [], 0

rng = np.random.RandomState(12345)
while nb < NBOX:
    t = TIMES[nb % len(TIMES)]
    corner = rng.randint(1, NGRID - B, size=3)            # grid-aligned, 1-based
    t0 = time.time()
    try:
        fld = cutout_field(t, corner)
        a, p, c = depletion_of_box(fld)
    except Exception as e:
        print(f"  box {nb} FAILED: {type(e).__name__}: {str(e)[:160]}")
        import traceback; traceback.print_exc(); break
    W2 += list(a); PP += list(p); DN += list(c); nb += 1
    np.savez(pool, w2=np.asarray(W2), P=np.asarray(PP), den=np.asarray(DN), nb=nb)
    m = np.asarray(PP); d = np.asarray(DN)
    print(f"  box {nb:2d} t={t} corner={corner.tolist()} ({time.time()-t0:5.1f}s) "
          f"+{len(a)} pts  N={len(W2)}  running_mean={m.sum()/(d.sum()+1e-30):.4f}", flush=True)

report(W2, PP, DN, f"{DATASET} CUTOUT FINAL (Re-high)")
print(f"\nsaved pool -> {pool}")
