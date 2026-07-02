#!/usr/bin/env python
r"""coarse_grain_run.py — the RESOLUTION-BIAS experiment (a-priori coarse-graining test) on isotropic32768.

WHY (verification pass 2026-07-01, paper punch-list item 5): the v0.5 cross-Re de-depletion claim rests on an
ASSERTED bias direction — "under-resolution pushes D(q)/D_mean toward 1" — which a skeptic argued is plausibly
BACKWARDS (filtered fields align MORE with the intermediate eigenvector -> ratio biased LOW -> the resolution
confound would run WITH the observed 0.87->0.93 trend). This experiment measures the bias sign and size
directly, with data already in hand: take isotropic32768 (k_max*eta~2.8), spectrally coarse-grain the SAME
velocity fields to the isotropic8192-equivalent resolution (k_max*eta~1.4), and re-measure the depletion
statistic PAIRED per box.

DESIGN (per box, all variants on the SAME assembled field -> paired differences, K boxes as blocks):
  V0  unfiltered, 4th-order FD on the native grid              (baseline; cross-check vs published pool)
  V1  sharp spectral low-pass at k_c = KCUT_FRAC_NYQ * k_Nyq, FD on the FINE grid   (spectral-content effect only)
  V2  same sharp filter, then downsample x2, FD on the COARSE grid (dx'=2dx)        (full i8192-emulation: content + stencil sampling)
  V3  Gaussian filter (transfer exp(-k^2 Delta^2 / 24), Delta = pi/k_c), downsample x2, FD coarse  (filter-shape robustness)

  Default KCUT_FRAC_NYQ = sqrt(2)/3 ~ 0.4714: Yeung PSDNS convention k_max = (sqrt2/3) N dk, so halving the
  resolved k_max*eta (2.8 -> 1.4) means k_c = k_max/2 = (sqrt2/3)*(N/2)*dk*... = 0.4714 * k_Nyquist.
  (If the 2/3-rule convention k_max=(2/3)k_Nyq is preferred, set KCUT_FRAC_NYQ=0.3333 — the bias SIGN
  question does not hinge on the convention; both halve the resolved range.)

EDGE HANDLING: cutouts are windows of a periodic global field, NOT periodic themselves; FFT filtering wraps.
We therefore assemble a REGION^3 volume (default 256^3 from 2x2x2 fetches of B_FETCH=128) and trim a MARGIN
(default 48 fine cells/side) AFTER filtering+FD, so every retained point is >= MARGIN cells from any window
edge (sharp-filter ringing decays ~1/(k_c d); at d=48, contamination is negligible for tail statistics; the
FD stencil halo (2 cells fine / 2 coarse=4 fine) is inside the margin).

STATISTICS (lessons from the verification pass are BAKED IN):
  - per-box tail re-selection (top-q of |w|^2 within the box) and per-box D_mean -> per-box ratio r_v(q)
  - paired differences Delta_b(q) = r_V(q) - r_V0(q) per box; report mean +/- t_{K-1} * SEM  (t, NOT z)
  - pooled estimates also reported; box count K pinned in the output
  - DECISION LINE: sign of Delta at top-1% / top-0.1%:
      Delta > 0  -> coarsening RAISES the ratio toward 1  -> paper's "conservative direction" CONFIRMED
                    (the i8192 point would read artificially HIGH; true cross-Re de-depletion is >= observed)
      Delta < 0  -> coarsening LOWERS the ratio           -> skeptic's reversal CONFIRMED
                    (the confound runs WITH the observed trend; cross-Re attribution is compromised)

RUN (in the SciServer container `_salucaTurbulence`, in the dir holding this file):
  export JHTDB_TOKEN=<real token>          # vault: JHTDB_TOKEN
  # quick probe first (~1 box, small region, validates fetch+filter+shapes end to end):
  NBOX=1 REGION=128 B_FETCH=64 MARGIN=24 python coarse_grain_run.py
  # production (resumable; ~8 fetches/box; kill+restart safe):
  nohup python coarse_grain_run.py > coarse_grain.log 2>&1 &
  tail -f coarse_grain.log
  # analysis only (no fetching), e.g. after downloading the npz pools locally:
  python coarse_grain_run.py --analyze

Durable: per-variant pools  pool_cgrain_<VARIANT>_isotropic32768.npz  (w2, P, den, nb, per-box offsets)
checkpointed after every box in the persistent dir; plus CGRAIN_RESULTS.md appended per box.
If B_FETCH=128 hangs (known intermittent gotcha on the 32768 cube), restart with B_FETCH=64 — resume skips
completed boxes; the region is assembled from (REGION/B_FETCH)^3 fetches either way.

LOCAL SELF-TEST (no JHTDB; synthetic periodic field, checks machinery + FKEEP=1 no-op equivalence):
  python coarse_grain_run.py --selftest
"""
import os, sys, time, json
import numpy as np

# ---------------- config ----------------
DATASET   = os.environ.get("DATASET", "isotropic32768")
REGION    = int(os.environ.get("REGION", "256"))       # assembled cube edge (fine cells)
B_FETCH   = int(os.environ.get("B_FETCH", "128"))      # per-getCutout edge (128 verified on 32768; 64 fallback)
NBOX      = int(os.environ.get("NBOX", "24"))          # number of assembled regions
MARGIN    = int(os.environ.get("MARGIN", "48"))        # fine cells trimmed per side after filter+FD
KCUT      = float(os.environ.get("KCUT_FRAC_NYQ", str(np.sqrt(2.0) / 3.0)))   # k_c / k_Nyquist
TOKEN     = os.environ.get("JHTDB_TOKEN", "edu.jhu.pha.turbulence.testing-201406")
SEED      = int(os.environ.get("SEED", "424242"))
QS        = (0.05, 0.01, 0.001)                        # top-5%, 1%, 0.1% (0.01% only printed, never claimed)

NGRID = {"isotropic8192": 8192, "isotropic32768": 32768}.get(DATASET, 32768)
TIMES = {"isotropic8192": [1, 2, 3, 4, 5], "isotropic32768": [1]}[DATASET]
LBOX = 2 * np.pi
DX = LBOX / NGRID

PERSIST = "/home/idies/workspace/Storage/saluca/persistent"
if not os.path.isdir(PERSIST):
    PERSIST = os.getcwd()
OUT = os.path.join(PERSIST, "depletion_out"); os.makedirs(OUT, exist_ok=True)
SCRATCH = "/home/idies/workspace/Temporary/saluca/scratch"
WORK = SCRATCH if os.path.isdir(SCRATCH) else os.path.join(OUT, "_giv"); os.makedirs(WORK, exist_ok=True)
REPORT = os.path.join(OUT, "CGRAIN_RESULTS.md")

VARIANTS = ("V0_unfiltered", "V1_sharp_fine", "V2_sharp_ds2", "V3_gauss_ds2")


# ---------------- numerics ----------------
def d4(f, axis, dx):
    """4th-order central derivative (np.roll wraps; caller must trim contaminated halo)."""
    return (-np.roll(f, -2, axis) + 8 * np.roll(f, -1, axis)
            - 8 * np.roll(f, 1, axis) + np.roll(f, 2, axis)) / (12 * dx)


def stats_of_field(field, dx, trim):
    """velocity (n,n,n,3) -> pointwise (w2, P, den) on the interior after trimming `trim` cells/side.
    Identical statistic to sciserver_cutout.depletion_of_box (Frobenius |S|=sqrt(2 SijSij))."""
    A = np.empty(field.shape[:3] + (3, 3))
    for i in range(3):
        for j in range(3):
            A[..., i, j] = d4(field[..., i], j, dx)
    s = (slice(trim, -trim),) * 3
    A = A[s]
    S = 0.5 * (A + np.swapaxes(A, -1, -2))
    w = np.stack([A[..., 2, 1] - A[..., 1, 2],
                  A[..., 0, 2] - A[..., 2, 0],
                  A[..., 1, 0] - A[..., 0, 1]], -1)
    w2 = (w * w).sum(-1)
    P = np.einsum("...i,...ij,...j->...", w, S, w)
    Smag = np.sqrt(2 * (S * S).sum((-1, -2)))
    return w2.ravel(), P.ravel(), (w2 * Smag).ravel()


def spectral_filter(field, frac_nyq, kind):
    """low-pass the (n,n,n,3) field: 'sharp' spherical cutoff at k_c=frac_nyq*k_Nyq, or 'gauss'
    (transfer exp(-k^2 Delta^2/24), Delta=pi/k_c — the standard Gaussian LES filter). FFT treats the
    window as periodic; caller trims MARGIN to discard wrap-contaminated edges."""
    n = field.shape[0]
    k1 = np.fft.fftfreq(n, d=1.0 / n)          # integer wavenumbers, Nyquist = n/2
    kx, ky, kz = np.meshgrid(k1, k1, k1[:n // 2 + 1], indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    kc = frac_nyq * (n / 2.0)
    if kind == "sharp":
        H = (k2 <= kc * kc).astype(np.float64)
    else:
        delta2 = (np.pi / kc) ** 2              # in grid-index units (dx=1 internally; ratio-invariant)
        H = np.exp(-k2 * delta2 / 24.0)
    out = np.empty_like(field)
    ax = (0, 1, 2)
    for c in range(3):
        out[..., c] = np.fft.irfftn(np.fft.rfftn(field[..., c], axes=ax) * H, s=field.shape[:3], axes=ax)
    return out


def per_box_ratios(w2, P, den, qs):
    """per-box D_mean and tail ratios r(q)=D(q)/D_mean with within-box tail re-selection.
    (sums in float64 regardless of pool storage dtype)"""
    dmean = P.sum(dtype=np.float64) / (den.sum(dtype=np.float64) + 1e-30)
    order = np.argsort(w2)
    out = {"D_mean": float(dmean)}
    for q in qs:
        k = max(1, int(q * len(w2)))
        idx = order[-k:]
        dq = P[idx].sum(dtype=np.float64) / (den[idx].sum(dtype=np.float64) + 1e-30)
        out[f"r_top{q*100:g}%"] = float(dq / dmean)
        out[f"n_top{q*100:g}%"] = int(k)
    return out


# ---------------- data access ----------------
_ds = None
def dataset():
    global _ds
    if _ds is None:
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
        _ds = (turb_dataset(dataset_title=DATASET, output_path=WORK, auth_token=TOKEN), getCutout, pkg)
        print(f"package={pkg} dataset={DATASET} N={NGRID}", flush=True)
    return _ds


def fetch_block(t, corner, b):
    """one getCutout of b^3 velocity at 1-based corner -> (b,b,b,3), axes (x,y,z,comp).
    (identical conventions to sciserver_cutout.cutout_field)"""
    ds, getCutout, _ = dataset()
    c = corner
    xyzt = np.array([[int(c[0]), int(c[0] + b - 1)],
                     [int(c[1]), int(c[1] + b - 1)],
                     [int(c[2]), int(c[2] + b - 1)],
                     [int(t), int(t)]], dtype=np.int64)
    strides = np.array([1, 1, 1, 1], dtype=np.int64)
    res = getCutout(ds, "velocity", xyzt, strides)
    dv = list(res.data_vars)
    da = res[dv[0]] if dv else res.to_array()
    arr = np.squeeze(np.asarray(da.values, dtype=np.float64))
    cax = [i for i, s in enumerate(arr.shape) if s == 3]
    arr = np.moveaxis(arr, cax[-1], -1)
    return arr.transpose(2, 1, 0, 3)            # reverse spatial (z,y,x)->(x,y,z)


def assemble_region(t, corner):
    """REGION^3 velocity assembled from (REGION/B_FETCH)^3 contiguous fetches."""
    m = REGION // B_FETCH
    assert m * B_FETCH == REGION, "REGION must be a multiple of B_FETCH"
    fld = np.empty((REGION, REGION, REGION, 3))
    for i in range(m):
        for j in range(m):
            for k in range(m):
                sub = corner + np.array([i, j, k]) * B_FETCH
                t0 = time.time()
                fld[i * B_FETCH:(i + 1) * B_FETCH,
                    j * B_FETCH:(j + 1) * B_FETCH,
                    k * B_FETCH:(k + 1) * B_FETCH] = fetch_block(t, sub, B_FETCH)
                print(f"    fetch {i}{j}{k} corner={sub.tolist()} ({time.time()-t0:.1f}s)", flush=True)
    return fld


# ---------------- variants ----------------
def compute_variants(fld):
    """all four variant statistics on one assembled region. Returns {variant: (w2,P,den)}."""
    out = {}
    # V0 unfiltered, fine grid
    out["V0_unfiltered"] = stats_of_field(fld, DX, MARGIN)
    # sharp filter once, reuse for V1+V2
    fs = spectral_filter(fld, KCUT, "sharp")
    out["V1_sharp_fine"] = stats_of_field(fs, DX, MARGIN)
    ds2 = fs[::2, ::2, ::2]                               # coarse grid, spacing 2*DX
    out["V2_sharp_ds2"] = stats_of_field(ds2, 2 * DX, MARGIN // 2)
    fg = spectral_filter(fld, KCUT, "gauss")
    out["V3_gauss_ds2"] = stats_of_field(fg[::2, ::2, ::2], 2 * DX, MARGIN // 2)
    return out


# ---------------- persistence ----------------
def pool_path(variant):
    return os.path.join(OUT, f"pool_cgrain_{variant}_{DATASET}.npz")


def load_pool(variant):
    """pools hold float32 point data (halves disk/RAM; the statistic is a float64-summed ratio,
    per-point 1e-7 relative storage error is irrelevant at the 1e-3 level we resolve).
    In-memory representation: lists of ndarray CHUNKS (one per box), concatenated on save."""
    p = pool_path(variant)
    if os.path.exists(p):
        z = np.load(p)
        return [z["w2"]], [z["P"]], [z["den"]], int(z["nb"]), list(z["offsets"])
    return [], [], [], 0, [0]


def save_pool(variant, W2, PP, DN, nb, offsets):
    np.savez(pool_path(variant),
             w2=np.concatenate(W2).astype(np.float32, copy=False),
             P=np.concatenate(PP).astype(np.float32, copy=False),
             den=np.concatenate(DN).astype(np.float32, copy=False),
             nb=nb, offsets=np.asarray(offsets))


def append_report(md):
    with open(REPORT, "a", encoding="utf-8") as f:
        f.write(md + "\n")


# ---------------- analysis ----------------
def analyze(pools=None, echo=print):
    """paired per-box analysis across variants. pools: {variant:(W2,PP,DN,nb,offsets)} or load from disk."""
    from math import sqrt
    try:
        from scipy import stats as _st
        tmult = lambda k: float(_st.t.ppf(0.975, k - 1)) if k > 1 else float("nan")
    except Exception:
        # t 97.5% quantiles fallback (df 1..30)
        T = [12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.36, 2.31, 2.26, 2.23, 2.20, 2.18, 2.16, 2.14,
             2.13, 2.12, 2.11, 2.10, 2.09, 2.09, 2.08, 2.07, 2.07, 2.06, 2.06, 2.06, 2.05, 2.05, 2.05, 2.04]
        tmult = lambda k: T[min(k - 2, len(T) - 1)] if k > 1 else float("nan")

    if pools is None:
        pools = {v: load_pool(v) for v in VARIANTS}
    nb = min(p[3] for p in pools.values())
    if nb < 1:
        echo("no boxes yet"); return None
    echo(f"\n==== PAIRED ANALYSIS  ({DATASET}, K={nb} boxes, KCUT_FRAC_NYQ={KCUT:.4f}, "
         f"REGION={REGION}, MARGIN={MARGIN}) ====")
    cat = {v: tuple(np.concatenate(pools[v][i]) for i in range(3)) for v in VARIANTS}
    perbox = {v: [] for v in VARIANTS}
    for v in VARIANTS:
        W2, PP, DN = cat[v]
        off = pools[v][4]
        for b in range(nb):
            s = slice(off[b], off[b + 1])
            perbox[v].append(per_box_ratios(W2[s], PP[s], DN[s], QS))
    # pooled + per-box summaries
    for v in VARIANTS:
        W2, PP, DN = cat[v]
        off = pools[v][4]
        pooled = per_box_ratios(W2[:off[nb]], PP[:off[nb]], DN[:off[nb]], QS)
        echo(f"\n{v}:  N={off[nb]}  ({nb} boxes)  pooled D_mean={pooled['D_mean']:.4f}")
        for q in QS:
            key = f"r_top{q*100:g}%"
            r = np.array([pb[key] for pb in perbox[v]])
            if nb > 1:
                tm = tmult(nb)
                echo(f"  {key:<12} pooled={pooled[key]:.4f}   box-mean={r.mean():.4f} "
                     f"+/- {tm*r.std(ddof=1)/sqrt(nb):.4f} (t_{nb-1} 95%)")
            else:
                echo(f"  {key:<12} pooled={pooled[key]:.4f}   box-mean={r.mean():.4f} (K=1, no interval)")
    # paired differences vs V0 — THE decision numbers
    echo("\n---- PAIRED DIFFERENCES (variant - V0, same boxes; positive = coarsening RAISES ratio toward 1) ----")
    verdicts = {}
    for v in VARIANTS[1:]:
        for q in QS[1:]:                                   # top-1%, top-0.1%
            key = f"r_top{q*100:g}%"
            d = np.array([perbox[v][b][key] - perbox["V0_unfiltered"][b][key] for b in range(nb)])
            tm = tmult(nb)
            m, hw = d.mean(), tm * d.std(ddof=1) / sqrt(nb)
            sig = "RESOLVED" if abs(m) > hw else "not resolved"
            echo(f"  {v:<14} Delta {key:<10} = {m:+.4f} +/- {hw:.4f}  [{sig}]")
            verdicts[(v, key)] = (m, hw)
    # decision line (primary = V2, the full i8192 emulation)
    for q in QS[1:]:
        key = f"r_top{q*100:g}%"
        m, hw = verdicts[("V2_sharp_ds2", key)]
        if not (np.isfinite(hw) and abs(m) > hw):    # nan hw (K=1) must land here, not in a direction call
            call = "INCONCLUSIVE at this K"
        elif m > 0:
            call = "paper's CONSERVATIVE direction confirmed (coarse reads LESS depleted)"
        else:
            call = "skeptic's REVERSAL confirmed (coarse reads MORE depleted; confound runs WITH the cross-Re trend)"
        echo(f"DECISION {key}: Delta_V2 = {m:+.4f} +/- {hw:.4f}  -> {call}")
    return verdicts


# ---------------- self-test (local, synthetic) ----------------
def selftest():
    """no JHTDB: synthetic band-limited periodic field; checks shapes, filter no-op at KCUT=1 (V1==V0 on a
    TRULY periodic field), monotone spectral truncation, and paired-analysis plumbing."""
    global REGION, MARGIN, KCUT
    REGION, MARGIN = 64, 12
    rng = np.random.RandomState(7)
    n = REGION
    k1 = np.fft.fftfreq(n, d=1.0 / n)
    kx, ky, kz = np.meshgrid(k1, k1, k1[:n // 2 + 1], indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    spec_shape = np.exp(-k2 / (2 * 6.0 ** 2))              # band-limited ~ k<18
    fld = np.empty((n, n, n, 3))
    for c in range(3):
        ph = rng.randn(n, n, n // 2 + 1) + 1j * rng.randn(n, n, n // 2 + 1)
        fld[..., c] = np.fft.irfftn(ph * spec_shape, s=(n, n, n))
        fld[..., c] /= fld[..., c].std()
    # 1) KCUT=1 sharp filter must be a no-op on a periodic field (V1 == V0 exactly)
    KCUT = 1.0
    a = stats_of_field(fld, DX, MARGIN)
    b = stats_of_field(spectral_filter(fld, 1.0, "sharp"), DX, MARGIN)
    # threshold 1e-4: at KCUT=1 the spherical mask still clips CORNER modes outside the inscribed
    # sphere (harmless in production, where KCUT<0.5 sits entirely inside the cube), and the FFT
    # round-trip error is amplified by the 4th-order derivative + squaring — 1e-5-ish is expected.
    err = max(np.abs(a[i] - b[i]).max() / (np.abs(a[i]).max() + 1e-30) for i in range(3))
    print(f"selftest 1 (KCUT=1 near-no-op): max rel err = {err:.2e}  {'PASS' if err < 1e-4 else 'FAIL'}")
    # 2) energy monotone under truncation
    e_full = (fld ** 2).sum()
    e_half = (spectral_filter(fld, 0.4714, "sharp") ** 2).sum()
    print(f"selftest 2 (truncation removes energy): {e_half/e_full:.3f} < 1  "
          f"{'PASS' if e_half < e_full else 'FAIL'}")
    # 3) full variant + paired-analysis plumbing on 2 synthetic 'boxes'
    KCUT = np.sqrt(2.0) / 3.0
    pools = {v: [[], [], [], 0, [0]] for v in VARIANTS}
    for _ in range(2):
        var = compute_variants(fld + 0.01 * rng.randn(*fld.shape))
        for v in VARIANTS:
            W2, PP, DN, nb, off = pools[v]
            w2, p, d = (a.astype(np.float32) for a in var[v])
            W2.append(w2); PP.append(p); DN.append(d)
            pools[v] = [W2, PP, DN, nb + 1, off + [off[-1] + len(w2)]]
    analyze({v: tuple(p) for v, p in pools.items()})
    print("selftest 3 (plumbing): PASS (analysis ran end-to-end)")


# ---------------- main ----------------
def main():
    if "--selftest" in sys.argv:
        selftest(); return
    if "--analyze" in sys.argv:
        analyze(); return
    pools = {v: list(load_pool(v)) for v in VARIANTS}
    nb_done = min(p[3] for p in pools.values())
    print(f"coarse-grain run: {DATASET} REGION={REGION} B_FETCH={B_FETCH} NBOX={NBOX} "
          f"MARGIN={MARGIN} KCUT={KCUT:.4f} seed={SEED}  (resumed {nb_done} boxes)", flush=True)
    if not os.path.exists(REPORT):
        append_report(f"# Coarse-graining resolution-bias experiment — {DATASET}\n\n"
                      f"_REGION={REGION} B_FETCH={B_FETCH} NBOX={NBOX} MARGIN={MARGIN} "
                      f"KCUT_FRAC_NYQ={KCUT:.4f} seed={SEED}_\n")
    rng = np.random.RandomState(SEED)
    # regenerate the deterministic corner sequence, skipping already-done boxes
    corners = [rng.randint(1, NGRID - REGION, size=3) for _ in range(NBOX)]
    for b in range(nb_done, NBOX):
        t = TIMES[b % len(TIMES)]
        corner = corners[b]
        print(f"box {b+1}/{NBOX} t={t} corner={corner.tolist()}", flush=True)
        t0 = time.time()
        fld = assemble_region(t, corner)
        var = compute_variants(fld)
        del fld
        for v in VARIANTS:
            W2, PP, DN, nb, off = pools[v]
            w2, p, d = (a.astype(np.float32) for a in var[v])
            W2.append(w2); PP.append(p); DN.append(d)
            pools[v] = [W2, PP, DN, nb + 1, off + [off[-1] + len(w2)]]
            save_pool(v, *pools[v])
        rb = {v: per_box_ratios(*map(np.asarray, var[v][:3]), QS) for v in VARIANTS}
        line = (f"box {b+1}: " + "  ".join(
            f"{v.split('_')[0]}[mean={rb[v]['D_mean']:.3f} r1%={rb[v]['r_top1%']:.3f} "
            f"r0.1%={rb[v]['r_top0.1%']:.3f}]" for v in VARIANTS) +
            f"  ({time.time()-t0:.0f}s)")
        print("  " + line, flush=True)
        append_report(line)
    verdicts = analyze({v: tuple(p) for v, p in pools.items()})
    # persist the analysis to the report too
    if verdicts:
        buf = []
        analyze({v: tuple(p) for v, p in pools.items()}, echo=lambda s: buf.append(str(s)))
        append_report("\n```\n" + "\n".join(buf) + "\n```\n")
    print(f"\nreport -> {REPORT}")


if __name__ == "__main__":
    main()
