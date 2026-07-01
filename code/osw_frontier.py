r"""osw_frontier.py — does De Gregorio (a=1) blow up for LOW-REGULARITY data? (Elgindi-Jeong frontier)

Two upgrades over osw_1d.py:

1. GROWTH-RATE / BKM CLASSIFIER. Instead of "did ||omega||_inf cross a threshold by time T" (which
   misclassifies late-forming singularities as regular), we read blow-up from the TREND of the logarithmic
   growth rate g(t) = d/dt log||omega||_inf. For a (T*-t)^{-1}-type singularity, g(t) -> infinity and
   1/g(t) -> 0 LINEARLY at the singular time T*. We fit 1/g(t) over the late window: a clean negative slope
   extrapolating to a finite T* (with a good fit) => BLOW-UP at T*; a flat/positive 1/g => REGULAR. We also
   report the BKM integral I(t)=int_0^t ||omega||_inf ds (blow-up <=> I diverges).

2. LOW-HOLDER INITIAL DATA. omega0(x) = sign(sin x) * |sin x|^alpha is C^alpha (a cusp at the stagnation
   points x=0,pi) but NOT C^1 for alpha<1; alpha=1 recovers smooth sin x. Smooth data (alpha=1) is known to
   be regularized by the De Gregorio advection (a=1); the question is whether ROUGH data (small alpha) blows
   up. We scan alpha and CHECK N-CONVERGENCE (a real cusp is resolved only down to the grid scale, so any
   blow-up must persist under refinement to be physical, not a Gibbs artifact).

  python osw_frontier.py                       # a=1 De Gregorio, scan alpha, with N-convergence
  python osw_frontier.py --a 1.0 --alphas 1.0,0.7,0.5,0.3 --N 8192
"""
from __future__ import annotations
import sys, os, json, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "osw_frontier_results.md")


def make_ops(N):
    k = np.fft.fftfreq(N, d=1.0 / N)
    ik = 1j * k
    H_mult = -1j * np.sign(k)
    nz = k != 0
    u_from_omega = np.zeros_like(k, dtype=complex)
    u_from_omega[nz] = H_mult[nz] / (ik[nz])
    kmax = np.max(np.abs(k))
    mask = (np.abs(k) <= (2.0 / 3.0) * kmax).astype(float)
    return ik, H_mult, u_from_omega, mask


def holder_ic(x, alpha):
    s = np.sin(x)
    w = np.sign(s) * np.abs(s) ** alpha          # C^alpha cusp at x=0,pi (alpha<1); alpha=1 -> sin x
    return w / max(np.max(np.abs(w)), 1e-12)


def rhs(omega, a, ik, H_mult, u_from_omega, mask):
    wh = np.fft.fft(omega) * mask
    Hw = np.real(np.fft.ifft(H_mult * wh))
    wx = np.real(np.fft.ifft(ik * wh))
    u = np.real(np.fft.ifft(u_from_omega * wh))
    d = Hw * omega - a * u * wx
    return np.real(np.fft.ifft(np.fft.fft(d) * mask))


def classify(ts, Ms, hit_cap=False):
    """Decide blow-up (and T*) vs regular. If the run drove ||omega|| into the cap, that IS blow-up
    (T* ~ last time). Otherwise read the trend of 1/g(t): g=d/dt log||omega||; for a (T*-t)^{-1}
    singularity 1/g -> 0 linearly, so a clean negative-slope extrapolation to a finite T* near/ahead of
    the window => impending blow-up; a flat/positive 1/g => regular."""
    ts = np.asarray(ts); Ms = np.asarray(Ms)
    peak = float(np.max(Ms))
    if hit_cap or not np.isfinite(Ms[-1]):
        return {"verdict": "BLOW-UP", "Tstar": float(ts[-1]), "inv_slope": None,
                "fit_resid": None, "peak": peak}
    if len(ts) < 8:
        return {"verdict": "inconclusive", "Tstar": None, "inv_slope": None,
                "fit_resid": None, "peak": peak}
    g = np.gradient(np.log(Ms), ts)
    j = max(2, int(0.6 * len(ts)))
    tw = ts[j:]; gw = np.clip(g[j:], 1e-6, None)
    inv = 1.0 / gw
    slope, b = np.polyfit(tw, inv, 1)
    resid = np.sqrt(np.mean((np.polyval([slope, b], tw) - inv) ** 2)) / (np.mean(inv) + 1e-12)
    Tstar = -b / slope if slope != 0 else np.inf
    horizon = ts[-1] + 4 * (ts[-1] - ts[0]) + 5
    blow = (slope < -0.2) and np.isfinite(Tstar) and (ts[0] < Tstar < horizon) \
        and resid < 0.35 and Ms[-1] > 3.0 and g[-1] > 0.3
    return {"verdict": "BLOW-UP" if blow else "regular", "Tstar": float(Tstar) if blow else None,
            "inv_slope": float(slope), "fit_resid": float(resid), "peak": peak}


def run(a, alpha, N=8192, T=25.0, blowcap=1e4, dt0=2e-4):
    ik, H_mult, u_from_omega, mask = make_ops(N)
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    omega = holder_ic(x, alpha)
    t = 0.0
    ts, Ms = [0.0], [float(np.max(np.abs(omega)))]
    bkm = 0.0
    next_rec = 0.1
    hit_cap = False
    while t < T:
        M = np.max(np.abs(omega))
        if not np.isfinite(M) or M > blowcap:
            ts.append(t); Ms.append(float(M)); hit_cap = True
            break
        dt = min(dt0, 0.05 / max(M, 1.0))
        if t + dt > T:
            dt = T - t
        k1 = rhs(omega, a, ik, H_mult, u_from_omega, mask)
        k2 = rhs(omega + 0.5 * dt * k1, a, ik, H_mult, u_from_omega, mask)
        k3 = rhs(omega + 0.5 * dt * k2, a, ik, H_mult, u_from_omega, mask)
        k4 = rhs(omega + dt * k3, a, ik, H_mult, u_from_omega, mask)
        omega = omega + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        bkm += dt * M
        t += dt
        if t >= next_rec:
            ts.append(t); Ms.append(float(np.max(np.abs(omega)))); next_rec += 0.1
    c = classify(ts, Ms, hit_cap=hit_cap)
    c.update({"a": a, "alpha": alpha, "N": N, "T_reached": t, "bkm_integral": float(bkm)})
    return c


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--alphas", default="1.0,0.7,0.5,0.3,0.2")
    ap.add_argument("--N", type=int, default=8192)
    ap.add_argument("--T", type=float, default=25.0)
    ap.add_argument("--converge", default="", help="comma N list to N-converge a single alpha (with --alphas one value)")
    a = ap.parse_args()
    alphas = [float(s) for s in a.alphas.split(",")]
    print(f"OSW frontier — a={a.a} (a=1=De Gregorio), low-Holder data |sin x|^alpha, N={a.N}, T={a.T}")
    out = []
    for al in alphas:
        r = run(a.a, al, N=a.N, T=a.T)
        Tstr = f"T*={r['Tstar']:.2f}" if r["verdict"] == "BLOW-UP" and r["Tstar"] else "-"
        sl = f"{r['inv_slope']:+.2f}" if r["inv_slope"] is not None else "cap"
        rd = f"{r['fit_resid']:.2f}" if r["fit_resid"] is not None else "-"
        print(f"  alpha={al:.2f}: {r['verdict']:8} {Tstr:10} peak={r['peak']:.2e} "
              f"BKM={r['bkm_integral']:.2f} 1/g-slope={sl} resid={rd}", flush=True)
        out.append(r)
    # N-convergence on any blow-up alpha
    conv = []
    if a.converge:
        Ns = [int(s) for s in a.converge.split(",")]
        al = alphas[0]
        print(f"\n  N-convergence at alpha={al}, a={a.a}:")
        for NN in Ns:
            r = run(a.a, al, N=NN, T=a.T)
            print(f"    N={NN:6}: {r['verdict']:8} peak={r['peak']:.2e} T*={r.get('Tstar')}", flush=True)
            conv.append(r)
    with open(RESULTS, "w", encoding="utf-8") as f:
        f.write(f"# OSW frontier — De Gregorio (a={a.a}) on low-Holder data\n\n")
        f.write("omega0 = sign(sin x)|sin x|^alpha (C^alpha cusp; alpha=1 = smooth sin x). "
                "Growth-rate/BKM classifier.\n\n")
        f.write("| alpha | verdict | T* | peak | BKM | 1/g slope | fit resid |\n|---|---|---|---|---|---|---|\n")
        for r in out:
            sl = f"{r['inv_slope']:+.2f}" if r['inv_slope'] is not None else "cap"
            rd = f"{r['fit_resid']:.2f}" if r['fit_resid'] is not None else "-"
            f.write(f"| {r['alpha']:.2f} | {r['verdict']} | {r['Tstar'] if r['Tstar'] else '-'} | "
                    f"{r['peak']:.2e} | {r['bkm_integral']:.2f} | {sl} | {rd} |\n")
        if conv:
            f.write(f"\n**N-convergence** (alpha={alphas[0]}): "
                    + "; ".join(f"N={r['N']}->{r['verdict']}(peak {r['peak']:.1e})" for r in conv) + "\n")
    print("\n" + open(RESULTS, encoding="utf-8").read())
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
