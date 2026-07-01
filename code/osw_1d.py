r"""osw_1d.py — the Okamoto-Sakajo-Wunsch 1D model family: does advection deplete blow-up?

The 1D models that isolate the SAME mechanism question as 3D NS regularity — whether transport/advection
geometrically depletes a nonlinear blow-up:

    omega_t + a * u * omega_x = u_x * omega ,    with   u_x = H[omega]   (Hilbert transform)

  a = 0  -> Constantin-Lax-Majda (1985): EXACT finite-time blow-up (no advection). Our validation case.
  a = 1  -> De Gregorio (1990): advection believed to REGULARIZE (no blow-up from smooth data on the line).
  scan a in between -> map the blow-up <-> regularity boundary = how much advection depletes the CLM blow-up.

This is the 1D shadow of the 3D "does vortex stretching (omega.grad u) get depleted by transport" question.
Pure numpy FFT spectral, 2/3 dealiasing, RK4 with CFL-adaptive dt, blow-up detected by ||omega||_inf.

  python osw_1d.py                       # validate (a=0 blows up, a=1 regular) + scan a
  python osw_1d.py --N 16384 --avals 0,0.25,0.5,0.75,1.0
"""
from __future__ import annotations
import sys, json, os, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "osw_results.md")


def make_ops(N):
    k = np.fft.fftfreq(N, d=1.0 / N)            # integer wavenumbers on [0,2pi)
    ik = 1j * k
    H_mult = -1j * np.sign(k)                    # Hilbert transform multiplier:  H[e^{ikx}] = -i sgn(k) e^{ikx}
    inv_k = np.zeros_like(k, dtype=float)
    nz = k != 0
    inv_k[nz] = 1.0 / k[nz]
    # u from u_x = H[omega]:  u_hat = (H omega)_hat / (ik) = (-i sgn k)/(ik) omega_hat = -(1/|k|) omega_hat
    u_from_omega = np.zeros_like(k, dtype=complex)
    u_from_omega[nz] = (H_mult[nz]) / (ik[nz])
    # 2/3 dealiasing mask
    kmax = np.max(np.abs(k))
    mask = (np.abs(k) <= (2.0 / 3.0) * kmax).astype(float)
    return k, ik, H_mult, u_from_omega, mask


def rhs(omega, a, ik, H_mult, u_from_omega, mask):
    wh = np.fft.fft(omega) * mask
    Hw = np.real(np.fft.ifft(H_mult * wh))      # u_x = H[omega]
    wx = np.real(np.fft.ifft(ik * wh))          # omega_x
    u = np.real(np.fft.ifft(u_from_omega * wh)) # velocity
    dwdt = Hw * omega - a * u * wx
    # project back through the dealias mask
    return np.real(np.fft.ifft(np.fft.fft(dwdt) * mask))


def initial_data(x, ic):
    if ic == "sin":
        return np.sin(x)
    if ic in ("tri", "cusp"):                   # C^0 triangle wave: corners (not C^1) = low-regularity
        return (2.0 / np.pi) * np.arcsin(np.sin(x))
    if ic == "rough":                           # sin + high-frequency roughness
        return np.sin(x) + 0.3 * np.sin(7 * x) + 0.2 * np.sin(13 * x)
    if ic == "ic2":                             # different smooth datum (universality check)
        return np.sin(x) + 0.5 * np.sin(2 * x)
    raise ValueError(ic)


def run_one(a, N=8192, T=8.0, blowup_thresh=2e4, dt0=2e-4, log=None, ic="sin"):
    k, ik, H_mult, u_from_omega, mask = make_ops(N)
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    omega = initial_data(x, ic).astype(float)
    omega = omega / max(np.max(np.abs(omega)), 1e-12)   # normalize peak to 1 for fair comparison
    t = 0.0
    w0 = np.max(np.abs(omega))
    peak = w0
    steps = 0
    while t < T:
        wmax = np.max(np.abs(omega))
        peak = max(peak, wmax)
        if wmax > blowup_thresh or not np.isfinite(wmax):
            return {"a": a, "blowup": True, "t_blow": t, "peak": float(wmax),
                    "amplification": float(wmax / w0), "steps": steps, "N": N}
        # CFL-ish adaptive dt: keep dt * max|omega| bounded (the nonlinearity timescale)
        dt = min(dt0, 0.05 / max(wmax, 1.0))
        if t + dt > T:
            dt = T - t
        k1 = rhs(omega, a, ik, H_mult, u_from_omega, mask)
        k2 = rhs(omega + 0.5 * dt * k1, a, ik, H_mult, u_from_omega, mask)
        k3 = rhs(omega + 0.5 * dt * k2, a, ik, H_mult, u_from_omega, mask)
        k4 = rhs(omega + dt * k3, a, ik, H_mult, u_from_omega, mask)
        omega = omega + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        steps += 1
        if log and steps % 20000 == 0:
            print(f"   a={a:+.3f} t={t:.4f} max|w|={np.max(np.abs(omega)):.3e}", flush=True)
    return {"a": a, "blowup": False, "t_blow": None, "peak": float(peak),
            "amplification": float(peak / w0), "steps": steps, "N": N}


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=8192)
    ap.add_argument("--T", type=float, default=8.0)
    ap.add_argument("--avals", default="0.0,0.25,0.5,0.75,1.0,1.25")
    ap.add_argument("--ic", default="sin", help="sin | tri/cusp | rough | ic2")
    args = ap.parse_args()
    a_list = [float(s) for s in args.avals.split(",")]
    global RESULTS
    RESULTS = os.path.join(HERE, f"osw_results_{args.ic}.md")
    print(f"OSW 1D model family — N={args.N}, T={args.T}, ic={args.ic}")
    print("a=0 (CLM) should BLOW UP; a=1 (De Gregorio) should stay REGULAR for smooth data.\n")
    out = []
    for a in a_list:
        r = run_one(a, N=args.N, T=args.T, log=True, ic=args.ic)
        tag = f"BLOW-UP @ t={r['t_blow']:.4f}" if r["blowup"] else f"regular (peak {r['peak']:.3f})"
        print(f"a={a:+.3f}:  {tag}  [amp x{r['amplification']:.2e}, {r['steps']} steps]", flush=True)
        out.append(r)
    # boundary estimate
    blow = [r["a"] for r in out if r["blowup"]]
    reg = [r["a"] for r in out if not r["blowup"]]
    with open(RESULTS, "w", encoding="utf-8") as f:
        f.write("# OSW 1D model family — advection vs blow-up\n\n")
        f.write(f"omega_t + a*u*omega_x = u_x*omega,  u_x=H[omega],  data=sin(x), N={args.N}, T={args.T}\n\n")
        f.write("| a | outcome | t_blow | peak |omega| | amplification |\n|---|---|---|---|---|\n")
        for r in out:
            f.write(f"| {r['a']:+.3f} | {'BLOW-UP' if r['blowup'] else 'regular'} | "
                    f"{('%.4f'%r['t_blow']) if r['blowup'] else '-'} | {r['peak']:.3e} | "
                    f"{r['amplification']:.2e} |\n")
        if blow and reg:
            f.write(f"\n**Blow-up <-> regularity boundary** lies in a in ({max(blow):.3f}, {min(reg):.3f}): "
                    f"advection strength a depletes the CLM blow-up between these.\n")
    print("\n" + open(RESULTS, encoding="utf-8").read())
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
