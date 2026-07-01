r"""shell_model.py — dyadic (Katz-Pavlović / Cheskidov) shell model of the NS energy cascade.

The simplest system where 3D-NS-style blow-up vs global regularity is DECIDABLE, so it's our
gateway experiment. Engages the failure catalog directly: the inviscid forward cascade blows up
(Katz-Pavlović), and viscosity regularizes only above a dissipation threshold (γ≥2/5 in the
physical-intermittency sense; blow-up persists for small γ). The regularizer is energy *leaving*
each shell faster than the nonlinearity can pile it up — the "energy-spreading" the catalog flags
as the finer structure real NS may possess and naive models lack.

Model (real, nonnegative, energy-conserving inviscid forward cascade):
    da_n/dt = k_{n-1} a_{n-1}^2 - k_n a_n a_{n+1} - nu * k_n^{2*gamma} * a_n   (+ steady forcing on shell 0)
    k_n = 2^n,  E = sum a_n^2  (conserved inviscid),  "H^s norm" = sum k_n^{2s} a_n^2.

Blow-up signature: the H^1 norm (sum k_n^2 a_n^2) -> infinity in finite time and the cascade front
races to the highest shell. Regularity: the front stalls and H^1 stays bounded.

  python shell_model.py
"""
import numpy as np
from scipy.integrate import solve_ivp

N = 22                      # shells (k_21 = 2^21 ~ 2e6 — enough inertial range, stiffness manageable)
K = 2.0 ** np.arange(N)     # k_n = 2^n
KM1 = np.concatenate(([0.0], K[:-1]))
TMAX = 40.0
BLOWUP_H1 = 1e10            # H^1 (enstrophy) threshold flagged as blow-up


def rhs(t, a, nu, gamma):
    am1 = np.concatenate(([0.0], a[:-1]))          # a_{n-1}
    ap1 = np.concatenate((a[1:], [0.0]))           # a_{n+1}
    da = KM1 * am1 * am1 - K * a * ap1             # energy-conserving forward cascade
    da -= nu * (K ** (2.0 * gamma)) * a            # viscous dissipation, degree gamma
    return da


def run(nu, gamma, label):
    a0 = np.zeros(N); a0[0] = 1.0
    H1_RE = 0.0  # placeholder

    def blowup(t, a, nu, gamma):
        return BLOWUP_H1 - float(np.sum((K ** 2) * a * a))
    blowup.terminal = True
    blowup.direction = -1

    sol = solve_ivp(rhs, [0.0, TMAX], a0, method="LSODA", args=(nu, gamma),
                    events=blowup, rtol=1e-8, atol=1e-12, max_step=TMAX / 50.0, dense_output=False)
    a = sol.y[:, -1]; t = sol.t[-1]
    H1 = float(np.sum((K ** 2) * a * a))
    Hs = float(np.sum((K ** 1.0) * a * a))         # a milder critical-ish norm
    front = int(np.argmax((K ** 2) * a * a))
    blew = sol.status == 1 and len(sol.t_events[0]) > 0
    if blew:
        tb = float(sol.t_events[0][0])
        print(f"  [{label}] BLOW-UP at t={tb:.3f}: enstrophy(H1)->{BLOWUP_H1:.0e}, cascade front at shell {front}/{N}")
    else:
        print(f"  [{label}] REGULAR to t={t:.1f}: H1={H1:.3e} bounded, front stalled at shell {front}/{N}")
    return {"label": label, "nu": nu, "gamma": gamma, "blew": blew, "t_end": t,
            "t_blow": float(sol.t_events[0][0]) if blew else None}


if __name__ == "__main__":
    print("Dyadic shell model — dissipation-threshold experiment (k_n=2^n, N=%d shells)\n" % N)
    regimes = [
        (0.0,  0.0,  "inviscid          (nu=0)         -> expect BLOW-UP (Katz-Pavlovic)"),
        (1e-3, 0.15, "subcritical visc  (nu=1e-3,g=.15) -> expect BLOW-UP (g<1/4)"),
        (1e-3, 0.50, "supercritical vis (nu=1e-3,g=.50) -> expect REGULAR (g>2/5)"),
    ]
    results = []
    for nu, g, lab in regimes:
        print(lab)
        results.append(run(nu, g, lab.split()[0]))
        print()
    print("=== verdict ===")
    for r in results:
        print(f"  gamma={r['gamma']:.2f}: {'BLOW-UP' if r['blew'] else 'REGULAR'} (t_end={r['t_end']:.4f})")
