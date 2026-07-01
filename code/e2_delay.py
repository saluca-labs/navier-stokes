r"""e2_delay.py (v2, amplitude form) — Tao-delay concentration vs energy-spreading, SAME viscosity.

E1 (standard dyadic, amplitude form) was REGULAR at nu=1e-3, gamma=0.5 — the forward cascade SPREADS
across shells and dissipation catches it. Catalog claim: that spreading is the regularizer; forcing
sequential near-complete transfer (Tao's delay) destroys it and forces blow-up. Test at the SAME viscosity.

Conservative gated dyadic (amplitude form, energy sum a_n^2 conserved inviscidly):
    da_n/dt = k_{n-1} a_{n-1}^2 G_{n-1}  -  k_n a_n a_{n+1} G_n  -  nu k_n^{2g} a_n
The transfer n->n+1 is the PAIR (+k_{n-1}a_{n-1}^2 into n) & (-k_{n-1}a_{n-1}a_n out of n-1); gating both
by the SAME G_{n-1} preserves the telescoping cancellation -> energy still conserved.
    SPREAD:      G_n = 1                              (standard; parallel transfers -> spreading)
    CONCENTRATE: G_n = a_n^2/(a_n^2 + a_{n-1}^2 + eps)  (n releases forward only after predecessor depletes
                                                        -> sequential single-front cascade, no spreading)
Blow-up = enstrophy sum(k_n^2 a_n^2) -> infinity.
"""
import numpy as np
from scipy.integrate import solve_ivp

N = 22
K = 2.0 ** np.arange(N)
KM1 = np.concatenate(([0.0], K[:-1]))
EPS = 1e-30
TMAX = 60.0
BLOWUP = 1e10


def rhs(t, a, nu, gamma, mode):
    am1 = np.concatenate(([0.0], a[:-1]))           # a_{n-1}
    ap1 = np.concatenate((a[1:], [0.0]))            # a_{n+1}
    if mode == "spread":
        G = np.ones(N)
    else:  # concentrate: gate n's forward release until its predecessor is depleted
        G = (a * a) / (a * a + am1 * am1 + EPS)
    Gm1 = np.concatenate(([1.0], G[:-1]))           # G_{n-1}
    inflow  = KM1 * am1 * am1 * Gm1                 # k_{n-1} a_{n-1}^2 G_{n-1}
    outflow = K * a * ap1 * G                       # k_n a_n a_{n+1} G_n
    return inflow - outflow - nu * (K ** (2.0 * gamma)) * a


def run(mode, nu, gamma, label):
    a0 = np.zeros(N); a0[0] = 1.0

    def blow(t, a, nu, gamma, mode):
        return BLOWUP - float(np.sum((K ** 2) * a * a))
    blow.terminal = True; blow.direction = -1

    sol = solve_ivp(rhs, [0.0, TMAX], a0, method="LSODA", args=(nu, gamma, mode),
                    events=blow, rtol=1e-9, atol=1e-13, max_step=TMAX / 200.0)
    a = sol.y[:, -1]
    enstr = float(np.sum((K ** 2) * a * a))
    w = (K ** 2) * (a * a); w = w / (w.sum() + 1e-300)
    idx = np.arange(N); mean = (idx * w).sum()
    spread = float(np.sqrt(((idx - mean) ** 2 * w).sum()))
    front = int(np.argmax((K ** 2) * a * a))
    blew = sol.status == 1 and len(sol.t_events[0]) > 0
    if blew:
        print(f"  [{label}] BLOW-UP at t={float(sol.t_events[0][0]):.3f}: front shell {front}/{N}, active-shell spread={spread:.2f}")
    else:
        print(f"  [{label}] REGULAR to t={sol.t[-1]:.1f}: enstrophy={enstr:.3e}, front shell {front}/{N}, active-shell spread={spread:.2f}")
    return {"label": label, "blew": blew, "spread": spread}


if __name__ == "__main__":
    print("E2 (amplitude form) — spreading vs Tao-delay concentration, SAME viscosity\n")
    r1 = run("spread",      1e-3, 0.50, "SPREAD  g=0.50 (E1 control)")
    r2 = run("concentrate", 1e-3, 0.50, "DELAY   g=0.50 (concentration)")
    print("\ncontrols:")
    run("concentrate", 0.0,  0.50, "DELAY   inviscid (expect blow-up)")
    run("spread",      0.0,  0.50, "SPREAD  inviscid (expect blow-up)")
    print("\n=== verdict ===")
    print(f"  spreading @ g=0.50:     {'BLOW-UP' if r1['blew'] else 'REGULAR'}  (active-shell spread {r1['spread']:.2f})")
    print(f"  concentration @ g=0.50: {'BLOW-UP' if r2['blew'] else 'REGULAR'}  (active-shell spread {r2['spread']:.2f})")
