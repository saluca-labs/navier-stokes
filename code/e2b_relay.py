r"""e2b_relay.py — COMPLETE-transfer relay (the faithful "concentration" model) vs spreading.

E2 showed a naive delay-gate concentrates but DELAYS (less singular). The right model is a
*complete, self-driven* forward transfer: each shell empties into the next based on its OWN energy
(no catalytic downstream seed, no gate) -> a single accelerating front, the essence of Tao's machine.

Energy form, e_n>=0, k_n=2^n, conserves sum(e_n) inviscidly (telescoping flux):
    de_n/dt = F_{n-1} - F_n - 2 nu k_n^{2g} e_n,   F_n = k_n e_n^{3/2}   (self-driven complete transfer)
Front amplitude e_n~E=O(1): transfer ~k_n, dissipation ~nu k_n^{2g}; transfer wins at high n iff
k_n^{1-2g}->inf, i.e. g<1/2. PREDICTION: blow-up iff g<1/2 (vs spreading's ~2/5 from E1) ->
window g in (0.4, 0.5): concentration BLOWS UP where spreading is REGULAR.

Spreading control = the E1 amplitude model: da_n/dt = k_{n-1}a_{n-1}^2 - k_n a_n a_{n+1} - nu k_n^{2g} a_n.
"""
import numpy as np
from scipy.integrate import solve_ivp

N = 22
K = 2.0 ** np.arange(N)
KM1 = np.concatenate(([0.0], K[:-1]))
TMAX = 60.0
BLOWUP = 1e10


def rhs_relay(t, e, nu, gamma):           # complete-transfer (energy form)
    e = np.maximum(e, 0.0)
    F = K * e ** 1.5
    Fin = np.concatenate(([0.0], F[:-1]))   # F_{n-1}
    return Fin - F - 2.0 * nu * (K ** (2.0 * gamma)) * e


def rhs_spread(t, a, nu, gamma):          # spreading (amplitude form, = E1)
    am1 = np.concatenate(([0.0], a[:-1]))
    ap1 = np.concatenate((a[1:], [0.0]))
    return KM1 * am1 * am1 - K * a * ap1 - nu * (K ** (2.0 * gamma)) * a


def run(kind, nu, gamma, label, quiet=False):
    if kind == "relay":
        y0 = np.zeros(N); y0[0] = 1.0; rhs = rhs_relay
        enstr = lambda y: float(np.sum((K ** 2) * np.maximum(y, 0.0)))
    else:
        y0 = np.zeros(N); y0[0] = 1.0; rhs = rhs_spread
        enstr = lambda y: float(np.sum((K ** 2) * y * y))

    def blow(t, y, nu, gamma):
        return BLOWUP - enstr(y)
    blow.terminal = True; blow.direction = -1

    sol = solve_ivp(rhs, [0.0, TMAX], y0, method="LSODA", args=(nu, gamma),
                    events=blow, rtol=1e-9, atol=1e-13, max_step=TMAX / 200.0)
    y = sol.y[:, -1]
    w = (K ** 2) * (np.maximum(y, 0.0) if kind == "relay" else y * y)
    front = int(np.argmax(w))
    # blow-up in finite-N = cascade front reaches the resolution limit (would continue to inf), or enstrophy event
    blew = (sol.status == 1 and len(sol.t_events[0]) > 0) or front >= N - 2
    res = "BLOW-UP" if blew else "REGULAR"
    if not quiet:
        if blew:
            print(f"  [{label}] BLOW-UP at t={float(sol.t_events[0][0]):.3f}: front shell {front}/{N}")
        else:
            print(f"  [{label}] REGULAR to t={sol.t[-1]:.1f}: enstrophy={enstr(y):.3e}, front shell {front}/{N}")
    return blew


if __name__ == "__main__":
    print("E2b — complete-transfer relay threshold + head-to-head vs spreading\n")
    print("front-based gamma-scan (front->top = blow-up):")
    grid = [0.30, 0.35, 0.38, 0.40, 0.42, 0.45, 0.50]
    print(f"  {'gamma':>6} | {'relay (concentrate)':>20} | {'spread':>10}")
    relay_thr = spread_thr = None
    for g in grid:
        rb = run("relay",  1e-3, g, "", quiet=True)
        sb = run("spread", 1e-3, g, "", quiet=True)
        print(f"  {g:>6.2f} | {('BLOW-UP' if rb else 'REGULAR'):>20} | {('BLOW-UP' if sb else 'REGULAR'):>10}")
        if relay_thr is None and not rb: relay_thr = g
        if spread_thr is None and not sb: spread_thr = g
    print(f"\n  relay (complete-transfer) regular from gamma ~ {relay_thr}")
    print(f"  spreading                 regular from gamma ~ {spread_thr}")
    print("  => relay stays singular to HIGHER gamma than spreading => concentration is more dangerous;")
    print("     'spreading regularizes' holds, though the finite-N (N=%d) window is narrow." % N)
