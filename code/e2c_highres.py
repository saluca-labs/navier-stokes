r"""e2c_highres.py — does the relay-vs-spread threshold gap open up at higher resolution?

E2b at N=22 found both thresholds ~0.42 (no gap). The heuristic says the complete-transfer relay's
threshold -> 1/2 only ASYMPTOTICALLY (k_n^{1-2g} effect needs large n). Push N up (more shells) and a
fine gamma grid; LSODA handles the stiff viscosity. Prediction: relay threshold climbs toward 0.5 as N
grows while spreading stays ~0.4 -> the gap (spreading-regularizes) becomes visible.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 34
K = 2.0 ** np.arange(N)
KM1 = np.concatenate(([0.0], K[:-1]))
TMAX = 80.0


def rhs_relay(t, e, nu, g):
    e = np.maximum(e, 0.0)
    F = K * e ** 1.5
    return np.concatenate(([0.0], F[:-1])) - F - 2.0 * nu * (K ** (2 * g)) * e


def rhs_spread(t, a, nu, g):
    am1 = np.concatenate(([0.0], a[:-1])); ap1 = np.concatenate((a[1:], [0.0]))
    return KM1 * am1 * am1 - K * a * ap1 - nu * (K ** (2 * g)) * a


def blows(kind, nu, g):
    y0 = np.zeros(N); y0[0] = 1.0
    rhs = rhs_relay if kind == "relay" else rhs_spread
    sol = solve_ivp(rhs, [0, TMAX], y0, method="LSODA", args=(nu, g),
                    rtol=1e-8, atol=1e-12, max_step=TMAX / 300.0)
    y = sol.y[:, -1]
    w = (K ** 2) * (np.maximum(y, 0.0) if kind == "relay" else y * y)
    front = int(np.argmax(w))
    return front >= N - 2          # cascade reached resolution limit = blow-up


if __name__ == "__main__":
    print(f"E2c — threshold gap vs resolution (N={N} shells, nu=1e-3)\n")
    grid = np.arange(0.36, 0.521, 0.02)
    print(f"  {'gamma':>6} | {'relay':>9} | {'spread':>9}")
    relay_thr = spread_thr = None
    for g in grid:
        rb = blows("relay", 1e-3, g); sb = blows("spread", 1e-3, g)
        print(f"  {g:>6.2f} | {('BLOW-UP' if rb else 'REGULAR'):>9} | {('BLOW-UP' if sb else 'REGULAR'):>9}")
        if relay_thr is None and not rb: relay_thr = g
        if spread_thr is None and not sb: spread_thr = g
    print(f"\n  relay regular from gamma~{relay_thr} ; spread regular from gamma~{spread_thr}")
    if relay_thr and spread_thr and relay_thr > spread_thr + 1e-9:
        print(f"  >>> GAP OPENS: concentration stays singular to higher gamma ({relay_thr}) than spreading "
              f"({spread_thr}) — spreading regularizes, as predicted.")
    else:
        print("  >>> still no clear gap at this N — push N higher (needs integrating-factor for N>~40).")
