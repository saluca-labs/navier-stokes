r"""d3_rescan.py — does geometric depletion strengthen or weaken with Reynolds number?

The crux measurement. For each Re, run Taylor-Green to past the enstrophy peak and report the
geometric-depletion factor  depl = <omega . S . omega> / <|omega|^2 |S|>  at the peak (and time-averaged
in the developed window). Resolution is matched (N ~ Re^{3/4}) to hold k_max*eta ~ 1 constant, and
k_max*eta is reported per run so under-resolution can't masquerade as physics.

  depl -> 0 as Re grows  => nonlinearity self-depletes => evidence the flow protects itself (regularity).
  depl -> O(1) as Re grows => de-depletion at high Re   => evidence toward singularity.

  python d3_rescan.py
"""
import math, sys, time
import torch

DEV = "cuda"
for i in range(torch.cuda.device_count()):
    if "5090" in torch.cuda.get_device_name(i):
        DEV = f"cuda:{i}"; break


def run_case(N, Re, T, dt):
    nu = 1.0 / Re
    x = torch.arange(N, device=DEV) * (2 * math.pi / N)
    X, Y, Z = torch.meshgrid(x, x, x, indexing="ij")
    u = torch.stack([torch.sin(X) * torch.cos(Y) * torch.cos(Z),
                     -torch.cos(X) * torch.sin(Y) * torch.cos(Z), torch.zeros_like(X)]).float()
    del X, Y, Z
    kx = torch.fft.fftfreq(N, d=1.0 / N, device=DEV).float()
    kz = torch.fft.rfftfreq(N, d=1.0 / N, device=DEV).float()
    KX = kx.view(N, 1, 1); KY = kx.view(1, N, 1); KZ = kz.view(1, 1, -1); KS = [KX, KY, KZ]
    K2 = KX ** 2 + KY ** 2 + KZ ** 2; K2c = K2.clone(); K2c[0, 0, 0] = 1.0
    kmax = N // 3
    mask = ((KX.abs() <= kmax) & (KY.abs() <= kmax) & (KZ.abs() <= kmax)).float()

    def fft(f):   return torch.fft.rfftn(f, dim=(1, 2, 3))
    def ifft(F):  return torch.fft.irfftn(F, dim=(1, 2, 3), s=(N, N, N))
    def ifft1(F): return torch.fft.irfftn(F, dim=(0, 1, 2), s=(N, N, N))
    def project(F):
        kdu = KX * F[0] + KY * F[1] + KZ * F[2]
        return torch.stack([F[i] - KS[i] * kdu / K2c for i in range(3)])
    def vhat(F):
        return torch.stack([1j * (KY * F[2] - KZ * F[1]), 1j * (KZ * F[0] - KX * F[2]),
                            1j * (KX * F[1] - KY * F[0])])
    def rhs(F):
        uu = ifft(F); ww = ifft(vhat(F))
        nl = torch.stack([uu[1] * ww[2] - uu[2] * ww[1], uu[2] * ww[0] - uu[0] * ww[2],
                          uu[0] * ww[1] - uu[1] * ww[0]])
        return project(fft(nl) * mask) - nu * K2 * F

    def diag(F):
        A = [[ifft1(1j * KS[j] * F[i]) for j in range(3)] for i in range(3)]
        wx = A[2][1] - A[1][2]; wy = A[0][2] - A[2][0]; wz = A[1][0] - A[0][1]
        w2 = wx * wx + wy * wy + wz * wz
        S = [[0.5 * (A[i][j] + A[j][i]) for j in range(3)] for i in range(3)]
        wv = [wx, wy, wz]
        Pf = sum(wv[i] * S[i][j] * wv[j] for i in range(3) for j in range(3))
        Sm = (2.0 * sum(S[i][j] * S[i][j] for i in range(3) for j in range(3))).sqrt()
        ens = float(0.5 * w2.mean()); wmax = float(w2.max().sqrt())
        depl = float(Pf.mean()) / (float((w2 * Sm).mean()) + 1e-30)
        return ens, wmax, depl

    Fh = project(fft(u)); del u
    nsteps = int(T / dt); rec = []
    for s in range(nsteps + 1):
        if s % 40 == 0 or s == nsteps:
            ens, wmax, depl = diag(Fh)
            rec.append((s * dt, ens, wmax, depl))
        if s == nsteps:
            break
        k1 = rhs(Fh); k2 = rhs(Fh + 0.5 * dt * k1); k3 = rhs(Fh + 0.5 * dt * k2); k4 = rhs(Fh + dt * k3)
        Fh = Fh + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    rec = torch.tensor(rec)
    ip = int(torch.argmax(rec[:, 1]))                # enstrophy peak index
    tpk, enspk, wmaxpk, deplpk = (float(rec[ip, j]) for j in range(4))
    win = (rec[:, 0] >= tpk - 1) & (rec[:, 0] <= tpk + 1)
    depl_mean = float(rec[win, 3].mean())
    # resolution: eps = nu<|w|^2> = 2 nu * ens(peak);  eta=(nu^3/eps)^.25 ; k_max*eta
    eps = 2 * nu * enspk
    eta = (nu ** 3 / eps) ** 0.25
    kme = kmax * eta
    return dict(N=N, Re=Re, tpk=tpk, enspk=enspk, deplpk=deplpk, depl_mean=depl_mean,
                wmaxpk=wmaxpk, kmax_eta=kme)


if __name__ == "__main__":
    print(f"D3 — geometric-depletion Re-scan on {torch.cuda.get_device_name(DEV)}\n")
    cases = [(128, 400, 10.0, 0.006), (192, 800, 10.0, 0.005), (256, 1600, 10.0, 0.005)]
    print(f"  {'Re':>6} {'N':>5} {'kmax*eta':>9} {'t_peak':>7} {'enstr_pk':>9} {'max|w|':>8} "
          f"{'depl@pk':>8} {'depl_avg':>9}")
    res = []
    for N, Re, T, dt in cases:
        t0 = time.time(); r = run_case(N, Re, T, dt); res.append(r)
        print(f"  {r['Re']:>6.0f} {r['N']:>5} {r['kmax_eta']:>9.2f} {r['tpk']:>7.2f} {r['enspk']:>9.3f} "
              f"{r['wmaxpk']:>8.1f} {r['deplpk']:>8.3f} {r['depl_mean']:>9.3f}   ({time.time()-t0:.0f}s)", flush=True)
    print("\n=== trend ===")
    for r in res:
        print(f"  Re={r['Re']:.0f}: geometric depletion (avg near peak) = {r['depl_mean']:.3f}")
    d0, d1 = res[0]['depl_mean'], res[-1]['depl_mean']
    print(f"\n  depletion {d0:.3f} (Re={res[0]['Re']:.0f}) -> {d1:.3f} (Re={res[-1]['Re']:.0f}): "
          + ("STRENGTHENS with Re (depl down) => flow self-regulates (toward regularity)."
             if d1 < d0 - 0.005 else
             "WEAKENS with Re (depl up) => de-depletion (toward singularity)."
             if d1 > d0 + 0.005 else "roughly FLAT in Re."))
    print("  (resolution honest: compare only rows with kmax*eta >~ 1.)")
