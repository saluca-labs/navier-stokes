r"""d4_extreme.py — does the EXTREME tail deplete or de-deplete with Re? (the singularity-relevant test)

D3 showed the MEAN geometric depletion strengthens with Re. But blow-up is about the rare points of
maximum vortex stretching. Here we measure depletion conditioned on the top-q quantile of |omega|:
    depl(q) = sum_{top-q |w|} (w.S.w) / sum_{top-q |w|} (|w|^2 |S|)
for q = 1%, 0.1%, 0.01%, at the enstrophy peak, vs Re. If the extreme-tail depletion WEAKENS (rises
toward 1) with Re while the mean falls, that's the Buaria-Pumir de-depletion-of-extremes signature -
the knife-edge that matters for regularity.

CAVEAT: extreme events are small-scale; they need k_max*eta >~ 1 to be trusted. Reported per run.

  python d4_extreme.py
"""
import math, sys, time
import torch

DEV = "cuda"
for i in range(torch.cuda.device_count()):
    if "5090" in torch.cuda.get_device_name(i):
        DEV = f"cuda:{i}"; break
QS = [0.01, 0.001, 0.0001]


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
        w2 = (wx * wx + wy * wy + wz * wz)
        S = [[0.5 * (A[i][j] + A[j][i]) for j in range(3)] for i in range(3)]
        wv = [wx, wy, wz]
        Pf = sum(wv[i] * S[i][j] * wv[j] for i in range(3) for j in range(3))
        Smag = (2.0 * sum(S[i][j] * S[i][j] for i in range(3) for j in range(3))).sqrt()
        den = w2 * Smag
        ens = float(0.5 * w2.mean()); wmax = float(w2.max().sqrt())
        mean_depl = float(Pf.mean()) / (float(den.mean()) + 1e-30)
        w2f = w2.flatten(); Pff = Pf.flatten(); df = den.flatten(); M = w2f.numel()
        ex = []
        for q in QS:
            k = max(1, int(q * M))
            idx = torch.topk(w2f, k).indices
            ex.append(float(Pff[idx].sum()) / (float(df[idx].sum()) + 1e-30))
        return ens, wmax, mean_depl, ex

    Fh = project(fft(u)); del u
    nsteps = int(T / dt); rec = []
    for s in range(nsteps + 1):
        if s % 40 == 0 or s == nsteps:
            ens, wmax, md, ex = diag(Fh)
            rec.append((s * dt, ens, wmax, md, *ex))
        if s == nsteps:
            break
        k1 = rhs(Fh); k2 = rhs(Fh + 0.5 * dt * k1); k3 = rhs(Fh + 0.5 * dt * k2); k4 = rhs(Fh + dt * k3)
        Fh = Fh + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    rec = torch.tensor(rec); ip = int(torch.argmax(rec[:, 1]))
    eps = 2 * nu * float(rec[ip, 1]); eta = (nu ** 3 / eps) ** 0.25
    return dict(N=N, Re=Re, tpk=float(rec[ip, 0]), enspk=float(rec[ip, 1]),
                mean=float(rec[ip, 3]), ex=[float(rec[ip, 4 + i]) for i in range(len(QS))],
                kme=kmax * eta)


if __name__ == "__main__":
    print(f"D4 — mean vs EXTREME-tail depletion vs Re on {torch.cuda.get_device_name(DEV)}\n")
    cases = [(128, 400, 10, 0.006), (192, 800, 10, 0.005), (256, 1600, 10, 0.005), (384, 3200, 10, 0.004)]
    print(f"  {'Re':>6} {'N':>5} {'kme':>5} {'mean':>6} {'top1%':>7} {'top0.1%':>8} {'top0.01%':>9}")
    res = []
    for N, Re, T, dt in cases:
        try:
            t0 = time.time(); r = run_case(N, Re, T, dt); res.append(r)
            print(f"  {r['Re']:>6.0f} {r['N']:>5} {r['kme']:>5.2f} {r['mean']:>6.3f} "
                  f"{r['ex'][0]:>7.3f} {r['ex'][1]:>8.3f} {r['ex'][2]:>9.3f}   ({time.time()-t0:.0f}s)", flush=True)
        except RuntimeError as e:
            print(f"  Re={Re} N={N}: FAILED ({str(e)[:60]})", flush=True); torch.cuda.empty_cache()
    print("\n=== trend (depletion factor; 1=no depletion, <<1=depleted) ===")
    for r in res:
        print(f"  Re={r['Re']:>5.0f} (kmax*eta={r['kme']:.2f}): mean={r['mean']:.3f}  "
              f"top1%={r['ex'][0]:.3f}  top0.1%={r['ex'][1]:.3f}  top0.01%={r['ex'][2]:.3f}")
    if len(res) >= 2:
        a, b = res[0], res[-1]
        print(f"\n  MEAN:       {a['mean']:.3f} -> {b['mean']:.3f}")
        print(f"  EXTREME(top0.01%): {a['ex'][2]:.3f} -> {b['ex'][2]:.3f}  "
              + ("(rises => extremes DE-deplete with Re — the singularity-relevant signature)"
                 if b['ex'][2] > a['ex'][2] + 0.01 else
                 "(falls => extremes also deplete more — toward regularity)"
                 if b['ex'][2] < a['ex'][2] - 0.01 else "(flat)"))
    print("  (trust extreme columns only where kmax*eta >~ 1.)")
