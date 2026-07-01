r"""dns_instrumented.py — instrumented GPU DNS of 3D incompressible Navier-Stokes.

Adds the structural diagnostics the failure catalog cares about, on the real equation:
  * enstrophy production P = <omega . S . omega>  (the vortex-stretching term that grows enstrophy;
    dOmega/dt = P - 2 nu <|grad omega|^2>). The supercriticality question is whether P can outrun
    dissipation; we track P and the production ratio P / (|omega|^2 |S|).
  * strain-vorticity ALIGNMENT: enstrophy-weighted <|cos(omega, e_i)|> for the three strain eigenvectors
    e1>=e2>=e3. The famous DNS fingerprint: omega aligns with the INTERMEDIATE eigenvector e2 (>0
    eigenvalue) -> self-stretching. Constantin-Fefferman: coherence of the vorticity direction in
    high-|omega| regions geometrically DEPLETES the nonlinearity (a candidate regularizing structure).
  * BKM integral  int_0^t max|omega| dt   (Beale-Kato-Majda: smooth <=> this stays finite).
  * energy spectrum E(k).

Taylor-Green vortex.  python dns_instrumented.py [N] [Re] [T]
"""
import math, sys, time
import numpy as np
import torch

N  = int(sys.argv[1])   if len(sys.argv) > 1 else 192
RE = float(sys.argv[2]) if len(sys.argv) > 2 else 1600.0
T  = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
NU = 1.0 / RE
DT = 0.005

dev = "cuda"
for i in range(torch.cuda.device_count()):
    if "5090" in torch.cuda.get_device_name(i):
        dev = f"cuda:{i}"; break
print(f"device {dev} | {N}^3 | Re={RE:.0f} (nu={NU:.2e}) | T={T} dt={DT}", flush=True)

x = torch.arange(N, device=dev) * (2 * math.pi / N)
X, Y, Z = torch.meshgrid(x, x, x, indexing="ij")
u = torch.stack([torch.sin(X) * torch.cos(Y) * torch.cos(Z),
                 -torch.cos(X) * torch.sin(Y) * torch.cos(Z),
                 torch.zeros_like(X)]).float()
del X, Y, Z

kx = torch.fft.fftfreq(N, d=1.0 / N, device=dev).float()
kz = torch.fft.rfftfreq(N, d=1.0 / N, device=dev).float()
KX = kx.view(N, 1, 1); KY = kx.view(1, N, 1); KZ = kz.view(1, 1, -1)
KS = [KX, KY, KZ]
K2 = KX ** 2 + KY ** 2 + KZ ** 2
K2c = K2.clone(); K2c[0, 0, 0] = 1.0
KMAG = K2.sqrt()
kmax = N // 3
mask = ((KX.abs() <= kmax) & (KY.abs() <= kmax) & (KZ.abs() <= kmax)).float()


def fft(f):  return torch.fft.rfftn(f, dim=(1, 2, 3))
def ifft(F): return torch.fft.irfftn(F, dim=(1, 2, 3), s=(N, N, N))
def ifft1(F): return torch.fft.irfftn(F, dim=(0, 1, 2), s=(N, N, N))   # single 3D field
def project(F):
    kdu = KX * F[0] + KY * F[1] + KZ * F[2]
    return torch.stack([F[i] - KS[i] * kdu / K2c for i in range(3)])
def vort_hat(F):
    return torch.stack([1j * (KY * F[2] - KZ * F[1]), 1j * (KZ * F[0] - KX * F[2]),
                        1j * (KX * F[1] - KY * F[0])])
def rhs(F):
    uu = ifft(F); ww = ifft(vort_hat(F))
    nl = torch.stack([uu[1] * ww[2] - uu[2] * ww[1], uu[2] * ww[0] - uu[0] * ww[2],
                      uu[0] * ww[1] - uu[1] * ww[0]])
    return project(fft(nl) * mask) - NU * K2 * F


def cheap_diag(F):
    uu = ifft(F)
    # velocity gradient A_ij = d u_i/d x_j
    A = [[ifft1(1j * KS[j] * F[i]) for j in range(3)] for i in range(3)]
    wx = A[2][1] - A[1][2]; wy = A[0][2] - A[2][0]; wz = A[1][0] - A[0][1]
    w2 = wx * wx + wy * wy + wz * wz
    S = [[0.5 * (A[i][j] + A[j][i]) for j in range(3)] for i in range(3)]
    # enstrophy production  P = sum_ij w_i S_ij w_j
    wv = [wx, wy, wz]
    Pfield = sum(wv[i] * S[i][j] * wv[j] for i in range(3) for j in range(3))
    Smag2 = 2.0 * sum(S[i][j] * S[i][j] for i in range(3) for j in range(3))
    E = float(0.5 * (uu ** 2).sum(0).mean())
    ens = float(0.5 * w2.mean())
    wmax = float(w2.max().sqrt())
    P = float(Pfield.mean())
    # geometric-depletion ratio in high-|w| regions: <P_field> / <|w|^2 |S|> (1 = no depletion, <<1 = depleted)
    denom = float((w2 * Smag2.sqrt()).mean()) + 1e-30
    depl = float(Pfield.mean()) / denom
    return E, ens, wmax, P, depl, (A, S, wv, w2)


def alignment(A, S, wv, w2):
    # enstrophy-weighted alignment of omega-direction with strain eigenvectors (CPU eigh; Blackwell
    # cusolver batched eigh is broken). Restrict to the high-enstrophy points (where it matters).
    Sm = torch.stack([torch.stack([S[i][j] for j in range(3)], -1) for i in range(3)], -2).reshape(-1, 3, 3)
    wn = torch.stack(wv, -1).reshape(-1, 3); w2f = w2.reshape(-1)
    K = 300000
    if Sm.shape[0] > K:
        idx = torch.topk(w2f, K).indices                       # the most intense vorticity
        Sm, wn, w2f = Sm[idx], wn[idx], w2f[idx]
    Sn = Sm.detach().cpu().numpy().astype(np.float64)
    wnp = wn.detach().cpu().numpy(); w2n = w2f.detach().cpu().numpy()
    wnp = wnp / (np.linalg.norm(wnp, axis=-1, keepdims=True) + 1e-12)
    ev, evec = np.linalg.eigh(Sn)                              # ascending eval; evec[:,:,k]
    cos = np.abs((wnp[:, :, None] * evec).sum(1))              # |cos(omega, e_k)|  (K,3)
    wgt = w2n / (w2n.sum() + 1e-30)
    al = [float((cos[:, k] * wgt).sum()) for k in range(3)]    # [<|cos e_min|>, <|cos e_mid|>, <|cos e_max|>]
    return al, None


def spectrum(F):
    e = 0.5 * (F.abs() ** 2).sum(0)                 # 0.5|u_hat|^2 per mode (one-sided rfft; factor ~ ok for shape)
    kbin = KMAG.round().long().flatten()
    Ek = torch.zeros(int(KMAG.max()) + 2, device=dev).scatter_add_(0, kbin, e.flatten())
    return Ek.cpu().numpy()


Fh = project(fft(u)); del u
nsteps = int(T / DT); bkm = 0.0; prev_wmax = None
snap_times = {2.0, 4.0, 6.0, 8.0, 9.0, 10.0}
print(f"{'t':>5} {'energy':>10} {'enstrophy':>10} {'P(stretch)':>11} {'max|w|':>8} {'BKM_int':>8} {'depl':>6}", flush=True)
t0 = time.time()
for s in range(nsteps + 1):
    t = s * DT
    if s % 20 == 0 or s == nsteps:
        E, ens, wmax, P, depl, grads = cheap_diag(Fh)
        if prev_wmax is not None:
            bkm += 0.5 * (wmax + prev_wmax) * (DT * 20)     # trapezoid BKM integral
        prev_wmax = wmax
        line = f"{t:5.2f} {E:10.4e} {ens:10.4e} {P:11.4e} {wmax:8.3f} {bkm:8.3f} {depl:6.3f}"
        if round(t, 3) in snap_times:
            al, _ = alignment(*grads)
            line += f"   align[e3min,e2mid,e1max]=[{al[0]:.2f},{al[1]:.2f},{al[2]:.2f}]"
        print(line, flush=True)
    if s == nsteps:
        break
    k1 = rhs(Fh); k2 = rhs(Fh + 0.5 * DT * k1); k3 = rhs(Fh + 0.5 * DT * k2); k4 = rhs(Fh + DT * k3)
    Fh = Fh + (DT / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

print(f"\n{nsteps} steps / {time.time()-t0:.0f}s.  Spectrum (final), E(k) at k=1,2,4,8,16,32:")
Ek = spectrum(Fh)
for k in [1, 2, 4, 8, 16, 32]:
    if k < len(Ek): print(f"   k={k:3d}: {Ek[k]:.3e}")
print("\nReading: P>0 = net vortex stretching (enstrophy growth); align e2mid highest = the self-stretching")
print("fingerprint; depl<<1 = geometric depletion active; BKM_int finite & leveling = staying regular.")
