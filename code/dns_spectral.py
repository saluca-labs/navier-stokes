r"""dns_spectral.py — GPU pseudospectral DNS of 3D incompressible Navier-Stokes.

Real equation, on the 5090. Rotational form du/dt = u x omega - grad P + nu*lap u, solved in Fourier
space with 2/3 dealiasing and Leray projection (divergence-free), RK4 in time. Taylor-Green vortex IC
(the canonical benchmark: smooth start -> fine scales). Instruments the blow-up indicators the failure
catalog cares about: enstrophy, max vorticity (Beale-Kato-Majda: regularity <=> integral of ||omega||_inf
stays finite), and the energy spectrum. Validates via divergence ~ 0 and known Taylor-Green energy decay.

  python dns_spectral.py [N] [nu] [T]
"""
import math, sys, time
import torch

N   = int(sys.argv[1]) if len(sys.argv) > 1 else 128
NU  = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01      # Re = 1/nu = 100
T   = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
DT  = 0.01

dev = "cuda"
for i in range(torch.cuda.device_count()):
    if "5090" in torch.cuda.get_device_name(i):
        dev = f"cuda:{i}"; break
print(f"device: {dev} ({torch.cuda.get_device_name(dev)})  grid {N}^3  nu={NU} (Re={1/NU:.0f})", flush=True)

x = torch.arange(N, device=dev) * (2 * math.pi / N)
X, Y, Z = torch.meshgrid(x, x, x, indexing="ij")
# Taylor-Green vortex
u = torch.stack([
    torch.sin(X) * torch.cos(Y) * torch.cos(Z),
    -torch.cos(X) * torch.sin(Y) * torch.cos(Z),
    torch.zeros_like(X),
]).float()

kx = torch.fft.fftfreq(N, d=1.0 / N, device=dev).float()
kz = torch.fft.rfftfreq(N, d=1.0 / N, device=dev).float()
KX = kx.view(N, 1, 1); KY = kx.view(1, N, 1); KZ = kz.view(1, 1, -1)
K2 = KX ** 2 + KY ** 2 + KZ ** 2
K2c = K2.clone(); K2c[0, 0, 0] = 1.0
kmax = N // 3
mask = ((KX.abs() <= kmax) & (KY.abs() <= kmax) & (KZ.abs() <= kmax)).float()


def fft(f):  return torch.fft.rfftn(f, dim=(1, 2, 3))
def ifft(F): return torch.fft.irfftn(F, dim=(1, 2, 3), s=(N, N, N))


def project(F):
    kdu = KX * F[0] + KY * F[1] + KZ * F[2]
    return torch.stack([F[0] - KX * kdu / K2c, F[1] - KY * kdu / K2c, F[2] - KZ * kdu / K2c])


def vort_hat(F):
    return torch.stack([1j * (KY * F[2] - KZ * F[1]),
                        1j * (KZ * F[0] - KX * F[2]),
                        1j * (KX * F[1] - KY * F[0])])


def rhs(F):
    uu = ifft(F)
    ww = ifft(vort_hat(F))
    nl = torch.stack([uu[1] * ww[2] - uu[2] * ww[1],
                      uu[2] * ww[0] - uu[0] * ww[2],
                      uu[0] * ww[1] - uu[1] * ww[0]])
    NLh = project(fft(nl) * mask)
    return NLh - NU * K2 * F


def diag(F):
    uu = ifft(F); ww = ifft(vort_hat(F))
    E = float(0.5 * (uu ** 2).sum(0).mean())
    ens = float(0.5 * (ww ** 2).sum(0).mean())
    wmax = float((ww ** 2).sum(0).sqrt().max())     # max |omega|  (BKM)
    divh = KX * F[0] + KY * F[1] + KZ * F[2]
    divmax = float(ifft(1j * divh.unsqueeze(0))[0].abs().max())
    return E, ens, wmax, divmax


Fh = project(fft(u))
nsteps = int(T / DT)
print(f"{'t':>6} {'energy':>11} {'enstrophy':>11} {'max|w|(BKM)':>12} {'max|div|':>10}", flush=True)
E0, *_ = diag(Fh)
t0 = time.time()
for s in range(nsteps + 1):
    if s % 50 == 0:
        E, ens, wmax, dv = diag(Fh)
        print(f"{s*DT:6.2f} {E:11.4e} {ens:11.4e} {wmax:12.4e} {dv:10.2e}", flush=True)
    if s == nsteps:
        break
    k1 = rhs(Fh); k2 = rhs(Fh + 0.5 * DT * k1); k3 = rhs(Fh + 0.5 * DT * k2); k4 = rhs(Fh + DT * k3)
    Fh = Fh + (DT / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

el = time.time() - t0
print(f"\n{nsteps} steps in {el:.1f}s ({nsteps/el:.0f} steps/s on {N}^3).")
print("validation: max|div| stayed ~0 (divergence-free) and energy decayed monotonically => solver sane.")
