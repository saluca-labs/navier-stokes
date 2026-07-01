# OSW 1D model family — advection vs blow-up

omega_t + a*u*omega_x = u_x*omega,  u_x=H[omega],  data=sin(x), N=4096, T=40.0

| a | outcome | t_blow | peak |omega| | amplification |
|---|---|---|---|---|
| +0.920 | BLOW-UP | 14.2226 | 2.019e+04 | 2.02e+04 |
| +0.940 | BLOW-UP | 29.8525 | 2.014e+04 | 2.01e+04 |
| +0.960 | regular | - | 2.104e+00 | 2.10e+00 |
| +0.980 | regular | - | 1.388e+00 | 1.39e+00 |

**Blow-up <-> regularity boundary** lies in a in (0.940, 0.960): advection strength a depletes the CLM blow-up between these.
