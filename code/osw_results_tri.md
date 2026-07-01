# OSW 1D model family — advection vs blow-up

omega_t + a*u*omega_x = u_x*omega,  u_x=H[omega],  data=sin(x), N=4096, T=20.0

| a | outcome | t_blow | peak |omega| | amplification |
|---|---|---|---|---|
| +0.000 | BLOW-UP | 2.7013 | 2.001e+04 | 2.00e+04 |
| +1.000 | regular | - | 1.000e+00 | 1.00e+00 |
| +1.250 | regular | - | 1.000e+00 | 1.00e+00 |

**Blow-up <-> regularity boundary** lies in a in (0.000, 1.000): advection strength a depletes the CLM blow-up between these.
