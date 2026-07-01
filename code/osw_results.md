# OSW 1D model family — advection vs blow-up

omega_t + a*u*omega_x = u_x*omega,  u_x=H[omega],  data=sin(x), N=4096, T=14.0

| a | outcome | t_blow | peak |omega| | amplification |
|---|---|---|---|---|
| +0.820 | BLOW-UP | 6.3769 | 2.007e+04 | 2.01e+04 |
| +0.860 | BLOW-UP | 7.8085 | 2.017e+04 | 2.02e+04 |
| +0.900 | BLOW-UP | 10.5006 | 2.020e+04 | 2.02e+04 |
| +0.940 | regular | - | 3.575e+00 | 3.57e+00 |
| +0.980 | regular | - | 1.335e+00 | 1.33e+00 |

**Blow-up <-> regularity boundary** lies in a in (0.900, 0.940): advection strength a depletes the CLM blow-up between these.
