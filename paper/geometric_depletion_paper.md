---
title: "Geometric depletion of enstrophy production in the extreme vorticity tail: a Reynolds-number study toward the Navier–Stokes regularity question"
author: "Cristian Ruvalcaba and the Agentic AI Research Team · Saluca Labs"
date: "Draft v0.6 — 2026-07-01"
geometry: margin=1in
fontsize: 11pt
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
---

*Status (v0.6, a correction-and-extension revision): the extreme-tail depletion ratio remains **resolved
below unity** at the top-$0.1\%$ and top-$1\%$ at both $R_\lambda\approx1300$ and $R_\lambda=2500$
(contiguous-block 95% $t$-intervals exclude 1; the result is $3.8$–$10.7$ standard errors below unity and
survives every correction in this revision). The v0.5 claim that the **cross-Reynolds weakening of tail
depletion was "statistically resolved at top-$1\%$" is withdrawn**: an independent adversarial verification
pass found (i) the block intervals used a normal rather than $t_{19}$ multiplier — under the correct
multiplier the intervals touch — and (ii) the asserted "conservative" direction of the resolution confound is
**wrong**: a new paired coarse-graining experiment (§6.1) measures the bias directly and finds it runs *with*
the observed trend, absorbing $\approx60\%$ of the top-$1\%$ difference. The corrected statement is a
**directionally Buaria–Pumir-consistent difference that is not statistically resolved at top-$1\%$ once the
measured resolution bias is applied** (suggestive at top-$0.1\%$). §7 is reframed from "corroboration" to
"illustration" and now credits Lushnikov–Silantyev–Siegel's prior periodic-domain threshold [LSS21], which
our scan reproduces. A full change log is given before the references.*

---

## Abstract

The global regularity of the three-dimensional incompressible Navier–Stokes equations turns on whether the
vortex-stretching term can drive the enstrophy to infinity in finite time. A classical obstruction to blow-up
is **geometric depletion** (Constantin–Fefferman): the nonlinear production of enstrophy, $\omega\cdot S\omega$,
is suppressed below its dimensional bound $|\omega|^2|S|$ because vorticity preferentially aligns with the
*intermediate* eigenvector of the strain tensor. Whether this depletion *persists* in the most intense
vorticity events — the events that would actually seed a singularity — and how it scales with Reynolds number,
is the quantity that matters. We build a validated pseudospectral DNS and a direct-data pipeline over the Johns
Hopkins Turbulence Database (JHTDB) to measure the normalized production ratio
$D \equiv \langle\omega\cdot S\omega\rangle/\langle|\omega|^2|S|\rangle$ as a function of vorticity percentile,
from the bulk into the extreme tail. We reproduce the canonical intermediate-eigenvector alignment and a mean
depletion $D\!\approx\!0.10$. At Taylor-microscale Reynolds number $R_\lambda\!\approx\!433$, the mean depletion
is firmly established ($D\!\approx\!0.12$) and the tail point estimates decrease into the extreme ($D_{\rm
tail}/D_{\rm mean}\!\to\!0.77$–$0.81$ in a $4.1\times10^5$-point sample, for which no resampled interval is
available); bootstrap intervals on a separate, smaller ($6\times10^4$-point) sample that we can fully resample
straddle unity at the top-$0.1\%$, so the extreme-tail ratio is not statistically resolved there. At higher
Reynolds number we now resolve it. Executing a high-statistics, per-snapshot, single-$R_\lambda$ measurement —
contiguous cutouts of JHTDB `isotropic8192` ($R_\lambda\!\approx\!1300$, snapshots 0–4 only, 340 cutout reads,
$\sim\!7.7\times10^7$ interior points) and `isotropic32768` (the exascale $32768^3$ field, $R_\lambda\!=\!2500$,
86 cutouts of $128^3$, $\sim\!1.6\times10^8$ points) — and bounding the tail with contiguous-block (batch-means)
95% $t$-intervals that respect within-cutout spatial correlation, we find the extreme-tail ratio is **resolved
below unity** at both Reynolds numbers: $D_{\rm tail}/D_{\rm mean}=0.91\,[0.89,0.93]$ (top-$1\%$) and
$0.87\,[0.82,0.92]$ (top-$0.1\%$) at $R_\lambda\!\approx\!1300$, and $0.94\,[0.93,0.96]$, $0.93\,[0.90,0.97]$ at
$R_\lambda\!=\!2500$ — i.e. the most intense vorticity remains *more* geometrically depleted than the bulk
(regularizing) at both. Comparing the two, the tail depletion is **directionally weaker at the higher Reynolds
number** (two-sample $p\!\approx\!0.003$ at top-$1\%$ before systematics) — the Buaria–Pumir direction — but a
new paired coarse-graining experiment (§6.1) **measures the resolution bias directly and finds it runs *with*
this trend**, absorbing $\approx60\%$ of the top-$1\%$ difference; after the measured correction the
cross-Reynolds difference is **not statistically resolved** at top-$1\%$ (suggestive at top-$0.1\%$). We
therefore report a two-point comparison whose attribution between $R_\lambda$ and small-scale resolution
($k_{\max}\eta\!\approx\!1.4$ vs $2.8$) remains open; the deepest tail (top-$0.01\%$ and below) is not
gradient-resolution-converged (§3.2). A stencil-order cross-check (on-grid 4th-order vs JHTDB 6th/8th-order/spline)
confirms the *stencil choice* contributes $<1\%$ through top-$0.1\%$; the coarse-graining experiment bounds the
*grid-resolution* systematic that stencil agreement cannot see. We do not claim priority on the
Reynolds-dependence of tail depletion — Buaria and co-workers established the direction at higher fidelity; our
contribution is an independent public-data measurement with explicit, measured uncertainty and resolution
systematics, plus the $R_\lambda\!=\!2500$ datum. We also report the confirmed $R_\lambda\!\approx\!433$ baseline
and a decidable 1D-model illustration of the depletion mechanism (§7).

---

## 1. Introduction

The Millennium Problem for 3D incompressible Navier–Stokes asks whether smooth, finite-energy initial data can
develop a singularity in finite time, or whether solutions remain smooth for all time. The Beale–Kato–Majda
criterion [BKM84] reduces this to the time-integrability of the vorticity maximum: a solution remains smooth so
long as $\int_0^T \|\omega(\cdot,t)\|_\infty\,dt < \infty$. The enstrophy $\Omega(t)=\tfrac12\int|\omega|^2$
evolves by

$$\frac{d\Omega}{dt} = \int \omega\cdot S\omega \;-\; \nu\int|\nabla\omega|^2,$$

(exact for unforced periodic flow; the forced, statistically steady flows analyzed below add a forcing-input
term to their steady balance, which does not affect the pointwise production statistic studied here), so
blow-up requires the **production** term $\int\omega\cdot S\omega$ to outrun dissipation. Dimensionally
$\omega\cdot S\omega \le |\omega|^2|S|$, and if this bound were saturated the enstrophy growth would be
supercritical with respect to the energy-controlled scaling — the heart of Tao's *supercriticality barrier*
[Tao16]. The reason 3D turbulence does not (apparently) blow up is that the bound is **not** saturated:
Constantin–Fefferman [CF93] showed that sufficient smoothness of the *vorticity direction*
$\xi=\omega/|\omega|$ in regions of high vorticity rules out singularities, and DNS [AKKG87, JWSR93] established
that vorticity aligns preferentially with the **intermediate** strain eigenvector, whose eigenvalue is often
small or positive but typically far below the most-extensional eigenvalue. The result is *geometric depletion*
of production.

The regularity-relevant question is not the *mean* depletion but its behaviour in the **extreme tail** of the
vorticity distribution — the events heuristically closest to a putative singularity, though we stress that
extremes of a forced statistically steady flow are not proto-singularities — and its **Reynolds-number
trend**. Recent high-$R_\lambda$ DNS by Buaria and co-workers
[BBP20, BPBY19, BP22] found that the strongest enstrophy production derives from the *most-extensive* strain
direction (driven by nonlocal strain), and reported a *mild* reduction of geometric depletion confined to the
deepest tail at the highest $R_\lambda$ — a subtle effect, not a wholesale "de-depletion." If such a reduction
grows with $R_\lambda$, the geometric protection might weaken precisely where it
matters. This paper measures the percentile-resolved depletion ratio $D_{\rm tail}/D_{\rm mean}$ across Reynolds
number, using both in-house DNS and direct analysis of the JHTDB datasets. We confirm the regularizing
(ratio $<1$, strengthening into the tail) behaviour at $R_\lambda\approx433$, and — with a high-statistics,
per-snapshot, single-$R_\lambda$ measurement on $\sim\!10^8$-point cutout samples (§6) — resolve the extreme-tail
ratio below unity at the top-$0.1\%$ and top-$1\%$ at both $R_\lambda\approx1300$ and $R_\lambda=2500$. Across
these two Reynolds numbers the tail depletion is directionally weaker at the higher $R_\lambda$ — the direction
reported by Buaria and co-workers — but §6.1's paired coarse-graining experiment shows the datasets'
resolution difference biases the comparison in the *same* direction, and after the measured correction the
difference is not statistically resolved at top-$1\%$. Because only two Reynolds numbers are sampled and they
differ in small-scale resolution, we treat this as a two-point comparison whose $R_\lambda$-vs-resolution
attribution is open; the deepest tail (top-$0.01\%$ and below) remains resolution-limited (§6). We do not claim
priority on the Reynolds trend itself [BBP20, BPBY19, BP22]; our contribution is the percentile-resolved
statistic with explicit uncertainty discipline, the measured resolution systematic, and the $R_\lambda=2500$
datum.

## 2. The geometric-depletion statistic

Let $A_{ij}=\partial_j u_i$ be the velocity-gradient tensor, $S=\tfrac12(A+A^\top)$ the strain-rate tensor, and
$\omega=\nabla\times u$ the vorticity. The pointwise enstrophy production is $P=\omega_i S_{ij}\omega_j$ and the
strain magnitude $|S|=\sqrt{2\,S_{ij}S_{ij}}$. We define the **production efficiency** (geometric-depletion
statistic)

$$D \;=\; \frac{\langle \omega\cdot S\omega\rangle}{\langle |\omega|^2\,|S|\rangle},$$

an enstrophy-weighted, dimensionless ratio. With the **Frobenius** strain norm $|S|=\sqrt{2S_{ij}S_{ij}}$ used
here, the pointwise efficiency $\omega\!\cdot\!S\omega/(|\omega|^2|S|)$ is bounded by $[-1/\sqrt3,\,+1/\sqrt3]$
(extremum at axisymmetric extension/compression with $\omega$ along the distinguished eigenvector); $D\to1/\sqrt3$
would mean perfect alignment with the most-extensional eigenvector (no depletion), and the observed
$D\approx0.10$–$0.12$ reflects the intermediate-eigenvector preference (strong depletion). **Normalization note (corrected in v0.6):**
because we use the Frobenius norm rather than the strain operator-norm (largest $|\lambda|$), our $D$ is smaller
than operator-norm-normalized efficiencies by the ratio $\sqrt{2\sum_i\lambda_i^2}/\max_i|\lambda_i|$, which for
a traceless strain is **hard-bounded in $[\sqrt3,\,2]$** ($\sqrt3$ at axisymmetric extension/compression, $2$ at
planar strain; $\approx\!1.8$ for the canonical $3\!:\!1\!:\!-4$ state). Our $D_{\rm mean}\approx0.12$ therefore
corresponds to at most $\approx\!0.24$ under an operator norm — consistent with the lower end of values quoted
elsewhere ($\sim\!0.2$; cf. Tsinober; values derivable from [BBP20] Table II), while higher quoted efficiencies
($\sim\!0.3$–$0.4$) reflect *different conditioning or statistics*, not merely a normalization change. The
diagnostic of interest is
the **conditional** version: restrict the averages to the top-$q$ fraction of points ranked by $|\omega|^2$, and
form the ratio $D(q)/D_{\rm mean}$. A ratio **below 1** means the extreme tail is *more* depleted than the bulk
(regularizing); a ratio **above 1** means extreme **de-depletion** (the singularity-relevant direction).

## 3. Methods

**3.1 In-house DNS.** A standard Fourier pseudospectral solver (rotational form, $2/3$ dealiasing, exact Leray
projection, RK4 time stepping) integrates the incompressible NS equations on a triply periodic box, GPU-resident
(PyTorch FFT, single RTX 5090). Validation is on the Taylor–Green vortex; resolution is controlled by
$k_{\max}\eta$.

**3.2 JHTDB direct-data pipeline.** For Reynolds numbers beyond single-GPU reach we analyze existing
state-of-the-art DNS in the Johns Hopkins Turbulence Database [JHTDB08] via the `giverny` interface inside
SciServer Compute, where the datasets are mounted locally. Two access modes: (i) JHTDB higher-order
velocity-gradient queries (`getData`, 4th–8th-order FD stencils with Lagrange interpolation) at points;
(ii) contiguous sub-volume *cutouts* with the
velocity-gradient computed on the grid by 4th-order central differences (interior, halo-trimmed). Datasets:
`isotropic1024coarse` ($R_\lambda\approx433$); `isotropic8192`, whose six snapshots are *not* a single Reynolds
number — snapshots 0–4 at $R_\lambda\approx1200$–$1300$ and snapshot 5 at $R_\lambda\approx610$, per the JHTDB
dataset documentation [JHTDB-i8192] and confirmed by the $\sim\!6\times$ enstrophy gap we measure between them
(§6); and `isotropic32768`, a single snapshot of the exascale $32768^3$ forced-isotropic DNS of Ravikumar, Yeung
and co-workers [Yeung32768], at $R_\lambda=2500$ with exceptionally fine small-scale resolution
$k_{\max}\eta\approx2.8$ (versus $k_{\max}\eta\approx1.4$ for `isotropic8192`) — i.e. roughly double the
Reynolds number *and* better resolved, a covariance we account for in §6.
**Differentiation order is not a tail-error source (verified).** The most intense velocity gradients live at
the highest resolved wavenumbers, where a low-order finite-difference stencil could in principle deviate from
the true derivative — making the extreme-tail statistic the part most exposed to *differentiation-order* error.
JHTDB serves pointwise gradients through finite-difference stencils (4th/6th/8th order; methods
`fd4noint`,`fd6noint`,`fd8noint`,`fd4lag4`,…) rather than a pointwise spectral derivative, so the relevant
validation is an **order-convergence check**. We performed it as a paired test at $\sim\!5\times10^5$ *identical*
grid points (300 boxes per dataset), comparing our on-grid 4th-order gradient against the database's 6th- and
8th-order finite-difference stencils and its `m2q8` spline (`fd4noint`,`fd6noint`,`fd8noint`,`m2q8`), on **both**
`isotropic8192` and `isotropic32768`. The per-point vorticity correlates at $0.9998$–$1.0000$, and the
conditional depletion ratios are flat across the order ladder through the **top-$0.1\%$**: on `isotropic8192`,
$D_{\rm mean}$ ratio FD/reference $=1.000$–$1.003$ and the top-$0.1\%$ ratio agrees to $<1\%$
($0.778$ vs $0.776$–$0.780$ across fd6/fd8/m2q8); on `isotropic32768` likewise ($D_{\rm mean}$ ratio $1.000$,
top-$0.1\%$ $0.670$ vs $0.673$–$0.676$). **Scope of this check (rescoped in v0.6):** it establishes that the
*stencil choice* contributes $<1\%$ through top-$0.1\%$ — on-grid 4th-order FD is as good as the database's
higher-order operators at that depth. It does **not** establish grid-resolution adequacy: all stencils share the
same grid and Nyquist, so they agree identically on whatever sub-grid gradient content the *grid itself*
attenuates (a common-mode error this test is structurally blind to). That grid-resolution systematic is instead
**measured directly** by the paired coarse-graining experiment of §6.1, which finds a real $\approx\!-0.02$
common-mode bias on the top-$1\%$ ratio at `isotropic8192`-equivalent resolution. (The small stencil-gate
samples here also carry a different *sampling design* — $\sim\!5\times10^5$ scattered points with a pooled tail
threshold vs §6's per-cutout thresholds on $10^8$ points — which shifts the conditional-ratio *magnitude*
systematically; the $0.778$ here vs $0.87$ in §6 for the same dataset reflects that estimand difference, not a
contradiction.) At the **top-$0.01\%$** the stencils scatter $5$–$9\%$ and the *sign* of the disagreement flips
between the two datasets — the signature of a noise-limited comparison at $\sim\!50$ tail points per method with
independently selected tails — so we do not claim the gradient is resolved at top-$0.01\%$ or below.

## 4. Results — validation and the moderate-Reynolds arc

**4.1 Alignment and mean depletion.** The in-house DNS reproduces the canonical PDF of the cosine between
$\omega$ and the strain eigenvectors, peaked on the **intermediate** eigenvector, and a mean depletion
$D_{\rm mean}\approx0.10$ across $\mathrm{Re}=400\text{–}3200$ (Taylor–Green). This recovers [AKKG87, JWSR93]
from scratch and validates the pipeline.

**4.2 Taylor–Green Reynolds arc.** In single-realization in-house DNS (one flow per Reynolds number,
$k_{\max}\eta\ge1.1$ at $512^3$), the extreme-tail depletion (top $0.01\%$ of $|\omega|$) trends *downward* with
Reynolds number through $\mathrm{Re}=3200$ ($D_{\rm tail}$: $0.22\to0.11\to0.09$), converging toward the mean
from above — i.e. at these moderate Reynolds numbers both the bulk and the resolved extremes trend toward *more*
regularity. We caution that, as the high-$R_\lambda$ audit of §6 makes vivid, single-realization extreme-tail
estimates carry large statistical uncertainty; we read this arc as *suggestive* of strengthening depletion at
moderate Re, not a converged trend — noting also that at fixed grid, higher Re means coarser $k_{\max}\eta$,
and §6.1's measured bias direction implies these single-realization tail values read systematically *low* at
the higher Re, which would exaggerate exactly this apparent strengthening. An apparent extreme de-depletion at $\mathrm{Re}=800$ was identified as a
low-Reynolds transient, and a
$512^3$ control confirmed it was **not** a resolution artifact.

## 5. Results — high-Reynolds isotropic baseline ($R_\lambda\approx433$)

Direct analysis of JHTDB `isotropic1024coarse` ($4.1\times10^5$ points) gives $D_{\rm mean}=0.117$ and tail
ratios whose point estimates decrease into the tail (statistical uncertainty quantified below):

| top-$q$ of $|\omega|$ | $D(q)$ | $D(q)/D_{\rm mean}$ | points |
|---|---|---|---|
| 10% | 0.113 | 0.96 | 41 369 |
| 1% | 0.105 | 0.89 | 4 136 |
| 0.1% | 0.095 | 0.81 | 413 |
| 0.05% | 0.090 | **0.77** | 206 |
| 0.01% | 0.091 | 0.78 | 41 |

These point estimates suggest the most intense vorticity is *more* depleted than the bulk — a regularizing
tendency — and several structural facts are robust: $24.2\%$ of points have *negative* production
($\omega\cdot S\omega<0$, local enstrophy destruction); the top $1\%$ of enstrophy carries only $\sim\!26\%$ of
total production; intermittency $\max|\omega|^2/\langle|\omega|^2\rangle\approx135$.

**We now quantify the tail uncertainty honestly, and it is large.** Bootstrapping ($2000$ resamples) the
archived independent $6\times10^4$-point sample of the same dataset gives a tight $D_{\rm mean}=0.121$ but tail
ratios with wide intervals that **straddle unity**: $0.974$ $[0.913,1.026]$ (top $5\%$), $0.977$ $[0.875,1.077]$
(top $1\%$), and $0.915$ $[0.638,1.174]$ (top $0.1\%$). The $4.1\times10^5$-point sample's point estimates fall
lower ($0.77$–$0.81$ in the tail) with $\sim\!2.6\times$ tighter intervals, but we cannot currently bootstrap it
(the source data are not in our present mount). The honest status at $R_\lambda\approx433$ is therefore: the
**mean depletion is firmly established**; the **tail shows a regularizing tendency** (point estimates $<1$,
deepening into the tail in the large sample); but the **extreme top-$0.1\%$ ratio is not statistically resolved**
even here — its confidence interval includes unity in the sample we can fully bootstrap. A definitive
extreme-tail error bar requires a full bootstrap of the large sample and, ultimately, the higher-statistics
program of §9. This is the same convergence wall reaching even the most accessible Reynolds number.

![Geometric-depletion ratio $D(q)/D_{\rm mean}$ versus vorticity percentile at $R_\lambda\approx433$ (JHTDB `isotropic1024coarse`, $4.1\times10^5$ points). The ratio stays **below** the de-depletion threshold (dashed line at 1.0) and *decreases* into the extreme tail — i.e. the most intense vorticity is the most geometrically depleted (regularizing). The high-statistics $R_\lambda\approx1300$ and $R_\lambda=2500$ curves are resolved in §6 (Fig. 2).](../experiments/depletion_ratio_vs_Re.png){width=82%}

## 6. The high-Reynolds tail, resolved: $R_\lambda\approx1300$ and $R_\lambda=2500$

The prior revision of this paper reported a *cautionary null* here: at $\sim\!5\times10^5$ points the
`isotropic8192` extreme tail would not converge, and a data-provenance audit showed the dataset must be treated
per-snapshot — it mixes $R_\lambda\approx1200$–$1300$ (snapshots 0–4) with a $R_\lambda\approx610$ snapshot 5
whose mean enstrophy is $\sim\!6\times$ lower, so global pooling is invalid (including snapshot 5 gave a
top-$0.1\%$ ratio $\approx0.57$; excluding it flipped the pooled value to $\approx1.2$ — a sign change from
dropping $1/6$ of the data, the signature of an artifact). Those cautions stand and dictate the present
protocol. We now execute the high-statistics measurement that resolves the tail.

**Protocol.** Per dataset we draw many contiguous velocity cutouts, compute the on-grid 4th-order
velocity gradient (interior, halo-trimmed), and pool the pointwise $(|\omega|^2,\,P,\,|\omega|^2|S|)$ — strictly
single-$R_\lambda$: `isotropic8192` uses **only** the $R_\lambda\approx1300$ snapshots (0–4); `isotropic32768`
is a single $R_\lambda=2500$ snapshot. Pinning the sample composition exactly (v0.6): `isotropic8192` pools
**340 cutout reads** (predominantly $64^3$, accumulated across resumed runs; $N=7.678\times10^7$ interior
points — v0.5 loosely described these as $128^3$); `isotropic32768` pools **86 cutouts of $128^3$**
($N=1.640\times10^8$). Because points inside one cutout are spatially correlated, the raw $N$ overstates the
independent sample; we therefore bound the tail ratio with a **contiguous-block (batch-means) estimator** —
splitting each pool into $K=20$ contiguous blocks, computing the conditional ratio *within* each block (tail
re-selected per block, $D_{\rm mean}$ recomputed per block, so the correlated numerator and denominator are
handled inside the ratio), and quoting the across-block mean $\pm\,t_{19}\cdot{\rm SEM}$ — rather than a point
bootstrap, which would give falsely tight intervals. Two disclosures (v0.6): (i) v0.5 called this a "block
resample" and used a normal ($1.96$) multiplier; it is a batch-means partition (no resampling), and with
$K=20$ blocks the correct multiplier is $t_{19}=2.093$ ($\approx\!7\%$ wider) — the table below is
$t$-corrected. (ii) The across-block mean of within-block ratios is a slightly different estimand from the
pooled ratio-of-sums (each block's top-$q$ threshold is shallower than the pool's); both are reported below and
the difference is within the quoted intervals. Block boundaries are cut at equal point counts and can split a
cutout; the effect is mildly anti-conservative and small at $K=20$. We also note the family of tests: below-unity
is assessed at 2 datasets $\times$ 2 depths (all far from marginal), the cross-Re difference at 2 depths.

**Result — the tail is resolved below unity at both Reynolds numbers.**

| | $D_{\rm mean}$ | $D(q)/D_{\rm mean}$, top-$1\%$ | $D(q)/D_{\rm mean}$, top-$0.1\%$ |
|---|---|---|---|
| `isotropic8192` ($R_\lambda\approx1300$) | $0.116$ | $0.91\;[0.89,\,0.93]$ (pooled $0.909$) | $0.87\;[0.82,\,0.92]$ (pooled $0.859$) |
| `isotropic32768` ($R_\lambda=2500$) | $0.119$ | $0.94\;[0.93,\,0.96]$ (pooled $0.955$) | $0.93\;[0.90,\,0.97]$ (pooled $0.938$) |

(Brackets are contiguous-block batch-means 95% $t_{19}$-intervals; centers are across-block means, with the
pooled ratio-of-sums in parentheses.) At both Reynolds numbers the conditional ratio is **resolved below 1**
at top-$1\%$ and top-$0.1\%$ — the intervals exclude unity by $3.8$–$10.7$ standard errors, robust to the
multiplier and to $K\in\{5,20,40\}$ — so the most intense vorticity is *more* geometrically depleted than the
bulk, the regularizing direction, established at $R_\lambda$ up to $2500$. The ratio decreases monotonically
into the tail at both. This is the measurement the previous revision could not make; the
$\sim\!10^2$–$10^3\times$ larger, provenance-clean sample is what converts the order-unity per-snapshot scatter
of the old data into stable, interval-bounded estimates. One honesty note on the interval itself: the
`isotropic32768` blocks are spatial sub-samples of a **single snapshot**, so the interval quantifies
within-realization sampling scatter under an ergodic assumption, not realization-to-realization variability
(`isotropic8192` has 5 independent snapshots); and the block SEM is not $K$-independent (it grows mildly with
finer partitions, indicating residual heterogeneity across cutouts), so the below-unity margins — which survive
even a doubled SEM — are trustworthy while marginal comparisons built on these intervals are not.

![Conditional geometric-depletion ratio $D(q)/D_{\rm mean}$ versus vorticity percentile for the two
high-Reynolds JHTDB datasets, with block-resampled 95% bands ($K=20$ contiguous blocks). Both curves lie
**below** the de-depletion threshold (dashed, $=1$) and decrease into the tail — the most intense vorticity is
the most geometrically depleted (regularizing) — at both $R_\lambda\approx1300$ and $R_\lambda=2500$. The
higher-Reynolds curve sits systematically **above** the lower (less depleted) — a difference that is
directionally consistent with Buaria–Pumir but, after the measured resolution correction of §6.1, not
statistically resolved at top-$1\%$ (see text).
Solid markers are the claimed depth (top-$0.1\%$ and shallower); faded dotted markers beyond it are
resolution-limited (§3.2) and not claimed.](depletion_crossRe.png){width=86%}

**The Reynolds comparison (reworked in v0.6).** Between the two datasets the tail ratio is higher at the
higher $R_\lambda$: $0.91\to0.94$ at top-$1\%$ and $0.87\to0.93$ at top-$0.1\%$. v0.5 called the top-$1\%$
difference "statistically resolved" on the grounds that the two intervals did not overlap; that argument was a
knife-edge artifact of the normal multiplier — under the correct $t_{19}$ intervals they touch
($[0.889,0.925]$ vs $[0.925,0.962]$) — and CI-overlap is in any case the wrong test. The proper two-sample
comparison of block means gives, **before systematics**, a difference of $+0.036\pm0.012$ ($p\approx0.003$) at
top-$1\%$ and $+0.061\pm0.029$ ($p\approx0.037$) at top-$0.1\%$. The decisive question is then systematic, not
statistical: the two datasets differ in small-scale resolution as well as $R_\lambda$, and §6.1 **measures**
that resolution bias directly — it is $-0.018\pm0.005$ at top-$1\%$, i.e. it runs *with* the observed
difference and accounts for $\approx60\%$ of it. After the measured correction the top-$1\%$ cross-Reynolds
difference is $\approx\!+0.012$, comparable to its combined uncertainty — **not resolved**; at top-$0.1\%$ a
suggestive $\approx\!+0.043$ remains, with the caveat that the deep-tail bias magnitude is filter-shape
sensitive (§6.1). The direction — milder geometric depletion of the extreme tail at higher Reynolds number — is
what Buaria, Pumir and co-workers report from independent, better-resolved DNS [BBP20, BPBY19, BP22]; our data
are consistent with that trend but do not independently establish it.

### 6.1 The resolution bias, measured: a paired coarse-graining experiment (new in v0.6)

v0.5 *argued* that the resolution confound was conservative — that under-resolution smooths tail/bulk contrast
and pushes the conditional ratio toward 1, so the coarser `isotropic8192` would read artificially *less*
depleted, against the observed direction. Adversarial review identified a competing, literature-grounded
mechanism with the opposite sign (under-resolution preferentially misses the small-scale extensional extreme
events; filtered fields align *more* with the intermediate eigenvector), so we measured the bias directly.

**Design.** 23 regions of $256^3$ (assembled from contiguous $128^3$ cutouts) of `isotropic32768`, each
measured four ways *on the same field*: **V0** unfiltered ($k_{\max}\eta\approx2.8$); **V1** sharp-spectrally
low-passed to $k_c\eta\approx1.4$ (the `isotropic8192`-equivalent spectral content), gradients on the fine
grid; **V2** the same filter followed by $2\times$ downsampling, gradients on the coarse grid — the full
emulation of the `isotropic8192` measurement (content *and* stencil sampling); **V3** a Gaussian filter at
matched cutoff plus downsampling (filter-shape/windowing control). A generous margin (48 cells) is trimmed
after filtering to kill window-wraparound; paired per-region differences are quoted as mean
$\pm\,t_{22}\cdot{\rm SEM}$.

**Result — the v0.5 direction argument is falsified.** Coarse-graining *lowers* the tail ratio:

| paired $\Delta$ (variant $-$ V0) | top-$1\%$ | top-$0.1\%$ |
|---|---|---|
| **V2 (full emulation)** | $\mathbf{-0.018\pm0.005}$ (negative in 22/23 regions) | $-0.017\pm0.016$ |
| V1 (spectral content only) | $-0.009\pm0.005$ | $-0.012\pm0.010$ |
| V3 (Gaussian control) | $-0.011\pm0.006$ | $+0.002\pm0.022$ |

Both components — spectral truncation (V1) and coarser stencil sampling (V2$-$V1) — contribute with the same
sign; V3's agreement at top-$1\%$ rules out window-ringing as the driver, while its null at top-$0.1\%$ marks
the deep-tail bias magnitude as filter-shape sensitive. The coarse field reads systematically **more**
depleted, so the resolution confound runs *with* the observed cross-Reynolds difference. Applying the measured
bias to the `isotropic8192` point (an a-priori transfer — the bias is measured at $R_\lambda=2500$ statistics)
corrects $0.910\to\approx0.928$ at top-$1\%$: still clearly below unity (the below-unity headline *gains* a
measured systematic bound), but the cross-Reynolds difference shrinks from $+0.036$ to $\approx\!+0.012$,
inside its combined uncertainty. Caveats: an a-priori filtered fine DNS is not a coarser DNS (no nonlinear
feedback at the truncated scales); the bias could itself be Reynolds-dependent; and as an independent check,
V0 on these 23 fresh regions reproduces the §6 production measurement ($D_{\rm mean}=0.118$, pooled top-$1\%$
$0.958$, top-$0.1\%$ $0.931$). Code and pools are released with this revision (`coarse_grain_run.py`,
`pool_cgrain_*`).

**Why we do not call this a Reynolds-number scaling law.** Three reasons, stated plainly. (1) **Two points.**
Two Reynolds numbers establish a difference, not a trend or a law. (2) **A confounded second axis, now
measured.** `isotropic32768` is not only higher-$R_\lambda$ but better resolved ($k_{\max}\eta\approx2.8$ vs
$1.4$), and §6.1 shows the resolution bias is real, runs in the *same* direction as the observed difference,
and accounts for the majority of it at top-$1\%$ — so with two datasets the $R_\lambda$-vs-resolution
attribution cannot be closed from our data. (3) **The deepest tail is not resolved.** §3.2 shows the gradient
stencils scatter $5$–$9\%$ at top-$0.01\%$ (sign-flipping between datasets), so we restrict every claim to
top-$0.1\%$ and shallower. The cross-check that the gate (point-gradient `getData` boxes) gave the *opposite*
cross-$R_\lambda$ ordering on a $\sim\!300\times$ smaller sample — a disagreement whose plausible mechanism
(the sampling-design and estimand differences of §3.2, compounded by the §6.1 bias) we can now name but have
not fully decomposed — is itself a reminder that the second-decimal comparison is at the edge of what this
data supports.

## 7. Model-problem illustration: advection depletes blow-up (1D OSW family)

The 3D measurement is statistical and indirect — forced isotropic turbulence is a proxy, not a blow-up
scenario (§10). To exhibit the *same* competition — nonlinear stretching versus depleting transport — in a
setting where finite-time blow-up is **decidable**, we turn to the one-dimensional Okamoto–Sakajo–Wunsch (OSW)
family [OSW08]:

$$\omega_t + a\,u\,\omega_x = u_x\,\omega, \qquad u_x = H[\omega],$$

with $H$ the Hilbert transform and $a$ the advection strength. At $a=0$ this is the Constantin–Lax–Majda model
[CLM85], which possesses an *exact* finite-time singularity; at $a=1$ it is the De Gregorio model [DG90]. The
rigorous picture, which we take care to state, is regularity- and domain-dependent: for **smooth** data on the
**circle** (periodic domain, our setting) the De Gregorio advection is *conjectured* (and numerically observed)
to suppress blow-up, but **low-regularity ($C^\alpha$) data blows up for *every* $a$** — Chen, Hou & Huang
[CHH21] prove finite-time, self-similar singularities for $C^\alpha$ data at all $a$ (on $\mathbb{R}$ and $S^1$),
and Elgindi & Jeong [EJ20] establish $C^\alpha$ blow-up and isolate the advection-vs-stretching competition.
(On the line $\mathbb{R}$, [CHH21] also construct *smooth* compactly-supported blow-up for $a=1$; the
smooth-data regularity conjecture is specific to the periodic domain.) The quantitative numerical landscape of
this family was mapped by **Lushnikov, Silantyev & Siegel [LSS21]** (uncited in v0.5 — a prior-art omission we
correct here): with spectrally accurate, variable-precision methods on the periodic domain with sin-type smooth
data they report collapse for $a<a_c=0.6890665...$, self-similar blow-up for $a_c<a\le0.95$, **inconclusive
behaviour for $0.95<a<1$**, and global existence for $a\ge1$. The coefficient $a$ thus tunes the strength of
the depleting transport for fixed (smooth) data. We use it as a *mechanistic analogue* of the 3D story — a
secondary transport term competing with the stretching nonlinearity — while stressing the analogy is loose: the
1D model has no strain tensor and no eigenvector geometry, so it illustrates that *transport can deplete
blow-up*, not the specific geometric alignment mechanism of §2.

We integrate the family by a Fourier pseudospectral method ($2/3$ dealiasing, CFL-adaptive RK4, $N=4096$) from
smooth data $\omega_0=\sin x$, declaring blow-up when $\|\omega\|_\infty$ amplifies by $2\times10^4$. The solver
reproduces the known endpoints: $a=0$ blows up at $t\approx2.01$ (the exact CLM solution for this datum gives
$T^*=2$ analytically, so the $0.3\%$ overshoot calibrates the threshold criterion), while $a=1$ stays flat —
as it must, since $\omega_0=\sin x$ is an **exact steady state** of the De Gregorio model ($u=-\sin x$ makes
$u\omega_x-u_x\omega\equiv0$); the flatness validates the code rather than evidencing generic $a=1$
regularization. Scanning $a$ and fitting the blow-up time to $t_{\rm blow}(a)=t_0+C\,(a^*-a)^{-p}$ reveals a
**depletion threshold** at which the singularity time *diverges*:

| $a$ | 0.0 | 0.4 | 0.6 | 0.7 | 0.8 | 0.86 | 0.90 | 0.92 | 0.94 | 0.96 |
|---|---|---|---|---|---|---|---|---|---|---|
| $t_{\rm blow}$ | 2.01 | 2.79 | 3.63 | 4.39 | 5.87 | 7.81 | 10.50 | 14.22 | 29.85 | — (regular to $T{=}40$) |

The blow-up time **diverges steeply** as $a$ increases. The blow-up side of our scan (through $a=0.94$, where
the singularity forms at $t\approx29.9$) **reproduces the published boundary of [LSS21]** (blow-up for
$a\le0.95$ on the circle with sin-type data); we present it as a consistency check of our solver against their
higher-precision result, not as a new threshold. At $a=0.96$ we observe no blow-up under integration to $T=40$
(amplification $\times2.1$, $2\times10^5$ steps) — but $0.96$ lies inside the band $0.95<a<1$ that [LSS21]
explicitly declare *inconclusive* with far stronger methods, and a fixed-time survival test is exactly the
instrument our own analysis warns against (below), so **we classify $a=0.96$ as inconclusive, not regular**;
the true smooth-data threshold may sit at $a=1$, the original conjecture of [OSW08], consistent with the
nonlinear stability of the $a=1$ equilibrium proved by Jia–Stewart–Šverák [JSS19]. A power-law fit
$t_{\rm blow}(a)=C\,(a^*-a)^{-p}$ over the blow-up points returns $a^*\approx0.95$, $p\approx0.5$ — the
divergence-law fit appears to be new for the periodic case (we found no published blow-up-time exponent in
[LSS21]), with the caveats that $p$ moves between $0.46$ and $0.53$ depending on the offset convention, the
smooth-data scan has no grid-refinement study at $N=4096$, and the $a=0.94$ trajectory is non-monotone
(rise–dip–re-runaway), the signature of a marginal case; the fit's $a^*$ moves $0.947\to0.951$ without that
point. We stress a methodological point that is *itself* the heart of the open De Gregorio regularity question:
because $t_{\rm blow}$ diverges so steeply near $a^*$, a fixed-time "survival" test **systematically
misclassifies late-forming singularities as regular** — extending the integration window successively converted
$a=0.90$, $0.92$, and $0.94$ from apparently-regular to blow-up. The threshold must therefore be read from the
divergence law (or a growth-rate criterion), not from survival to a finite $T$ — which is precisely why we do
not read our own $a=0.96$ survival as regularity. Consistent with the rigorous result of Chen–Hou–Huang [CHH21]
that $C^\alpha$ data blows up for *all* $a$, a separate low-regularity probe (cusp data $\sim\!|\sin x|^\alpha$)
did exhibit blow-up at $a=1$, but its singular time collapsed toward zero under grid refinement — indicating
the raw point-cusp is too singular to provide a clean finite-time test, and that a faithful reproduction of
[CHH21] requires their specially-structured data rather than a generic cusp. What this section contributes is
therefore an *illustration with a quantitative anchor*: strengthening the depleting transport continuously
delays and then removes the finite-time blow-up, with the blow-up-time divergence law made explicit. The open
3D question (§5–6) is whether the geometric depletion survives into the extreme tail at high Reynolds number
(§6) — a statistical question the 1D model cannot answer.

## 8. Discussion

At $R_\lambda\approx433$ the *mean* depletion is firmly established and the tail point estimates indicate a
regularizing tendency (ratio $\approx0.77$–$0.81$ in the large sample, decreasing into the tail), though the
extreme top-$0.1\%$ ratio there is statistically marginal — the resampleable bootstrap straddles unity (§5).
At the higher Reynolds numbers, where the §6 high-statistics measurement provides $10^2$–$10^3\times$ more
points, the picture sharpens into a resolved result: at $R_\lambda\approx1300$ and $R_\lambda=2500$ the
extreme-tail ratio is **below unity with batch-means $t$-intervals that exclude 1** at top-$1\%$ and
top-$0.1\%$ — the geometric depletion of the most intense vorticity is real and persists to $R_\lambda=2500$,
now with its dominant systematic (grid resolution, §6.1) measured rather than argued. Comparing the two, the
tail depletion is directionally weaker at the higher Reynolds number — the mild de-depletion direction reported
by Buaria–Pumir — but the §6.1 experiment shows the resolution confound pushes the same way and accounts for
most of the top-$1\%$ difference, so **our two datasets do not independently resolve a Reynolds trend**. If
such erosion of the geometric protection with $R_\lambda$ is confirmed (as the better-resolved Buaria–Pumir
data indicate), it matters in the direction the supercriticality barrier identifies as decisive: the
protection would weaken in the events heuristically closest to a putative singularity — though we repeat that
forced steady-state extremes are not proto-singularities. We are careful not to over-read any of it — two
Reynolds numbers are a comparison, not a law; resolution co-varies with $R_\lambda$ and its bias is now known
to be non-conservative; and the deepest tail (top-$0.01\%$) is not gradient-resolved (§3.2, §6). And none of
this bears on blow-up directly: these are statistically steady forced flows, not a putative self-similar
collapse; the measurement quantifies a geometric *trend*, not a singularity.

## 9. Completed measurements and what remains

The high-statistics program of the prior revision has been executed (a standard JHTDB authorization token lifts
the 4096-point testing cap):

1. **High-statistics $R_\lambda\approx1300$ tail — done.** Per-snapshot ($R_\lambda\approx1300$ only) contiguous
   $128^3$ cutouts of `isotropic8192`, $7.7\times10^7$ interior points; tail resolved below unity at top-$1\%$
   and top-$0.1\%$ with block-resampled intervals (§6). As anticipated, the effective sample is far below the
   raw count (cutout points are spatially correlated), which is why we bound the tail by block (not point)
   resampling; the top-$0.01\%$–$0.001\%$ are *not* claimed.
2. **The exascale DNS — done.** Same diagnostic on `isotropic32768` ($R_\lambda=2500$, $k_{\max}\eta\approx2.8$),
   $1.6\times10^8$ points: a second, higher-Reynolds point on the $D_{\rm tail}/D_{\rm mean}$ curve (§6).
3. **Differentiation-order cross-check — done.** The full stencil ladder (`fd4`/`fd6`/`fd8`/`m2q8`) at
   $\sim\!5\times10^5$ identical grid points on both datasets confirms the *stencil* contribution is $<1\%$
   through top-$0.1\%$ ($D_{\rm mean}$ ratio $1.00$), with $5$–$9\%$ sign-flipping scatter at top-$0.01\%$
   (§3.2; note the rescoped interpretation there — stencil convergence does not certify grid resolution).
4. **The resolution-bias measurement — done (new in v0.6).** The paired coarse-graining experiment of §6.1
   measures the grid-resolution systematic directly: $-0.018\pm0.005$ on the top-$1\%$ ratio at
   `isotropic8192`-equivalent resolution, falsifying v0.5's conservative-direction argument and demoting the
   cross-Reynolds claim.

What remains, and would convert the two-point comparison into a trend:

5. **A third Reynolds number** between $433$ and $2500$ at matched resolution, to test monotonicity and begin a
   genuine $D_{\rm tail}/D_{\rm mean}(R_\lambda)$ curve, and ideally a resolution-matched pair to break the
   $R_\lambda$–$k_{\max}\eta$ degeneracy of §6 (the §6.1 bias measurement makes a resolution-matched
   comparison constructible by filtering the finer member of any pair).
6. **The deepest tail.** Resolving top-$0.01\%$–$0.001\%$ — where any *crossing* of 1 (true de-depletion) would
   live — needs both more effective-independent samples ($\mathcal{O}(10^8$–$10^9)$) and a finer gradient than
   on-grid 4th-order FD provides at that depth.
7. **Conditioning.** Resolve $D$ jointly on $|\omega|$ and $|S|$, and on the local strain-state ($R$–$Q$
   invariants), to identify *which* extreme events de-deplete.

## 10. Limitations

(i) Forced statistically-steady isotropic turbulence is a proxy for the regularity question, not a blow-up
scenario. (ii) The high-Reynolds results resolve the top-$1\%$ and top-$0.1\%$ but **not** the
top-$0.01\%$–$0.001\%$, where any true *crossing* of unity would live: the gradient stencils scatter $5$–$9\%$
there (§3.2) and the effective independent sample is insufficient (§6, §9). (iii) The Reynolds comparison rests
on **two** datasets that differ in small-scale resolution as well as $R_\lambda$ ($k_{\max}\eta\approx1.4$ vs
$2.8$); the §6.1 experiment measures the resolution bias and finds it runs *with* the observed direction and
comparable in size, so we claim a directionally-consistent two-point difference whose Reynolds attribution is
open — not a resolved trend, and not a scaling law. (iv) The batch-means intervals use $K=20$ contiguous
blocks; the block SEM is not $K$-independent (it *grows* mildly with finer partitions, reflecting cutout-scale
heterogeneity — a larger $K$ does not simply tighten it, correcting a v0.5 misstatement), the `isotropic32768`
blocks sub-sample a single snapshot (an ergodic assumption stands in for realization replication), and true
independent replication is $\approx\!5$ snapshots at $R_\lambda\approx1300$ and $1$ at $R_\lambda=2500$ — the
below-unity margins survive even a doubled SEM, but marginal inferences on these intervals do not. (v) JHTDB
and in-house Taylor–Green differ in forcing/decay; absolute $D$ values are not directly comparable across
them, so we rely on the *bulk-vs-tail ratio* as the within-dataset discriminator. (vi) The §6.1 bias
measurement is an a-priori filtering test at $R_\lambda=2500$, transferred to the $R_\lambda\approx1300$ point;
it is not a substitute for a genuinely resolution-matched dataset pair.

## 11. Conclusion

We have built and validated an apparatus that measures the percentile-resolved geometric depletion of enstrophy
production across Reynolds number, reproduced the canonical alignment, and firmly established the mean depletion
($D\approx0.12$). At $R_\lambda\approx433$ the tail point estimates indicate a regularizing tendency (ratio
$\to0.77$–$0.81$), though the resampleable bootstrap straddles unity at the top-$0.1\%$. At
$R_\lambda\approx1300$ and $R_\lambda=2500$, a high-statistics, per-snapshot, single-$R_\lambda$ measurement
($7.7\times10^7$ and $1.6\times10^8$ interior points) **resolves the extreme tail below unity** at top-$1\%$ and
top-$0.1\%$ — batch-means $t$-intervals exclude 1 by $3.8$–$10.7$ standard errors — establishing that the most
intense vorticity is more geometrically depleted than the bulk, and that this regularizing depletion persists
to $R_\lambda=2500$, now with the grid-resolution systematic measured (§6.1) rather than argued. Across the two
Reynolds numbers the tail depletion is directionally weaker at higher $R_\lambda$ (the Buaria–Pumir direction;
two-sample $p\approx0.003$ at top-$1\%$ before systematics), but the measured resolution bias runs the same way
and absorbs the majority of the difference — so we report a two-point comparison whose Reynolds attribution is
open, not a resolved trend, and certainly not a scaling law; the deepest tail (top-$0.01\%$), where a true
crossing of unity would live, is not gradient-resolved. As an independent, *decidable* illustration of the
depletion *mechanism*, §7 exhibits the transport-depletes-stretching competition in the one-dimensional
Okamoto–Sakajo–Wunsch model family, reproducing the blow-up/regularity landscape of [LSS21] and adding an
explicit blow-up-time divergence law. The contributions of this work are thus the validated apparatus, the
resolved extreme-tail depletion at $R_\lambda\le2500$ with its measured resolution systematic, the honestly
demoted Reynolds comparison, the $R_\lambda\approx433$ baseline, the 1D-model illustration, and a reproducible,
self-audited account — including its own published errata — that holds each claim to the depth the data
support.

## 12. Methodology note: an AI-augmented, consumer-hardware research apparatus

This work was produced by a small human-led, agentic-AI research process, and we report the methodology in the
same spirit of disclosure as the results, because *how* it was done is part of what it demonstrates.

**The apparatus.** A human investigator (C.R.) set direction, posed questions, and served as the accountable
validation layer; an agentic large-language-model system did the bulk of the literature mapping, code
authoring, numerical experimentation, and drafting, under continuous human review. Two ingredients are worth
naming. First, a **Topic-Keyed Hash Routing (TKHR)** memory substrate — a tamper-evident, hash-chained index
that mirrors a soul/agent-memory pattern — was used to organize the relevant literature and the catalogue of
"what does not work" for the Navier–Stokes regularity problem, and to surface candidate *cross-subdomain
connections* (e.g. linking the 1D model-problem mechanism of §7 to the 3D geometric-depletion picture). Second,
a capacity-aware LLM router fanned ideation and an adversarial literature-scrutiny pass across a mix of
local and free/flat-rate cloud models. The proposed experiments were then implemented and run as ordinary,
inspectable code.

**Consumer hardware.** The entire effort ran on commodity equipment: a single workstation (one RTX 5090 and one
RTX 5070), the public Johns Hopkins Turbulence Database accessed through a free SciServer Compute container, and
free or flat-rate LLM tiers. The in-house pseudospectral DNS, the JHTDB cutout pipeline, the 1D model-problem
solver (§7), and the literature scrutiny were all feasible at this scale. What was *not* feasible at this scale
is equally informative: the definitive extreme-tail measurement (§6, §9) exceeded the statistics our access tier
could provide, and we say so rather than paper over it.

**The integrity layer is the point.** During this study the AI system produced, at three separate stages,
confident and internally coherent conclusions about the high-$R_\lambda$ extreme tail that did not survive
scrutiny: (1) an apparent Reynolds trend that was an artifact of pooling mixed-$R_\lambda$ snapshots; (2) a
"converged" tail that was rare-event under-sampling at $\sim\!5\times10^5$ points; and (3), in this revision, a
small-sample point-gradient ("gate") comparison that put the cross-$R_\lambda$ ordering the *wrong way*, reversed
only when the $10^8$-point cutout measurement and a block-resampled error bar were brought in. A fourth catch was
external-source-shaped: a web search confidently asserted the two datasets were at *comparable* Reynolds number,
which would have voided the entire comparison — the primary simulation paper [Yeung32768] shows
`isotropic32768` is at $R_\lambda=2500$, twice `isotropic8192`. Each was caught — by a human-prompted provenance
audit, an error-bar discipline that respects sample correlation, and verification against primary sources rather
than an LLM's paraphrase — before it could stand. **This revision adds three more catches, found
post-publication by re-running the audit with a newer-generation frontier model placed in an adversarial
verifier configuration** (independent agents briefed to refute each claim and to verify every reference against
its primary record): (5) the v0.5 reference for the exascale source simulation carried a confabulated author
list ([Yeung32768], corrected here); (6) three prior literature-scrutiny passes missed directly relevant prior
art on the §7 threshold ([LSS21], credited here); and (7) v0.5's asserted "conservative" resolution-bias
direction — plausible, internally coherent, and load-bearing for the cross-Reynolds claim — was *reversed* by
the direct measurement of §6.1 within hours of being challenged (along with a normal-vs-Student-$t$ interval
slip that had made a marginal comparison look resolved). The pattern is instructive in both directions:
verification depth scaled with the capability of the model doing the verifying, and yet the decisive
instrument was not the model but the *measurement* it demanded. We present this as the central methodological
result: an AI research collaborator moves fast and finds real structure, but *validation* — deciding what is
true enough to claim, and to what depth — must sit with an accountable human and an apparatus that distrusts
its own first answer, and that treats its own published claims as standing targets. The final status of this
paper (a firmly-established mean depletion; a tail depletion *resolved* below unity at $R_\lambda\le2500$ with
a measured resolution systematic; a Reynolds comparison honestly demoted to an open attribution; an unresolved
deepest tail; and published errata) is itself the product of that layer working.

**On time-to-findings, and on disciplinary reach.** This work was carried out over a span of days, by a human
investigator with a background in engineering and security — *not* in turbulence or fluid dynamics. From a
standing start the apparatus reproduced canonical vortex-stretching results, stood up a working GPU DNS,
ingested and analysed state-of-the-art database turbulence, located the sharp threshold of a 1D model problem,
and produced this manuscript with its literature scrutiny and self-audit — on commodity hardware. The point is
not that the result is large (it is deliberately modest and hedged); it is that an agentic AI research partner
can compress the *time to a checkable, honestly-bounded research artifact*, and can lower the disciplinary
barrier to entering an unfamiliar field, provided an accountable human holds the validation layer. This mirrors
the broader 2024–2026 evidence that agentic systems accelerate the *workflow* of science — literature synthesis,
hypothesis generation, code authoring, and closed-loop analysis — with reported compressions from months to days
and hours to minutes in some pipelines, and with experimentally validated outputs in others (e.g. multi-agent
"AI co-scientist" systems) [GoogleAICoSci25]. Our contribution to that conversation is methodological and
cautionary in equal measure: the *speed* is real, and so is the failure mode (confident, well-formed, wrong
intermediate conclusions); the two are reconciled only by a disciplined human-in-the-loop validation layer of
exactly the kind we document here. We disclose AI use explicitly, in keeping with emerging journal norms on
AI-assisted research and authorship [NatureAI23].

## Author contributions and AI transparency statement

C.R. (Cristian Ruvalcaba) conceived and directed the study, made all substantive judgments, validated every
claim, and is solely accountable for the contents. The **Agentic AI Research Team** (a large-language-model
system operated by Saluca Labs) performed literature search and synthesis, authored the numerical codes
(`dns_*`, `osw_*`, `sciserver_cutout.py`, analysis and bootstrap scripts), executed the computational
experiments, drafted the manuscript, and ran the adversarial literature scrutiny — all under human direction and
review. Listing the AI system as a co-contributor reflects the reality of how the work was done; it does not
transfer scientific responsibility, which remains with the human author. No claim in this paper should be
treated as established solely because an AI produced it; the human-validated, reproducible artifacts (code,
data, and this disclosure) are the basis for trust.

## Data and code availability

DNS solver, JHTDB pipeline (`sciserver_cutout.py`), the coarse-graining experiment of §6.1
(`coarse_grain_run.py`), the 1D OSW model solver (`osw_1d.py`, `osw_frontier.py`), analysis and batch-means
scripts, and figure code are archived with this work and released under the MIT license
at **https://github.com/saluca-labs/navier-stokes**. JHTDB datasets are public
[JHTDB08]; the geometric-depletion statistic, the per-snapshot audit of §6, and the §6.1 bias measurement are
reproducible from them.

## Changes in v0.6 (2026-07-01)

This revision follows an independent adversarial verification pass (ten agents: five verifying every reference
against its primary record, five briefed to refute each derived claim, with code audit) and one new
measurement. **Corrections:** the cross-Reynolds "statistically resolved at top-$1\%$" claim is withdrawn
(normal-vs-$t_{19}$ interval slip; proper two-sample test reported instead; §6); the "conservative resolution
bias" argument is falsified by direct measurement and replaced by the §6.1 experiment; §3.2 is rescoped to
stencil-order convergence only; §7 is reframed from "corroboration" to "illustration", credits the previously
uncited [LSS21] whose periodic threshold our scan reproduces, and downgrades the $a=0.96$ classification to
inconclusive; the [Yeung32768] author list and the [CHH21] DOI are corrected; the Frobenius-vs-operator-norm
factor is corrected to the bounded range $[\sqrt3,2]$; the §6 sample composition is pinned exactly.
**Additions:** the §6.1 paired coarse-graining experiment (the measured resolution bias, $-0.018\pm0.005$ on
the top-$1\%$ ratio); $t$-corrected intervals throughout; estimand and single-snapshot disclosures; three new
entries in the §12 integrity ledger. The central below-unity result is unchanged and strengthened by the
measured systematic.

## References

- [Leray34] J. Leray, *Sur le mouvement d'un liquide visqueux emplissant l'espace*, Acta Math. **63** (1934) 193–248. doi:10.1007/BF02547354
- [BKM84] J.T. Beale, T. Kato, A. Majda, *Remarks on the breakdown of smooth solutions for the 3-D Euler equations*, Comm. Math. Phys. **94** (1984) 61–66. doi:10.1007/BF01212349
- [CF93] P. Constantin, C. Fefferman, *Direction of vorticity and the problem of global regularity for the Navier–Stokes equations*, Indiana Univ. Math. J. **42** (1993) 775–789.
- [AKKG87] W.T. Ashurst, A.R. Kerstein, R.M. Kerr, C.H. Gibson, *Alignment of vorticity and scalar gradient with strain rate in simulated Navier–Stokes turbulence*, Phys. Fluids **30** (1987) 2343–2353. doi:10.1063/1.866513
- [JWSR93] J. Jiménez, A.A. Wray, P.G. Saffman, R.S. Rogallo, *The structure of intense vorticity in isotropic turbulence*, J. Fluid Mech. **255** (1993) 65–90. doi:10.1017/S0022112093002393
- [JHTDB08] Y. Li, E. Perlman, M. Wan, Y. Yang, C. Meneveau, R. Burns, S. Chen, A. Szalay, G. Eyink, *A public turbulence database cluster and applications to study Lagrangian evolution of velocity increments in turbulence*, J. Turbulence **9** (2008) N31. doi:10.1080/14685240802376389. arXiv:0804.1703
- [JHTDB-i8192] Johns Hopkins Turbulence Database, *Forced isotropic turbulence dataset (`isotropic8192`)*, dataset documentation, turbulence.pha.jhu.edu/Isotropic8192.aspx (accessed 2026-06-27): $R_\lambda\approx1200$–$1300$ for snapshots 0–4, $\approx610$ for snapshot 5.
- [JHTDB-i32768] Johns Hopkins Turbulence Database, *Isotropic turbulence dataset (`isotropic32768`)* — a single snapshot of the $32768^3$ exascale forced-isotropic DNS, $R_\lambda=2500$, $k_{\max}\eta\approx2.8$; dataset documentation, turbulence.idies.jhu.edu/datasets/homogeneousTurbulence/isotropic32768 (accessed 2026-06-30).
- [Yeung32768] P.K. Yeung, K. Ravikumar, R. Uma-Vaideswaran, D.L. Dotson, K.R. Sreenivasan, S.B. Pope, C. Meneveau, S. Nichols, *Small-scale properties from exascale computations of turbulence on a $32768^3$ periodic cube*, J. Fluid Mech. **1019** (2025) R2. doi:10.1017/jfm.2025.10493. Forced isotropic turbulence at $R_\lambda=2500$ with dissipation-range resolution $k_{\max}\eta\approx2.8$; the source simulation for the JHTDB `isotropic32768` snapshot. *(v0.5's entry for this reference carried an incorrect author list and year; corrected here.)*
- [BPBY19] D. Buaria, A. Pumir, E. Bodenschatz, P.K. Yeung, *Extreme velocity gradients in turbulent flows*, New J. Phys. **21** (2019) 043004. doi:10.1088/1367-2630/ab0756. arXiv:1901.09989
- [BP22] D. Buaria, A. Pumir, *Vorticity-strain rate dynamics and the smallest scales of turbulence*, Phys. Rev. Lett. **128** (2022) 094501. doi:10.1103/PhysRevLett.128.094501. arXiv:2202.01364
- [BBP20] D. Buaria, E. Bodenschatz, A. Pumir, *Vortex stretching and enstrophy production in high Reynolds number turbulence*, Phys. Rev. Fluids **5** (2020) 104602. doi:10.1103/PhysRevFluids.5.104602. arXiv:2006.01312
- [Tao16] T. Tao, *Finite time blowup for an averaged three-dimensional Navier–Stokes equation*, J. Amer. Math. Soc. **29** (2016) 601–674. arXiv:1402.0290
- [CLM85] P. Constantin, P.D. Lax, A. Majda, *A simple one-dimensional model for the three-dimensional vorticity equation*, Comm. Pure Appl. Math. **38** (1985) 715–724. doi:10.1002/cpa.3160380605
- [DG90] S. De Gregorio, *On a one-dimensional model for the three-dimensional vorticity equation*, J. Stat. Phys. **59** (1990) 1251–1263. doi:10.1007/BF01334750
- [OSW08] H. Okamoto, T. Sakajo, M. Wunsch, *On a generalization of the Constantin–Lax–Majda equation*, Nonlinearity **21** (2008) 2447–2461. doi:10.1088/0951-7715/21/10/013
- [CHH21] J. Chen, T.Y. Hou, D. Huang, *On the finite time blowup of the De Gregorio model for the 3D Euler equations*, Comm. Pure Appl. Math. **74** (2021) 1282–1350. doi:10.1002/cpa.21991. arXiv:1905.06387. *(v0.5 listed an incorrect DOI; corrected here.)*
- [LSS21] P.M. Lushnikov, D.A. Silantyev, M. Siegel, *Collapse versus blow-up and global existence in the generalized Constantin–Lax–Majda equation*, J. Nonlinear Sci. **31** (2021) 82. doi:10.1007/s00332-021-09737-x. arXiv:2010.01201
- [JSS19] H. Jia, S. Stewart, V. Šverák, *On the De Gregorio modification of the Constantin–Lax–Majda model*, Arch. Ration. Mech. Anal. **231** (2019) 1269–1304. doi:10.1007/s00205-018-1298-1
- [EJ20] T.M. Elgindi, I.-J. Jeong, *On the effects of advection and vortex stretching*, Arch. Ration. Mech. Anal. **235** (2020) 1763–1817. doi:10.1007/s00205-019-01455-9. arXiv:1701.04050
- [GoogleAICoSci25] J. Gottweis et al., *Towards an AI co-scientist* (multi-agent, Gemini-based system for hypothesis generation; reports experimentally validated hypotheses, e.g. in-vitro-validated drug-repurposing candidates), 2025. arXiv:2502.18864.
- [NatureAI23] Nature Editorial, *Tools such as ChatGPT threaten transparent science; here are our ground rules for their use*, Nature **613** (2023) 612. doi:10.1038/d41586-023-00191-1.
