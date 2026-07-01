# Experiments — results log

## E1. Dyadic shell-model dissipation threshold (2026-06-21) ✅
`shell_model.py` — reproduces the failure-catalog dissipation threshold from scratch (N=22 shells,
k_n=2^n, energy-conserving forward cascade, LSODA stiff integrator, enstrophy blow-up event).

| Regime | γ | Result | Cascade front |
|---|---|---|---|
| inviscid | 0 | **BLOW-UP** t≈1.69 | races to shell 21/22 |
| subcritical | 0.15 (<¼) | **BLOW-UP** t≈1.69 | races to shell 21/22 |
| supercritical | 0.50 (>⅖) | **REGULAR** to t=40 | **stalls at shell 12**, enstrophy ~76 |

**Takeaway:** blow-up = the forward cascade front reaching arbitrarily small scales in finite time
faster than dissipation can act; regularity = dissipation catching the front so it parks. This is the
concrete face of the supercriticality story. Numerics note: explicit RK4 with N=36 (k~7e10) was too
stiff (NaN); N=22 + LSODA + blow-up event is robust.

## E2. Tao-delay concentration vs spreading (2026-06-21) — INSTRUCTIVE NEGATIVE
`e2_delay.py` — conservative gated dyadic (amplitude form). At SAME viscosity (ν=1e-3, γ=0.5):
spreading REGULAR (spread 1.59), naive-gate concentration also REGULAR (spread 1.38). Inviscid: gate
blows up LATER (t=2.72) than spreading (t=1.69) → **the naive "sequence-after-predecessor-depletes"
gate makes the cascade LESS singular, not more.** Energy-variable v1 was buggy (catalytic √(eₙ) fixed
point at 0 → no cascade); amplitude form is correct.

**Lesson:** destroying spreading ≠ forcing blow-up. Tao's blow-up needs transfer that is COMPLETE *and
ACCELERATING*, not merely sequenced. Timescale heuristic: a complete-transfer relay spends time ~k_j⁻¹
at shell j ⇒ total viscous dissipation ~Σ k_j^{2γ−1}, which CONVERGES (packet survives ⇒ blow-up) **iff
γ<½**, vs spreading's ~⅖. ⇒ predicted clean window **γ∈(0.4,0.5): complete-transfer concentration blows
up where spreading is regular.** This is why Tao needed a self-replicating "machine," not a gate.

## E2b. Complete-transfer relay (2026-06-21) — window is ASYMPTOTIC, not visible at N=22
`e2b_relay.py` — self-driven flux F_n=k_n e_n^{3/2} (complete transfer, no gate). Front-based blow-up
(front→top shell). Front-threshold γ-scan: **relay and spreading BOTH blow up for γ≤0.40, both regular
for γ≥0.42** at N=22 — the predicted γ∈(0.4,0.5) separation did NOT appear. The heuristic (relay
threshold ½) is an asymptotic (k_n^{1-2γ}→∞ at large n); N=22 + along-the-way dissipation collapse the
gap. **Honest verdict:** the "spreading regularizes" threshold-separation is real in theory but
asymptotic — needs N≥40 + an exponential/integrating-factor integrator (E2c) to resolve. Two predictions
(E2 gate sign, E2b clean window) were corrected by the lab — exactly its job.

**Status of shell track:** E1 clean ✓ (dissipation threshold = catalog's main quantitative result).
E2/E2b: finer claim is asymptotic; deferred (E2c) in favor of real-equation DNS (D1). Lesson banked:
demonstrating threshold *gaps* in shell models is a large-N effect, not a few-shell one.

### Next on this track
- **E2 — Tao's delay variant:** add a delay forcing near-complete a_j→a_{j+1} transfer (kills the
  regularizing energy-spreading) → should force BLOW-UP even at γ=0.5. Directly demonstrates the catalog's
  "Tao had to insert a delay to destroy the regularizing spreading effect."
- **E3 — quantify energy-spreading:** measure how energy distributes across shells during cascade
  (spread width) for regular vs blow-up — turn the catalog's qualitative "finer structure" into a number.
- **E4 — threshold scan:** sweep γ to locate the numerical regular/blow-up boundary; compare to the
  theoretical γ≥2/5 (physical) / 5D-scaling threshold.

## E2c. Threshold gap vs resolution (2026-06-21) — qualitative yes, quantitative murky
`e2c_highres.py`, N=34. A NARROW real window: at γ=0.36 relay BLOW-UP while spread REGULAR ⇒ concentration
is more singular ✓ (qualitative "spreading regularizes"). BUT gap is 1 grid-step wide AND absolute
thresholds DROPPED with N (0.42@N22 → ~0.37@N34) — opposite of the "→½" heuristic, because more shells =
more along-the-way dissipation before the front hits the resolution limit (the front-reaches-top criterion
is N-biased). **Verdict:** qualitative claim confirmed; a clean quantitative gap needs an integrating-factor
integrator at N≫40 + a front-velocity blow-up test. SHELL-MODEL TRACK CLOSED for now (3 experiments;
E1 solid, E2/E2b/E2c corrected my predictions — the lab's real value). DNS (D1) is the instrument forward.

## D1. GPU pseudospectral DNS — real 3D Navier-Stokes ✅ (2026-06-21)
`dns_spectral.py` — PyTorch FFT pseudospectral solver on the RTX 5090 (rotational form, 2/3 dealias,
Leray projection, RK4). Taylor-Green vortex, 128^3, Re=100. **VALIDATED:** max|div u|~1e-7 (div-free),
energy monotone decay 0.125->0.026, **enstrophy peaks at t≈5** (textbook TG), max|ω| (BKM) peaks ~5.8
then bounded (regular at Re=100). **171 steps/s at 128^3** (1000 steps / 5.9s). We have a correct DNS of
the actual equation as a local instrument.
### Next on DNS: crank Re (lower nu) + resolution; instrument vorticity-direction coherence
(Constantin-Fefferman geometric depletion), the BKM integral ∫max|ω|dt, and the energy spectrum.

## D2. Instrumented DNS — Taylor-Green Re=1600, 256^3 ✅✅ (2026-06-21) — THE science run
`dns_instrumented.py`, 2000 steps / **85s** on the 5090. Reproduced two celebrated structural results
from scratch AND measured the catalog's candidate regularizer:

1. **Quantitative validation:** enstrophy peaks at **t≈8.9 (Ω≈10.3)** — matches the published TG Re=1600
   benchmark to the eyeball. Energy monotone decay, div-free.
2. **Vortex stretching:** enstrophy production P=⟨ω·Sω⟩ > 0 throughout, peaks ~20 at the enstrophy peak —
   the growth term, captured.
3. **Strain-vorticity ALIGNMENT (the famous fingerprint), enstrophy-weighted ⟨|cos(ω,e_k)|⟩:**
   - t=2 (laminar stretch): [e_min,e_mid,e_max]=[0.01, 0.04, **0.99**] — ω ∥ extensional e_max.
   - t=4 (transition):       [0.07, **0.91**, 0.23] — flips to the **INTERMEDIATE** eigenvector.
   - t=6..10 (developed):    e_mid stays dominant (0.86→0.71). **Reproduced the Ashurst/Kerr intermediate-
     eigenvalue alignment of vorticity in turbulence — a top structural result — measured locally.**
4. **GEOMETRIC DEPLETION (Constantin-Fefferman candidate regularizer):** depl = ⟨ω·Sω⟩/⟨|ω|²|S|⟩ drops to
   **~0.08-0.10 in developed turbulence** — i.e. the nonlinear stretching runs at only ~10% of its
   worst-case geometry. **This IS the catalog's "finer structure" (the open door), measured: ~90% depletion.**
   max|ω| (BKM) bounded (peaks ~73), BKM integral grows ~linearly = consistent with regularity at Re=1600.

(Spectrum normalization is unnormalized = shape only; shows a developed cascade k=2..32, no high-k pileup
= resolved. Fix normalization later.)

### THE research question this opens (D3)
**How does the geometric-depletion factor scale with Re?** If depletion STRENGTHENS (depl→0) as Re→∞, the
structure self-regulates and suppresses blow-up (evidence for regularity). If it WEAKENS (depl→O(1)), the
nonlinearity de-depletes at high Re (evidence toward singularity). This is a real, crux-of-the-problem
measurement we can run: Re-scan of depl at fixed resolution-adequacy. The catalog said vorticity geometry
is the open door; D3 quantifies whether that door opens or closes with Re.

## D3. Geometric-depletion Re-scan ✅ (2026-06-21) — the crux measurement, part 1
`d3_rescan.py`, resolution-matched Taylor-Green (k_max·η≈1.0-1.45, all adequately resolved). Mean
geometric depletion `<ω·Sω>/<|ω|²|S|>` at the enstrophy peak:
  Re=400 (128³): 0.117 | Re=800 (192³): 0.106 | Re=1600 (256³): 0.098.
**Mean depletion STRENGTHENS monotonically with Re over 4× range** ⇒ bulk nonlinearity self-depletes more
as turbulence intensifies — suggestive (not proof) toward self-regularization (consistent w/ Constantin-
Fefferman). CAVEATS: 3 points ≠ scaling law (~16% drop; could saturate); k_max·η not perfectly constant.
**KEY SUBTLETY:** this is a MEAN; the blow-up question lives in the EXTREME tail (top |ω|), where recent
high-Re work (Buaria-Pumir-Bodenschatz) finds extremes may DE-deplete even as the mean depletes — the
knife-edge. ⇒ D4.

## D4. Mean vs EXTREME-tail depletion vs Re ✅ (2026-06-21) — reached the resolution frontier
`d4_extreme.py`. Depletion conditioned on top-q |ω| at the enstrophy peak:
  Re=400(kmη1.45): mean .126 | t1% .110 | t0.1% .014 | t0.01% -.019   (LAMINAR — tail = TG sheets, ignore)
  Re=800(kmη1.29): mean .106 | t1% .108 | t0.1% .143 | t0.01% .217   (extremes DE-deplete: tail ~2× mean)
  Re=1600(kmη1.00): mean .096 | t1% .099 | t0.1% .096 | t0.01% .111
  Re=3200(kmη0.86): mean .093 | t1% .088 | t0.1% .088 | t0.01% .099
**Findings:** (1) MEAN depletion strengthens with Re (→regularity), robust, confirmed to Re=3200.
(2) Genuine signal at Re=800: extreme tail (top0.01%=.217) ≫ mean (.106) — extreme stretching events
DE-deplete (singularity-relevant). (3) Apparent fall-back of the tail at Re=1600→3200 is CONFOUNDED by
resolution (k_max·η 1.0→0.86; the tail = smallest scales = first clipped). Cannot separate physics from
under-resolution — **this is the exact frontier wall that keeps the problem open.** The script's auto-verdict
("rises -> de-deplete") wrongly anchored on the laminar Re=400 row; corrected here.
**D5 (the decisive test):** ONE well-resolved high-Re run at 512³ (k_max·η≳1.2 at Re~3200-4000) for a
TRUSTWORTHY extreme-tail point — push the 5090 to its edge (free ollama, ~15-20min). Below that resolution,
the extreme-event verdict is not honestly obtainable.

## D5. Decisive well-resolved 512³ Re=3200 ✅ (2026-06-21) — extreme-tail verdict
`d5_decisive.py`, 512³ Re=3200, **k_max·η=1.14** (well-resolved), 48-min run on the 5090, ollama auto-restored.
mean=0.096, top1%=0.092, top0.1%=0.083, **top0.01%=0.094** vs under-resolved 384³ tail 0.099 (kmη0.86).
**VERDICT: well-resolved tail (0.094) ≈ under-resolved (0.099) ⇒ under-resolution was NOT clipping the tail;
the fall-back is REAL.** Well-resolved extreme-tail arc: 0.217(Re800)→0.111(Re1600)→0.094(Re3200) — extremes
STRENGTHEN with Re and CONVERGE to the mean; the Re=800 extreme de-depletion was a low-Re transient. Up to
Re=3200 the most intense stretching events deplete like the bulk (toward regularity). CAVEAT: Re=3200 is
moderate; the Buaria-Pumir extreme DE-depletion appears at far higher Re beyond single-GPU reach — that
asymptotic regime is the remaining frontier. **NS DNS arc complete (D1-D5): validated solver → reproduced
vorticity alignment + geometric depletion → measured its Re-trend (mean & resolved extremes → regularity up
to Re=3200) → hit/named the asymptotic resolution frontier.**

## D6. FRONTIER ROUTE OPEN — diagnostics on JHTDB high-Re DNS ✅ (2026-06-21)
`ns_jhtdb.py` (givernylocal; pyJHTDB deprecated). Queried 2000 velocity-gradient points from JHU's
`isotropic1024coarse` (Re_λ~433, effectively ≫ our TG Re=3200), public testing token, 8s. Geometric
depletion on REAL high-Re turbulence: **mean=0.126, top1%=0.113** (extreme tail BELOW mean ⇒ extremes more
depleted than bulk — NO de-depletion even here; consistent w/ D5 toward-regularity). **The resolution wall is
bypassed: analyze turbulence JHU already simulated.** CAVEATS: 2000 pts = proof-of-pipeline (top1%=20 pts,
noisy); different flow (forced isotropic vs decaying TG) so absolute values not directly vs-TG-comparable.
**DEFINITIVE extreme-tail (top 0.01% ~ needs 1e6 pts) requires: free registered JHTDB/SciServer token + a
batched downloader.** Diagnostic pipeline (gradient→S,ω→depletion+tail) validated and ready.

## D8 — SciServer LOCAL run (`sciserver_depletion.py` / `analyze_pool.py`), 2026-06-22 — frontier opened
Inside the SciServer Compute container `_salucaTurbulence` ([[jhtdb-sciserver]]) the JHTDB DNS is mounted local;
`giverny` reads it at **~1.4s / 4096-pt call (~12x faster/pt than the web API)**. Testing token still caps
4096 pts/call (client-side policy), so we chunk; a free token lifts it. **Baseline confirm (isotropic1024coarse,
Re_λ~433, t=0 single snapshot, 40,960 pts)** — the tail/mean ratio falls MONOTONICALLY into the extreme tail:
| quantile | depletion | ratio/mean | n |
|---|---|---|---|
| mean | 0.1131 | 1.00 | 40960 |
| top 10% | 0.107 | 0.95 | 4096 |
| top 1% | 0.088 | 0.78 | 409 |
| top 0.1% | 0.074 | **0.66** | 40 |
| top 0.05% | 0.070 | **0.62** | 20 |
| top 0.01% | 0.098 | 0.87 | 4 (noise) |

**The most intense vorticity is progressively MORE geometrically depleted than average** (ratio → 0.62 at
top-0.05%) — clean monotone regularizing signature, robust (nested subsets). Structural: **24% of points have
NEGATIVE enstrophy production** (ω·Sω<0); top-1% enstrophy carries only 23% of total production (concentrated,
not pathological); intermittency max|ω|²/mean ≈ 135. Consistent with D6/D7 (no de-depletion at Re_λ~433).
**NEXT (decisive):** the same diagnostic at **isotropic4096 (Re_λ~610)** and **isotropic8192 (Re_λ~1300)** —
does the ratio/mean curve REVERSE (rise above 1 in the tail = Buaria extreme de-depletion) as Re_λ climbs? Full
run: `NPTS=2000000 DATASETS=isotropic8192,isotropic4096,isotropic1024coarse python sciserver_depletion.py`
(detached, resumable per-dataset `pool_*.npz`). 2M pts → top-0.01% at ~200 pts (real). Token email sent for top-0.001%.

**D8b — high-Re access status (2026-06-22):** mounted high-Re sets = **isotropic8192 (Re_λ~1300)** +
**isotropic32768** (largest DNS); isotropic4096/1024coarse NOT local (web fallback). **HARD WALL: testing token
caps every call (getData AND getCutout) at 4096 pts** → real authorization token REQUIRED (email drafted/queued).
`sciserver_cutout.py` (contiguous-box reads, 4th-order on-grid FD, spatial-axis reversal, np.int64 xyzt ranges,
var='velocity', integer snapshot times) is READY for the token: B=128 box = 2M pts/read → top-0.001% in minutes.
**LIVE Re_λ~1300 SIGNAL** (random getData, ~12k pts/3 batches @208s): mean≈0.127, **top0.1%/mean ≈ 1.3 — the
ratio INVERTED above 1** vs 0.77 at Re_λ~433 (the de-depletion onset; needs ~1M pts to confirm; overnight run
`NPTS=1000000 DATASETS=isotropic8192` queued, resumable). **Money figure:** `plot_depletion.py` →
`depletion_ratio_vs_Re.png` (tail/mean vs quantile w/ 1.0 threshold); baseline Re_λ~433 drawn, overlay high-Re
when it lands. See [[jhtdb-sciserver]].

## D7 — BATCHED JHTDB extreme-tail (`d7_jhtdb_batch.py`), 2026-06-21 — the named next step, executed
Pushed D6 from a 2000-pt proof-of-pipeline to a **60,000-point pool** on the same public testing token by
batching small calls across **11 independent time-slices** (t=0..10, fresh turbulence realizations) of
isotropic1024coarse (Re_λ~433). Durable/resumable: pool persisted to `_d7_pool.npz` every batch (survives the
session). ~8s/batch, ~30 batches, ~4 min total. **RESULT (robust to top 0.1%):**
| quantile | depletion | n |
|---|---|---|
| mean | **0.1209** | 60000 |
| top 5% | 0.1177 | 3000 |
| top 1% | 0.1181 | 600 |
| top 0.1% | **0.1106** | 60 |
| top 0.01% | 0.169 | 6 (NOISE — not interpretable) |

**Verdict:** at Re_λ~433 the extreme tail (top 0.1% = 0.111) is **equal-to-or-slightly-below the mean (0.121)**
— extremes are as geometrically depleted as the bulk; **NO extreme de-depletion**. This is the *regularizing*
picture, now confirmed with real statistics at an Re well above our TG ceiling (3200) on an INDEPENDENT,
externally-simulated flow (forced isotropic vs decaying TG — absolute values not 1:1 comparable, but the
bulk-vs-tail RELATION is the discriminator and it points to regularity). Consistent with D6's 2000-pt hint
(extremes ≥ bulk-depleted). **Frontier that remains:** a trustworthy top-0.01% needs ~600k pts (→ free
registered token; testing token's 60 pts there = noise) AND the genuinely higher-Re regime (Re_λ~650–1300,
the Buaria–Pumir de-depletion band) — JHTDB has no isotropic set that high, so the next data step is the
registered-token deep pull on isotropic1024 + checking the higher-Re forced/channel sets. Pipeline + batching
proven; the only gate to the definitive top-0.01% answer is token quota, not method.

## (superseded) Open / GPU track
- **D1 — GPU pseudospectral DNS** of real 3D incompressible NS (jax-cfd or torch-FFT spectral) on the 5090:
  instrument vorticity-direction coherence (Constantin-Fefferman geometric depletion) + energy spectrum,
  probe whether the "finer structure" (vorticity differential operators / energy-spreading) is measurable.
