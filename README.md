# Geometric depletion of enstrophy production in the extreme vorticity tail

**A Reynolds-number study toward the Navier–Stokes regularity question.**

Reproducibility package for the paper:

> Cristian Ruvalcaba and the Agentic AI Research Team · Saluca Labs.
> *Geometric depletion of enstrophy production in the extreme vorticity tail:
> a Reynolds-number study toward the Navier–Stokes regularity question.* (v0.5, 2026)
> Zenodo: [10.5281/zenodo.21084339](https://doi.org/10.5281/zenodo.21084339)

The paper itself (Markdown + PDF) is in [`paper/`](paper/).

---

## What this is

The global regularity of the 3D incompressible Navier–Stokes equations turns on whether
vortex stretching can drive the enstrophy to infinity in finite time. A classical
obstruction is **geometric depletion** (Constantin–Fefferman): the nonlinear enstrophy
production `ω·Sω` is suppressed below its dimensional bound `|ω|²|S|` because vorticity
preferentially aligns with the *intermediate* eigenvector of the strain tensor.

This work measures the depletion ratio

```
D = ⟨ω·Sω⟩ / ⟨|ω|² |S|⟩
```

**as a function of vorticity percentile** — i.e. in the most intense events, the ones that
would actually seed a singularity — and asks **how it scales with Reynolds number**.

**Main result (v0.5).** Extreme-tail depletion is *resolved below unity* (regularizing) at
the top-1% and top-0.1% at both `R_λ ≈ 1300` (`isotropic8192`, 77M points) and `R_λ = 2500`
(`isotropic32768` exascale, 164M points); block-resampled 95% intervals exclude 1. A **mild
de-depletion with Reynolds number** is statistically resolved at top-1%. Reported as a
two-point comparison, **not** a scaling law; the deepest tail (top-0.01%) is resolution-limited.
Consistent with Buaria–Pumir.

## Repository layout

```
paper/    the paper (.md + .pdf) and the cross-Reynolds figure
code/     DNS solver, JHTDB pipeline, 1D model, analysis, bootstrap, figure code
```

### Key scripts (`code/`)

| File | Role |
|------|------|
| `dns_spectral.py`, `dns_instrumented.py` | In-house pseudospectral DNS (validation + canonical alignment) |
| `fd_vs_spectral.py` | Finite-difference vs. spectral gradient sanity check |
| `sciserver_cutout.py`, `sciserver_depletion.py` | Production JHTDB pipeline (SciServer cutouts → depletion statistic) |
| `ns_jhtdb.py`, `d7_jhtdb_batch.py` | Early JHTDB pipeline / batched testing-token accumulation |
| `master_ns_run.py` | Orchestrator for the cross-Reynolds production runs |
| `osw_1d.py`, `osw_frontier.py` | 1D OSW model solver (finite-time blow-up scan of the depletion mechanism) |
| `shell_model.py` | Shell-model corroboration |
| `analyze_pool.py`, `block_uncertainty.py` | Analysis and block-resample uncertainty (not point bootstrap) |
| `plot_crossRe.py`, `plot_depletion.py` | Figure code |
| `d3_rescan.py`, `d4_extreme.py`, `d5_decisive.py`, `e2*.py` | Exploratory percentile / extreme-tail scans (research trail) |
| `*_results.md`, `RESULTS.md` | Run provenance |

## Data

The turbulence fields are the public **Johns Hopkins Turbulence Database** (JHTDB), datasets
`isotropic8192` and `isotropic32768`, accessed via `giverny` / SciServer. The geometric-depletion
statistic and the per-snapshot audit are reproducible directly from them — no private data is required.

- The scripts default to the **public** JHTDB testing token `edu.jhu.pha.turbulence.testing-201406`,
  which is rate-limited to small queries (fine for the DNS validation and small cutouts).
- The high-statistics production runs need a **free registered JHTDB token**. Set it via
  `export JHTDB_TOKEN=...` (a placeholder `YOUR_REGISTERED_JHTDB_TOKEN` appears in the docs/help text).
- Register at <https://turbulence.pha.jhu.edu/>.

### Reproduction notes / gotchas

- Production cutouts run on SciServer; `getCutout` with `B=128` hangs on the `isotropic32768`
  cube — use `B=64`.
- The per-box pool is checkpointed as an `.npz` that is rewritten whole each box: **never analyze a
  pool `.npz` while a run is in progress** (it will be a torn/corrupt read). Read the result from the
  `MASTER_RESULTS.md` / `CUTOUT FINAL` log line instead.

## Provenance & method

This is a human-validated, AI-accelerated research artifact. An agentic large-language-model system
(operated by Saluca Labs) performed literature mapping, code authoring, and closed-loop analysis;
every quantitative claim rests on the human-validated, reproducible artifacts in this repository. No
result is treated as established solely because an AI produced it. See the paper's *Provenance* and
*Data and code availability* sections for the full disclosure.

## License

[MIT](LICENSE) © 2026 Cristian Ruvalcaba / Saluca Labs.
