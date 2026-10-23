# AGENTS.md

Instructions for AI coding agents working in this repository.

## What this is

The reproducibility package for the paper "Geometric depletion of enstrophy production in the
extreme vorticity tail" (v0.6). Research scripts that produced the paper's numbers, not a library.

## Read in this order

1. [`README.md`](README.md): the v0.6 result, its 2026-10-22 correction note, and the data path diagram.
2. [`code/sciserver_cutout.py`](code/sciserver_cutout.py): the production measurement (`depletion_of_box`).
3. [`code/block_uncertainty.py`](code/block_uncertainty.py): how the 95% intervals are computed.
4. [`code/coarse_grain_run.py`](code/coarse_grain_run.py): the v0.6 resolution-bias experiment.
5. [`paper/geometric_depletion_paper.md`](paper/geometric_depletion_paper.md): methods and change log.

## What you can and cannot see

- You can read and run every script in `code/`. `shell_model.py`, `osw_1d.py`, `osw_frontier.py` and
  `e2*.py` run on a CPU.
- You cannot reproduce the production numbers from this repository alone. The per-point pools
  (`.npz`) are not committed, and regenerating them needs JHTDB access through `giverny` (SciServer)
  or `givernylocal`, with a registered `JHTDB_TOKEN` for full statistics.
- The DNS scripts (`dns_*.py`, `d3_rescan.py`, `d4_extreme.py`, `d5_decisive.py`) need PyTorch with
  CUDA and have no CPU fallback.
- `experiments/` is referenced by the paper and by `d5_decisive.py` but is not in the repository.

## The only numbers to trust

- The main result (tail ratios, intervals, the resolution-bias correction): the v0.6 paper,
  [`paper/geometric_depletion_paper.md`](paper/geometric_depletion_paper.md), sections 6 and 6.1.
  The README's "Main result (v0.6)" restates it.
- The shell-model verdict: [`code/RESULTS.md`](code/RESULTS.md) experiment E1, reproduced by
  `python code/shell_model.py` (captured in `docs/img/shell-model-threshold-2026-10-22.out`).
- Not to trust: any v0.5 statement that the cross-Reynolds difference is "statistically resolved"
  (withdrawn in v0.6), the title string in `code/plot_crossRe.py`, and the Fig. 2 bands, which use
  `1.96` rather than the t multiplier (see the README's reproduction notes).

## Setup

```bash
python -m pip install numpy scipy matplotlib
```

Python 3.12. The DNS scripts also need PyTorch with CUDA; the JHTDB scripts need `giverny`
(SciServer) or `givernylocal`. On Windows set `PYTHONIOENCODING=utf-8`.

## Test

There is no test suite. The CPU check, from the repository root:

```bash
python code/shell_model.py
```

Passing looks like the verdict block: `gamma=0.00: BLOW-UP`, `gamma=0.15: BLOW-UP`,
`gamma=0.50: REGULAR (t_end=40.0000)`.

## Conventions

- numpy, scipy, matplotlib (and torch / giverny where already used); do not add dependencies.
- Paths are relative to the script or come from environment variables; never hardcode a machine path.
- The JHTDB token comes from `JHTDB_TOKEN`; only the public testing token may appear in code.
- No em-dashes in docs.

## Do not

- Edit `paper/` or the archived `code/*_results.md` / `RESULTS.md` except as a dated correction.
  `osw_1d.py` and `osw_frontier.py` overwrite their results files; run copies elsewhere.
- Commit `.npz` pools, logs, or a registered JHTDB token.
- Commit secrets, hostnames, internal paths or personal emails.

## Done means

The CPU check above passes, any script you touched still compiles (`python -m py_compile`), and if
README.md, the paper, RESULTS.md or anything under `docs/` changed, `llms-full.txt` is regenerated in
the same change.
