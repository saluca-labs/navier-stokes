# Documentation map

Every document in this repository, grouped by what you want to know. One line each on what it
answers.

## Start here

- [README](../README.md): what the study measures, the v0.6 result and its correction, install,
  the CPU quickstart, and the production data path diagram.
- [AGENTS.md](../AGENTS.md): orientation for an AI agent working in the repository.

## The result, in depth

- [Paper, Markdown](../paper/geometric_depletion_paper.md): the v0.6 paper, methods, tables and
  change log. The figure it draws from code is [`paper/depletion_crossRe.png`](../paper/depletion_crossRe.png).
- [Paper, PDF](../paper/geometric_depletion_paper.pdf): the typeset v0.6.
- Archived record: concept DOI [10.5281/zenodo.20970675](https://doi.org/10.5281/zenodo.20970675),
  which resolves to the latest version.

## Run logs (archived output, do not edit except as a dated correction)

- [code/RESULTS.md](../code/RESULTS.md): dated log of the exploratory experiments (E1 shell model,
  E2 delay and relay variants, D1 to D5 GPU DNS scans).
- [code/osw_results_sin.md](../code/osw_results_sin.md),
  [code/osw_results_tri.md](../code/osw_results_tri.md) and
  [code/osw_results.md](../code/osw_results.md): archived 1D OSW model scans (outcome, blow-up time
  and peak vorticity per advection strength `a`). The current `osw_1d.py` writes
  `osw_results_<ic>.md` for its `--ic` value (default `sin`); `osw_results.md` predates that naming,
  and the header of `osw_results_tri.md` records `data=sin(x)`, so read each file's own header for
  the initial condition it used.
- [code/osw_frontier_results.md](../code/osw_frontier_results.md): low-Hölder data for De Gregorio,
  written by `osw_frontier.py`.

## Reproduce or build on it

- [README, Usage](../README.md#usage): each script with a real command and whether it was run.
- [README, Data](../README.md#data): JHTDB access, the token, and the reproduction gotchas.
- [README, Architecture](../README.md#architecture): which scripts need JHTDB, a CUDA GPU, or only
  a CPU, and what is not in this repository.
- [CITATION.cff](../CITATION.cff): how to cite.

## Images

- [docs/img/shell-model-threshold-2026-10-22.png](img/shell-model-threshold-2026-10-22.png): captured
  output of `python code/shell_model.py`; the command and text are in the matching `.txt`.
