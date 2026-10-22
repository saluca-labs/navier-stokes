# AGENTS.md

Instructions for AI coding agents working in this repository.

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
README.md, the paper or RESULTS.md changed, `llms-full.txt` is regenerated in the same change.
