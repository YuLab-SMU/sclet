# Basilisk / Python environment notes

Moved out of `AGENTS.md` on 2026-10-21 — these are facts and pins, not
instructions a reader needs on every task. Read this file when you touch
`R/basilisk.R`, a Python-backed action, or a Python env.

Every pin below was paid for by a failure. Do not "clean up" a version number
because it looks arbitrary.

## Environment declaration

- **`packages` are pip-installed; `channels` is a no-op** (basilisk ≥ 1.24).
  Non-PyPI or conda-only packages go in `pip = c(...)`. Assert pins by
  introspecting `env@packages` rather than trusting the declaration.
- sclet's basilisk envs live under
  `~/.cache/R/basilisk/<ver>/<pkg>/<pkg-ver>/<envname>/` (declared in
  `R/basilisk.R`; seven envs, see the pins there). They are built lazily — a
  missing directory means "not built yet", not "broken".
- **The path is keyed by the package version**, so any DESCRIPTION bump forces
  a rebuild. Observed 2026-10-20: the 1.0.1 → 1.0.2 bump made the CellRank env
  disappear and cost ~10 minutes to rebuild. Not a network failure.
- `obtainEnvironmentPath()` **CREATES** the env. Prefer reading the path
  directly when you only need to inspect it.
- Velocity tests execute in `velociraptor/1.22.0/env/`, not sclet's own env.

## Driving Python from R

- **Errors surface wrapped**: read the inner `error:` of
  `checkForRemoteErrors`; the full traceback is in
  `reticulate::py_last_error()`.
- **`full.activation = TRUE`** for conda envs, or libstdc++ ordering gives
  GLIBCXX errors.
- **Callbacks must be package-scope.** They are serialized to the worker;
  closures capturing S4 objects break.
- `basiliskRun` may be called repeatedly on one proc — do not assume a fresh
  interpreter.
- **`convert = FALSE`** when threading Python objects through multi-step
  libraries; call `py_to_r()` only on the final result.
- **Shims in the shared wrapper**, each guarded by `hasattr`: pandas
  `get_values` → `to_numpy`; `np.object` → `object`. Import and env-setup
  failures are wrapped as typed errors (`sclet_glibcxx_not_found`,
  `sclet_env_setup_failed`) carrying repair text, so a reader is told what to
  do rather than shown a traceback.

## pySCENIC pins

`prune2df(rnkdbs, modules, motif_annotations, ...)` with
`client_or_address = "custom_multiprocessing"`. The full pinned set:

```
setuptools==80.9.0
arboreto==0.1.6
dask==2024.2.1
distributed==2024.2.1
dask-expr==0.5.3
pyarrow==20.0.0
```

The dask family is pinned as one unit: mixing versions across
`dask`/`distributed`/`dask-expr` produces import errors that name none of them.

## Throwaway envs on this machine

Find an existing env before concluding a dependency is missing. `uv` creates
throwaway venvs:

```bash
uv venv --python 3.12
uv pip install <pkg>          # no activation
./.venv/bin/python script.py  # use the interpreter directly
```

Hardlink warnings under `/tmp` are harmless. `uv` needs no proxy override; if a
download times out, apply the machine's global proxy rule
(`127.0.0.1:7897`) rather than concluding the package is unavailable.

## Version bumps do not require a rebuild

The env path is keyed by package version, and basilisk treats a missing
directory as "not built yet" — so bumping DESCRIPTION *looks* like it demands
a fresh ~10-minute build. It does not, when the bump carries no dependency
change.

What was done for 1.0.2 -> 1.1.0 (2026-10-21):

```bash
cd ~/.cache/R/basilisk/1.24.0/sclet
ln -s 1.0.2 1.1.0          # one symlink, no copy
```

Verified end to end: `obtainEnvironmentPath()` resolves through it,
`dir.exists()` is TRUE, a basilisk worker launches, and cellrank 2.0.7 /
scvelo 0.3.4 import inside it (`sys.prefix` reports the 1.0.2 realpath, which is
correct — the symlink is a name, not a copy).

Two things that did NOT work and should not be repeated:

- **`cp -a` is the wrong tool here.** It copied 7 GB, and because it does not
  preserve hardlinks across the tree the copy was a full second copy, not a
  deduplicated one. Deleted immediately.
- **Editing the copied env's baked-in paths is a trap.** `pyvenv.cfg` and
  `bin/activate` carry absolute paths, and so does the shebang of every console
  script in `bin/` — patching them is a long list of edits with no benefit over
  a symlink. Do not go down that route.

Only rebuild when the pins actually change (`R/basilisk.R`), not when only the
package version does.
