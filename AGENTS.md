# AGENTS.md — sclet workspace layout

This is the single entry point for the whole `sclet` workspace. The workspace
contains **separate, nested `git` checkouts** of *different branches* of the
`YuLab-SMU/sclet` repository. The rule is simple:

> **Each directory is its own checkout of its own branch.**
> Entering a directory is equivalent to being on that branch. Do not assume a
> file you see in one directory is the same version as the one in another.

## Directory ↔ branch mapping

| Directory (`<workspace>/`)  | Git branch | What it contains                                                    | Remote |
|-----------------------------|------------|---------------------------------------------------------------------|--------|
| `.` (this top directory)    | `devel`    | **`sclet` R package** source (`R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `tests/`) | `git@github.com:YuLab-SMU/sclet.git` |
| `bookdown/`                 | `docs`     | Bookdown **source** (`*.Rmd`, `_bookdown.yml`, `index.Rmd`, `references.bib`, ...) | `git@github.com:YuLab-SMU/sclet.git` |
| `bookdown/docs/`            | `gh-pages` | **Published site** (built `.html`, `search_index.json`, `libs/`, `figures/`, `.nojekyll`) | `git@github.com:YuLab-SMU/sclet.git` |

`bookdown` is already git-ignored at this top level, so the nested checkouts do
not pollute the R package repository's index.

## What goes where

- **R package work** (`sclet` the package): this top-level checkout, branch `devel`.
- **Edits to the book**: change files under `bookdown/` (the `docs` branch).
  Add/edit `*.Rmd`, `_bookdown.yml`, `references.bib`, figures, etc. here.
- **Rebuilding / deploying the site**: rebuild the book, and the generated
  `docs/` output is the `bookdown/docs/` checkout (the `gh-pages` branch).
  Publish `gh-pages` to update the live book site.

## Pitfalls to avoid

- Do **not** mix `.Rmd` source and built `.html` output. They live in different
  branches and different directories for a reason.
- Each directory is an **independent repository**. Moving/renaming files across
  them does not move them across branches — each has its own `.git`.
- A change committed in `bookdown/` is **not** visible in `bookdown/docs/` (and
  vice-versa) until the book is actually rebuilt and the `gh-pages` branch is
  re-pushed.
- Never run a command in `bookdown/docs/` expecting it to affect the `docs`
  source branch, and never run bookdown builds from the wrong directory.
- When you `cd` into a directory, keep in mind it is now a **different branch** —
  this top directory is `devel`, not `docs`.

## R package development notes

Lessons captured from fixing real issues, so future maintenance avoids the same
traps. These apply to the `devel` checkout (`R/`, `man/`, `NAMESPACE`, ...).

- **Sparse assays break reshape/plot code.** The `logcounts` assay is a sparse
  `dgCMatrix` by default for large datasets. `t()` on a sparse matrix dispatches to
  `t.default()` and fails with `argument is not a matrix`; `as.data.frame()` on a
  still-sparse transpose can also fail. When reshaping an assay, materialize *just
  the requested features* with `as.matrix()` before transposing:
  `as.data.frame(t(as.matrix(assay(sce, assay_name)[features, , drop = FALSE])))`.
  Never transpose the full assay.
- **Keep zero-expression cells in plot/fit data.** `geom_smooth` and GAM fit
  `Expression ~ pseudotime`; cells with 0 expression anchor the low end of the
  curve. Dropping them — e.g. by going straight from sparse `summary()` triplets,
  which omit structural zeros — biases the curve upward (measured ~2x in a test
  case). Densify (or otherwise include) the zeros; otherwise the fitted curves are
  wrong.
- **Prefer positional alignment over key joins for speed.** In an SCE, assay columns
  and `colData` rows are guaranteed to be in the same order, so `cbind` +
  `reshape2::melt` is both correct and fast. Reserve `merge(by = "cell")` for real
  name-based joins — it is roughly 20x slower than the positional `cbind` + `melt`
  pipeline on large cell counts.

### Basilisk / Python backend

- **Basilisk `packages` are pip-installed and `channels` is ignored (basilisk ≥ 1.24).**
  `setupBasiliskEnv` only calls `virtualenv_create`: it rewrites `pkg=ver` →
  `pkg==ver`, merges `packages` + `pip`, and pip-installs everything; conda/mamba is
  never invoked and the `channels` argument is a no-op. Any package that isn't on PyPI,
  or whose version doesn't exist on a conda channel, must go in the `pip = c(...)` field
  of `BasiliskEnvironment` (confirmed: `pyscenic` is bioconda-only, `loompy=3.0.7` is PyPI-only).
- **A basilisk error is raised by `clusterCall`; the real exception is wrapped inside.**
  `basiliskStart` runs the callback in a fork/PSOCK cluster worker, so a Python failure
  surfaces as `checkForRemoteErrors(lapply(cl, recvResult))` wrapping the actual error.
  Read the inner `error:` message, and get the full Python traceback with `reticulate::py_last_error()`.
- **Use `full.activation = TRUE` for conda-backed envs; keep callbacks package-scope and
  call `basiliskRun` more than once on the same proc.** Partial activation doesn't source
  the env's own activation, so a bundled newer `libstdc++.so.6` isn't first in the worker's
  runtime search path; a module then loads the older system libstdc++ first and (glibc loads
  one `libstdc++.so.6` per process) every C++ extension binds to it → `GLIBCXX_/CXXABI_... not found`.
  `full.activation = TRUE` sources it and is a no-op on the pip/venv backend (which bundles no
  libstdc++). Because `basiliskRun(proc, fun)` serializes `fun` to the worker, use package-scope
  functions (closures capturing S4 objects break), and you can call `basiliskRun` repeatedly on
  the same `proc` to run a setup/patch before the real callback.
- **pandas ≥ 2.2 removed `get_values()` from `Series`/`Categorical`/`Index`** (error:
  `'Categorical' object has no attribute 'get_values'`). If a dependency still calls it, restore
  it as an alias for `to_numpy()` once in the shared basilisk wrapper (covers every backend):
  `for _c in (pd.Series, pd.Categorical, pd.Index): _c.get_values = _c.to_numpy`, guarded by
  `if not hasattr(_c, 'get_values')`.
- **Wrap Python import failures once, in the shared wrapper, and raise an actionable condition.**
  Detect `GLIBCXX_/CXXABI_... not found` (and similar) in `sclet_basilisk_run` and rethrow as a
  typed `sclet_glibcxx_not_found` error that states it is a host/toolchain issue and lists the two
  real options (newer system C++ runtime, or a conda-backed env) instead of a raw `ImportError`.
- **Environment setup failures deserve their own typed error.** basilisk deletes a
  partially-created env on setup failure, so it cannot be inspected in place and the useful
  pip/build traceback only shows in the R/Jupyter output. Wrap `obtainEnvironmentPath` +
  `basiliskStart` in `sclet_basilisk_run` and rethrow as `sclet_env_setup_failed` pointing at the
  external-Python repair (`Sys.setenv(BASILISK_CUSTOM_PYTHON_3_10 = ..., BASILISK_EXTERNAL_DIR = ...,
  BASILISK_NO_PYENV = '1')` + a manually created venv) instead of a generic message.
- **Use `convert = FALSE` when threading Python objects through a multi-step library.** For
  pySCENIC/pySCENIC-like pipelines you pass *lists of Python objects* (modules, ranking-database
  objects) between steps; default `reticulate` auto-conversion mangles them. Import with
  `convert = FALSE`, build lists with `import_builtins(convert = FALSE)$list(...)`, and only
  `reticulate::py_to_r()` the final result.
- **numpy ≥ 1.24 removed `np.object`; pySCENIC 0.12.1 still references it** (`AttributeError:
  module 'numpy' has no attribute 'object'`). Restore it as an alias for the builtin `object` in the
  shared wrapper (run *before* the callback imports pyscenic), guarded like the pandas `get_values`
  shim: `if not hasattr(np, 'object'): np.object = object`.

### Python environments on this machine (where they live)

Several ways of having Python here are all legitimate — **do not force one**. The point of
this section is that you can *find* an existing environment instead of concluding a
dependency is missing. What exists on this machine:

- **`uv` is available** and is the handiest way to create a throwaway venv (no conda/mamba
  is on `PATH`). Verified 2026-10-09:

  | | |
  |---|---|
  | `uv` binary | `~/.local/bin/uv` (0.12.4) |
  | `uv python dir` | `~/.local/share/uv/python` (includes cpython-3.12.13) |
  | `uv cache dir` | `~/.cache/uv` |
  | `uv tool dir` | `~/.local/share/uv/tools` |

  There is no `~/.config/uv/uv.toml`, no `~/.pip/pip.conf` and no `PIP_INDEX_URL`, so installs
  go straight to PyPI. If one hangs or times out, apply the proxy rule from the global
  instructions (`http://127.0.0.1:7897`) rather than concluding the package is unavailable.

  ```bash
  uv venv --python 3.12          # creates ./.venv from a uv-managed interpreter
  uv pip install <pkg>           # activation is unnecessary: use ./.venv/bin/python
  ```

  Sibling projects happen to keep a per-project `.venv` (`uv venv`'s default), so that layout
  will look familiar — but a plain `python -m venv` is equally fine, and whatever already
  exists should be reused rather than replaced. Seeing `warning: Failed to hardlink files;
  falling back to full copy` just means the venv is on a different filesystem from
  `~/.cache/uv` (e.g. under `/tmp`); it is harmless, and `UV_LINK_MODE=copy` silences it.
  Reproduced here: a venv under `/tmp` warns, one under `$HOME` does not.

- **Environment-variable overrides** for pointing tools at an existing interpreter:
  `BASILISK_CUSTOM_PYTHON_3_10`, `BASILISK_EXTERNAL_DIR`, `BASILISK_NO_PYENV=1` (see the
  environment-setup bullet in the Basilisk section above).

**basilisk — what sclet itself uses.** sclet's own Python dependencies are declared as
`BasiliskEnvironment` objects in `R/basilisk.R` and materialised under:

```
~/.cache/R/basilisk/<basilisk-version>/<pkgname>/<pkg-version>/<envname>/
# e.g. ~/.cache/R/basilisk/1.24.0/sclet/1.0.1/sclet_cellrank_env/
```

(`velociraptor` nests one level deeper: `.../velociraptor/1.22.0/env/`.)

The seven declared envs and their pins:

| envname | python | notable pins |
|---|---|---|
| `sclet_scenic_env` | 3.10 | scanpy 1.10.1, anndata 0.10.7, pySCENIC 0.12.1 + the SCENIC pin set in the section below |
| `sclet_scvi_env` | 3.12 | scanpy 1.10.1, anndata 0.11.4, scvi-tools 1.4.3 |
| `sclet_cellrank_env` | 3.12 | scanpy 1.10.1, anndata 0.11.4, scvelo 0.3.4, cellrank 2.0.7 |
| `sclet_regvelo_env` | 3.10 | scanpy 1.10.3, scvelo 0.3.2, scvi-tools 1.1.6, pytorch 2.5.1, regvelo 0.4.2 (pip) |
| `sclet_cell2location_env` | 3.10 | scvi-tools 1.2.2, cell2location 0.1.5 |
| `sclet_celloracle_env` | 3.10 | celloracle 0.20.0 |
| `sclet_cellphonedb_env` | 3.10 | cellphonedb 5.0.0 (pip) |

**Envs are built lazily on first use — a missing directory means "not built yet", not
"broken".** Only two Python envs existed on this machine as of 2026-10-09:

- `~/.cache/R/basilisk/1.24.0/velociraptor/1.22.0/env/` — scvelo 0.3.3, scanpy 1.11.3,
  anndata 0.11.4, numpy 2.2.6, pandas 2.3.1, scipy 1.15.3. **This is the env the velocity
  tests actually execute in**, via `velociraptor::scvelo()`, not one of sclet's own.
- `~/.cache/R/basilisk/1.24.0/sclet/1.0.1/sclet_cellrank_env/` — scvelo 0.3.4, cellrank 2.0.7,
  scanpy 1.10.1, anndata 0.11.4 (matches sclet's own pins).

The path is keyed by the R package **version**, so a `DESCRIPTION` version bump
invalidates the cache: the directory above says `sclet/1.0.1` while `DESCRIPTION` is now
1.0.2, so basilisk will look under `sclet/1.0.2/` and rebuild from scratch. Expect a long
first run after a version bump; do not misread it as a network failure.

To resolve an env's path from R — note that `obtainEnvironmentPath()` **creates** the env
when absent, which can trigger a long install, so prefer reading the on-disk path above:

```r
basilisk::obtainEnvironmentPath(sclet:::.env("sclet_cellrank_env"))
```

The legacy `inst/regvelo-smoke-test/setup_pyenv_regvelo.sh` builds a **pyenv** venv aimed at a
GPU cluster (Torch CUDA wheels) and is not the usual path here — but a pyenv- or
`python -m venv`-created environment is perfectly acceptable; if one already exists and has
what you need, use it rather than rebuilding with a different tool.

### SCENIC / pySCENIC

- **`pyscenic.prune.prune2df()` takes `(rnkdbs, modules, motif_annotations_fname, ...)`, NOT `(dbs, motif_f, adj)`.** The old sclet call passed the raw `.feather` paths and the GRNBoost2
  `adj` in the slots pySCENIC expects for the ranking databases and the co-expression modules — it
  either crashes or silently mis-prunes. The correct sequence is: `modules =
  list(modules_from_adjacencies(adj, ex_df))`, wrap each `.feather` as
  `ctxcore.rnkdb.FeatherRankingDatabase(fname = path, name = basename)`, then
  `prune2df(db_objects, modules, motif_annotations_path, ...)`, then `df2regulons()`.
- **Motif pruning and AUCell are not multi-step Dask consumers the same way.** Pass
  `client_or_address = "custom_multiprocessing"` to `prune2df` to avoid spinning up a Dask
  scheduler/`Nanny` (which newer dask breaks with `Nanny failed to start` / `diagnostics_port`);
  `pyscenic.aucell.aucell()` takes `num_workers`/`seed` directly, and `arboreto.grnboost2()` takes
  `client_or_address` (a `dask.distributed.Client`) rather than `num_workers`.
- **Pin a full, tested dependency set for SCENIC, not just the top-level packages.**
  pySCENIC 0.12.1 + current numpy/pandas need: `setuptools==80.9.0` (keeps `pkg_resources`
  importable; removed in setuptools ≥ 81), `arboreto==0.1.6`, `dask==2024.2.1`,
  `distributed==2024.2.1`, `dask-expr==0.5.3`, `pyarrow==20.0.0` (needed for the parquet
  checkpoints used by `RunSCENIC(save_intermediate = TRUE)`).
- **A `BasiliskEnvironment` merges `packages` + `pip` into its single `@packages` slot** (for a pip
  backend everything is rewritten to `pkg==ver` and pip-installed). Introspect `env@packages` to
  assert the pinned set in a test.

### Testing & doc maintenance

- **Basilisk-backed functions are tested by mocking the shared wrapper, not by building an env.**
  You do not need a Python env to unit-test `RunSCENIC`/`RunCellRank`/...: use
  `testthat::with_mocked_bindings(..., .package = "sclet")` to replace `sclet_basilisk_run` with a
  stub that returns a small matrix, then assert the altExp/state contract and the args the stub
  captured. `pkgload::load_all` defaults to `export_all = TRUE`; if you pass `export_all = FALSE`,
  testthat cannot see internal functions and every call to `sclet_*` fails with
  `could not find function "sclet_..."` — that is a harness misconfiguration, not a real bug.
- **Don't let roxygen2 regenerate the whole `NAMESPACE` for a small change.** A newer roxygen2
  reformats `NAMESPACE` and can *drop* imports that were only present in the hand-maintained file
  (this broke 30+ tests in unrelated files). When you only add params/comments, restore `NAMESPACE`
  and let the one `man/*.Rd` be the only regenerated artifact; leave `DESCRIPTION` and unrelated
  `man/*.Rd` alone.

- **Run package checks through the repository Makefile.** Do not invoke `R CMD check` directly for routine validation. Use `make check`; it runs `devtools::check()` with the package's local conventions. Use `make check2` only when an installable source tarball check is specifically needed, and use `make clean` to remove generated `*.Rcheck/` output. Before a broad check, run the focused `testthat` filter for the changed module.
- **Treat Makefile targets as the maintenance contract.** The `all` target runs documentation, checks, and cleanup; `rd` regenerates documentation and may touch unrelated Rd/NAMESPACE files, so review and revert unrelated changes before delivery. If check output reports package metadata problems, fix the metadata or environment rather than bypassing the Makefile with a direct check command.


```bash
# Confirm which branch you are on before editing
git -C .                 branch --show-current   # → devel
git -C bookdown          branch --show-current   # → docs
git -C bookdown/docs     branch --show-current   # → gh-pages

# All checkouts point at the same repository remote
git -C .                 remote get-url origin
git -C bookdown          remote get-url origin
git -C bookdown/docs     remote get-url origin
```
