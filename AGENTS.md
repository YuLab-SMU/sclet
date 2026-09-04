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

## Quick reference

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
