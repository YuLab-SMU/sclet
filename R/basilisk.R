#' @importFrom basilisk BasiliskEnvironment
# SCENIC (pySCENIC) environment.
#
# pySCENIC 0.12.1 has known incompatibilities with newer transitive
# dependencies, so this environment pins a full, tested lock instead of only the
# top-level packages (#27):
#   * setuptools==80.9.0  - keeps `pkg_resources` importable (setuptools >= 81
#     removed it, breaking `import pkg_resources` in several pySCENIC deps).
#   * dask/distributed 2024.2.1 + dask-expr 0.5.3 - the scheduler pair that
#     avoids the `Nanny failed to start` / `diagnostics_port` regressions seen
#     with newer Dask.
#   * arboreto 0.1.6 - GRNBoost2 backend, pinned.
#   * pyarrow 20.0.0 - required for the intermediate feather/parquet
#     checkpointing used by `RunSCENIC(save_intermediate = TRUE)`.
#   * numpy 1.26.4 - kept for the rest of the stack; `np.object` is restored at
#     runtime by `sclet_numpy_compat()` because pySCENIC 0.12.1 still references
#     it (removed in numpy >= 1.24).
# `pyscenic` and `loompy` are PyPI-only, so they are installed via `pip` on both
# the conda and the pip/venv basilisk backends (#23).
sclet_scenic_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_scenic_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2"
    ),
    pip = c(
        "setuptools==80.9.0",
        "pyscenic==0.12.1",
        "loompy==3.0.7",
        "arboreto==0.1.6",
        "dask==2024.2.1",
        "distributed==2024.2.1",
        "dask-expr==0.5.3",
        "pyarrow==20.0.0"
    )
)

sclet_scvi_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_scvi_env",
    pkgname = "sclet",
    packages = c(
        "python=3.12",
        "scanpy=1.10.1",
        "anndata=0.11.4",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "scvi-tools=1.4.3"
    )
)

sclet_cellrank_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_cellrank_env",
    pkgname = "sclet",
    packages = c(
        "python=3.12",
        "scanpy=1.10.1",
        "anndata=0.11.4",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "scvelo=0.3.4",
        "cellrank=2.0.7"
    )
)

sclet_regvelo_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_regvelo_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.3",
        "anndata=0.10.8",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "scvelo=0.3.2",
        "scvi-tools=1.1.6",
        "pytorch=2.5.1"
    ),
    pip = c(
        "regvelo==0.4.2"
    )
)

sclet_cell2location_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_cell2location_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "scvi-tools=1.2.2",
        "cell2location=0.1.5"
    )
)

sclet_celloracle_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_celloracle_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "celloracle=0.20.0"
    )
)

sclet_cellphonedb_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_cellphonedb_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "anndata=0.10.7",
        "pandas=2.2.2"
    ),
    pip = c(
        "cellphonedb==5.0.0"
    )
)

# Resolve a basilisk environment object by its short name. This is what the
# public troubleshooting guidance uses, e.g.
#   basilisk::obtainEnvironmentPath(sclet:::.env("sclet_cellrank_env"))
# It returns the package-namespace BasiliskEnvironment object for `name`.
.env <- function(name) {
    envs <- list(
        sclet_scenic_env = sclet_scenic_env,
        sclet_scvi_env = sclet_scvi_env,
        sclet_cellrank_env = sclet_cellrank_env,
        sclet_regvelo_env = sclet_regvelo_env,
        sclet_cell2location_env = sclet_cell2location_env,
        sclet_celloracle_env = sclet_celloracle_env,
        sclet_cellphonedb_env = sclet_cellphonedb_env
    )
    if (!name %in% names(envs)) {
        cli::cli_abort(
            c(
                "Unknown sclet basilisk environment: {.val {name}}.",
                i = "Valid names are: {.val {names(envs)}}."
            ),
            class = "sclet_unknown_env"
        )
    }
    envs[[name]]
}

# Locate the directory holding the shared C/C++ libraries bundled inside a
# basilisk-managed conda environment. Some layouts use `lib`, others `lib64`.
sclet_env_libdir <- function(envpath) {
    cands <- c(file.path(envpath, "lib"), file.path(envpath, "lib64"))
    found <- cands[dir.exists(cands)]
    if (length(found)) {
        found[[1]]
    } else {
        NULL
    }
}

# Prepend the environment's own library directory to LD_LIBRARY_PATH so that
# compiled Python extensions (e.g. optree, scipy, pytorch, c-extensions built
# against a recent GCC) resolve a libstdc++/libgcc that is new enough instead of
# falling back to an outdated system one (e.g. `GLIBCXX_3.4.31 not found`).
# This only helps when the environment actually bundles a newer libstdc++:
# conda-backed basilisk environments ship `libstdc++-ng`, whereas the pip/venv
# backend leaves the env `lib`/`lib64` directory empty, so on that backend this
# is a no-op and the host runtime is used. The preceding variable is restored by
# the caller after the basilisk process has finished.
sclet_prepend_env_lib <- function(env) {
    envpath <- basilisk::obtainEnvironmentPath(env)
    libdir <- sclet_env_libdir(envpath)
    if (is.null(libdir)) {
        return(invisible(NULL))
    }
    old <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
    if (nzchar(old)) {
        Sys.setenv(LD_LIBRARY_PATH = paste(libdir, old, sep = ":"))
    } else {
        Sys.setenv(LD_LIBRARY_PATH = libdir)
    }
    invisible(libdir)
}

# Turn a cryptic `GLIBCXX_... not found` / `CXXABI_... not found` /
# `libstdc++.so.6` import failure into an actionable error. The root cause is a
# host/toolchain mismatch (the C++ runtime the process actually uses is too old
# for a compiled Python extension). Detect both symbol families, since modern
# wheels can demand a newer `CXXABI_` version even when `GLIBCXX_` is satisfied.
sclet_glibcxx_error <- function(env, e) {
    msg <- conditionMessage(e)
    is_glibcxx <- grepl("(GLIBCXX|CXXABI)_[0-9]+\\.[0-9]+\\.[0-9]+", msg) &&
        grepl("not found", msg)
    if (!is_glibcxx) {
        stop(e)
    }
    envname <- env@envname
    cli::cli_abort(
        c(
            "The {.val {envname}} Python environment could not load a compiled C++ extension: the C++ runtime in use is too old.",
            "!" = "A Python extension was built against a newer {.file libstdc++} than is currently being resolved (e.g. {.val GLIBCXX_3.4.31} or {.val CXXABI_1.3.15} not found in {.file libstdc++.so.6}).",
            "i" = "This is a host/toolchain limitation, not a sclet bug, and it cannot be fixed from inside R.",
            "i" = "Either (a) provide a newer system C++ runtime (e.g. a newer {.pkg libstdc++6} / {.pkg gcc}; on Ubuntu that means >= 22.04 with GCC >= 12), or (b) run the environment on a conda-backed basilisk backend so the environment's own newer {.file libstdc++} is used. See the sclet NEWS for the minimum host requirement.",
            "i" = "Original error: {msg}"
        ),
        class = "sclet_glibcxx_not_found"
    )
}

# pandas >= 2.2 dropped the legacy `get_values()` API on Series/Categorical/Index
# (e.g. `'Categorical' object has no attribute 'get_values'`). Some dependency in the
# sclet Python stack still calls it, so restore it as an alias for the modern
# `to_numpy()` in the launched Python process. This runs before the user callback in
# `sclet_basilisk_run()`, is a no-op where `get_values` already exists (older pandas),
# and only affects the current process. It is package-scope (not a closure) so it
# serializes to the basilisk worker without capturing any S4 object.
sclet_pandas_compat <- function() {
    reticulate::py_run_string(paste0(
        "import pandas as pd\n",
        "for _c in (pd.Series, pd.Categorical, pd.Index):\n",
        "    if not hasattr(_c, 'get_values'):\n",
        "        _c.get_values = _c.to_numpy\n"
    ))
}

# numpy >= 1.24 removed the deprecated `np.object` alias, but pySCENIC 0.12.1
# still references it (e.g. `pyscenic/transform.py`), raising
# `AttributeError: module 'numpy' has no attribute 'object'`. Restore it as an
# alias for the Python builtin `object` in the launched process before any
# sclet callback imports pySCENIC. This runs in `sclet_basilisk_run()` ahead of
# the user callback; it is a no-op on numpy versions that still ship `np.object`
# and only affects the current process.
sclet_numpy_compat <- function() {
    reticulate::py_run_string(paste0(
        "try:\n",
        "    import numpy as np\n",
        "    if not hasattr(np, 'object'):\n",
        "        np.object = object\n",
        "except Exception:\n",
        "    pass\n"
    ))
}

# Turn a failure to create/activate an environment into an actionable error.
# basilisk deletes a partially-created virtualenv/conda env on setup failure, so
# a broken setup cannot be inspected in place — this points the user at the
# documented external-Python repair path instead of a generic R-side message.
sclet_env_setup_error <- function(env, e) {
    msg <- conditionMessage(e)
    is_glibcxx <- grepl("(GLIBCXX|CXXABI)_[0-9]+\\.[0-9]+\\.[0-9]+", msg) &&
        grepl("not found", msg)
    if (is_glibcxx) {
        return(sclet_glibcxx_error(env, e))
    }
    envname <- env@envname
    cli::cli_abort(
        c(
            "The {.val {envname}} Python environment could not be created or started.",
            "!" = "sclet cannot reliably create/repair a basilisk environment from inside R: on failure basilisk deletes the partially-created environment, and the useful pip/build traceback only appears in the R/Jupyter session output rather than in the R-side error.",
            "i" = "Original error: {msg}",
            "i" = "To inspect or repair the environment, build it manually and point basilisk at it (see the sclet NEWS 'RunSCENIC' entry for the documented steps). The typical setup is:",
            " " = "Sys.setenv(BASILISK_CUSTOM_PYTHON_3_10 = '<your-python-3.10>', BASILISK_EXTERNAL_DIR = '<your-env-root>', BASILISK_NO_PYENV = '1')",
            "i" = "If the message instead concerns a missing C++ symbol (GLIBCXX_/CXXABI_... not found), that is a host-toolchain limitation; see the sclet host-requirement guidance."
        ),
        class = "sclet_env_setup_failed"
    )
}

# Start a basilisk environment with its bundled libraries taking precedence
# over system ones, run `fun`, and tear the environment down again. This is the
# single entry point used for every sclet basilisk call.
#
# `full.activation = TRUE` makes basilisk run the environment's own activation
# (for a conda env that sources its `activate` script). This is what actually
# puts the environment's bundled `lib` directory first in the runtime search
# path of the worker that imports the Python extensions, so e.g. a conda env
# that ships a newer `libstdc++.so.6` is used instead of an outdated system one
# (the `GLIBCXX_/CXXABI_... not found` errors). It is a no-op for the default
# pip/venv backend beyond what the explicit LD_LIBRARY_PATH prepending below
# already does, so it is safe to leave on for every backend.
sclet_basilisk_run <- function(
    env,
    fun,
    ...,
    full.activation = TRUE,
    fork = basilisk::getBasiliskFork(),
    shared = basilisk::getBasiliskShared()
) {
    old_ld <- Sys.getenv("LD_LIBRARY_PATH", unset = NA)
    tryCatch(
        sclet_prepend_env_lib(env),
        error = function(e) sclet_env_setup_error(env, e)
    )
    on.exit({
        if (is.na(old_ld)) {
            Sys.unsetenv("LD_LIBRARY_PATH")
        } else {
            Sys.setenv(LD_LIBRARY_PATH = old_ld)
        }
    }, add = TRUE)

    proc <- tryCatch(
        basilisk::basiliskStart(
            env,
            full.activation = full.activation,
            fork = fork,
            shared = shared
        ),
        error = function(e) sclet_env_setup_error(env, e)
    )
    on.exit(basilisk::basiliskStop(proc), add = TRUE)

    # Restore the pandas `get_values()` API and the numpy `np.object` alias in the
    # worker before any callback runs, so dependencies that still rely on them
    # work on the pinned pandas/numpy versions (no-op otherwise).
    basilisk::basiliskRun(proc, sclet_pandas_compat)
    basilisk::basiliskRun(proc, sclet_numpy_compat)

    tryCatch(
        basilisk::basiliskRun(proc, fun = fun, ...),
        error = function(e) sclet_glibcxx_error(env, e)
    )
}
