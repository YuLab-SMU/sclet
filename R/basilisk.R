#' @importFrom basilisk BasiliskEnvironment
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
        "pyscenic==0.12.1",
        "loompy==3.0.7"
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
    sclet_prepend_env_lib(env)
    on.exit({
        if (is.na(old_ld)) {
            Sys.unsetenv("LD_LIBRARY_PATH")
        } else {
            Sys.setenv(LD_LIBRARY_PATH = old_ld)
        }
    }, add = TRUE)

    proc <- basilisk::basiliskStart(
        env,
        full.activation = full.activation,
        fork = fork,
        shared = shared
    )
    on.exit(basilisk::basiliskStop(proc), add = TRUE)

    tryCatch(
        basilisk::basiliskRun(proc, fun = fun, ...),
        error = function(e) sclet_glibcxx_error(env, e)
    )
}
