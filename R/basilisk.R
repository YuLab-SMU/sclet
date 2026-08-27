#' @importFrom basilisk BasiliskEnvironment
sclet_scenic_env <- basilisk::BasiliskEnvironment(
    envname = "sclet_scenic_env",
    pkgname = "sclet",
    packages = c(
        "python=3.10",
        "scanpy=1.10.1",
        "anndata=0.10.7",
        "numpy=1.26.4",
        "pandas=2.2.2",
        "pyscenic=0.12.1",
        "loompy=3.0.7"
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
# against a recent GCC) resolve a libstdc++/libgcc that is new enough instead
# of falling back to an outdated system one (e.g. `GLIBCXX_3.4.31 not found` on
# Ubuntu 20.04, whose libstdc++ only reaches GLIBCXX_3.4.29). The variable is
# restored by the caller after the basilisk process has finished.
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

# Start a basilisk environment with its bundled libraries taking precedence
# over system ones, run `fun`, and tear the environment down again. This is the
# single entry point used for every sclet basilisk call.
sclet_basilisk_run <- function(
    env,
    fun,
    ...,
    full.activation = NA,
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

    basilisk::basiliskRun(proc, fun = fun, ...)
}
