# Install sclet's version-free environment lookup into basilisk.
#
# `basilisk::basiliskStart()` resolves the environment directory itself, via
# `obtainEnvironmentPath()` in basilisk's namespace, and then hands that path to
# `useBasiliskEnv()`. Every path into a basilisk environment therefore goes
# through that one function: `basiliskStart()`, `basiliskRun()`, and
# `getBasiliskPython()`. There is no argument through which a caller can supply
# the directory.
#
# So the lookup is redirected once, here. For sclet's own environments it
# returns [sclet_env_path()], which is keyed on the environment specification
# instead of `DESCRIPTION`'s `Version`; every other package's environment is
# passed through untouched. See the comment above `sclet_env_spec_fingerprint()`
# for why the version must not reach the path.
#
# The rebinding has to modify basilisk's namespace, which is why it happens at
# load time rather than on first use.
sclet_env_shim_names <- function() {
    c(
        "sclet_scenic_env",
        "sclet_scvi_env",
        "sclet_cellrank_env",
        "sclet_regvelo_env",
        "sclet_cell2location_env",
        "sclet_celloracle_env",
        "sclet_cellphonedb_env",
        "sclet_gpu_env"
    )
}

# TRUE once the shim is live and points at sclet's own environments. Asserted by
# tests so the wiring cannot be silently lost; it never touches Python.
sclet_env_shim_installed <- function() {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        return(FALSE)
    }
    ns <- asNamespace("basilisk")
    if (!exists("obtainEnvironmentPath", envir = ns, inherits = FALSE)) {
        return(FALSE)
    }
    current <- get("obtainEnvironmentPath", envir = ns)
    isTRUE(attr(current, "sclet_shim")) && identical(
        current(sclet_scvi_env),
        sclet_env_path("sclet_scvi_env")
    )
}

sclet_install_env_shim <- function() {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        # basilisk backs only the Python-side actions; with it absent there is no
        # environment lookup to redirect.
        return(invisible(FALSE))
    }
    ns <- asNamespace("basilisk")
    if (!exists("obtainEnvironmentPath", envir = ns, inherits = FALSE)) {
        return(invisible(FALSE))
    }
    current <- get("obtainEnvironmentPath", envir = ns)
    if (isTRUE(attr(current, "sclet_shim"))) {
        return(invisible(TRUE))
    }
    original <- current
    handled <- sclet_env_shim_names()
    replacement <- function(env) {
        envname <- tryCatch(env@envname, error = function(e) NULL)
        if (is.character(envname) && length(envname) == 1L &&
            envname %in% handled) {
            return(sclet_env_path(envname))
        }
        original(env)
    }
    attr(replacement, "sclet_shim") <- TRUE
    unlockBinding("obtainEnvironmentPath", ns)
    assign("obtainEnvironmentPath", replacement, envir = ns)
    lockBinding("obtainEnvironmentPath", ns)
    invisible(TRUE)
}

.onLoad <- function(libname, pkgname) {
    sclet_install_env_shim()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(yulab.utils::yulab_msg(pkgname))
  packageStartupMessage(
    "AI-assisted analysis: scletAI(object, \"cluster the cells\") ",
    "(needs the optional 'aisdk' package and a provider key); ",
    "?scletAI, ?RunAIAnalysis, ?AskAI"
  )
}
