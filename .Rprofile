# Render artefacts stay inside the ignored render directory.
#
# `opts_chunk$set(cache = TRUE)` in index.Rmd makes knitr write one cache
# directory per chunk, named `<label>_<hash>`. The directory is resolved
# against the working directory *at knit time*, and the default is
# `sclet_cache/` at the book root -- which is also the directory a reader
# does `setwd()` into before knitting a single chapter. A local `make gitbook`
# left ~1 GB of them there (measured 2026-10-21: 992 MB, 493 files, oldest
# 2026-08-30). Pointing cache.path at the render output keeps the book root
# clean and lets `make clean` remove the whole thing in one move.
#
# CI does not depend on this: every run is a fresh runner, so the cache is
# never reused across builds, and the workflow renders to an explicit
# `output_dir` (`$GITHUB_WORKSPACE/site`) that this path follows.
#
# IMPORTANT -- this file does not run the user's own ~/.Rprofile.
#
# R reads exactly one user profile: the first `.Rprofile` found by walking up
# from the working directory. Adding this file to `bookdown/` therefore SHADOWS
# `~/.Rprofile` for every R process started in this directory, including a
# plain `Rscript -e 'library(sclet)'` while writing a chapter. On the machine
# where this file was added, the user's library lives in `~/Rlib` and is put on
# `.libPaths()` by `~/.Rprofile`; without it `library(sclet)` fails with
# "there is no package called 'sclet'" -- an error that looks like a broken
# package and is actually a shadowed profile. `Rscript` does honour
# `R_PROFILE_USER`, so R has an escape hatch, but a `setwd()`-ing researcher
# does not.
#
# So the body is conditional and does nothing unless `R_BOOKDOWN_RENDER` is
# set, which only the Makefile's render targets do. A reader who lands here
# still loses `~/.Rprofile` (that part cannot be avoided while this file
# exists), so it re-sources the home profile when it is not already running as
# part of a render. That keeps `library(sclet)` working from this directory.
local({
    home_profile <- path.expand("~/.Rprofile")
    rendering <- nzchar(Sys.getenv("R_BOOKDOWN_RENDER"))

    # Re-run the user's profile first, in BOTH modes: their library paths (and
    # their options) must be in place before `knitr` can even be loaded, and a
    # render needs them for exactly the same reason a bare `Rscript` does.
    if (identical(Sys.getenv("R_PROFILE_USER", unset = ""), "") &&
        file.exists(home_profile)) {
        try(sys.source(home_profile, envir = globalenv()), silent = TRUE)
    }

    if (rendering) {
        # Redirect the per-chunk knitr cache out of the book root.
        try(knitr::opts_chunk$set(cache.path = "docs/sclet_cache/"), silent = TRUE)
    }
})
