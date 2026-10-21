## The Python environment path must not follow the version of the R package.
##
## `basilisk::obtainEnvironmentPath()` builds
## `<external>/<basilisk-version>/<pkgname>/<packageVersion(pkgname)>/<envname>`,
## so bumping DESCRIPTION's Version used to move every environment: the next
## call reinstalled a multi-GB stack, and `dir.expiry` eventually deleted the
## previous one. sclet keys the directory on the environment *specification*
## instead. These tests assert the property directly, without ever building or
## importing an environment -- they must stay cheap and offline.

test_that("an environment path carries no package version", {
    skip_if_not_installed("basilisk")
    for (envname in sclet_env_shim_names()) {
        path <- sclet_env_path(envname)
        expect_type(path, "character")
        version <- as.character(utils::packageVersion("sclet"))
        # The whole version string must be absent, and so must its major.minor
        # prefix -- a path keyed on "1.1" would break on the same schedule.
        expect_false(grepl(version, path, fixed = TRUE))
        expect_false(grepl(paste0("/", strsplit(version, ".", fixed = TRUE)[[1L]][1L],
            "\\.", strsplit(version, ".", fixed = TRUE)[[1L]][2L]), path))
        expect_true(grepl(envname, basename(path), fixed = TRUE))
    }
})

test_that("the path is a function of the environment specification", {
    skip_if_not_installed("basilisk")
    env <- sclet_scvi_env
    path <- sclet_env_path("sclet_scvi_env")

    # Same spec, same path. The environment may already exist on this machine
    # (built by an earlier run), so the assertion is on stability and on not
    # *creating* anything, not on absence: `create = FALSE` must leave the
    # filesystem untouched, which is checked by comparing a directory that
    # cannot be created by a lookup, the env root itself.
    root <- dirname(path)
    before <- if (dir.exists(root)) sort(list.files(root)) else character()
    expect_identical(sclet_env_path("sclet_scvi_env"), path)
    after <- if (dir.exists(root)) sort(list.files(root)) else character()
    expect_identical(before, after)

    # A changed pin produces a different directory. Build the altered
    # environment object by hand rather than by editing the package source.
    changed <- basilisk::BasiliskEnvironment(
        envname = env@envname,
        pkgname = env@pkgname,
        packages = c(env@packages, "different-spec-sentinel=0.0.1")
    )
    expect_false(identical(
        sclet_env_spec_fingerprint(env),
        sclet_env_spec_fingerprint(changed)
    ))
    # And that different spec resolves to a path outside the existing tree.
    expect_false(identical(
        sclet_env_path("sclet_scvi_env"),
        file.path(dirname(path), sclet_env_spec_fingerprint(changed))
    ))
})

test_that("the basilisk lookup is shimmed once and passes others through", {
    skip_if_not_installed("basilisk")
    expect_true(sclet_env_shim_installed())

    ns <- asNamespace("basilisk")
    lookup <- get("obtainEnvironmentPath", envir = ns)
    expect_true(isTRUE(attr(lookup, "sclet_shim")))

    # sclet's own environments resolve to the version-free path ...
    expect_identical(
        lookup(sclet_scvi_env),
        sclet_env_path("sclet_scvi_env")
    )

    # ... and installing the shim twice is a no-op rather than a double wrap.
    before <- get("obtainEnvironmentPath", envir = ns)
    sclet_install_env_shim()
    expect_identical(
        attr(get("obtainEnvironmentPath", envir = ns), "sclet_shim"),
        attr(before, "sclet_shim")
    )
})

test_that("the sandboxed python runner clears the inherited loader path", {
    skip_if_not_installed("basilisk")
    # Regression: a child python inherited R's LD_LIBRARY_PATH and resolved the
    # system libpython, failing with `_ctypes ... undefined symbol:
    # _PyErr_SetLocaleString`. `sclet_env_run()` must pass a cleared
    # environment. Asserted on the constructed argument rather than by running
    # python, so the test needs no interpreter.
    expect_true(is.function(sclet_env_run))
    body_text <- paste(deparse(body(sclet_env_run)), collapse = "\n")
    expect_match(body_text, "LD_LIBRARY_PATH=", fixed = TRUE)
})
