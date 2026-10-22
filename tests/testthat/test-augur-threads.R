test_that("sclet_augur_threads stays within what parallel will allow", {
    # Augur hands its `n_threads` straight to mclapply(mc.cores = ...), and
    # `parallel` refuses to spawn more than two processes when
    # _R_CHECK_LIMIT_CORES_ is set -- which R CMD check does. Augur's own
    # default is 4, so sclet has to ask for fewer or the package's test
    # suite fails inside R CMD check while passing under devtools::test().
    # Measured 2026-10-22.
    local({
        old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
        on.exit({
            if (is.na(old)) {
                Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
            } else {
                Sys.setenv("_R_CHECK_LIMIT_CORES_" = old)
            }
        }, add = TRUE)
        Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
        expect_lte(sclet_augur_threads(), 2L)
        expect_equal(sclet_augur_threads(1L), 1L)
    })
    local({
        old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
        on.exit({
            if (is.na(old)) {
                Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
            } else {
                Sys.setenv("_R_CHECK_LIMIT_CORES_" = old)
            }
        }, add = TRUE)
        Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
        expect_gte(sclet_augur_threads(), 1L)
        expect_lte(sclet_augur_threads(), 4L)
        expect_equal(sclet_augur_threads(3L), 3L)
    })
})

test_that("the threads helper never asks for zero or a negative count", {
    old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
    on.exit({
        if (is.na(old)) {
            Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
        } else {
            Sys.setenv("_R_CHECK_LIMIT_CORES_" = old)
        }
    }, add = TRUE)
    Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    expect_gte(sclet_augur_threads(), 1L)
})
