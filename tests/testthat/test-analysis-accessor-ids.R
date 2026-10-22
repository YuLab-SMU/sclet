test_that("a nonexistent id returns nothing, instead of another record", {
    # The legacy fallback used to fire whenever the requested record was
    # missing, so asking for an id that does not exist returned a DIFFERENT
    # record, and `has_<type>(sce, id = "typo")` reported TRUE. Measured
    # 2026-10-22.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_analysis(sce, "velocity", list(id = "velocity", method = "velociraptor::scvelo"))
    sce <- sclet_set_analysis(sce, "trajectory", list(id = "trajectory", method = "slingshot"))
    sce <- sclet_set_analysis(sce, "milo", list(id = "milo", method = "miloR"))

    expect_null(get_velocity(sce, id = "typo"))
    expect_false(has_velocity(sce, id = "typo"))
    expect_null(get_trajectory(sce, id = "typo"))
    expect_false(has_trajectory(sce, id = "typo"))
    expect_null(get_milo(sce, id = "typo"))
    expect_false(has_milo(sce, id = "typo"))
})

test_that("the legacy record is still reachable without an id", {
    # The fallback exists because older objects kept their single analysis
    # record in the `analyses` slot. Fixing the wrong-id case must not take
    # that away.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_analysis(sce, "velocity", list(id = "velocity", method = "velociraptor::scvelo"))
    expect_equal(get_velocity(sce)$method, "velociraptor::scvelo")
    expect_true(has_velocity(sce))
})

test_that("the legacy record is reachable by its own id", {
    # The legacy slot does not index by id, so a caller who read the id off the
    # record and passed it back must still get the record.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_analysis(sce, "velocity", list(id = "velocity", method = "velociraptor::scvelo"))
    expect_equal(get_velocity(sce, id = "velocity")$method, "velociraptor::scvelo")
})

test_that("a state record asked for by id is the one returned", {
    # With two records on the object, `id` has to select, and it has to select
    # the one named -- not the active one, and not the legacy one.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_analysis(sce, "velocity", list(id = "legacy", method = "legacy-method"))
    sce <- sclet_set_analysis_state(
        object = sce, type = "velocity", id = "run_a", method = "method-a",
        summary = list(which = "a")
    )
    sce <- sclet_set_analysis_state(
        object = sce, type = "velocity", id = "run_b", method = "method-b",
        summary = list(which = "b")
    )
    expect_equal(get_velocity(sce, id = "run_a")$method, "method-a")
    expect_equal(get_velocity(sce, id = "run_b")$method, "method-b")
    # The active one is the last written, so asking for the other one proves the
    # selection is by id and not by "whatever is active".
    expect_equal(get_velocity(sce)$method, "method-b")
    # And an id that exists in neither place is still nothing.
    expect_null(get_velocity(sce, id = "run_c"))
})

test_that("sclet_legacy_analysis returns NULL when the caller named another id", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_analysis(sce, "milo", list(id = "milo", method = "miloR"))
    expect_null(sclet_legacy_analysis(sce, "milo", "milo", "other", "other"))
    expect_false(is.null(sclet_legacy_analysis(sce, "milo", "milo", NULL, "milo")))
    expect_false(is.null(sclet_legacy_analysis(sce, "milo", "milo", "milo", "milo")))
    # No legacy record at all -> NULL regardless.
    bare <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    expect_null(sclet_legacy_analysis(bare, "milo", "milo", NULL, NULL))
})
