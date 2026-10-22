test_that("get_annotation answers only about the annotation", {
    # Six analyses registered a record under `type = "annotation"`: SingleR
    # annotation, the CNV burden (two paths), the four-state malignant call,
    # the iterative malignant classification, and the phenotype call. The
    # bucket has ONE "active" slot, so the last one written took it, and
    # `get_annotation()` -- the question "what is this object annotated as?" --
    # returned a CNV run. Measured 2026-10-22:
    #   annotation -> SingleR; then CNV -> "inferCNV per-cell copy-number burden".
    #
    # Each analysis now records under its own type: cnv / malignancy /
    # phenotype. The annotation bucket keeps the annotation.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_analysis_state(
        object = sce, type = "annotation", id = "singler_run",
        method = "SingleR", summary = list(n_labels = 3)
    )
    expect_equal(get_annotation(sce)$method, "SingleR")

    sce <- sclet_set_state_record(
        object = sce, type = "cnv", id = "cnv1",
        value = list(method = "inferCNV per-cell copy-number burden")
    )
    expect_equal(get_annotation(sce)$method, "SingleR")

    sce <- sclet_set_state_record(
        object = sce, type = "malignancy", id = "malignant",
        value = list(method = "four-state malignant-cell call from CNV + signature evidence")
    )
    expect_equal(get_annotation(sce)$method, "SingleR")

    sce <- sclet_set_state_record(
        object = sce, type = "phenotype", id = "phenotype_call",
        value = list(method = "phenotype call from independent lines of evidence")
    )
    expect_equal(get_annotation(sce)$method, "SingleR")
})

test_that("a malignant call without an annotation does not masquerade as one", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_state_record(
        object = sce, type = "malignancy", id = "malignant",
        value = list(method = "four-state malignant-cell call")
    )
    expect_null(get_annotation(sce))
    expect_equal(sclet_get_state_record(sce, "malignancy", "malignant")$method,
                 "four-state malignant-cell call")
})

test_that("the new state types are in the vocabulary", {
    # Without these, `sclet_set_state_record()` refuses the type and the
    # analyses would have had to stay in the annotation bucket.
    types <- sclet_state_types()
    expect_true(all(c("cnv", "malignancy", "phenotype") %in% types))
    expect_true("annotation" %in% types)
})

test_that("each analysis's record keeps its own summary shape", {
    # The reason the records could not simply be dropped: their `summary` is
    # what a reader queries, and the `analyses` slot holds the flat result
    # instead. This pins that the summary survives the move.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_state_record(
        object = sce, type = "phenotype", id = "phenotype_call",
        value = list(
            method = "phenotype call from independent lines of evidence",
            summary = list(independent_providers = c("yes_a", "yes_b"),
                           independent_lines = 1L)
        )
    )
    rec <- sclet_get_state_record(sce, "phenotype", "phenotype_call")
    expect_equal(rec$summary$independent_lines, 1L)
    expect_setequal(rec$summary$independent_providers, c("yes_a", "yes_b"))
})

test_that("get_malignancy reads the analyses slot the markers also write", {
    # `ClassifyMalignantEpithelium()` writes `metadata$malignancy[[name]]` and
    # a state record; `CallMalignantCells()` writes the state record. The
    # accessor reads the analyses slot, which is where the older path put its
    # result -- so both are queryable, by their own channel.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- sclet_set_analysis(sce, "malignancy", list(
        run_a = list(id = "run_a", method = "call-a"),
        run_b = list(id = "run_b", method = "call-b")
    ))
    expect_equal(get_malignancy(sce, id = "run_a")$method, "call-a")
    expect_equal(get_malignancy(sce)$method, "call-b")
    expect_null(get_malignancy(sce, id = "run_c"))
    expect_null(get_malignancy(scuttle::mockSCE(ncells = 20, ngenes = 50)))
})
