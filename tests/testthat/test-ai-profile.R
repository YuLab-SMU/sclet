test_that("GetAIProfile reports missing design metadata without values", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 3, ncol = 2))
    )

    profile <- sclet:::GetAIProfile(sce)

    expect_equal(
        names(profile)[seq_len(10L)],
        c(
            "schema_version", "dataset", "design", "qc", "structure",
            "analysis_state", "diagnostics", "capabilities", "cost_estimates",
            "privacy"
        )
    )
    expect_equal(profile$schema_version, "1.0")
    expect_equal(profile$dataset$n_cells, 2L)
    expect_equal(profile$dataset$n_genes, 3L)
    expect_equal(profile$design$metadata$status, "missing")
    expect_equal(profile$design$sample$status, "missing")
    expect_false(profile$privacy$complete_matrix_included)
    expect_false(profile$privacy$metadata_values_included)
})

test_that("GetAIProfile exposes metadata design candidates as unknown", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(seq_len(12), nrow = 3, ncol = 4)),
        colData = S4Vectors::DataFrame(
            sample_id = c("sample-a", "sample-a", "sample-b", "sample-b"),
            batch = c("run-1", "run-1", "run-2", "run-2"),
            condition = c("control", "control", "treated", "treated"),
            subject_id = c("subject-a", "subject-a", "subject-b", "subject-b"),
            secret_value = rep("do-not-return", 4)
        )
    )

    profile <- sclet:::GetAIProfile(sce)

    expect_equal(profile$design$sample$status, "candidate")
    expect_equal(profile$design$sample$semantic, "unknown")
    expect_equal(profile$design$sample$column, "sample_id")
    expect_equal(profile$design$batch$column, "batch")
    expect_equal(profile$design$condition$column, "condition")
    expect_equal(profile$design$subject$column, "subject_id")
    expect_true(profile$design$sample$confirmation_required)
    expect_equal(profile$diagnostics$batch$status, "candidate")
    expect_equal(profile$diagnostics$batch$semantic, "unknown")
})

test_that("GetAIProfile fingerprint is stable across repeated calls", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 5, ncol = 3))
    )

    first <- sclet:::GetAIProfile(sce)
    second <- sclet:::GetAIProfile(sce)

    expect_match(first$fingerprint, "^sclet-ai-profile-[0-9a-f]{8}$")
    expect_identical(first$fingerprint, second$fingerprint)
    expect_identical(first, second)
})

test_that("sparse assay summaries stay within the configured bound", {
    sparse_counts <- Matrix::sparseMatrix(
        i = c(1L, 2L, 3L),
        j = c(1L, 2L, 3L),
        x = c(2, 4, 6),
        dims = c(3L, 3L)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = sparse_counts)
    )

    available <- sclet:::GetAIProfile(sce, max_sparse_elements = 20L)
    skipped <- sclet:::GetAIProfile(sce, max_sparse_elements = 3L)

    expect_true(available$dataset$assays$counts$sparse)
    expect_equal(available$dataset$assays$counts$bounded_summary$status, "available")
    expect_equal(available$dataset$assays$counts$bounded_summary$nnzero, 3)
    expect_equal(skipped$dataset$assays$counts$bounded_summary$status, "skipped")
    expect_equal(skipped$dataset$assays$counts$bounded_summary$reason, "size_limit")
    expect_false(available$dataset$assays$counts$values_included)
})

test_that("profile output does not contain matrix or metadata values", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(c(101, 202, 303, 404), nrow = 2)),
        colData = S4Vectors::DataFrame(
            sample_id = c("private-cell-a", "private-cell-b"),
            note = c("private-note-a", "private-note-b")
        )
    )

    profile <- sclet:::GetAIProfile(sce)
    serialized_profile <- paste(capture.output(str(profile)), collapse = " ")

    expect_false(grepl("private-cell|private-note|101|202|303|404", serialized_profile))
    expect_false(profile$privacy$assay_values_included)
    expect_false(profile$privacy$metadata_values_included)
    expect_false(profile$capabilities$complete_matrix_transfer)
})
