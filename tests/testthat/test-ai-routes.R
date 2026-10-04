test_that("RunIntegrationRoutes returns clarification_required when design is missing", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:80, nrow = 20, ncol = 4))
    )
    r <- RunIntegrationRoutes(sce)
    expect_equal(r$status, "clarification_required")
    expect_true(length(r$questions) > 0)
    expect_false(r$execution$performed)
})

test_that("RunIntegrationRoutes returns clarification_required when design$batch is NULL", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:80, nrow = 20, ncol = 4))
    )
    r <- RunIntegrationRoutes(sce, design = list(condition = "x"))
    expect_equal(r$status, "clarification_required")
    expect_false(r$execution$performed)
})

test_that("RunIntegrationRoutes returns clarification_required when batch is not a colData column", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:80, nrow = 20, ncol = 4))
    )
    r <- RunIntegrationRoutes(sce, design = list(batch = "not_a_column"))
    expect_equal(r$status, "clarification_required")
    expect_false(r$execution$performed)
})

test_that("RunIntegrationRoutes confirm=ask in non-interactive returns dry_run without modifying object", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:80, nrow = 20, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$sample_id <- rep(c("group_1", "group_2"), each = 2)
    before <- sce
    r <- RunIntegrationRoutes(sce, design = list(batch = "sample_id"),
        routes = c("raw", "fastMNN"), confirm = "ask")
    expect_equal(r$status, "dry_run")
    expect_false(r$execution$performed)
    expect_equal(r$execution$n_routes, 0L)
    before_state <- sclet_get_state(before)
    after_state <- sclet_get_state(r$object)
    expect_identical(before_state$states$records, after_state$states$records)
})

test_that("RunIntegrationRoutes confirm=never returns cancelled", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:80, nrow = 20, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$sample_id <- rep(c("group_1", "group_2"), each = 2)
    r <- RunIntegrationRoutes(sce, design = list(batch = "sample_id"), confirm = "never")
    expect_equal(r$status, "cancelled")
    expect_false(r$execution$performed)
})

test_that("RunIntegrationRoutes confirm=always with raw baseline registers metrics", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(80, 10), nrow = 20, ncol = 4))
    )
    rownames(sce) <- paste0("g", seq_len(20))
    colnames(sce) <- paste0("c", seq_len(4))
    SummarizedExperiment::colData(sce)$sample_id <- rep(c("group_1", "group_2"), each = 2)
    SummarizedExperiment::colData(sce)$condition <- rep(c("group_1", "group_2"), each = 2)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 10)
    sce <- RunPCA(sce, ncomponents = 2)

    r <- RunIntegrationRoutes(sce, design = list(batch = "sample_id", condition = "condition"),
        routes = "raw", confirm = "always")
    expect_equal(r$status, "completed")
    expect_true(r$execution$performed)
    expect_identical(r$report$recommendation, NULL)

    raw_state <- sclet_get_state_record(r$object, "integration", "raw")
    expect_false(is.null(raw_state))
    metrics <- raw_state$summary$metrics
    expect_gte(length(metrics), 4L)
    metric_names <- vapply(metrics, function(m) m$name, character(1))
    expect_true(all(c("batch_mixing", "biological_preservation", "cluster_stability", "runtime_sec") %in% metric_names))
    stab <- Filter(function(m) identical(m$name, "cluster_stability"), metrics)[[1]]
    expect_identical(stab$uncertainty$status, "not_available")
    expect_true(grepl("stability_resampling_skipped", stab$uncertainty$reason))
})

test_that("RunIntegrationRoutes with mocked fastMNN preserves original counts", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(80, 10), nrow = 20, ncol = 4))
    )
    rownames(sce) <- paste0("g", seq_len(20))
    colnames(sce) <- paste0("c", seq_len(4))
    SummarizedExperiment::colData(sce)$sample_id <- rep(c("group_1", "group_2"), each = 2)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 10)
    sce <- RunPCA(sce, ncomponents = 2)
    original_counts <- SummarizedExperiment::assay(sce, "counts")

    fake_integration <- function(object, ...) {
        args <- list(...)
        id <- args$name %||% "fastmnn"
        SingleCellExperiment::reducedDim(object, "correctedPCA") <- matrix(
            seq_len(ncol(object) * 2L), ncol = 2)
        object <- sclet_set_analysis_state(object, "integration", id,
            method = "mocked_integration", inputs = list(batch = args$batch),
            summary = list(), artifacts = list(), active = FALSE)
        object
    }
    r <- testthat::with_mocked_bindings({
        RunIntegrationRoutes(sce, design = list(batch = "sample_id"),
            routes = c("raw", "fastMNN"), confirm = "always")
    }, .package = "sclet", RunIntegration = fake_integration)

    expect_equal(r$status, "completed")
    expect_true(r$execution$performed)
    final_counts <- SummarizedExperiment::assay(r$object, "counts")
    expect_identical(as.matrix(original_counts), as.matrix(final_counts))
    recs <- sclet_get_state(r$object)$states$records$integration
    expect_true("raw" %in% names(recs))
    expect_true("fastmnn" %in% names(recs))
})

test_that("RunIntegrationRoutes metrics have at least name/value/direction or not_available", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(80, 10), nrow = 20, ncol = 4))
    )
    rownames(sce) <- paste0("g", seq_len(20))
    colnames(sce) <- paste0("c", seq_len(4))
    SummarizedExperiment::colData(sce)$sample_id <- rep(c("group_1", "group_2"), each = 2)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 10)
    sce <- RunPCA(sce, ncomponents = 2)
    r <- RunIntegrationRoutes(sce, design = list(batch = "sample_id"),
        routes = "raw", confirm = "always")
    rec <- sclet_get_state_record(r$object, "integration", "raw")
    for (m in rec$summary$metrics) {
        expect_true("name" %in% names(m))
        expect_true("direction" %in% names(m))
        has_value <- "value" %in% names(m)
        has_status <- !is.null(m$status) || !is.null(m$uncertainty$status)
        expect_true(has_value || has_status)
    }
})

test_that("no automatic best route is produced by metrics conflict fixture", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(80, 10), nrow = 20, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$sample_id <- rep(c("group_1", "group_2"), each = 2)
    # Register two fake integration routes with conflicting metrics
    sce <- sclet_set_analysis_state(sce, "integration", "route_A", "mocked",
        inputs = list(batch = "sample_id"),
        summary = list(
            route = "A",
            metrics = list(
                list(name = "batch_mixing", value = 0.9, direction = "higher_is_better"),
                list(name = "biological_preservation", value = 0.4, direction = "higher_is_better")
            )
        ), active = FALSE)
    sce <- sclet_set_analysis_state(sce, "integration", "route_B", "mocked",
        inputs = list(batch = "sample_id"),
        summary = list(
            route = "B",
            metrics = list(
                list(name = "batch_mixing", value = 0.6, direction = "higher_is_better"),
                list(name = "biological_preservation", value = 0.85, direction = "higher_is_better")
            )
        ), active = FALSE)
    sce <- sclet_set_analysis(sce, "route_A", list(id = "route_A", type = "integration", status = "completed"))
    sce <- sclet_set_analysis(sce, "route_B", list(id = "route_B", type = "integration", status = "completed"))
    r <- RunIntegrationRoutes(sce, design = list(batch = "sample_id"),
        routes = c("raw"), confirm = "never")
    # The fixture just ensures the conflict routes exist; CompareAIAnalyses handles tradeoff detection in T3.
    # For T2, we only assert RunIntegrationRoutes never sets a $report$best_route nor a non-NULL recommendation.
    expect_null(r$report$best_route)
    expect_null(r$report$recommendation)
})
