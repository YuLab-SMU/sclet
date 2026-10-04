test_that("CompareAIAnalyses reports no routes without recorded integration analyses", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 2L)))
    result <- sclet:::CompareAIAnalyses(sce)
    expect_equal(result$status, "not_available")
    expect_false(result$execution$performed)
})

test_that("CompareAIAnalyses compares recorded routes without executing", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(
        sce, type = "integration", id = "raw_baseline", method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 1, direction = "higher_is_better",
                uncertainty = list(status = "reported")),
            list(name = "runtime_sec", value = 1, direction = "lower_is_better",
                uncertainty = list(status = "reported"))
        ))
    )
    sce <- sclet:::sclet_set_analysis_state(
        sce, type = "integration", id = "harmony_1", method = "harmony",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.8, direction = "higher_is_better"),
            list(name = "biological_preservation", value = 0.9, direction = "higher_is_better"),
            list(name = "runtime_sec", value = 2, direction = "lower_is_better")
        ))
    )
    result <- sclet:::CompareAIAnalyses(sce, criterion = "biological_preservation")
    expect_equal(result$status, "available")
    expect_equal(result$baseline, "raw_baseline")
    expect_length(result$routes, 2L)
    expect_true(length(result$metrics) >= 3L)
    expect_null(result$recommendation)
    expect_false(result$execution$performed)
    expect_false(grepl("counts|1 1 1", paste(utils::capture.output(str(result)), collapse = " ")))
})

test_that("CompareAIAnalyses orders routes by criterion", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "raw_baseline",
        method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 0.5, direction = "higher_is_better"),
            list(name = "batch_mixing", value = 0.3, direction = "higher_is_better")
        )))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "route_A",
        method = "mocked_A",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 0.95, direction = "higher_is_better"),
            list(name = "batch_mixing", value = 0.2, direction = "higher_is_better")
        )))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "route_B",
        method = "mocked_B",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 0.7, direction = "higher_is_better"),
            list(name = "batch_mixing", value = 0.4, direction = "higher_is_better")
        )))
    r <- CompareAIAnalyses(sce, criterion = "biological_preservation")
    ord <- names(r$routes)
    # raw_baseline can be first or last; the others should have route_A before route_B
    non_baseline <- setdiff(ord, "raw_baseline")
    expect_equal(which(non_baseline == "route_A"), 1L)
})

test_that("CompareAIAnalyses detects batch_vs_bio tradeoff and refuses recommendation", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "raw_baseline",
        method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.4, direction = "higher_is_better"),
            list(name = "biological_preservation", value = 0.8, direction = "higher_is_better")
        )))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "harmony_conflict",
        method = "harmony",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.85, direction = "higher_is_better"),
            list(name = "biological_preservation", value = 0.6, direction = "higher_is_better")
        )))
    r <- CompareAIAnalyses(sce, criterion = "biological_preservation")
    tradeoffs <- r$tradeoffs
    expect_true(is.list(tradeoffs) && (length(tradeoffs) > 0L || !is.null(tradeoffs$status)))
    has_conflict <- any(vapply(tradeoffs, function(t) identical(t$type, "batch_vs_bio"), logical(1L)))
    expect_true(has_conflict)
    expect_null(r$recommendation)
    expect_false("best_route" %in% names(r))
})

test_that("CompareAIAnalyses normalizes missing metrics to not_available", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "raw_baseline",
        method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.3, direction = "higher_is_better")
        )))
    r <- CompareAIAnalyses(sce, criterion = "biological_preservation")
    baseline_metrics_key <- grep("^raw_baseline::", names(r$metrics), value = TRUE)
    metric_names <- sub("^raw_baseline::", "", baseline_metrics_key)
    expect_true("cluster_stability" %in% metric_names)
    stab_key <- "raw_baseline::cluster_stability"
    m <- r$metrics[[stab_key]]
    expect_identical(m$uncertainty$status, "not_available")
    expect_true(is.na(m$value))
})
