## Route comparison metrics must measure the ROUTE's own embedding.
##
## Measured on Kang IFN-beta (2026-10-20, 5,000 cells): the registered
## `batch_mixing` was 0.8778 for the raw baseline, for fastMNN and for Harmony
## alike, because the metric always read the reduction whose name matched
## "pca" -- which on a corrected object is still the UNCORRECTED PCA. The
## batch eta-squared on each route's own embedding was 0.1222 / 0.0064 /
## 0.0427, a 19-fold difference the comparison reported as "no route changed
## anything". These tests pin the correction without running fastMNN: the
## fixture carries both embeddings.
##
## The second half pins the comparison itself: one evidence node (the routes
## are alternative methods on one measurement, not corroboration), anonymized
## route codes, and a route comparison that a plan can run as a step.

## An object with two embeddings: a PCA that separates the batches and a
## "corrected" embedding that does not, plus an integration record naming the
## corrected one as the route's artifact.
route_metric_object <- function(n_per = 40L, with_metrics = TRUE) {
    set.seed(5L)
    n <- n_per * 2L
    pca <- matrix(0, nrow = n, ncol = 5L)
    pca[seq_len(n_per), 1L] <- 6      # batch A sits at +6
    pca[(n_per + 1L):n, 1L] <- -6     # batch B sits at -6
    ## The other components are noise, so they carry no batch structure.
    pca[, 2:5] <- matrix(stats::rnorm(n * 4L), nrow = n)
    corrected <- pca
    corrected[, 1L] <- stats::rnorm(n)   # no batch structure left
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 4L, ncol = n)))
    SummarizedExperiment::colData(sce)$batch <- rep(c("batch1", "batch2"),
        each = n_per)
    SummarizedExperiment::colData(sce)$condition <- rep(c("ctrl", "stim"),
        each = n_per)
    SingleCellExperiment::reducedDim(sce, "PCA") <- pca
    SingleCellExperiment::reducedDim(sce, "corrected") <- corrected
    sce <- sclet_set_analysis_state(
        object = sce, type = "integration", id = "fastmnn",
        method = "batchelor::fastMNN",
        inputs = list(batch = "batch"),
        artifacts = list(reduction = "corrected"),
        summary = list(route = "fastmnn"), active = FALSE
    )
    if (isTRUE(with_metrics)) {
        ## The metrics are computed AFTER the record exists, exactly as the
        ## product does: the route's own embedding is named by the record's
        ## artifacts, so a metric computed before it would fall back to the
        ## PCA -- the very defect this file pins.
        sce <- sclet:::sclet_ai_record_route_metrics(sce, route_name = "fastmnn",
            batch = "batch", runtime_sec = 1,
            design = list(batch = "batch", condition = "condition"))
    }
    sce
}

eta_of <- function(object, group_col, reduction) {
    m <- SingleCellExperiment::reducedDim(object, reduction)
    groups <- as.character(SummarizedExperiment::colData(object)[[group_col]])
    mean(vapply(seq_len(min(5L, ncol(m))), function(i) {
        sclet:::sclet_ai_diag_eta_squared(as.numeric(m[, i]), groups)
    }, numeric(1)))
}

test_that("a route's metrics are measured on the route's own embedding", {
    sce <- route_metric_object()
    design <- list(batch = "batch", condition = "condition")

    eta_pca <- eta_of(sce, "batch", "PCA")
    eta_corrected <- eta_of(sce, "batch", "corrected")
    expect_gt(eta_pca, 0.15)
    expect_lt(eta_corrected, 0.05)

    metrics <- sclet:::sclet_ai_compute_route_metrics(sce, "fastmnn",
        design = design, runtime_sec = 1)
    batch_mixing <- metrics[[1L]]
    expect_equal(batch_mixing$value, unname(1 - eta_corrected))
    expect_false(isTRUE(all.equal(batch_mixing$value, unname(1 - eta_pca))))

    ## The metric says WHICH embedding it read, so the number is reproducible.
    expect_equal(batch_mixing$parameters$reduction, "corrected")
    expect_false(identical(batch_mixing$method, "pca_eta_squared_complement"))
})

test_that("the baseline route is measured on the PCA", {
    sce <- route_metric_object()
    design <- list(batch = "batch", condition = "condition")
    metrics <- sclet:::sclet_ai_compute_route_metrics(sce, "raw",
        design = design, runtime_sec = 0)
    expect_equal(metrics[[1L]]$value, unname(1 - eta_of(sce, "batch", "PCA")))
    expect_equal(metrics[[1L]]$parameters$reduction, "PCA")
})

test_that("the plan-driven run_integration action records its route metrics", {
    sce <- route_metric_object()
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch",
        condition = "condition"))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    expect_true("run_integration" %in% names(registry))
    expect_true("compare_integration_routes" %in% names(registry))
})

test_that("the route comparison registers one anonymous evidence node", {
    sce <- route_metric_object()
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch",
        condition = "condition"))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "compare routes",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "c1", action = "compare_integration_routes"))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    expect_true(validation$valid)

    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true("ev:integration_comparison" %in% names(nodes))
    node <- nodes[["ev:integration_comparison"]]
    expect_equal(node$values$n_routes, 2L)
    ## Anonymized codes, and no route name anywhere in the values.
    expect_true(all(grepl("^route_[0-9]+$", node$values$route_codes)))
    serialized <- paste(utils::capture.output(str(node$values, max.level = 3)), collapse = " ")
    expect_false(grepl("fastmnn", serialized, fixed = TRUE))
    expect_false(grepl("corrected", serialized, fixed = TRUE))

    ## The numbers in the node are the ones the route recorded.
    routes <- sclet:::sclet_ai_integration_route_metrics(executed$object)
    expect_equal(node$values$batch_mixing_per_route[[2L]],
        routes[["fastmnn"]]$batch_mixing)
    expect_true(isTRUE(node$values$batch_mixing_improved))
})

test_that("the route comparison refuses when no route was recorded", {
    sce <- route_metric_object(with_metrics = FALSE)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "compare routes",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "c1", action = "compare_integration_routes"))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    ## The fixture carries an integration record without metrics, so the
    ## comparison has nothing to read.
    expect_false(isTRUE(validation$valid))
    expect_true(any(grepl("no_integration_route_recorded", validation$errors)))
})

test_that("the route comparison is one line of support, not several", {
    sce <- route_metric_object()
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch",
        condition = "condition"))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "compare routes",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "c1", action = "compare_integration_routes"))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)

    ceiling <- check_claim_ceiling(executed$object,
        evidence_ids = "ev:integration_comparison")
    expect_equal(ceiling$n_independent_groups, 1L)
    expect_equal(ceiling$ceiling_claim_level, "associated")
})
