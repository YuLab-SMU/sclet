test_that("analysis_story is empty on an SCE with zero recorded analyses", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L))
    )
    ledger <- GetAnalysisLedger(sce)
    story <- ledger$analysis_story
    expect_equal(story$timeline, list())
    expect_equal(story$n_steps, 0L)
    expect_equal(story$unordered_steps, character())
    expect_equal(story$user_decisions, list())
    expect_equal(story$design_confirmations, list())
    expect_equal(story$evidence_gaps, list())
    expect_equal(story$conflicts, list())
    expect_false(story$raw_values_included)
})

test_that("analysis_story timeline is ordered by created_at and n_steps matches", {
    set.seed(11L)
    counts <- matrix(rpois(50L * 40L, lambda = 5L), nrow = 50L, ncol = 40L,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("cell_", seq_len(40L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- RunPCA(sce, ncomponents = 5L)
    sce <- FindNeighbors(sce, dims = seq_len(3L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 0.8)
    ledger <- GetAnalysisLedger(sce)
    story <- ledger$analysis_story
    expect_true(story$n_steps >= 4L)
    expect_equal(story$n_steps, length(story$timeline))
    times <- vapply(story$timeline, function(x) {
        if (is.null(x$created_at)) NA_real_ else as.numeric(as.POSIXct(x$created_at))
    }, numeric(1L))
    times_with <- times[!is.na(times)]
    if (length(times_with) > 1L) {
        expect_true(all(diff(times_with) >= 0))
    }
})

test_that("analysis_story conflicts reports multiple rare_cells records", {
    skip_if_not_installed("BiocNeighbors")
    set.seed(11L)
    counts <- matrix(rpois(50L * 60L, lambda = 5L), nrow = 50L, ncol = 60L,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("cell_", seq_len(60L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    SingleCellExperiment::reducedDim(sce, "PCA") <- SingleCellExperiment::reducedDim(sce, "PCA")
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    plan1 <- sclet:::new_sclet_ai_plan(
        task = "rare_cell_review",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "rare1", action = "run_rare_cell_detection",
            params = list(name = "rare_a", rare_threshold = 1000L, k = 5L, dims = c(1L, 2L, 3L))))
    )
    val1 <- ValidateAIPlan(plan1, object = sce, registry = registry)
    result1 <- ExecuteAIPlan(sce, plan1, registry, validation = val1,
        dry_run = FALSE, confirmation = val1$confirmation_token)
    sce <- result1$object
    plan2 <- sclet:::new_sclet_ai_plan(
        task = "rare_cell_review",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "rare2", action = "run_rare_cell_detection",
            params = list(name = "rare_b", rare_threshold = 1000L, k = 5L, dims = c(1L, 2L, 3L))))
    )
    val2 <- ValidateAIPlan(plan2, object = sce, registry = registry)
    result2 <- ExecuteAIPlan(sce, plan2, registry, validation = val2,
        dry_run = FALSE, confirmation = val2$confirmation_token)
    sce <- result2$object
    ledger <- GetAnalysisLedger(sce)
    story <- ledger$analysis_story
    rare_conflicts <- Filter(function(x) identical(x$type, "rare_cells"), story$conflicts)
    expect_true(length(rare_conflicts) >= 1L)
    expect_true(length(rare_conflicts[[1L]]$keys) >= 2L)
})

test_that("analysis_story unordered_steps captures records without created_at", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L))
    )
    sce <- sclet:::sclet_set_analysis_state(
        sce,
        type = "preprocess",
        id = "no_time_record",
        method = "test_method",
        summary = list(value = "test")
    )
    state <- sclet:::sclet_get_state(sce)
    records_key <- names(state$states$records$preprocess)[1L]
    state$states$records$preprocess[[records_key]]$created_at <- NULL
    state$states$records$preprocess[[records_key]]$summary$created_at <- NULL
    sce <- sclet:::sclet_set_state(sce, state)
    ledger <- GetAnalysisLedger(sce)
    story <- ledger$analysis_story
    expect_true("no_time_record" %in% story$unordered_steps ||
        any(grepl("no_time_record", story$unordered_steps)))
})

test_that("analysis_story user_decisions captures user_decision evidence nodes", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    updated <- sclet:::sclet_ai_record_clarification_response(
        sce,
        question_id = "design_batch",
        answer = "batch",
        blocked_action = "run_integration"
    )
    ledger <- GetAnalysisLedger(updated)
    story <- ledger$analysis_story
    expect_true(length(story$user_decisions) >= 1L)
    decision <- story$user_decisions[[1L]]
    expect_true(is.character(decision$id))
    expect_true(is.character(decision$summary) || is.null(decision$summary))
    story_str <- utils::capture.output(str(story))
    expect_false(any(grepl("donor", story_str, ignore.case = TRUE)))
})

test_that("analysis_story design_confirmations includes design but not sensitive fields", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(
            batch_col = rep(c("donor_1", "donor_2"), each = 4L),
            condition = rep(c("ctrl", "stim"), 4L)
        )
    )
    confirmed <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch_col"))
    ledger <- GetAnalysisLedger(confirmed)
    story <- ledger$analysis_story
    expect_true(length(story$design_confirmations) >= 1L)
    dc <- story$design_confirmations[[1L]]
    expect_equal(dc$design, list(batch = "batch_col"))
    expect_true(is.character(dc$id))
    expect_false("design_value_key" %in% names(dc))
    expect_false("object_fingerprint" %in% names(dc))
    story_str <- utils::capture.output(str(story))
    expect_false(any(grepl("donor_1", story_str, fixed = TRUE)))
    expect_false(any(grepl("donor_2", story_str, fixed = TRUE)))
})

test_that("analysis_story evidence_gaps reports missing doublet evidence when rare_cells present", {
    skip_if_not_installed("BiocNeighbors")
    set.seed(11L)
    counts <- matrix(rpois(50L * 60L, lambda = 5L), nrow = 50L, ncol = 60L,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("cell_", seq_len(60L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "rare_cell_review",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "rare1", action = "run_rare_cell_detection",
            params = list(name = "rare_test", rare_threshold = 1000L, k = 5L, dims = c(1L, 2L, 3L))))
    )
    val <- ValidateAIPlan(plan, object = sce, registry = registry)
    result <- ExecuteAIPlan(sce, plan, registry, validation = val,
        dry_run = FALSE, confirmation = val$confirmation_token)
    sce <- result$object
    ledger <- GetAnalysisLedger(sce)
    story <- ledger$analysis_story
    expect_true(length(story$evidence_gaps) >= 1L)
    gap <- story$evidence_gaps[[1L]]
    expect_equal(gap$analysis_type, "rare_cells")
    expect_equal(gap$gap, "doublet_evidence_available")
    expect_true(is.character(gap$reason))
    expect_true(nzchar(gap$reason))
})

test_that("sclet_ai_blocked_actions returns NA for execute_analysis$blocked", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L))
    )
    ledger <- GetAnalysisLedger(sce)
    blocked <- ledger$blocked_actions$execute_analysis$blocked
    expect_true(is.na(blocked))
    expect_false(identical(blocked, TRUE))
    expect_false(identical(blocked, FALSE))
})

test_that("sclet_ai_capabilities returns NA for execution", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L))
    )
    ledger <- GetAnalysisLedger(sce)
    exec <- ledger$capabilities$execution
    expect_true(is.na(exec))
    expect_false(identical(exec, TRUE))
    expect_false(identical(exec, FALSE))
})

test_that("analysis_story never contains raw expression values", {
    set.seed(11L)
    counts <- matrix(rpois(50L * 40L, lambda = 5L), nrow = 50L, ncol = 40L,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("cell_", seq_len(40L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- RunPCA(sce, ncomponents = 5L)
    ledger <- GetAnalysisLedger(sce)
    story <- ledger$analysis_story
    expect_false(story$raw_values_included)
    story_str <- utils::capture.output(str(story))
    expect_false(any(grepl("gene_1", story_str, fixed = TRUE)))
    expect_false(any(grepl("cell_1", story_str, fixed = TRUE)))
})

test_that("analysis_story participates in fingerprint changes", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L))
    )
    fp1 <- GetAnalysisLedger(sce)$fingerprint
    sce <- sclet:::sclet_set_analysis_state(
        sce,
        type = "preprocess",
        id = "new_record",
        method = "test_method",
        summary = list(value = "test")
    )
    fp2 <- GetAnalysisLedger(sce)$fingerprint
    expect_false(identical(fp1, fp2))
})

test_that("GetAnalysisLedger is read-only and does not mutate the object", {
    set.seed(11L)
    counts <- matrix(rpois(50L * 40L, lambda = 5L), nrow = 50L, ncol = 40L,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("cell_", seq_len(40L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- RunPCA(sce, ncomponents = 5L)
    sce_before <- sce
    ledger <- GetAnalysisLedger(sce)
    expect_identical(sce, sce_before)
})
