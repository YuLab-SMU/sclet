## Design-aware differential expression.
##
## Kang IFN-beta real-data run (2026-10-18) asked for "stimulated vs control
## within donors"; no action could express it and a cluster-marker table was
## returned instead. `run_de_test(group=, covariates=)` now contrasts a colData
## column, and adjusts for donor by pseudobulk + limma instead of claiming a
## paired analysis the cell-level test cannot provide.

paired_de_object <- function(ng = 60L, per = 20L, fold = 4) {
    set.seed(42)
    cells <- expand.grid(donor = c("d1", "d2"), condition = c("ctrl", "stim"),
        rep = seq_len(per))
    counts <- matrix(stats::rnbinom(ng * nrow(cells), mu = 5, size = 2),
        nrow = ng, ncol = nrow(cells))
    rownames(counts) <- paste0("g", seq_len(ng))
    colnames(counts) <- paste0("c", seq_len(ncol(counts)))
    planted <- paste0("g", 1:6)
    up <- cells$condition == "stim"
    counts[match(planted, rownames(counts)), up] <-
        counts[match(planted, rownames(counts)), up] * fold
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    SummarizedExperiment::colData(sce)$condition <- cells$condition
    SummarizedExperiment::colData(sce)$donor <- cells$donor
    sce
}

condition_de_plan <- function(...) {
    sclet:::new_sclet_ai_plan(
        task = "condition de",
        actions = list(list(id = "d1", action = "run_de_test", params = list(...)))
    )
}

test_that("run_de_test contrasts a condition column with donor adjustment", {
    skip_if_not_installed("limma")
    sce <- paired_de_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- condition_de_plan(group = "condition", ident.1 = "stim",
        ident.2 = "ctrl", covariates = "donor", name = "cond_de")
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    expect_true(validation$valid)

    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    rec <- sclet:::sclet_get_state_record(executed$object, "detest", "cond_de")
    expect_true(grepl("condition", rec$inputs$design, fixed = TRUE))
    expect_true(grepl("donor", rec$inputs$design, fixed = TRUE))
    expect_true(grepl("pseudobulk", rec$inputs$adjustment, ignore.case = TRUE))
    table <- rec$artifacts$result
    expect_true(all(c("gene", "cluster", "avg_log2FC", "pval", "padj",
        "pct.1", "pct.2") %in% colnames(table)))
    expect_true(sum(paste0("g", 1:6) %in% utils::head(table$gene, 12L)) >= 4L)

    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true("ev:detest_condition_stim_vs_ctrl" %in% names(nodes))
})

test_that("a condition contrast without both groups is refused", {
    sce <- paired_de_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- condition_de_plan(group = "condition", ident.1 = "stim")
    err <- tryCatch(ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        error = function(e) e)
    expect_s3_class(err, "sclet_ai_invalid_plan")
    expect_true(grepl("contrast_missing", conditionMessage(err), fixed = TRUE))
})

test_that("an unknown group level or covariate is refused", {
    sce <- paired_de_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    bad_level <- condition_de_plan(group = "condition", ident.1 = "treated",
        ident.2 = "ctrl")
    err <- tryCatch(ValidateAIPlan(bad_level, object = sce, registry = registry,
        strict = TRUE), error = function(e) e)
    expect_true(grepl("contrast_invalid", conditionMessage(err), fixed = TRUE))

    bad_cov <- condition_de_plan(group = "condition", ident.1 = "stim",
        ident.2 = "ctrl", covariates = "patient")
    err2 <- tryCatch(ValidateAIPlan(bad_cov, object = sce, registry = registry,
        strict = TRUE), error = function(e) e)
    expect_true(grepl("covariates_missing", conditionMessage(err2), fixed = TRUE))
})

test_that("donors measured in one condition only are excluded, not fatal", {
    skip_if_not_installed("limma")
    ## Full Kang has 9 control-only donors. A paired contrast must use the
    ## paired donors and say how many combinations it dropped, instead of
    ## failing the whole step (measured 2026-10-18).
    sce <- paired_de_object()
    extra <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rnbinom(60L * 40L, mu = 5, size = 2), nrow = 60L,
            ncol = 40L, dimnames = list(paste0("g", seq_len(60L)),
                paste0("x", seq_len(40L))))))
    SummarizedExperiment::colData(extra)$condition <- "ctrl"
    SummarizedExperiment::colData(extra)$donor <- "d3"
    combined <- cbind(sce, extra)

    result <- RunDEtest(combined, group = "condition", ident.1 = "stim",
        ident.2 = "ctrl", covariates = "donor", name = "paired_only")
    rec <- sclet:::sclet_get_state_record(result, "detest", "paired_only")
    expect_equal(rec$inputs$n_pseudobulk_samples, 4L)
    expect_equal(rec$inputs$excluded_covariate_combinations, "d3")
    expect_gt(rec$inputs$n_cells_used, 0L)
    expect_gt(nrow(rec$artifacts$result), 0L)
})

test_that("a contrast with no paired donor at all is still refused", {
    skip_if_not_installed("limma")
    sce <- paired_de_object()
    SummarizedExperiment::colData(sce)$donor <- ifelse(
        SummarizedExperiment::colData(sce)$condition == "stim", "d2", "d1")
    expect_error(
        RunDEtest(sce, group = "condition", ident.1 = "stim", ident.2 = "ctrl",
            covariates = "donor", name = "no_pair"),
        "not identifiable|not identifiable otherwise")
})
