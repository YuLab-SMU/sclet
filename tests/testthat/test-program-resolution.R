library(testthat)
library(sclet)
library(SingleCellExperiment)

# Issue #28: source = "auto" was passed straight to the resolver, which never
# handled it, so plot_program_dotplot(source = "auto") failed with the opaque
# "No program activity data could be resolved." message even when SCENIC results
# were present. These tests lock in the fix at the shared resolver.

make_geneset_sce <- function() {
    counts <- matrix(rpois(200, 5), nrow = 10)
    rownames(counts) <- paste0("gene", 1:10)
    colnames(counts) <- paste0("cell", 1:20)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    score <- runif(20)
    names(score) <- colnames(counts)
    SummarizedExperiment::colData(sce)[["Score_ARID3A(+)"]] <- score
    sce <- sclet:::sclet_set_state_record(
        sce,
        type = "geneset_scoring",
        id = "Score",
        active = TRUE,
        value = list(
            method = "UCell",
            artifacts = list(storage = "colData", score_columns = c("Score_ARID3A(+)"))
        )
    )
    list(sce = sce, score = score)
}

make_scenic_sce <- function() {
    counts <- matrix(rpois(200, 5), nrow = 10)
    rownames(counts) <- paste0("gene", 1:10)
    colnames(counts) <- paste0("cell", 1:20)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    regulons <- c("ARID3A(+)", "ATF3(+)", "ATF4(+)", "GRHL1(+)")
    auc <- matrix(
        runif(length(regulons) * ncol(counts)),
        nrow = length(regulons),
        dimnames = list(regulons, colnames(counts))
    )
    SingleCellExperiment::altExp(sce, "SCENIC_AUC") <-
        SingleCellExperiment::SingleCellExperiment(assays = list(AUC = auc))
    sce <- sclet:::sclet_set_state_record(
        sce,
        type = "scenic",
        id = "scenic_main",
        active = TRUE,
        value = list(
            method = "pySCENIC",
            artifacts = list(altExp = "SCENIC_AUC", assay = "AUC")
        )
    )
    list(sce = sce, auc = auc)
}

test_that("get_program(source='auto') resolves gene set scoring colData when no SCENIC result exists", {
    skip_if_not_installed("SingleCellExperiment")
    fx <- make_geneset_sce()
    act <- get_program(fx$sce, "ARID3A(+)", source = "auto")
    expect_length(act, 20L)
    expect_equal(act, as.numeric(fx$score), ignore_attr = TRUE)
})

test_that("get_program(source='auto') resolves SCENIC AUC when gene set scoring is absent", {
    skip_if_not_installed("SingleCellExperiment")
    fx <- make_scenic_sce()
    act <- get_program(fx$sce, "ARID3A(+)", source = "auto")
    expect_length(act, 20L)
    expect_equal(act, as.numeric(fx$auc["ARID3A(+)", ]), ignore_attr = TRUE)
})

test_that("get_program(source='auto') reports every attempted source when nothing can resolve", {
    skip_if_not_installed("SingleCellExperiment")
    counts <- matrix(rpois(20, 5), nrow = 5)
    rownames(counts) <- paste0("gene", 1:5)
    colnames(counts) <- paste0("cell", 1:4)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    expect_error(
        get_program(sce, "ARID3A(+)", source = "auto"),
        "could not be resolved from any available source"
    )
    expect_error(
        get_program(sce, "ARID3A(+)", source = "auto"),
        "geneset_scoring: No gene set scoring results found"
    )
    expect_error(
        get_program(sce, "ARID3A(+)", source = "auto"),
        "scenic: No SCENIC results found"
    )
})

test_that("plot_program_dotplot(source='auto') uses SCENIC regulon activity (issue #28 scenario)", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("ggplot2")
    fx <- make_scenic_sce()
    SingleCellExperiment::colLabels(fx$sce) <- rep(c("A", "B"), each = 10)

    p <- plot_program_dotplot(
        fx$sce,
        programs = c("ARID3A(+)", "ATF3(+)"),
        source = "auto",
        group.by = "colLabels"
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_program_dotplot surfaces per-program resolution failures instead of a blank error", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("ggplot2")
    counts <- matrix(rpois(20, 5), nrow = 5)
    rownames(counts) <- paste0("gene", 1:5)
    colnames(counts) <- paste0("cell", 1:4)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SingleCellExperiment::colLabels(sce) <- rep(c("A", "B"), each = 2)

    err <- tryCatch(
        plot_program_dotplot(sce, programs = c("ARID3A(+)"), source = "auto", group.by = "colLabels"),
        error = function(e) e
    )
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "could not be resolved from any available source")
})

test_that("plot_program_heatmap and has_program honor source='auto'", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("ggplot2")
    fx <- make_scenic_sce()
    SingleCellExperiment::colLabels(fx$sce) <- rep(c("A", "B"), each = 10)

    hp <- plot_program_heatmap(fx$sce, programs = c("ARID3A(+)", "ATF3(+)"), source = "auto", group.by = "colLabels")
    expect_s3_class(hp, "ggplot")

    expect_true(has_program(fx$sce, "ARID3A(+)", source = "auto"))
    expect_false(has_program(fx$sce, "NOT_A_REGULON", source = "auto"))
})
