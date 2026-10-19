# Tests for the integration reduction as a referenceable plan output.
#
# The measured failure this guards (2026-10-19, real-data Kang benchmark):
# the model wrote the standard advanced pipeline --
# normalize -> HVG -> scale -> PCA -> confirm_design -> run_integration ->
# find_neighbors(reduction = "${step_6.output.reduction}") -> find_clusters ->
# run_annotation -> run_de_test(group, covariates) -- and validation refused it
# with
#   prerequisites not met for action step_7:
#   reduction is not available: ${step_6.output.reduction}
#
# The reference could not project because `run_integration`'s output_schema has
# no `reduction` key: the name it registers depends on the method
# (fastMNN -> "corrected", Harmony -> "HARMONY", scVI -> "scVI"), so it cannot
# be a static value like run_pca's `reduction = "PCA"`. The template therefore
# reached find_neighbors' prerequisite as a literal string.
#
# That made the single most standard advanced workflow unplannable: integrate,
# then cluster on the integrated embedding. It cost one wasted billed run
# (record 29 of results/realdata_kang_ifnb_bench.jsonl) to see it.

integration_reference_fixture <- function(seed = 1L) {
    set.seed(seed)
    n <- 60L
    g <- 80L
    mu <- matrix(runif(g, 1, 5), nrow = g, ncol = n)
    mu[1:20, 1:30] <- mu[1:20, 1:30] * 6
    counts <- matrix(rpois(g * n, as.vector(mu)), nrow = g)
    rownames(counts) <- paste0("g", seq_len(g))
    colnames(counts) <- paste0("c", seq_len(n))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    cd <- SummarizedExperiment::colData(sce)
    cd$batch <- rep(c("b1", "b2"), length.out = n)
    cd$condition <- rep(c("ctrl", "stim"), each = n / 2)
    cd$donor <- rep(paste0("d", 1:6), length.out = n)
    SummarizedExperiment::colData(sce) <- cd
    sce
}

# The chain from the failing plan, with the reference as the model wrote it.
integration_reference_plan <- function(method = "Harmony") {
    new_sclet_ai_plan(
        task = "integration_reference",
        actions = list(
            list(id = "s1", action = "normalize_data", params = list()),
            list(id = "s2", action = "run_pca", params = list(ncomponents = 10L),
                depends_on = "s1"),
            list(id = "s3", action = "confirm_design",
                params = list(batch = "batch", condition = "condition",
                    individual = "donor"), depends_on = "s2"),
            list(id = "s4", action = "run_integration",
                params = list(batch = "batch", method = method, name = "harmony"),
                depends_on = "s3"),
            list(id = "s5", action = "find_neighbors",
                params = list(dims = as.integer(1:10),
                    reduction = "${s4.output.reduction}", k = 10L),
                depends_on = "s4")
        )
    )
}

test_that("a plan can reference the reduction run_integration registers", {
    sce <- integration_reference_fixture()
    sce <- NormalizeData(sce)
    sce <- RunPCA(sce, ncomponents = 10L)
    sce <- ConfirmAIDesignSemantics(sce,
        design = list(batch = "batch", condition = "condition",
            subject = "donor"))
    reg <- AIDefaultExecutionRegistry(sce, include = "all")

    plan <- integration_reference_plan("Harmony")
    validation <- ValidateAIPlan(plan, registry = reg, object = sce)

    # Validation accepts it. That is the whole defect: before the fix this same
    # plan was refused with "reduction is not available:
    # ${s4.output.reduction}".
    expect_true(isTRUE(validation$valid))
})

test_that("the declared name matches what each method really registers", {
    # The mapping is the contract the projection relies on, so it is pinned
    # against the methods the action offers rather than restated.
    expect_identical(sclet:::sclet_ai_integration_reduction("fastMNN"), "corrected")
    expect_identical(sclet:::sclet_ai_integration_reduction("Harmony"), "HARMONY")
    expect_identical(sclet:::sclet_ai_integration_reduction("scVI"), "scVI")
    # Harmony is the method that can actually run here; confirm the name on the
    # object rather than trusting the switch above.
    skip_if_not_installed("harmony")
    sce <- integration_reference_fixture()
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 40L)
    sce <- RunPCA(sce, ncomponents = 10L)
    integrated <- RunIntegration(sce, method = "Harmony", batch = "batch",
        name = "harmony")
    expect_true("HARMONY" %in%
        SingleCellExperiment::reducedDimNames(integrated))
})

test_that("a static reduction declaration is left alone", {
    # run_pca declares `reduction = "PCA"`; the integration special case must
    # not touch plans that never integrate.
    sce <- integration_reference_fixture()
    sce <- NormalizeData(sce)
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "pca_only", actions = list(
        list(id = "p1", action = "run_pca", params = list(ncomponents = 10L))))
    validation <- ValidateAIPlan(plan, registry = reg, object = sce)
    expect_true(isTRUE(validation$valid))
    expect_identical(
        as.character(reg$run_integration$output_schema$reduction %||% NA_character_),
        NA_character_)
})
