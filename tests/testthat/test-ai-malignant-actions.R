## Malignant-cell identification as AI actions.
##
## The merged feature (36d2156) exposed RunCNVInference()/CallMalignantCells(),
## but no action reached them. These tests pin the two actions, their
## prerequisites (exogenous reference, two-line call) and the goal mapping.

malignant_test_object <- function() {
    set.seed(1)
    counts <- matrix(stats::rnbinom(50L * 20L, mu = 5, size = 2),
        nrow = 50L, ncol = 20L,
        dimnames = list(paste0("g", seq_len(50L)), paste0("c", seq_len(20L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    SummarizedExperiment::colData(sce)$reference <- rep(c(TRUE, FALSE), each = 10L)
    SummarizedExperiment::colData(sce)$malignant_score <- stats::rnorm(20L)
    sce
}

test_that("the malignant group is registered and mapped from goals", {
    sce <- malignant_test_object()
    registry <- AIDefaultExecutionRegistry(sce, include = "malignant")
    expect_true(all(c("run_cnv_inference", "call_malignant_cells") %in% names(registry)))
    expect_false("run_cnv_inference" %in% names(AIDefaultExecutionRegistry(sce, include = "read")))

    coverage <- sclet:::sclet_ai_goal_coverage
    plan <- sclet:::new_sclet_ai_plan(task = "t",
        actions = list(list(id = "s1", action = "call_malignant_cells", params = list())))
    result <- coverage(plan, goal = "Identify malignant cells in this tumour sample.")
    expect_true("call_malignant_cells" %in% result$required)
    cnv_only <- coverage(sclet:::new_sclet_ai_plan(task = "t", actions = list()),
        goal = "Infer copy-number burden from the normal cells.")
    expect_true("run_cnv_inference" %in% cnv_only$required)
})

test_that("run_cnv_inference refuses a missing or endogenous reference", {
    sce <- malignant_test_object()
    action <- AIDefaultExecutionRegistry(sce, include = "malignant")$run_cnv_inference
    ## inferCNV is optional: without it the action says so instead of failing
    ## later inside CreateInfercnvObject.
    missing_pkg <- action$prerequisites(sce, list(reference_col = "reference"))
    if (!requireNamespace("inferCNV", quietly = TRUE)) {
        expect_true(grepl("optional_package_missing", missing_pkg, fixed = TRUE))
    }
    ## No reference at all is refused; the message says why clusters are not one.
    if (requireNamespace("inferCNV", quietly = TRUE)) {
        reason <- action$prerequisites(sce, list())
        expect_true(grepl("reference_missing", reason, fixed = TRUE))
        expect_true(grepl("Do not use clusters", reason, fixed = TRUE))
        small <- action$prerequisites(sce, list(reference_ids = "c1"))
        expect_true(grepl("reference_too_small", small, fixed = TRUE))
    }
})

test_that("call_malignant_cells requires a CNV record and a signature line", {
    sce <- malignant_test_object()
    action <- AIDefaultExecutionRegistry(sce, include = "malignant")$call_malignant_cells

    ## No CNV record anywhere: refused, and the reason names the two-line rule.
    none <- action$prerequisites(sce, list(signature_score = "malignant_score"))
    expect_true(grepl("cnv_record_missing", none, fixed = TRUE))
    expect_true(grepl("TWO-line", none, fixed = TRUE))

    ## A CNV step earlier in the same plan satisfies validation.
    planned <- action$prerequisites(sce,
        list(cnv_id = "cnv", signature_score = "malignant_score"),
        planned = list(cnv_ids = "cnv"))
    expect_true(isTRUE(planned))

    ## Without a signature line there is no second observation.
    no_sig <- action$prerequisites(sce, list(cnv_id = "cnv"),
        planned = list(cnv_ids = "cnv"))
    expect_true(grepl("signature_missing", no_sig, fixed = TRUE))

    ## A signature_score naming a missing column is refused.
    bad_col <- action$prerequisites(sce, list(cnv_id = "cnv", signature_score = "nope"),
        planned = list(cnv_ids = "cnv"))
    expect_true(grepl("signature_score_missing", bad_col, fixed = TRUE))
})

test_that("a two-step malignant plan validates on an unconfirmed object", {
    sce <- malignant_test_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "malignant"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "malignant",
        actions = list(
            list(id = "cnv1", action = "run_cnv_inference",
                params = list(reference_col = "reference", name = "cnv",
                    assay = "counts")),
            list(id = "call1", action = "call_malignant_cells",
                params = list(cnv_id = "cnv", signature_score = "malignant_score",
                    name = "malignant_call"))
        ))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    ## If inferCNV is absent the first step is unavailable for that reason and
    ## the call step must still not complain about a missing CNV record.
    expect_false(any(grepl("cnv_record_missing", validation$errors)))
})
