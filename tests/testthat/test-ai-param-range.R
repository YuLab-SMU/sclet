## `dims`-style integer vectors accept the R range spelling a real model writes.
##
## Measured need: `dims="1:30"` appears as the first-attempt error in
## ai_plan_repair records 25-26, run_annotation_bench records 1 and 4, and
## fate_chain_bench -- one wasted repair round per run. The fix is narrow: only
## a range string, only for vector-typed parameters, with a size guard.

test_that("only a range string is parsed, and only into a vector", {
    parse <- sclet:::sclet_ai_range_vector
    expect_equal(parse("1:10"), 1:10)
    expect_equal(parse("1:10"), seq.int(1L, 10L))
    expect_equal(parse("c(1:10)"), 1:10)
    expect_equal(parse(" 3 : 7 "), 3:7)
    expect_equal(parse("-2:2"), -2:2)

    ## Not ranges: a scalar string, a selector, prose, a single number, NA.
    expect_null(parse("50"))
    expect_null(parse("hvg"))
    expect_null(parse("1:30 dims"))
    expect_null(parse(NA_character_))
    expect_null(parse(c("1:10", "2:12")))
    ## Guarded: a typo must not expand into a large allocation.
    expect_null(parse("1:10000000"))
})

test_that("vector-typed parameters accept a range while scalars still refuse a string", {
    valid <- sclet:::sclet_ai_param_type_valid
    expect_true(valid("1:10", "integer_vector"))
    expect_true(valid("c(1:10)", "integer_vector"))
    expect_true(valid("1:10", "numeric_vector"))
    ## A scalar is a scalar: the old failure mode ("50" for ncomponents) stays.
    expect_false(valid("50", "integer"))
    expect_false(valid("50", "numeric"))
    ## And a non-range string is still not a vector.
    expect_false(valid("hvg", "integer_vector"))
})

test_that("a plan with dims as a range validates and executes", {
    set.seed(1L)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(40L * 30L, 5), nrow = 40L, ncol = 30L,
            dimnames = list(paste0("gene_", seq_len(40L)), paste0("c", seq_len(30L)))))
    )
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20L)
    sce <- RunPCA(sce, ncomponents = 8L)
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- sclet:::new_sclet_ai_plan(
        task = "neighbors from a range spelling",
        actions = list(list(id = "n1", action = "find_neighbors",
            params = list(dims = "1:5", reduction = "PCA")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, validation$plan, registry,
        validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")
})

test_that("the range spelling does not make a scalar parameter accept a string", {
    set.seed(1L)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(30L * 20L, 5), nrow = 30L, ncol = 20L,
            dimnames = list(paste0("gene_", seq_len(30L)), paste0("c", seq_len(20L)))))
    )
    sce <- NormalizeData(sce)
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- sclet:::new_sclet_ai_plan(
        task = "scalar string must stay invalid",
        actions = list(list(id = "p1", action = "run_pca",
            params = list(ncomponents = "20")))
    )
    err <- tryCatch(ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        error = function(e) e)
    expect_s3_class(err, "sclet_ai_invalid_plan")
    expect_true(grepl("ncomponents", conditionMessage(err), fixed = TRUE))
})
