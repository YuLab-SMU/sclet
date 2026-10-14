# Regression: run_integration must execute, and must not forward the inert
# `dims` parameter to RunIntegration().
#
# Found offline (2026-10-14) by the integration benchmark's dry run before any
# provider was involved: the handler forwarded `dims = dims` to
# RunIntegration(), which has no such argument, so every call raised
# `unused argument (dims = ...)` -- the action could not succeed on any plan,
# with any parameters, regardless of model behaviour. The fix removed the
# forwarding; this file pins both halves: execution now completes, and a plan
# that DOES carry `dims` (the parameter is still declared in `input_schema`)
# also completes rather than dying on the unused argument.

# Minimal two-batch object: raw counts, a real batch column, no design
# confirmation (the benchmark fixture covers the confirmation prerequisite).
integration_regression_object <- function(seed = 1L) {
    set.seed(seed)
    g <- 60L
    n_per <- 60L
    n <- n_per * 2L
    mu <- matrix(8, nrow = g, ncol = n)
    # a real shift in genes 1:20 for batch B so the correction has something
    # to see; the rest is background
    mu[1:20, (n_per + 1L):n] <- 14
    counts <- matrix(rpois(g * n, as.vector(mu)), nrow = g)
    rownames(counts) <- paste0("g", seq_len(g))
    colnames(counts) <- paste0("c", seq_len(n))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce$batch <- rep(c("A", "B"), each = n_per)
    sce
}

two_batch_pipeline <- function(integration_params) {
    list(
        list(id = "s1", action = "normalize_data", params = list(), depends_on = character()),
        list(id = "s2", action = "find_variable_features",
            params = list(nfeatures = 30L), depends_on = "s1"),
        list(id = "s3", action = "run_integration",
            params = integration_params, depends_on = "s2")
    )
}

test_that("run_integration executes end to end and registers its reduction", {
    sce <- integration_regression_object()
    # the action's hard prerequisite: a confirmed batch design
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "reg", actions = two_batch_pipeline(list(
        method = "fastMNN", batch = "batch")))
    v <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_true(v$valid, info = paste(v$errors, collapse = "; "))

    ex <- ExecuteAIPlan(sce, plan, reg, validation = v, dry_run = FALSE,
        confirmation = v$confirmation_token, record = FALSE)
    integ <- Filter(function(r) identical(r$action, "run_integration"), ex$results)
    expect_identical(integ[[1L]]$status, "completed")
    # fastMNN registers the reduction named after the method
    expect_true("corrected" %in% SingleCellExperiment::reducedDimNames(ex$object))
    # integration must not drop or add cells
    expect_equal(ncol(ex$object), ncol(sce))
})

test_that("a plan carrying dims still executes (the parameter is inert, not fatal)", {
    # The exact pre-fix failure: `dims` in the params reached RunIntegration()
    # as an unused argument. A plan may still carry it -- it is declared in
    # input_schema -- so it must be accepted and simply ignored.
    sce <- integration_regression_object()
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "reg_dims", actions = two_batch_pipeline(list(
        method = "fastMNN", batch = "batch", dims = as.integer(1:20))))
    v <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_true(v$valid, info = paste(v$errors, collapse = "; "))

    ex <- ExecuteAIPlan(sce, plan, reg, validation = v, dry_run = FALSE,
        confirmation = v$confirmation_token, record = FALSE)
    integ <- Filter(function(r) identical(r$action, "run_integration"), ex$results)
    expect_identical(integ[[1L]]$status, "completed")
    expect_true("corrected" %in% SingleCellExperiment::reducedDimNames(ex$object))
})

test_that("the handler does not forward dims to RunIntegration", {
    # Pin the mechanism, not just the outcome: the forwarding line was the bug.
    src <- paste(readLines("R/ai-execution.R"), collapse = "\n")
    handler <- substr(src,
        regexpr("run_integration = AIAction", src),
        regexpr("run_integration = AIAction", src) + 900L)
    expect_false(grepl("dims = dims", handler, fixed = TRUE))
    # and RunIntegration genuinely has no such formal argument
    fmls <- names(formals(sclet::RunIntegration))
    expect_false("dims" %in% fmls)
})
