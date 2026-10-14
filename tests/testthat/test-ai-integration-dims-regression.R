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

# --- Deferred projection of planned-only reductions (2026-10-14) ------------
#
# A step that registers a METHOD-DEPENDENT reduction (run_integration:
# "corrected" / "HARMONY") declared it in `allowed_reductions`, which nothing
# read, so a downstream step referencing it failed validation even though
# execution really registers it. Fixing the propagation then exposed two
# further crashes: find_neighbors' dims prerequisite read
# `planned$reduction_dims[[reduction]]` (throws on a missing name) and then
# `reducedDim(object, reduction)` (throws while the reduction is only planned).
# The three tests below pin all three layers.

test_that("allowed_reductions propagates into the planned view", {
    planned <- sclet:::sclet_ai_plan_capabilities(
        SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 3L, 3L))))
    reg_env <- environment(sclet:::sclet_ai_apply_planned_output)
    descriptor <- list(
        name = "run_integration",
        output_schema = list(allowed_reductions = c("corrected", "HARMONY")))
    after <- sclet:::sclet_ai_apply_planned_output(planned, descriptor, list())
    expect_true("corrected" %in% after$reductions)
    expect_true("HARMONY" %in% after$reductions)
    # required_reductions still works alongside
    d2 <- list(name = "run_pca", output_schema = list(required_reductions = "PCA"))
    expect_true("PCA" %in% sclet:::sclet_ai_apply_planned_output(after, d2, list())$reductions)
})

test_that("find_neighbors tolerates a dims bound it cannot know at planning time", {
    # The exact pre-fix crash: reduction known only via `planned`, absent from
    # `planned$reduction_dims` AND from the object. Prerequisites must return
    # TRUE (deferred), not raise.
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 6L, 6L)))
    prereq <- AIDefaultExecutionRegistry(sce, include = "all")$find_neighbors$prerequisites
    planned <- list(reductions = c("PCA", "HARMONY"), reduction_dims = list(PCA = 30L),
        active_reduction = "PCA", assays = "logcounts")
    out <- tryCatch(
        prereq(sce, list(reduction = "HARMONY", dims = as.integer(1:10), k = 20L),
            planned = planned),
        error = function(e) e)
    # Pre-fix this raised `invalid subscript` from `planned$reduction_dims[[nm]]`
    # and then the reducedDim() crash; post-fix it defers with TRUE.
    expect_true(isTRUE(out), info = paste(format(out), collapse = " | "))
    # A dim count that IS known still bounds correctly.
    planned2 <- list(reductions = c("PCA"), reduction_dims = list(PCA = 5L),
        active_reduction = "PCA", assays = "logcounts")
    out2 <- prereq(sce, list(reduction = "PCA", dims = as.integer(1:10), k = 20L),
        planned = planned2)
    expect_true(is.character(out2))
    expect_true(grepl("dims", out2))
})

test_that("a plan referencing the planned integration reduction validates and executes", {
    # End to end, fastMNN (registers `corrected`, the measured name): before
    # the fix this failed at validation with "reduction is not available:
    # corrected"; after the propagation fix it failed deeper with the
    # `invalid subscript` crash. Both layers must hold now.
    sce <- integration_regression_object()
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "deferred", actions = list(
        list(id = "s1", action = "normalize_data", params = list(), depends_on = character()),
        list(id = "s2", action = "find_variable_features",
            params = list(nfeatures = 30L), depends_on = "s1"),
        list(id = "s3", action = "run_integration",
            params = list(method = "fastMNN", batch = "batch"), depends_on = "s2"),
        list(id = "s4", action = "find_neighbors",
            params = list(reduction = "corrected", dims = as.integer(1:5), k = 10L),
            depends_on = "s3")))
    v <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_true(v$valid, info = paste(v$errors, collapse = "; "))
    ex <- ExecuteAIPlan(sce, plan, reg, validation = v, dry_run = FALSE,
        confirmation = v$confirmation_token, record = FALSE)
    expect_identical(ex$status, "completed")
    expect_true("corrected" %in% SingleCellExperiment::reducedDimNames(ex$object))
})
