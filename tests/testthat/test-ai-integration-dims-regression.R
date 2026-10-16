# Regression: run_integration must execute without an inert `dims` argument.
# RunIntegration() chooses dimensions internally; accepting `dims` in the AI
# schema silently discarded user intent, so the parameter is now rejected at
# plan validation rather than accepted and ignored.

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

test_that("a plan cannot silently request unsupported integration dims", {
    sce <- integration_regression_object()
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "reg_dims", actions = two_batch_pipeline(list(
        method = "fastMNN", batch = "batch", dims = as.integer(1:20))))
    v <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_false(v$valid)
    expect_true(any(grepl("unknown parameter.*dims", v$errors)))
})

test_that("integration dimensions are not advertised as a supported AI parameter", {
    reg <- AIDefaultExecutionRegistry(
        SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 3L, 3L))),
        include = "all"
    )
    expect_false("dims" %in% names(reg$run_integration$input_schema))
    expect_false("dims" %in% names(formals(sclet::RunIntegration)))
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

test_that("a projected multi-valued reduction reference gives a message, not a crash", {
    # The 4th layer, found in the live L5 record (remaining_errors = "the
    # condition has length > 1"): a plan referencing run_integration's
    # THREE-valued allowed_reductions projects to a length-3 vector, and the
    # unguarded `%in%` test then raised inside prerequisites. Validation's
    # tryCatch had turned the raise into the opaque
    # "prerequisite check failed: <conditionMessage>"; the guard converts it
    # back into an actionable "prerequisites not met" message the repair loop
    # can act on (the shape sclet_ai_selector_problem set for features).
    sce <- integration_regression_object()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "vec", actions = two_batch_pipeline(list(
        method = "fastMNN", batch = "batch")))
    # simulate what projection hands to prerequisites AFTER resolving
    # ${s2.output.allowed_reductions}: a length-3 character
    projected <- c("corrected", "HARMONY", "scVI")
    for (nm in c("find_neighbors", "run_umap")) {
        prereq <- reg[[nm]]$prerequisites
        out <- tryCatch(
            prereq(sce,
                params <- list(reduction = projected, dims = as.integer(1:5),
                    k = 10L),
                planned = list(reductions = c("PCA", projected),
                    reduction_dims = list(PCA = 30L), assays = "logcounts")),
            error = function(e) e)
        expect_true(is.character(out) && !inherits(out, "error"),
            info = nm)
        expect_true(grepl("single reduction", out), info = nm)
    }
    # and the plan-level path through ValidateAIPlan still produces a
    # "prerequisites not met" entry rather than a "prerequisite check failed"
    # crash entry.
    plan2 <- new_sclet_ai_plan(task = "vec2", actions = list(
        list(id = "s1", action = "normalize_data", params = list(), depends_on = character()),
        list(id = "s2", action = "find_variable_features",
            params = list(nfeatures = 30L), depends_on = "s1"),
        list(id = "s2b", action = "run_integration",
            params = list(method = "fastMNN", batch = "batch"), depends_on = "s2"),
        list(id = "s3", action = "find_neighbors",
            params = list(reduction = "${s2b.output.allowed_reductions}",
                dims = as.integer(1:5), k = 10L), depends_on = "s2b")))
    v <- ValidateAIPlan(plan2, object = ConfirmAIDesignSemantics(
        sce, design = list(batch = "batch")), registry = reg, strict = FALSE)
    expect_false(v$valid)
    msg <- paste(v$errors, collapse = " | ")
    expect_true(grepl("prerequisites not met", msg, fixed = TRUE))
    expect_false(grepl("prerequisite check failed", msg, fixed = TRUE))
})
