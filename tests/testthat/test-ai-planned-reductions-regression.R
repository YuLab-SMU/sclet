# Regression: prerequisites must consult planned$reductions, not only the
# object. A plan that runs run_pca first and then references "PCA" (or UMAP
# for run_trajectory) is VALID, but validation happens before any step
# executes -- the reduction exists only in the planned view. Two prerequisites
# read reducedDimNames(object) alone, so every such plan failed with
# "reduction is not available" and the model could not recover by rephrasing.
# Found offline by the doublet/rare benchmark dry run (2026-10-14) before any
# billed call; run_trajectory's branch is the same pattern, fixed together but
# NOT end-to-end measured in this wave (monocle3 not installed here).

toy_counts <- function(g = 30L, n = 20L, seed = 1L) {
    set.seed(seed)
    m <- matrix(rpois(g * n, 8), nrow = g)
    rownames(m) <- paste0("g", seq_len(g))
    colnames(m) <- paste0("c", seq_len(n))
    SingleCellExperiment::SingleCellExperiment(list(counts = m))
}

test_that("run_rare_cell_detection validates after a planned run_pca (not only an existing one)", {
    sce <- toy_counts()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "rare", actions = list(
        list(id = "s1", action = "normalize_data", params = list(),
            depends_on = character()),
        list(id = "s2", action = "find_variable_features",
            params = list(nfeatures = 20L), depends_on = "s1"),
        list(id = "s3", action = "run_pca", params = list(ncomponents = 10L),
            depends_on = "s2"),
        list(id = "s4", action = "run_rare_cell_detection",
            params = list(name = "rareq"), depends_on = "s3")))
    v <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_true(v$valid, info = paste(v$errors, collapse = "; "))
})

test_that("run_rare_cell_detection still rejects when NO reduction exists anywhere", {
    # The relaxation must not become a rubber stamp: with no planned PCA and no
    # object reduction the message is still required, or execution dies later.
    prereq <- AIDefaultExecutionRegistry(toy_counts(), include = "all")$
        run_rare_cell_detection$prerequisites
    out <- prereq(toy_counts(), list(reduction = "PCA"),
        planned = list(reductions = character()))
    expect_true(is.character(out))
    expect_true(grepl("reduction_missing", out))
})

test_that("run_trajectory validates a reduction an upstream step will register", {
    # Same mechanism; measured end-to-end NOT done (monocle3 absent), but the
    # prerequisite decision itself is offline and testable.
    sce <- toy_counts()
    cd <- SummarizedExperiment::colData(sce)
    cd$cluster <- rep(c("c1", "c2"), length.out = ncol(sce))
    SummarizedExperiment::colData(sce) <- cd
    sce <- sclet:::sclet_set_active_ident(sce, "cluster")
    prereq <- AIDefaultExecutionRegistry(sce, include = "all")$
        run_trajectory$prerequisites
    planned <- list(reductions = c("PCA", "UMAP"), assays = "logcounts")
    # the prerequisite demands an explicit group column and a real label; give
    # it both so the reduction check (the layer under test) is what runs.
    out <- prereq(sce, list(start_cluster = "c1", group = "cluster",
        reduction = "UMAP"), planned = planned)
    expect_true(isTRUE(out),
        info = "UMAP known only via the planned view must defer, not reject")
    # and still rejects a reduction nothing promised
    out2 <- prereq(sce, list(start_cluster = "c1", group = "cluster",
        reduction = "nope"), planned = planned)
    expect_true(is.character(out2))
})

test_that("run_rare_cell_detection executes end to end after run_pca", {
    # Validation accepting the promise is only half the contract -- execution
    # must also complete. Toy object: no rare population may qualify, which is
    # itself the honest result (the action records evidence nodes CONDITIONALLY,
    # so zero nodes with zero labels is a pass, not a silent failure).
    sce <- toy_counts()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(task = "rare_exec", actions = list(
        list(id = "s1", action = "normalize_data", params = list(),
            depends_on = character()),
        list(id = "s2", action = "find_variable_features",
            params = list(nfeatures = 20L), depends_on = "s1"),
        list(id = "s3", action = "run_pca", params = list(ncomponents = 5L),
            depends_on = "s2"),
        list(id = "s4", action = "run_rare_cell_detection",
            params = list(name = "rareq", dims = as.integer(1:5)),
            depends_on = "s3")))
    v <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    skip_if_not(isTRUE(v$valid))
    ex <- ExecuteAIPlan(sce, plan, reg, validation = v, dry_run = FALSE,
        confirmation = v$confirmation_token, record = FALSE)
    rare <- Filter(function(r) identical(r$action, "run_rare_cell_detection"),
        ex$results)
    expect_length(rare, 1L)
    expect_identical(rare[[1L]]$status, "completed")
    expect_equal(ncol(ex$object), ncol(sce))  # labelling never removes cells
})
