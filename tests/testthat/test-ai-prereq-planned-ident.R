# Prerequisite checks must account for state an earlier step in the same plan
# will create.
#
# run_de_test's prerequisite read ActiveIdent(object) and ignored
# planned$active_ident. Since find_clusters declares
# output_schema$active_ident = "colLabels", a plan of
# cluster -> find markers was rejected at validation time on any object that
# did not *already* have an identity -- i.e. on raw counts, which is the
# ordinary case for the request "cluster the cells and find markers for each
# cluster". The plan could not be proposed at all.
#
# The identity a later step will see is not knowable before the producing step
# runs: neither its group sizes nor its label set. So the check is deferred,
# and the producing step's output contract plus the execution-time check remain
# the authority.

de_ident_test_sce <- function(n_cells = 8L, n_genes = 12L, groups = NULL) {
    counts <- matrix(
        seq_len(n_genes * n_cells) %% 7L,
        nrow = n_genes, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(n_genes)), paste0("cell_", seq_len(n_cells)))
    )
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    if (!is.null(groups)) {
        SummarizedExperiment::colData(sce)$grp <- rep(groups, length.out = n_cells)
        ActiveIdent(sce) <- "grp"
    }
    sce
}

test_that("run_de_test defers when an earlier step will set the identity", {
    raw <- de_ident_test_sce()
    expect_null(ActiveIdent(raw))
    descriptor <- AIDefaultExecutionRegistry(raw, include = "all")[["run_de_test"]]

    # find_clusters declares active_ident = "colLabels".
    check <- descriptor$prerequisites(raw, list(), list(active_ident = "colLabels"))
    expect_true(isTRUE(check))
})

test_that("run_de_test still refuses when no earlier step will change the identity", {
    one_group <- de_ident_test_sce(groups = "a")
    descriptor <- AIDefaultExecutionRegistry(one_group, include = "all")[["run_de_test"]]

    # planned identity equals the current one, so the check is meaningful.
    check <- descriptor$prerequisites(one_group, list(), list(active_ident = "grp"))
    expect_true(is.character(check))
    expect_true(grepl("at least two distinct identity groups", check))

    two_groups <- de_ident_test_sce(groups = c("a", "b"))
    descriptor2 <- AIDefaultExecutionRegistry(two_groups, include = "all")[["run_de_test"]]
    expect_true(isTRUE(
        descriptor2$prerequisites(two_groups, list(), list(active_ident = "grp"))
    ))
})

test_that("the deferral also covers re-clustering an already-clustered object", {
    # Here the object HAS an identity, but the plan replaces it. Checking the
    # object's current groups would validate the wrong identity.
    clustered <- de_ident_test_sce(groups = c("a", "b"))
    descriptor <- AIDefaultExecutionRegistry(clustered, include = "all")[["run_de_test"]]
    check <- descriptor$prerequisites(clustered, list(), list(active_ident = "colLabels"))
    expect_true(isTRUE(check))
})

de_pipeline_sce <- function(n_per_group = 20L, n_genes = 40L, seed = 1L) {
    set.seed(seed)
    n <- n_per_group * 2L
    mu <- matrix(stats::runif(n_genes, 1, 5), nrow = n_genes, ncol = n)
    mu[1:10, 1:n_per_group] <- mu[1:10, 1:n_per_group] * 6
    mu[11:20, (n_per_group + 1L):n] <- mu[11:20, (n_per_group + 1L):n] * 6
    counts <- matrix(stats::rpois(n_genes * n, lambda = as.vector(mu)), nrow = n_genes)
    dimnames(counts) <- list(paste0("gene_", seq_len(n_genes)), paste0("cell_", seq_len(n)))
    SingleCellExperiment::SingleCellExperiment(list(counts = counts))
}

test_that("cluster-then-markers validates and runs from a raw counts object", {
    sce <- de_pipeline_sce()
    expect_null(ActiveIdent(sce))

    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(
            list(id = "s1", action = "normalize_data", params = list()),
            list(id = "s2", action = "find_variable_features", params = list(nfeatures = 20L)),
            list(id = "s3", action = "scale_data", params = list()),
            list(id = "s4", action = "run_pca", params = list(ncomponents = 10L)),
            list(id = "s5", action = "find_neighbors", params = list(dims = as.integer(1:10)), depends_on = "s4"),
            list(id = "s6", action = "find_clusters", params = list()),
            list(id = "s7", action = "run_de_test", params = list())
        )
    )

    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    expect_true(validation$valid)

    execution <- RunAIPlan(
        object = sce, plan = plan, registry = registry,
        dry_run = FALSE, confirm = TRUE
    )
    expect_equal(execution$status, "completed")
    # The deferral is only safe because execution proves the identity really
    # exists and has more than one group.
    idents <- as.character(Idents(execution$object))
    expect_gte(length(unique(idents)), 2L)
    expect_true(all(vapply(execution$results, function(r) identical(r$status, "completed"), logical(1))))
})
