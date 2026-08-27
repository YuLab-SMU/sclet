make_workflow_sce <- function(n_per_cluster = 30) {
    centers <- list(
        A = c(-2, -1.5),
        B = c(0, 0),
        C = c(2, 1.5)
    )
    clusters <- rep(names(centers), each = n_per_cluster)
    set.seed(20260827)
    coords <- do.call(rbind, lapply(clusters, function(cl) {
        centers[[cl]] + rnorm(2, sd = 0.3)
    }))
    counts <- matrix(
        rpois(length(coords) * 10, lambda = 2),
        nrow = 20,
        ncol = length(clusters)
    )
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
    rownames(counts) <- paste0("gene", seq_len(nrow(counts)))

    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    pca <- coords
    rownames(pca) <- colnames(counts)
    colnames(pca) <- paste0("PC", seq_len(ncol(pca)))
    SingleCellExperiment::reducedDim(sce, "PCA") <- pca
    SummarizedExperiment::colData(sce)$cluster <- factor(clusters)
    sce
}

test_that("RunTrajectoryWorkflow runs the slingshot mainline and records state", {
    skip_if_not_installed("slingshot")

    sce <- make_workflow_sce()
    sce <- RunTrajectoryWorkflow(
        sce,
        group = "cluster",
        reduction = "PCA",
        run_diffusion_map = FALSE,
        run_slingshot = TRUE,
        run_velocity = "never",
        run_fate = "never",
        name = "trajectory_workflow"
    )

    expect_s4_class(sce, "SingleCellExperiment")
    expect_true(any(grepl("^slingPseudotime", names(SummarizedExperiment::colData(sce)))))

    workflow <- sclet_get_analysis(sce, "trajectory_workflow")
    expect_identical(workflow$id, "trajectory_workflow")
    expect_identical(workflow$inputs$group, "cluster")
    expect_identical(workflow$inputs$reduction, "PCA")
    expect_true(workflow$artifacts$trajectory == "trajectory_workflow_slingshot")
    expect_null(workflow$artifacts$velocity)
    expect_null(workflow$artifacts$fate)
    expect_true(isTRUE(workflow$summary$has_slingshot))
    expect_false(isTRUE(workflow$summary$has_velocity))
    expect_false(isTRUE(workflow$summary$has_fate))

    traj_id <- workflow$artifacts$trajectory
    expect_true(has_trajectory(sce, id = traj_id))
    traj <- get_trajectory(sce, id = traj_id)
    expect_identical(traj$method, "slingshot")
    expect_true(traj$summary$n_lineages >= 1)
    expect_true(all(grepl("^slingPseudotime", traj$artifacts$pseudotime_cols)))

    commands <- CommandLog(sce)
    expect_true("RunTrajectoryWorkflow" %in% commands$command)
})

test_that("RunTrajectoryWorkflow computes diffusion map when enabled", {
    skip_if_not_installed("slingshot")
    skip_if_not_installed("destiny")

    sce <- make_workflow_sce()
    sce <- RunTrajectoryWorkflow(
        sce,
        group = "cluster",
        reduction = "PCA",
        run_diffusion_map = TRUE,
        diffusion_components = 2,
        run_slingshot = FALSE,
        run_velocity = "never",
        run_fate = "never"
    )

    expect_true("DM" %in% SingleCellExperiment::reducedDimNames(sce))
    dm <- SingleCellExperiment::reducedDim(sce, "DM")
    expect_equal(ncol(dm), 2)
})

test_that("RunTrajectoryWorkflow auto mode skips velocity without spliced assays", {
    skip_if_not_installed("slingshot")

    sce <- make_workflow_sce()
    sce <- RunTrajectoryWorkflow(
        sce,
        group = "cluster",
        reduction = "PCA",
        run_diffusion_map = FALSE,
        run_slingshot = FALSE,
        run_velocity = "auto",
        run_fate = "never"
    )

    workflow <- sclet_get_analysis(sce, "trajectory_workflow")
    expect_false(isTRUE(workflow$summary$has_velocity))
    expect_false(has_velocity(sce))
})

test_that("RunTrajectoryWorkflow validates velocity and fate prerequisites", {
    sce <- make_workflow_sce()

    expect_error(
        RunTrajectoryWorkflow(
            sce,
            group = "cluster",
            reduction = "PCA",
            run_diffusion_map = FALSE,
            run_slingshot = FALSE,
            run_velocity = "always"
        ),
        regexp = "spliced"
    )

    expect_error(
        RunTrajectoryWorkflow(
            sce,
            group = "cluster",
            reduction = "PCA",
            run_diffusion_map = FALSE,
            run_slingshot = FALSE,
            run_velocity = "never",
            run_fate = "always"
        ),
        regexp = "Velocity results are required"
    )
})

test_that("RunTrajectoryWorkflow validates group and reduction inputs", {
    sce <- make_workflow_sce()

    expect_error(
        RunTrajectoryWorkflow(sce, group = "missing_cluster"),
        regexp = "not found in colData"
    )

    bare <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 5, ncol = 6))
    )
    SummarizedExperiment::colData(bare)$cluster <- factor(rep(c("A", "B"), 3))
    expect_error(
        RunTrajectoryWorkflow(bare, group = "cluster"),
        regexp = "No active reduction found"
    )

    no_ident <- make_workflow_sce()
    expect_error(
        RunTrajectoryWorkflow(no_ident, group = NULL),
        regexp = "No active identity found"
    )
})
