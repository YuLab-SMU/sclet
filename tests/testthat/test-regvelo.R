test_that("RegVelo prior accepts matrix and edge table inputs", {
    genes <- paste0("gene", 1:5)
    W <- matrix(
        c(1, 0, 0, 2, 3, 0),
        nrow = 3,
        dimnames = list(c("gene1", "gene2", "missing"), c("gene3", "gene4"))
    )

    prior <- sclet_prepare_regvelo_prior(W, genes)
    expect_equal(prior$genes, c("gene1", "gene2"))
    expect_equal(prior$regulators, c("gene3", "gene4"))
    expect_equal(dim(prior$W), c(2L, 2L))

    edge_prior <- sclet_prepare_regvelo_prior(
        data.frame(
            target = c("gene1", "gene2", "missing"),
            regulator = c("gene3", "gene4", "gene3"),
            weight = c(0.5, 1.5, 9)
        ),
        genes
    )
    expect_equal(edge_prior$genes, c("gene1", "gene2"))
    expect_equal(edge_prior$regulators, c("gene3", "gene4"))
    expect_equal(edge_prior$W["gene1", "gene3"], 0.5)
    expect_equal(edge_prior$W["gene2", "gene4"], 1.5)
})

test_that("plot_velocity_latent_time uses RegVelo latent-time artifacts", {
    skip_if_not_installed("ggplot2")

    counts <- matrix(rpois(50, 5), nrow = 10)
    rownames(counts) <- paste0("gene", 1:10)
    colnames(counts) <- paste0("cell", 1:5)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SingleCellExperiment::reducedDim(sce, "UMAP") <- matrix(rnorm(10), ncol = 2)
    sce$regvelo_latent_time <- seq(0, 1, length.out = ncol(sce))

    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            artifacts = list(colData = "regvelo_latent_time"),
            inputs = list(reduction = "UMAP")
        )
    )

    p <- plot_velocity_latent_time(sce, id = "regvelo")
    expect_s3_class(p, "ggplot")
})

test_that("VelocityMagnitude summarizes stored RegVelo velocity assays", {
    counts <- matrix(0, nrow = 3, ncol = 2)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:2)
    velocity <- matrix(
        c(3, 4, 0, 0, 0, 12),
        nrow = 3,
        dimnames = dimnames(counts)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts, regvelo_velocity = velocity)
    )
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(rnorm(4), ncol = 2)
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            inputs = list(reduction = "PCA"),
            artifacts = list(velocity_assay = "regvelo_velocity")
        )
    )

    expect_equal(VelocityMagnitude(sce, id = "regvelo", margin = "cell"), c(cell1 = 5, cell2 = 12))
    expect_equal(VelocityMagnitude(sce, id = "regvelo", margin = "gene"), c(gene1 = 3, gene2 = 4, gene3 = 12))

    sce2 <- VelocityMagnitude(sce, id = "regvelo", margin = "cell", name = "velocity_norm", store = TRUE)
    expect_true("velocity_norm" %in% colnames(SummarizedExperiment::colData(sce2)))
})

test_that("plot_velocity_magnitude uses stored velocity assay", {
    skip_if_not_installed("ggplot2")

    counts <- matrix(0, nrow = 3, ncol = 2)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:2)
    velocity <- matrix(
        c(3, 4, 0, 0, 0, 12),
        nrow = 3,
        dimnames = dimnames(counts)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts, regvelo_velocity = velocity)
    )
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(rnorm(4), ncol = 2)
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            inputs = list(reduction = "PCA"),
            artifacts = list(velocity_assay = "regvelo_velocity")
        )
    )

    p <- plot_velocity_magnitude(sce, id = "regvelo")
    expect_s3_class(p, "ggplot")
})

test_that("TopVelocityGenes ranks genes from stored velocity assay", {
    counts <- matrix(0, nrow = 3, ncol = 2)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:2)
    velocity <- matrix(
        c(3, 4, 0, 0, 0, 12),
        nrow = 3,
        dimnames = dimnames(counts)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts, regvelo_velocity = velocity)
    )
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            artifacts = list(velocity_assay = "regvelo_velocity")
        )
    )

    top <- TopVelocityGenes(sce, id = "regvelo", n = 2)
    expect_equal(top$gene, c("gene3", "gene2"))
    expect_equal(top$velocity_magnitude, c(12, 4))
})

test_that("plot_top_velocity_genes returns ggplot", {
    skip_if_not_installed("ggplot2")

    counts <- matrix(0, nrow = 3, ncol = 2)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:2)
    velocity <- matrix(
        c(3, 4, 0, 0, 0, 12),
        nrow = 3,
        dimnames = dimnames(counts)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts, regvelo_velocity = velocity)
    )
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            artifacts = list(velocity_assay = "regvelo_velocity")
        )
    )

    p <- plot_top_velocity_genes(sce, id = "regvelo")
    expect_s3_class(p, "ggplot")
})

test_that("CellRank input preparation uses RegVelo velocity assay from object", {
    counts <- matrix(0, nrow = 3, ncol = 2)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:2)
    spliced <- counts + 1
    unspliced <- counts + 2
    velocity <- matrix(
        c(3, 4, 0, 0, 0, 12),
        nrow = 3,
        dimnames = dimnames(counts)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = counts,
            spliced = spliced,
            unspliced = unspliced,
            regvelo_velocity = velocity
        )
    )
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            artifacts = list(velocity_assay = "regvelo_velocity")
        )
    )

    input <- sclet_prepare_cellrank_velocity_input(
        sce = sce,
        vel_state = get_velocity(sce, id = "regvelo"),
        velocity_id = "regvelo"
    )

    expect_equal(input$velocity_source, "object")
    expect_equal(input$velocity_assay, "regvelo_velocity")
    expect_equal(input$cell_names, colnames(sce))
    expect_equal(input$gene_names, rownames(sce))
    expect_equal(dim(input$velocity), c(ncol(sce), nrow(sce)))
    expect_equal(as.matrix(input$velocity), t(velocity))
})

test_that("CellRank input preparation falls back to velocity result object", {
    counts <- matrix(0, nrow = 3, ncol = 2)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:2)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    velocity <- matrix(
        c(3, 4, 0, 0, 0, 12),
        nrow = 3,
        dimnames = dimnames(counts)
    )
    vel_res <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            velocity = velocity,
            spliced = counts + 1,
            unspliced = counts + 2
        )
    )
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "scvelo",
        active = TRUE,
        value = list(
            method = "velociraptor::scvelo",
            results = vel_res
        )
    )

    input <- sclet_prepare_cellrank_velocity_input(
        sce = sce,
        vel_state = get_velocity(sce, id = "scvelo"),
        velocity_id = "scvelo"
    )

    expect_equal(input$velocity_source, "results")
    expect_equal(input$velocity_assay, "velocity")
    expect_equal(input$cell_names, colnames(vel_res))
    expect_equal(input$gene_names, rownames(vel_res))
    expect_equal(as.matrix(input$velocity), t(velocity))
})

test_that("RunCellRank and RunCellFate expose reticulate backend arguments", {
    expect_true("backend" %in% names(formals(RunCellRank)))
    expect_true("python" %in% names(formals(RunCellRank)))
    expect_true("velocity_id" %in% names(formals(RunCellRank)))
    expect_true("backend" %in% names(formals(RunCellFate)))
    expect_true("python" %in% names(formals(RunCellFate)))
    expect_true("velocity_id" %in% names(formals(RunCellFate)))
})

test_that("CellRankSummary and VelocityFateCorrelation diagnose single-lineage results", {
    counts <- matrix(0, nrow = 3, ncol = 4)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:4)
    velocity <- matrix(
        c(3, 4, 0, 0, 0, 12, 1, 2, 3, 4, 5, 6),
        nrow = 3,
        dimnames = dimnames(counts)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts, regvelo_velocity = velocity)
    )
    sce$cellrank_terminal_states <- c("0", NA, "", "0")
    sce$cellrank_fate_X0 <- rep(1, 4)
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            artifacts = list(velocity_assay = "regvelo_velocity")
        )
    )
    sce <- sclet_set_state_record(
        sce,
        type = "trajectory",
        id = "cellrank",
        active = TRUE,
        value = list(
            method = "CellRank",
            inputs = list(velocity_id = "regvelo"),
            artifacts = list(
                terminal_state_col = "cellrank_terminal_states",
                fate_probability_cols = "cellrank_fate_X0",
                fate_probability_names = "0",
                has_lineage_drivers = FALSE
            )
        )
    )

    summary <- CellRankSummary(sce, id = "cellrank")
    expect_equal(summary$terminal_states$n_cells, c(2L, 2L))
    expect_equal(summary$fate_probabilities$min, 1)
    expect_equal(summary$fate_probabilities$max, 1)

    cor <- VelocityFateCorrelation(sce, velocity_id = "regvelo", trajectory_id = "cellrank")
    expect_equal(cor$fate_probability_col, "cellrank_fate_X0")
    expect_true(is.na(cor$estimate))
    expect_true(is.na(cor$p_value))
})

test_that("plot_velocity_fate_correlation returns ggplot", {
    skip_if_not_installed("ggplot2")

    counts <- matrix(0, nrow = 3, ncol = 4)
    rownames(counts) <- paste0("gene", 1:3)
    colnames(counts) <- paste0("cell", 1:4)
    velocity <- matrix(seq_len(12), nrow = 3, dimnames = dimnames(counts))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts, regvelo_velocity = velocity)
    )
    sce$cellrank_terminal_states <- c("0", "0", NA, "")
    sce$cellrank_fate_X0 <- seq(0.1, 0.9, length.out = 4)
    sce <- sclet_set_state_record(
        sce,
        type = "velocity",
        id = "regvelo",
        active = TRUE,
        value = list(
            method = "regvelo::REGVELOVI",
            artifacts = list(velocity_assay = "regvelo_velocity")
        )
    )
    sce <- sclet_set_state_record(
        sce,
        type = "trajectory",
        id = "cellrank",
        active = TRUE,
        value = list(
            method = "CellRank",
            artifacts = list(
                terminal_state_col = "cellrank_terminal_states",
                fate_probability_cols = "cellrank_fate_X0",
                fate_probability_names = "0"
            )
        )
    )

    p <- plot_velocity_fate_correlation(sce, velocity_id = "regvelo", trajectory_id = "cellrank")
    expect_s3_class(p, "ggplot")
})
