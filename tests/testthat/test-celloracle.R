test_that("RunInSilicoPerturbation validates inputs before calling basilisk", {
    counts <- matrix(rpois(20, 5), nrow = 4)
    rownames(counts) <- paste0("gene", 1:4)
    colnames(counts) <- paste0("cell", 1:5)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(rnorm(10), ncol = 2)

    expect_error(
        RunInSilicoPerturbation(sce, target_gene = "missing_gene", base_grn_path = "grn.csv"),
        "not found in sce"
    )
    expect_error(
        RunInSilicoPerturbation(sce, target_gene = "gene1", base_grn_path = "nonexistent_grn.csv"),
        "not found"
    )
})

test_that("get_perturbation reads perturbation state records", {
    counts <- matrix(0, nrow = 2, ncol = 3)
    rownames(counts) <- c("geneA", "geneB")
    colnames(counts) <- paste0("cell", 1:3)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SingleCellExperiment::reducedDim(sce, "co_shift_geneA") <- matrix(rnorm(6), ncol = 2)

    sce <- sclet_set_state_record(
        sce,
        type = "perturbation",
        id = "celloracle_geneA",
        active = TRUE,
        value = list(
            method = "celloracle::Oracle",
            inputs = list(reduction = "PCA"),
            artifacts = list(shift_reduction = "co_shift_geneA"),
            params = list(target_gene = "geneA", perturbation_value = 0),
            shift_reduction = "co_shift_geneA"
        )
    )

    rec <- get_perturbation(sce)
    expect_equal(rec$id, "celloracle_geneA")
    expect_equal(rec$params$target_gene, "geneA")

    expect_true(has_perturbation(sce))
    expect_equal(get_perturbation(sce, element = "shift_reduction"), "co_shift_geneA")
    expect_null(get_perturbation(sce, id = "celloracle_geneB"))
})

test_that("get_perturbation keeps the most recent record active", {
    counts <- matrix(0, nrow = 2, ncol = 3)
    rownames(counts) <- c("geneA", "geneB")
    colnames(counts) <- paste0("cell", 1:3)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    for (gene in c("geneA", "geneB")) {
        reduction <- paste0("co_shift_", gene)
        SingleCellExperiment::reducedDim(sce, reduction) <- matrix(rnorm(6), ncol = 2)
        sce <- sclet_set_state_record(
            sce,
            type = "perturbation",
            id = paste0("celloracle_", gene),
            active = TRUE,
            value = list(
                method = "celloracle::Oracle",
                params = list(target_gene = gene),
                artifacts = list(shift_reduction = reduction)
            )
        )
    }

    expect_equal(get_perturbation(sce)$id, "celloracle_geneB")
    expect_equal(get_perturbation(sce, id = "celloracle_geneA")$params$target_gene, "geneA")
})

test_that("get_perturbation falls back to legacy celloracle analysis", {
    counts <- matrix(0, nrow = 2, ncol = 3)
    rownames(counts) <- c("geneA", "geneB")
    colnames(counts) <- paste0("cell", 1:3)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    S4Vectors::metadata(sce)$sclet$analyses$celloracle <- list(
        method = "CellOracle",
        target_gene = "geneA",
        shift_reduction = "co_shift_geneA"
    )

    rec <- get_perturbation(sce)
    expect_equal(rec$method, "CellOracle")
    expect_equal(rec$target_gene, "geneA")
    expect_true(has_perturbation(sce))
})
