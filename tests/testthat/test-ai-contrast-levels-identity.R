test_that("the clustering is not offered as a contrast column on the standard path", {
    # The red line: ActiveIdent is the result the analysis is meant to
    # DISCOVER, so it goes out only as anonymised group codes. `identity_levels`
    # keeps that promise; `contrast_levels` used to break it.
    #
    # `FindClusters()` sets the identity to the special string "colLabels",
    # while SingleCellExperiment stores the labels in a colData column called
    # `label`. Excluding by name compared "label" against "colLabels", which is
    # FALSE, so the exclusion never fired and the clustering was published
    # verbatim -- and a cluster-vs-cluster contrast naming the real labels
    # validated as TRUE. Measured 2026-10-22.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 60, ngenes = 200)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 50)
    sce <- RunPCA(sce, ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce, resolution = 1.0)

    block <- sclet_ai_build_contrast_levels(sce)
    offered <- vapply(block$columns, `[[`, character(1), "column")

    expect_false("label" %in% offered)
    expect_true("label" %in% block$withheld)
    # And the labels themselves appear nowhere in the block.
    flat <- paste(unlist(block), collapse = " ")
    expect_false(grepl("colLabels", flat, fixed = TRUE))
})

test_that("the identity exclusion does not over-withhold a real design column", {
    # The other half of the rule in AGENTS.md: design levels must stay
    # available, because a planner that cannot see them has to invent them,
    # which is what made one of twelve same-tree trials die with
    # `contrast_invalid`. A column whose partition merely has the same NUMBER of
    # levels as the clustering is a design column, not a copy of it.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 60, ngenes = 200)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 50)
    sce <- RunPCA(sce, ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce, resolution = 1.0)

    block <- sclet_ai_build_contrast_levels(sce)
    offered <- vapply(block$columns, `[[`, character(1), "column")
    expect_true("Treatment" %in% offered)
    expect_true("Mutation_Status" %in% offered)
})

test_that("a renamed copy of the clustering is withheld too", {
    # A `cell_type` derived from the clusters, or a `seurat_clusters` copy, is
    # the same disclosure under a different name. The test is the PARTITION,
    # not the label text.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 60, ngenes = 200)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 50)
    sce <- RunPCA(sce, ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce, resolution = 1.0)

    clusters <- as.character(SingleCellExperiment::colLabels(sce))
    cd <- SummarizedExperiment::colData(sce)
    cd$cell_type <- paste0("CLUSTER_", clusters, "_SIGNATURE")
    SummarizedExperiment::colData(sce) <- cd

    block <- sclet_ai_build_contrast_levels(sce)
    offered <- vapply(block$columns, `[[`, character(1), "column")
    expect_false("cell_type" %in% offered)
    expect_true("cell_type" %in% block$withheld)
    expect_false(grepl("CLUSTER_", paste(unlist(block), collapse = " "), fixed = TRUE))
})

test_that("sclet_ai_same_partition compares partitions, not labels", {
    # Same grouping, different labels -> the same partition.
    expect_true(sclet_ai_same_partition(c("1", "2", "1", "2"), c("T", "B", "T", "B")))
    # Same grouping written in a different order -> still the same partition.
    expect_true(sclet_ai_same_partition(c("a", "b", "a", "b"), c("2", "1", "2", "1")))
    # A different grouping with the same number of levels -> NOT the same.
    expect_false(sclet_ai_same_partition(c("a", "a", "b", "b"), c("1", "2", "1", "2")))
    # Different level counts -> NOT the same.
    expect_false(sclet_ai_same_partition(c("a", "b", "c"), c("1", "2", "2")))
    # A single group is not a partition worth excluding.
    expect_false(sclet_ai_same_partition(c("a", "a"), c("1", "1")))
    # Missing values are dropped pairwise.
    expect_true(sclet_ai_same_partition(c("a", NA, "b"), c("1", "1", "2")))
})

test_that("the withheld columns are named, so nothing vanishes silently", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 60, ngenes = 200)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 50)
    sce <- RunPCA(sce, ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce, resolution = 1.0)
    # Precondition: the clustering has to have more than one level, or there is
    # nothing to disclose and `label` is dropped by the `n_levels < 2` rule
    # rather than withheld. A fixture that clusters into one group cannot test
    # this, which is why the first version of this test failed.
    skip_if(length(unique(as.character(SingleCellExperiment::colLabels(sce)))) < 2L,
            "the mock clustering produced a single group")
    block <- sclet_ai_build_contrast_levels(sce)
    expect_true(length(block$withheld) >= 1L)
    expect_true(all(nzchar(block$withheld)))
})
