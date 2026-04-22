test_that("FindVariableFeatures stores HVG state under metadata sclet", {
    skip_if_not_installed("scuttle")

    set.seed(1)
    sce <- scuttle::mockSCE(ncells = 30, ngenes = 100)
    sce <- scuttle::logNormCounts(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")

    state <- S4Vectors::metadata(sce)$sclet

    expect_false(is.null(state))
    expect_equal(state$features$hvg$method, "scran")
    expect_equal(state$features$hvg$n, 20)
    expect_true(all(c("mean", "total", "tech", "bio", "p.value", "FDR") %in% state$features$hvg$rowData_cols))
})

test_that("FindNeighbors stores graph under metadata sclet and FindClusters can reuse it", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(2)
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 120)
    sce <- scuttle::logNormCounts(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30, method = "scran")
    sce <- runPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 10)
    sce <- FindNeighbors(sce, dims = 1:5)

    state <- S4Vectors::metadata(sce)$sclet

    expect_false(is.null(state$graphs$knn_graph$object))
    expect_equal(state$graphs$knn_graph$params$k, 10)

    sce <- FindClusters(sce)
    expect_false(is.null(SingleCellExperiment::colLabels(sce)))
})

test_that("subset SingleCellExperiment no longer depends on Barcode", {
    skip_if_not_installed("scuttle")

    sce <- scuttle::mockSCE(ncells = 12, ngenes = 50)
    sce$sample_group <- rep(c("A", "B"), each = 6)
    if ("Barcode" %in% colnames(SummarizedExperiment::colData(sce))) {
        SummarizedExperiment::colData(sce)$Barcode <- NULL
    }

    sub <- subset(sce, sample_group == "A")

    expect_equal(ncol(sub), 6)
    expect_true(all(sub$sample_group == "A"))
})

test_that("VariableFeatures remains compatible with legacy metadata fields", {
    skip_if_not_installed("scuttle")

    set.seed(3)
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 80)
    sce <- scuttle::logNormCounts(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 15, method = "scran")

    md <- S4Vectors::metadata(sce)
    md$sclet <- NULL
    md$hvgmethod <- "scran"
    md$nVariableFeatures <- 15
    md$hvgcols <- c("mean", "total", "tech", "bio", "p.value", "FDR")
    S4Vectors::metadata(sce) <- md

    expect_no_error({
        hvgs <- VariableFeatures(sce)
    })
    expect_equal(length(hvgs), 15)
})

test_that("DefaultReduction tracks active reduction and can be set explicitly", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(rnorm(20), ncol = 5)
    SingleCellExperiment::reducedDim(sce, "UMAP") <- matrix(rnorm(8), ncol = 2)

    expect_equal(DefaultReduction(sce), "UMAP")

    sce <- `DefaultReduction<-`(sce, "PCA")
    expect_equal(DefaultReduction(sce), "PCA")
    expect_equal(S4Vectors::metadata(sce)$sclet$active$reduction, "PCA")
})

test_that("DefaultAssay and DefaultGraph can be managed through active state", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(
            counts = matrix(1, nrow = 5, ncol = 4),
            logcounts = matrix(2, nrow = 5, ncol = 4),
            scaled = matrix(3, nrow = 5, ncol = 4)
        )
    )

    expect_equal(DefaultAssay(sce), "scaled")

    sce <- `DefaultAssay<-`(sce, "logcounts")
    expect_equal(DefaultAssay(sce), "logcounts")

    fake_graph <- structure(list(), class = "igraph")
    sce <- sclet_set_graph(sce, fake_graph, name = "custom_graph", active = FALSE)
    expect_equal(DefaultGraph(sce), "custom_graph")

    sce <- `DefaultGraph<-`(sce, "custom_graph")
    expect_equal(DefaultGraph(sce), "custom_graph")
})

test_that("ActiveIdent can switch between colLabels and colData fields", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 6))
    )
    SummarizedExperiment::colData(sce)$celltype <- rep(c("T", "B"), each = 3)
    SingleCellExperiment::colLabels(sce) <- factor(rep(c("1", "2"), each = 3))

    expect_equal(ActiveIdent(sce), "colLabels")

    sce <- `ActiveIdent<-`(sce, "celltype")
    expect_equal(ActiveIdent(sce), "celltype")
    expect_equal(unname(Idents(sce)), rep(c("T", "B"), each = 3))
})

test_that("CommandLog records preprocessing and dimred steps", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")

    set.seed(4)
    sce <- scuttle::mockSCE(ncells = 24, ngenes = 80)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")
    sce <- runPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)

    log_df <- CommandLog(sce)
    log_df_details <- CommandLog(sce, details = TRUE)
    raw_log <- S4Vectors::metadata(sce)$sclet$commands

    expect_true(nrow(log_df) >= 3)
    expect_equal(tail(log_df$command, 3), c("NormalizeData", "FindVariableFeatures", "RunPCA"))
    expect_true(all(c("params_summary", "outputs_summary") %in% colnames(log_df)))
    expect_true(all(c("params", "outputs") %in% colnames(log_df_details)))
    expect_type(log_df_details$params[[1]], "list")
    expect_type(log_df_details$outputs[[1]], "list")
    expect_equal(DefaultReduction(sce), "PCA")
    expect_equal(tail(vapply(raw_log, function(x) x$command, character(1)), 1), "RunPCA")
})

test_that("main pipeline updates active assay, graph and ident", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(5)
    sce <- scuttle::mockSCE(ncells = 30, ngenes = 100)
    sce <- NormalizeData(sce)
    expect_equal(DefaultAssay(sce), "logcounts")

    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")
    sce <- runPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    expect_equal(DefaultGraph(sce), "knn_graph")

    sce <- FindClusters(sce)
    expect_equal(ActiveIdent(sce), "colLabels")

    sce <- ScaleData(sce, features = rownames(sce)[1:10], assay = "logcounts")
    expect_equal(DefaultAssay(sce), "scaled")
})

test_that("RunUMAP can infer source reduction when dims are omitted", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")

    set.seed(6)
    sce <- scuttle::mockSCE(ncells = 30, ngenes = 100)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")
    sce <- runPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)
    sce <- RunUMAP(sce)

    expect_true("UMAP" %in% SingleCellExperiment::reducedDimNames(sce))
    expect_equal(DefaultReduction(sce), "UMAP")
})

test_that("FindMarkers prefers logcounts when active assay is scaled", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("presto")

    set.seed(7)
    sce <- scuttle::mockSCE(ncells = 24, ngenes = 60)
    sce$cluster <- rep(c("A", "B"), each = 12)
    SingleCellExperiment::colLabels(sce) <- factor(sce$cluster)
    sce <- NormalizeData(sce)
    sce <- ScaleData(sce, features = rownames(sce)[1:20], assay = "logcounts")

    expect_equal(DefaultAssay(sce), "scaled")
    expect_no_error({
        res <- FindMarkers(sce, ident.1 = "A")
    })
    expect_true(is.data.frame(res))
})

test_that("analysis getters return stored records and support legacy fallback", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )

    sce <- sclet_set_analysis(sce, "milo", list(da_results = data.frame(id = 1:2)))
    expect_equal(nrow(get_milo(sce, "da_results")), 2)

    sce <- sclet_set_analysis(sce, "cellchat", list(method = "CellChat", object = list(net = TRUE)))
    expect_true(has_cellchat(sce))
    expect_equal(get_cellchat(sce, "method"), "CellChat")
    expect_true(isTRUE(get_cellchat(sce, "object")$net))

    sce <- sclet_set_analysis(sce, "supercell", list(object = list(membership = 1:4)))
    expect_equal(get_supercell(sce, "object")$membership, 1:4)

    md <- S4Vectors::metadata(sce)
    md$sclet$analyses$trajectory <- NULL
    md$sclet$analyses$milo <- NULL
    md$sclet$analyses$supercell <- NULL
    md$slingshot_info <- structure(list(dummy = TRUE), class = "mock_slingshot")
    md$milo_results <- list(da_results = data.frame(id = 3:4))
    md$SuperCell <- list(membership = 5:8)
    S4Vectors::metadata(sce) <- md

    expect_s3_class(get_trajectory(sce)$dataset, "mock_slingshot")
    expect_true(has_trajectory(sce))
    expect_true(has_milo(sce))
    expect_true(has_supercell(sce))
    expect_equal(nrow(get_milo(sce, "da_results")), 2)
    expect_equal(get_supercell(sce, "object")$membership, 5:8)
})

test_that("get_hvg, get_batch and get_graph expose unified records", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")

    set.seed(8)
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 80)
    sce <- scuttle::logNormCounts(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 12, method = "scran")
    hvg <- get_hvg(sce)
    expect_true(has_hvg(sce))
    expect_equal(hvg$method, "scran")
    expect_equal(hvg$n, 12)
    expect_true(all(c("mean", "total", "tech", "bio", "p.value", "FDR") %in% colnames(hvg$rowData)))

    sce <- runPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    graph_info <- get_graph(sce)
    expect_true(has_graph(sce))
    expect_equal(graph_info$name, "knn_graph")
    expect_equal(graph_info$params$k, 10)
    expect_false(is.null(graph_info$object))

    sce <- sclet_set_analysis(sce, "batch", list(method = "batchelor::batchCorrect", hvg_n = 12))
    batch_info <- get_batch(sce)
    expect_true(has_batch(sce))
    expect_equal(batch_info$method, "batchelor::batchCorrect")
    expect_equal(batch_info$hvg_n, 12)
})

test_that("has_* accessors return FALSE when results are absent", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )

    expect_false(has_trajectory(sce))
    expect_false(has_cellchat(sce))
    expect_false(has_milo(sce))
    expect_false(has_supercell(sce))
    expect_false(has_batch(sce))
    expect_false(has_hvg(sce))
    expect_false(has_graph(sce))
})
