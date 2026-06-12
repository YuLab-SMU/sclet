test_that("FindVariableFeatures stores HVG state under metadata sclet", {
    skip_if_not_installed("scuttle")

    set.seed(1)
    sce <- scuttle::mockSCE(ncells = 30, ngenes = 100)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")

    state <- S4Vectors::metadata(sce)$sclet

    expect_false(is.null(state))
    expect_equal(state$features$hvg$method, "scran")
    expect_equal(state$features$hvg$n, 20)
    expect_true(all(c("mean", "total", "tech", "bio", "p.value", "FDR") %in% state$features$hvg$rowData_cols))
})

test_that("NormalizeData supports explicit non-count assays", {
    skip_if_not_installed("scuttle")

    set.seed(11)
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 60)
    spliced <- SummarizedExperiment::assay(sce, "counts")
    SummarizedExperiment::assays(sce) <- S4Vectors::SimpleList(spliced = spliced)

    sce <- NormalizeData(sce, assay = "spliced")

    expect_true("logcounts" %in% SummarizedExperiment::assayNames(sce))
    expect_equal(DefaultAssay(sce), "logcounts")

    state <- S4Vectors::metadata(sce)$sclet
    expect_equal(state$layers$spliced$assay, "spliced")
    expect_equal(state$commands[[length(state$commands)]]$params$assay, "spliced")
})

test_that("FindNeighbors stores graph under metadata sclet and FindClusters can reuse it", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(2)
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 120)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 10)
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
    sce <- NormalizeData(sce)
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

test_that("FindVariableFeatures supports optional scrapper backend", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scrapper")

    set.seed(31)
    sce <- scuttle::mockSCE(ncells = 24, ngenes = 90)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scrapper")

    expect_equal(get_hvg(sce, "method"), "scrapper")
    expect_true("hvg" %in% get_hvg(sce, "rowData_cols"))
    expect_true("residuals" %in% get_hvg(sce, "rowData_cols"))
    expect_true(any(SummarizedExperiment::rowData(sce)$hvg))

    hvgs <- VariableFeatures(sce, method = "scrapper")
    expect_true(length(hvgs) > 0)
    expect_true(all(hvgs %in% rownames(sce)))
})

test_that("VariableFeatures handles scrapper HVG tables without residuals", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    rownames(sce) <- paste0("g", 1:4)
    SummarizedExperiment::rowData(sce)$hvg <- c(TRUE, FALSE, TRUE, FALSE)
    sce <- sclet_set_hvg_state(
        sce,
        nfeatures = 2,
        method = "scrapper",
        hvgcols = "hvg"
    )

    hvgs <- VariableFeatures(sce, method = "scrapper")
    expect_equal(hvgs, c("g1", "g3"))
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

test_that("DefaultAssay, DefaultLayer and DefaultGraph can be managed through active state", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(
            counts = matrix(1, nrow = 5, ncol = 4),
            logcounts = matrix(2, nrow = 5, ncol = 4),
            scaled = matrix(3, nrow = 5, ncol = 4)
        )
    )

    expect_equal(DefaultAssay(sce), "scaled")
    expect_equal(DefaultLayer(sce), "scaled")
    expect_equal(Layers(sce), c("counts", "logcounts", "scaled"))
    expect_equal(unique(as.vector(LayerData(sce, "counts"))), 1)

    sce <- `DefaultAssay<-`(sce, "logcounts")
    expect_equal(DefaultAssay(sce), "logcounts")

    sce <- `DefaultLayer<-`(sce, "logcounts")
    expect_equal(DefaultLayer(sce), "logcounts")

    fake_graph <- structure(list(), class = "igraph")
    sce <- sclet_set_graph(sce, fake_graph, name = "custom_graph", active = FALSE)
    expect_equal(DefaultGraph(sce), "custom_graph")

    sce <- `DefaultGraph<-`(sce, "custom_graph")
    expect_equal(DefaultGraph(sce), "custom_graph")
})

test_that("legacy knn_graph is resolved through unified helper path", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    S4Vectors::metadata(sce)$knn_graph <- structure(list(), class = "igraph")

    expect_equal(DefaultGraph(sce), "knn_graph")
    expect_true(has_graph(sce, "knn_graph"))

    sce <- `DefaultGraph<-`(sce, "knn_graph")
    expect_equal(DefaultGraph(sce), "knn_graph")
    expect_equal(get_graph(sce, "knn_graph", "name"), "knn_graph")
    expect_false(is.null(get_graph(sce, "knn_graph", "object")))
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
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)

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
    expect_equal(DefaultLayer(sce), "logcounts")

    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)
    expect_equal(sclet_get_active_state(sce, "reduction"), "pca")
    expect_equal(sclet_get_state_record(sce, "reduction")$artifacts$reduction, "PCA")
    sce <- FindNeighbors(sce, dims = 1:5)
    expect_equal(DefaultGraph(sce), "knn_graph")
    expect_equal(sclet_get_active_state(sce, "graph"), "knn_graph")
    expect_equal(sclet_get_state_record(sce, "graph")$method, "bluster::makeSNNGraph")

    sce <- FindClusters(sce)
    expect_equal(ActiveIdent(sce), "colLabels")
    expect_equal(sclet_get_active_state(sce, "clustering"), "louvain_clusters")
    expect_equal(sclet_get_state_record(sce, "clustering")$artifacts$ident, "colLabels")

    sce <- ScaleData(sce, features = rownames(sce)[1:10], assay = "logcounts")
    expect_equal(DefaultAssay(sce), "scaled")
    expect_equal(DefaultLayer(sce), "scaled")
    expect_equal(sclet_get_layer(sce, "scaled")$assay, "scaled")
})

test_that("NormalizeData registers preprocess state", {
    skip_if_not_installed("scuttle")

    set.seed(9)
    sce <- scuttle::mockSCE(ncells = 18, ngenes = 50)
    sce <- NormalizeData(sce)

    expect_equal(sclet_get_active_state(sce, "preprocess"), "normalize_logcounts")
    preprocess_state <- sclet_get_state_record(sce, "preprocess")
    expect_equal(preprocess_state$method, "scuttle::calculateCPM + log1p")
    expect_equal(preprocess_state$inputs$assay, "counts")
    expect_equal(preprocess_state$artifacts$assay, "logcounts")
    expect_equal(sclet_get_active_layer(sce), "logcounts")
    expect_equal(sclet_get_layer(sce, "logcounts")$role, "normalized")
})

test_that("RunUMAP can infer source reduction when dims are omitted", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")

    set.seed(6)
    sce <- scuttle::mockSCE(ncells = 30, ngenes = 100)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)
    sce <- RunUMAP(sce)

    expect_true("UMAP" %in% SingleCellExperiment::reducedDimNames(sce))
    expect_equal(DefaultReduction(sce), "UMAP")
    expect_equal(sclet_get_state_record(sce, "reduction", "umap")$inputs$layer, "logcounts")
})

test_that("RunPCA resolves assay from corrected layer", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")

    set.seed(11)
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 60)
    SummarizedExperiment::assay(sce, "reconstructed") <- SummarizedExperiment::assay(sce, "counts")
    sce <- sclet_set_layer(sce, "corrected", assay = "reconstructed", role = "corrected")

    sce <- RunPCA(sce, layer = "corrected", ncomponents = 5)
    pca_state <- sclet_get_state_record(sce, "reduction", "pca")
    expect_equal(pca_state$inputs$layer, "corrected")
    expect_equal(pca_state$inputs$assay, "reconstructed")

    sce <- RunUMAP(sce)
    umap_state <- sclet_get_state_record(sce, "reduction", "umap")
    expect_equal(umap_state$inputs$reduction, "PCA")
    expect_equal(umap_state$inputs$layer, "corrected")
})

test_that("downstream states inherit integration provenance from corrected layer", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(12)
    sce <- scuttle::mockSCE(ncells = 30, ngenes = 80)
    SummarizedExperiment::assay(sce, "reconstructed") <- SummarizedExperiment::assay(sce, "counts")
    sce <- sclet_set_layer(sce, "corrected", assay = "reconstructed", role = "corrected")
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "integration",
        id = "batchcorrect",
        method = "batchelor::batchCorrect",
        inputs = list(layer = "corrected", assay = "reconstructed"),
        artifacts = list(layer = "corrected", assay = "reconstructed")
    )

    sce <- RunPCA(sce, layer = "corrected", ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce)
    sce <- RunUMAP(sce)

    pca_state <- sclet_get_state_record(sce, "reduction", "pca")
    graph_state <- sclet_get_state_record(sce, "graph", "knn_graph")
    cluster_state <- sclet_get_state_record(sce, "clustering", "louvain_clusters")
    umap_state <- sclet_get_state_record(sce, "reduction", "umap")

    expect_equal(pca_state$inputs$integration$id, "batchcorrect")
    expect_equal(graph_state$inputs$integration$id, "batchcorrect")
    expect_equal(cluster_state$inputs$integration$id, "batchcorrect")
    expect_equal(umap_state$inputs$integration$id, "batchcorrect")
    expect_equal(graph_state$inputs$layer, "corrected")
    expect_equal(cluster_state$inputs$reduction, "PCA")
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

test_that("FindMarkers handles delayed marker assay for presto", {
    skip_if_not_installed("DelayedArray")
    skip_if_not_installed("scuttle")
    skip_if_not_installed("presto")

    set.seed(8)
    sce <- scuttle::mockSCE(ncells = 24, ngenes = 60)
    sce$cluster <- rep(c("A", "B"), each = 12)
    SingleCellExperiment::colLabels(sce) <- factor(sce$cluster)
    sce <- NormalizeData(sce)
    SummarizedExperiment::assay(sce, "logcounts") <- DelayedArray::DelayedArray(
        SummarizedExperiment::assay(sce, "logcounts")
    )

    expect_no_error({
        res <- FindMarkers(sce, ident.1 = "A")
    })
    expect_true(is.data.frame(res))
})

test_that("expression source resolver prefers active layer and falls back from scaled", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(
            counts = matrix(1, nrow = 5, ncol = 4),
            logcounts = matrix(2, nrow = 5, ncol = 4),
            scaled = matrix(3, nrow = 5, ncol = 4),
            reconstructed = matrix(4, nrow = 5, ncol = 4)
        )
    )
    sce <- sclet_set_layer(sce, "counts", assay = "counts", role = "counts", active = FALSE)
    sce <- sclet_set_layer(sce, "logcounts", assay = "logcounts", role = "normalized", active = FALSE)
    sce <- sclet_set_layer(sce, "scaled", assay = "scaled", role = "scaled")
    sce <- sclet_set_layer(sce, "corrected", assay = "reconstructed", role = "corrected", active = FALSE)

    marker_source <- sclet_resolve_expression_source(
        object = sce,
        prefer_nonscaled = TRUE,
        context = "marker testing"
    )
    expect_equal(marker_source$layer, "logcounts")
    expect_equal(marker_source$assay, "logcounts")

    corrected_source <- sclet_resolve_expression_source(object = sce, layer = "corrected")
    expect_equal(corrected_source$layer, "corrected")
    expect_equal(corrected_source$assay, "reconstructed")

    expect_equal(resolve_marker_assay(sce, layer = "scaled"), "logcounts")
    expect_equal(resolve_marker_assay(sce, layer = "corrected"), "reconstructed")
})

test_that("analysis-state registry stores typed records and active selections", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )

    sce <- sclet_set_analysis_state(
        object = sce,
        type = "annotation",
        id = "singler_main",
        method = "SingleR",
        inputs = list(assay = "logcounts"),
        artifacts = list(labels_col = "SingleR_labels")
    )
    sce <- sclet_set_state_record(
        object = sce,
        type = "annotation",
        id = "singler_alt",
        active = FALSE,
        value = list(
            method = "SingleR",
            artifacts = list(labels_col = "SingleR_labels_alt")
        )
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "trajectory",
        id = "slingshot_main",
        method = "slingshot",
        inputs = list(reduction = "UMAP", group = "label"),
        artifacts = list(analysis_key = "trajectory")
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "communication",
        id = "cellchat_main",
        method = "CellChat",
        inputs = list(assay = "logcounts", group = "label"),
        artifacts = list(analysis_key = "cellchat"),
        active = FALSE
    )

    expect_equal(sclet_get_active_state(sce, "annotation"), "singler_main")
    expect_equal(sclet_get_state_record(sce, "annotation")$artifacts$labels_col, "SingleR_labels")
    expect_equal(length(sclet_get_state_records(sce, "annotation")), 2)

    trajectory_state <- sclet_get_state_record(sce, "trajectory")
    expect_equal(trajectory_state$method, "slingshot")
    expect_equal(trajectory_state$artifacts$analysis_key, "trajectory")
    expect_s3_class(trajectory_state$created_at, "POSIXct")
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "trajectory",
        id = "slingshot_alt",
        method = "slingshot",
        inputs = list(reduction = "PCA"),
        artifacts = list(analysis_key = "trajectory"),
        active = FALSE
    )
    expect_true(has_trajectory(sce))
    expect_true(has_trajectory(sce, id = "slingshot_alt"))
    expect_equal(get_trajectory(sce, id = "slingshot_alt")$inputs$reduction, "PCA")

    expect_null(sclet_get_active_state(sce, "communication"))
    sce <- sclet_set_active_state(sce, "communication", "cellchat_main")
    expect_equal(sclet_get_active_state(sce, "communication"), "cellchat_main")
    expect_equal(sclet_get_state_record(sce, "communication")$method, "CellChat")
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "communication",
        id = "cellchat_alt",
        method = "CellChat",
        inputs = list(assay = "counts", group = "sample"),
        artifacts = list(analysis_key = "cellchat"),
        active = FALSE
    )
    expect_true(has_cellchat(sce, id = "cellchat_alt"))
    expect_equal(get_cellchat(sce, id = "cellchat_alt")$inputs$group, "sample")
})

test_that("analysis getters return stored records and prefer current contract", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )

    sce <- sclet_set_analysis(sce, "milo", list(da_results = data.frame(id = 1:2)))
    expect_equal(nrow(get_milo(sce, "da_results")), 2)

    sce <- sclet_set_analysis(sce, "cellchat", list(method = "CellChat", object = list(net = TRUE)))
    expect_true(has_cellchat(sce))
    expect_equal(get_cellchat(sce, "method"), "CellChat")
    expect_true(isTRUE(get_cellchat(sce, "object")$net))
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "communication",
        id = "cellchat",
        method = "CellChat",
        inputs = list(layer = "logcounts", group = "label"),
        artifacts = list(analysis_key = "cellchat")
    )
    expect_equal(get_cellchat(sce, "inputs")$layer, "logcounts")

    sce <- sclet_set_analysis_state(
        object = sce,
        type = "milo",
        id = "milo_main",
        method = "miloR",
        inputs = list(reduction = "PCA", sample_col = "sample"),
        artifacts = list(
            analysis_key = "milo",
            da_results = data.frame(id = 11:12),
            design_df = data.frame(sample = c("s1", "s2")),
            formula = "~ 0 + condition",
            contrasts = "B - A"
        )
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "milo",
        id = "milo_alt",
        method = "miloR",
        inputs = list(reduction = "PCA", sample_col = "sample2"),
        artifacts = list(
            analysis_key = "milo",
            da_results = data.frame(id = 21:22)
        ),
        active = FALSE
    )
    expect_true(has_milo(sce))
    expect_true(has_milo(sce, id = "milo_alt"))
    expect_equal(get_milo(sce, "inputs")$sample_col, "sample")
    expect_equal(nrow(get_milo(sce, id = "milo_alt", element = "artifacts")$da_results), 2)

    sce <- sclet_set_analysis(sce, "supercell", list(object = list(membership = 1:4)))
    expect_equal(get_supercell(sce, "object")$membership, 1:4)
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "aggregation",
        id = "supercell",
        method = "SuperCell::SCimplify",
        inputs = list(layer = "logcounts", parent_n_cells = 4),
        artifacts = list(analysis_key = "supercell", size_col = "size")
    )
    expect_true(has_supercell(sce))
    expect_equal(get_supercell(sce, "inputs")$layer, "logcounts")
    expect_equal(get_supercell(sce, "artifacts")$analysis_key, "supercell")

    sce <- sclet_set_analysis_state(
        object = sce,
        type = "annotation",
        id = "singler",
        method = "SingleR",
        inputs = list(layer = "logcounts"),
        artifacts = list(labels_col = "SingleR_labels")
    )
    expect_true(has_annotation(sce))
    expect_equal(get_annotation(sce, "method"), "SingleR")
    expect_equal(get_annotation(sce, "inputs")$layer, "logcounts")
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "mapping",
        id = "singler",
        method = "SingleR",
        inputs = list(mode = "label_transfer", layer = "logcounts"),
        artifacts = list(labels_col = "SingleR_labels", mapping_type = "reference_mapping")
    )
    expect_true(has_mapping(sce))
    expect_equal(get_mapping(sce, "method"), "SingleR")
    expect_equal(get_mapping(sce, "inputs")$mode, "label_transfer")

    md <- S4Vectors::metadata(sce)
    md$sclet$analyses$trajectory <- NULL
    md$sclet$analyses$milo <- NULL
    md$sclet$analyses$supercell <- NULL
    md$sclet$states$records$trajectory <- NULL
    md$sclet$states$records$milo <- NULL
    md$sclet$states$records$annotation <- NULL
    md$sclet$states$records$mapping <- NULL
    md$sclet$states$records$aggregation <- NULL
    md$sclet$states$active$trajectory <- NULL
    md$sclet$states$active$milo <- NULL
    md$sclet$states$active$annotation <- NULL
    md$sclet$states$active$mapping <- NULL
    md$sclet$states$active$aggregation <- NULL
    md$slingshot_info <- structure(list(dummy = TRUE), class = "mock_slingshot")
    md$milo_results <- list(da_results = data.frame(id = 3:4))
    md$SuperCell <- list(membership = 5:8)
    md$batch_correction <- list(method = "batchelor::batchCorrect", hvg_n = 9)
    SummarizedExperiment::colData(sce)$SingleR_labels <- rep(c("A", "B"), length.out = ncol(sce))
    S4Vectors::metadata(sce) <- md

    expect_false(has_trajectory(sce))
    expect_false(has_milo(sce))
    expect_false(has_supercell(sce))
    expect_false(has_batch(sce))
    expect_null(get_trajectory(sce))
    expect_null(get_milo(sce))
    expect_null(get_supercell(sce))
    expect_null(get_batch(sce))
    expect_true(has_annotation(sce))
    expect_true(has_mapping(sce))
    expect_equal(get_annotation(sce, "method"), "SingleR")
    expect_equal(get_mapping(sce, "artifacts")$mapping_type, "reference_mapping")
})

test_that("RunSingleR records annotation provenance from corrected layer", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("SingleR")

    set.seed(13)
    sce <- scuttle::mockSCE(ncells = 18, ngenes = 50)
    sce <- NormalizeData(sce)
    SummarizedExperiment::assay(sce, "reconstructed") <- SummarizedExperiment::assay(sce, "logcounts")
    sce <- sclet_set_layer(sce, "corrected", assay = "reconstructed", role = "corrected")
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "integration",
        id = "batchcorrect",
        method = "batchelor::batchCorrect",
        inputs = list(layer = "corrected", assay = "reconstructed"),
        artifacts = list(layer = "corrected", assay = "reconstructed")
    )

    ref <- scuttle::mockSCE(ncells = 12, ngenes = 50)
    ref <- NormalizeData(ref)
    rownames(ref) <- rownames(sce)
    SummarizedExperiment::colData(ref)$label.main <- rep(c("T", "B"), each = 6)

    sce <- RunSingleR(sce, ref = ref, layer = "corrected")
    annotation_state <- get_annotation(sce)
    mapping_state <- get_mapping(sce, id = "singler")
    expect_equal(annotation_state$method, "SingleR")
    expect_equal(annotation_state$inputs$layer, "corrected")
    expect_equal(annotation_state$inputs$integration$id, "batchcorrect")
    expect_equal(annotation_state$artifacts$labels_col, "SingleR_labels")
    expect_equal(mapping_state$inputs$mode, "label_transfer")
    expect_equal(mapping_state$inputs$integration$id, "batchcorrect")
    expect_equal(mapping_state$artifacts$mapping_type, "reference_mapping")
    expect_true("SingleR_labels" %in% colnames(SummarizedExperiment::colData(sce)))
})

test_that("get_analysis_context summarizes active analysis records", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    skip_if_not_installed("SingleR")

    set.seed(14)
    sce <- scuttle::mockSCE(ncells = 18, ngenes = 50)
    sce <- NormalizeData(sce)
    SummarizedExperiment::assay(sce, "reconstructed") <- SummarizedExperiment::assay(sce, "logcounts")
    sce <- sclet_set_layer(sce, "corrected", assay = "reconstructed", role = "corrected")
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "integration",
        id = "batchcorrect",
        method = "batchelor::batchCorrect",
        inputs = list(layer = "corrected", assay = "reconstructed"),
        artifacts = list(layer = "corrected", assay = "reconstructed")
    )
    sce <- RunPCA(sce, layer = "corrected", ncomponents = 5)
    sce <- sclet_set_graph(sce, structure(list(), class = "igraph"), name = "knn_graph")
    SingleCellExperiment::colLabels(sce) <- factor(rep(c("a", "b"), each = 9))
    sce <- sclet_set_active_ident(sce, "colLabels")

    ref <- scuttle::mockSCE(ncells = 12, ngenes = 50)
    ref <- NormalizeData(ref)
    rownames(ref) <- rownames(sce)
    SummarizedExperiment::colData(ref)$label.main <- rep(c("T", "B"), each = 6)
    sce <- RunSingleR(sce, ref = ref, layer = "corrected")

    context <- get_analysis_context(sce)
    expect_equal(context$active$layer, "corrected")
    expect_equal(context$active$reduction, "PCA")
    expect_equal(context$active$graph, "knn_graph")
    expect_equal(context$active$ident, "colLabels")
    expect_equal(context$records$integration$method, "batchelor::batchCorrect")
    expect_equal(context$records$annotation$method, "SingleR")
    expect_equal(context$records$mapping$inputs$mode, "label_transfer")

    status <- Status(sce, details = TRUE)
    expect_equal(status$active$layer, "corrected")
    expect_equal(status$active$reduction, "PCA")
    expect_equal(status$last_command, "RunSingleR")
    expect_equal(status$n_commands, nrow(status$commands))
    expect_true("corrected" %in% status$available$layers)
    expect_true("integration" %in% status$available$analyses)
    expect_equal(status$records$annotation$method, "SingleR")
})

test_that("Status handles objects without recorded analysis history", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )

    status <- Status(sce)
    expect_equal(status$active$assay, "counts")
    expect_null(status$last_command)
    expect_equal(status$n_commands, 0)
    expect_equal(status$available$assays, "counts")
})

test_that("get_hvg, get_batch and get_graph expose unified records", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    skip_if_not_installed("scran")

    set.seed(8)
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 80)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 12, method = "scran")
    hvg <- get_hvg(sce)
    expect_true(has_hvg(sce))
    expect_equal(hvg$method, "scran")
    expect_equal(hvg$n, 12)
    expect_true(all(c("mean", "total", "tech", "bio", "p.value", "FDR") %in% colnames(hvg$rowData)))

    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    graph_info <- get_graph(sce)
    expect_true(has_graph(sce))
    expect_equal(graph_info$name, "knn_graph")
    expect_equal(graph_info$params$k, 10)
    expect_false(is.null(graph_info$object))

    sce <- sclet_set_analysis(sce, "batch", list(method = "batchelor::batchCorrect", hvg_n = 12))
    batch_info <- get_batch(sce)
    expect_true(has_batch(sce))
    expect_equal(batch_info$type, "batch")
    expect_equal(batch_info$id, "batchcorrect")
    expect_equal(batch_info$method, "batchelor::batchCorrect")
    expect_equal(batch_info$hvg_n, 12)
    expect_true(has_integration(sce))
    expect_equal(get_integration(sce)$type, "integration")
    expect_equal(get_integration(sce, "method"), "batchelor::batchCorrect")
})

test_that("get_integration can access non-active merge provenance records", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "integration",
        id = "merged_inputs",
        method = "sce_merge",
        inputs = list(objects = c("a", "b")),
        summary = list(n_inputs = 2),
        active = FALSE
    )

    expect_true(has_integration(sce, id = "merged_inputs"))
    expect_equal(get_integration(sce, id = "merged_inputs", element = "method"), "sce_merge")
    expect_equal(get_integration(sce, id = "merged_inputs")$summary$n_inputs, 2)
    expect_null(get_integration(sce))
})


test_that("has_* accessors return FALSE when results are absent", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )

    expect_false(has_trajectory(sce))
    expect_false(has_cellchat(sce))
    expect_false(has_milo(sce))
    expect_false(has_supercell(sce))
    expect_false(has_mapping(sce))
    expect_false(has_batch(sce))
    expect_false(has_hvg(sce))
    expect_false(has_graph(sce))
    expect_false(has_detest(sce))
    expect_false(has_enrichment(sce))
})

test_that("RunDEtest and RunEnrichment states are recorded correctly", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("presto")
    
    set.seed(42)
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- NormalizeData(sce)
    SingleCellExperiment::colLabels(sce) <- factor(rep(c("A", "B"), each = 10))
    sce <- sclet_set_active_ident(sce, "colLabels")
    
    sce <- RunDEtest(sce, ident.1 = "A", ident.2 = "B")
    
    expect_true(has_detest(sce))
    detest_state <- get_detest(sce)
    expect_equal(detest_state$method, "FindMarkers")
    expect_equal(detest_state$inputs$ident.1, "A")
    expect_true(is.data.frame(detest_state$artifacts$result))
    
    # We mock the result to ensure there is at least one significant marker
    detest_state$artifacts$result$pval <- 0.01
    detest_state$artifacts$result$avg_log2FC <- 1.0
    detest_state$artifacts$result$gene <- rep(c("TP53", "TNF", "EGFR", "IL6", "AKT1"), length.out = nrow(detest_state$artifacts$result))
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "detest",
        id = detest_state$id,
        method = detest_state$method,
        inputs = detest_state$inputs,
        artifacts = detest_state$artifacts,
        active = TRUE
    )
    
    # Let's just check if it runs without error if orgDb is missing
    # Actually, we can skip enrichment if org.Hs.eg.db is missing
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        sce <- RunEnrichment(sce, db = "GO", orgDb = "org.Hs.eg.db", keyType = "SYMBOL")
        expect_true(has_enrichment(sce))
        enrich_state <- get_enrichment(sce)
        expect_equal(enrich_state$method, "clusterProfiler::compareCluster")
        expect_equal(enrich_state$inputs$detest_id, detest_state$id)
    }
})
