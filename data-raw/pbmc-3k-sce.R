## Build a cached PBMC KNN reference mapping result for bookdown.
## Output: ../data/knn-refmap-pbmc.rds
##
## Keep this example lightweight by downsampling the PBMC object before the
## reference/query split. The rendered book still shows the full code pattern,
## but the cache refresh runs on a smaller demo subset so it finishes faster.

library(sclet)
library(TENxPBMCData)
library(scuttle)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)

set.seed(1)

pbmc <- TENxPBMCData::TENxPBMCData("pbmc3k")
rownames(pbmc) <- scuttle::uniquifyFeatureNames(rownames(pbmc), SummarizedExperiment::rowData(pbmc)$Symbol_TENx)
pbmc <- SingleCellExperiment(
    assays = list(counts = SummarizedExperiment::assay(pbmc, "counts")),
    rowData = S4Vectors::DataFrame(SummarizedExperiment::rowData(pbmc)),
    colData = S4Vectors::DataFrame(SummarizedExperiment::colData(pbmc))
)

pbmc <- QCMetrics(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, "^MT-")
pbmc <- subset_feature(pbmc, mincell = 3, peek = FALSE)
keep <- colData(pbmc)$nFeature_RNA > 200 &
    colData(pbmc)$nFeature_RNA < 2500 &
    colData(pbmc)$percent.mt < 5
pbmc <- pbmc[, keep, drop = FALSE]

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, method = "scran")
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, subset_row = VariableFeatures(pbmc), layer = "scaled")
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
cluster_label <- SingleCellExperiment::colLabels(pbmc)
if (is.null(cluster_label)) {
    cluster_label <- factor(seq_len(ncol(pbmc)) %% 5)
}
SingleCellExperiment::colLabels(pbmc) <- cluster_label

if (ncol(pbmc) > 400) {
    pbmc <- pbmc[, sample(seq_len(ncol(pbmc)), 400), drop = FALSE]
}

ref_n <- min(200, max(1, floor(ncol(pbmc) / 2)))
ref_idx <- sample(seq_len(ncol(pbmc)), ref_n)
ref_sce <- pbmc[, ref_idx]
query_sce <- pbmc[, -ref_idx]
ref_sce$label <- as.character(SingleCellExperiment::colLabels(ref_sce))

query_sce <- RunReferenceMapping(
    object = query_sce,
    ref = ref_sce,
    labels = "label",
    method = "KNN",
    layer = "logcounts",
    k = 5,
    name = "knn_demo"
)

cache_object <- query_sce
S4Vectors::metadata(cache_object)$sclet_bookdown_cache <- list(
    source = "TENxPBMCData::pbmc3k",
    ref_cells = ncol(ref_sce),
    query_cells = ncol(query_sce),
    created_at = Sys.time()
)

cache_object$knn_reference_labels <- ref_sce$label

dir.create("../data", showWarnings = FALSE)
saveRDS(cache_object, file = "../data/knn-refmap-pbmc.rds", compress = "xz")
