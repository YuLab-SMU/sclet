## Build a compact processed PBMC 3k SingleCellExperiment example for bookdown.
## Output: ../data/pbmc-3k-sce.rds
##
## Keep this object small enough for GitHub by:
## - preserving sparse counts instead of densifying with as.matrix()
## - downsampling cells for documentation examples
## - removing the dense scaled assay after PCA is computed
## - saving with xz compression

library(sclet)
library(scuttle)
library(TENxPBMCData)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)

set.seed(20241001)

pbmc <- TENxPBMCData::TENxPBMCData("pbmc3k")
rownames(pbmc) <- scuttle::uniquifyFeatureNames(rownames(pbmc), rowData(pbmc)$Symbol_TENx)
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

if (ncol(pbmc) > 800) {
    pbmc <- pbmc[, sample(seq_len(ncol(pbmc)), 800), drop = FALSE]
}

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, method = "scran")
var_features <- VariableFeatures(pbmc)
if (length(var_features) > 2000) {
    var_features <- var_features[seq_len(2000)]
}

pbmc <- ScaleData(pbmc)
pca <- stats::prcomp(
    t(as.matrix(SummarizedExperiment::assay(pbmc, "scaled")[var_features, , drop = FALSE])),
    center = FALSE,
    scale. = FALSE,
    rank. = 50
)
SingleCellExperiment::reducedDim(pbmc, "PCA") <- pca$x
pbmc <- sclet:::sclet_set_active_reduction(pbmc, "PCA")
SummarizedExperiment::assays(pbmc)$scaled <- NULL

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.88)
pbmc$ident <- as.character(SingleCellExperiment::colLabels(pbmc))
pbmc <- RunUMAP(pbmc, 1:10)

SummarizedExperiment::assays(pbmc) <- S4Vectors::SimpleList(lapply(
    SummarizedExperiment::assays(pbmc),
    function(x) Matrix::Matrix(as.matrix(x), sparse = TRUE)
))

S4Vectors::metadata(pbmc)$sclet_bookdown_cache <- list(
    source = "TENxPBMCData::pbmc3k",
    cells = ncol(pbmc),
    genes = nrow(pbmc),
    created_at = Sys.time()
)

saveRDS(pbmc, file = "../data/pbmc-3k-sce.rds", compress = "xz")
