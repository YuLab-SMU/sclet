## Build a processed PBMC 3k SingleCellExperiment example for bookdown.
## Output: ../data/pbmc-3k-sce.rds

library(sclet)
library(scuttle)
library(TENxPBMCData)
library(SingleCellExperiment)
library(SummarizedExperiment)

set.seed(20241001)

pbmc <- TENxPBMCData::TENxPBMCData("pbmc3k")
rownames(pbmc) <- scuttle::uniquifyFeatureNames(rownames(pbmc), rowData(pbmc)$Symbol_TENx)
pbmc <- SingleCellExperiment(
    assays = list(counts = as.matrix(SummarizedExperiment::assay(pbmc, "counts"))),
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
pca <- stats::prcomp(
    t(as.matrix(SummarizedExperiment::assay(pbmc, "scaled"))),
    center = FALSE,
    scale. = FALSE,
    rank. = 50
)
SingleCellExperiment::reducedDim(pbmc, "PCA") <- pca$x
pbmc <- sclet:::sclet_set_active_reduction(pbmc, "PCA")
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.88)
pbmc$ident <- as.character(SingleCellExperiment::colLabels(pbmc))
pbmc <- RunUMAP(pbmc, 1:10)

saveRDS(pbmc, file = "../data/pbmc-3k-sce.rds")
