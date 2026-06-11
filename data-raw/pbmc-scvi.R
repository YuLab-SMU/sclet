## Build a cached PBMC scVI integration result for bookdown.
## Output: ../data/rm_scvi_pbmc.rds
##
## The bookdown chapter only needs the scVI reduction and sclet integration
## metadata. Save a lightweight in-memory SCE instead of serializing the full
## HDF5-backed TENxPBMCData object.

library(sclet)
library(TENxPBMCData)
library(scuttle)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)

set.seed(20241001)

process_pbmc <- function(dataset) {
    sce <- TENxPBMCData::TENxPBMCData(dataset = dataset)
    rownames(sce) <- scuttle::uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol_TENx)
    sce <- QCMetrics(sce)
    sce[["percent.mt"]] <- PercentageFeatureSet(sce, "^MT-")
    keep <- colData(sce)$nFeature_RNA > 200 &
        colData(sce)$nFeature_RNA < 2500 &
        colData(sce)$percent.mt < 5
    sce <- sce[, keep, drop = FALSE]
    if (ncol(sce) > 500) {
        sce <- sce[, sample(seq_len(ncol(sce)), 500), drop = FALSE]
    }
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce)
    sce
}

pbmc3k <- process_pbmc("pbmc3k")
pbmc4k <- process_pbmc("pbmc4k")

combined_pbmc <- sce_merge(list(pbmc3k = pbmc3k, pbmc4k = pbmc4k))
rm_scvi_pbmc <- RunIntegration(combined_pbmc, method = "scVI", batch = "batch")

placeholder <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = c(1, ncol(rm_scvi_pbmc)),
    dimnames = list("placeholder", colnames(rm_scvi_pbmc))
)
light_scvi <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = placeholder),
    colData = S4Vectors::DataFrame(SummarizedExperiment::colData(rm_scvi_pbmc))
)
SingleCellExperiment::reducedDims(light_scvi) <- SingleCellExperiment::reducedDims(rm_scvi_pbmc)
S4Vectors::metadata(light_scvi) <- S4Vectors::metadata(rm_scvi_pbmc)
S4Vectors::metadata(light_scvi)$sclet_bookdown_cache <- list(
    source = "TENxPBMCData::pbmc3k + pbmc4k",
    cells = ncol(light_scvi),
    created_at = Sys.time()
)

dir.create("../data", showWarnings = FALSE)
saveRDS(light_scvi, "../data/rm_scvi_pbmc.rds", compress = "xz")
