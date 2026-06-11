## Build a cached PBMC scVI integration result for bookdown.
## Output: ../data/rm_scvi_pbmc.rds
##
## The bookdown chapter uses this cached object instead of training scVI during render.

library(sclet)
library(TENxPBMCData)
library(scuttle)

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
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce)
    sce
}

pbmc3k <- process_pbmc("pbmc3k")
pbmc4k <- process_pbmc("pbmc4k")

combined_pbmc <- sce_merge(list(pbmc3k = pbmc3k, pbmc4k = pbmc4k))
rm_scvi_pbmc <- RunIntegration(combined_pbmc, method = "scVI", batch = "batch")

S4Vectors::metadata(rm_scvi_pbmc)$sclet_bookdown_cache <- list(
    source = "TENxPBMCData::pbmc3k + pbmc4k",
    created_at = Sys.time()
)

dir.create("../data", showWarnings = FALSE)
saveRDS(rm_scvi_pbmc, "../data/rm_scvi_pbmc.rds")
