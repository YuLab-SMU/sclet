## Prepare a small SpatialExperiment for spatial chapter demos.
## Source: STexampleData::Visium_humanDLPFC() (Bioconductor)
## Output: ../data/visium-dlpfc-sub.rds
##
## The human DLPFC Visium sample is subset for fast documentation rendering.
## SVP colocalization/niche results are precomputed here so bookdown only reads
## the cached object.

library(sclet)
library(STexampleData)
library(SpatialExperiment)
library(scuttle)

set.seed(20241001)

spe <- STexampleData::Visium_humanDLPFC()

keep_spots <- sample(seq_len(ncol(spe)), 500)
spe <- spe[, keep_spots]

keep_genes <- sample(seq_len(nrow(spe)), 500)
spe <- spe[keep_genes, ]

spe <- spe[, Matrix::colSums(SingleCellExperiment::counts(spe)) > 0]
spe <- scuttle::logNormCounts(spe)

expr <- as.matrix(SummarizedExperiment::assay(spe, "logcounts"))
top_features <- rownames(spe)[order(apply(expr, 1, var), decreasing = TRUE)[1:5]]
for (f in top_features) {
    SummarizedExperiment::colData(spe)[[f]] <- expr[f, ]
}

if (requireNamespace("SVP", quietly = TRUE)) {
    spe <- RunSpatialColocalization(
        spe,
        features = top_features,
        method = "SVP",
        name = "coloc_dlpfc"
    )
    spe <- RunSpatialNiche(
        spe,
        features = top_features[1],
        method = "SVP",
        name = "niche_dlpfc"
    )
}

S4Vectors::metadata(spe)$sclet_bookdown_cache <- list(
    source = "STexampleData::Visium_humanDLPFC",
    top_features = top_features,
    created_at = Sys.time()
)

dir.create("../data", showWarnings = FALSE)
saveRDS(spe, "../data/visium-dlpfc-sub.rds")
