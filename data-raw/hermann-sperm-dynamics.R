## Build cached velocity and CellRank results for bookdown.
## Input:  ../data/hermann-sperm-velo.rds
## Output: ../data/hermann-sperm-dynamics.rds
##
## This script intentionally lives outside the main bookdown build because it
## provisions Python-backed basilisk environments through velociraptor/CellRank.

library(sclet)
library(S4Vectors)

input <- "../data/hermann-sperm-velo.rds"
if (!file.exists(input)) {
    source("hermann-sperm-velo.R")
}

sce <- readRDS(input)
top_genes <- S4Vectors::metadata(sce)$sclet_bookdown_cache$top_genes
if (is.null(top_genes)) {
    top_genes <- rownames(sce)
}

sce <- RunVelocity(sce, subset_row = top_genes, assay.X = "spliced", use.dimred = "PCA")
sce <- RunCellFate(sce, reduction = "PCA", cluster_key = "cluster", name = "cellrank")

dir.create("../data", showWarnings = FALSE)
saveRDS(sce, "../data/hermann-sperm-dynamics.rds")
