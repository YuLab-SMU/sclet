## Prepare a small spliced/unspliced SCE for velocity & fate demos.
## Sources: scRNAseq::HermannSpermatogenesisData() (Bioconductor)
## Output: ../data/hermann-sperm-velo.rds
##
## The full dataset (2325 cells, ~12GB memory) is downsampled to 300 cells
## and 2000 genes for fast CI execution.

library(scRNAseq)
library(scuttle)
library(scran)
library(SingleCellExperiment)

set.seed(20241001)

sce <- HermannSpermatogenesisData()
sce <- sce[, 1:300]

sce <- logNormCounts(sce, assay.type = "spliced")
dec <- modelGeneVar(sce)
top_genes <- getTopHVGs(dec, n = 2000)
sce <- sce[top_genes, ]

dir.create("../data", showWarnings = FALSE)
saveRDS(sce, "../data/hermann-sperm-velo.rds")
