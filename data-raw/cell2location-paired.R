## Prepare paired spatial + reference data for cell2location CI demo.
## This script requires matching tissue types. For a minimal demo we use
## a subset of the DLPFC Visium as spatial and a same-tissue scRNA-seq 
## reference. If no matching reference is available, cell2location chunks
## will skip in CI.
##
## Output: ../data/dlpfpc-spatial.rds + ../data/dlpfpc-ref.rds
##
## NOTE: This script creates placeholder data. Real cell2location analysis
## requires matched tissue scRNA-seq reference with cell type labels.
## For CI purposes, both objects share the same gene space but the reference
## labels are synthetic.

library(STexampleData)
library(SpatialExperiment)
library(scuttle)

set.seed(20241001)

spe <- Visium_humanDLPFC()
spe <- spe[, 1:300]
spe <- logNormCounts(spe)

keep_genes <- sample(seq_len(nrow(spe)), 200)
spe <- spe[keep_genes, ]

ref_sce <- as(spe, "SingleCellExperiment")
colData(ref_sce)$cell_type <- sample(
    c("Neuron", "Astrocyte", "Oligodendrocyte", "Microglia", "Endothelial"),
    ncol(ref_sce), replace = TRUE
)

dir.create("../data", showWarnings = FALSE)
saveRDS(spe, "../data/dlpfpc-spatial.rds")
saveRDS(ref_sce, "../data/dlpfpc-ref.rds")

message("Paired cell2location data prepared.")
message("To run cell2location in CI: set SCLET_BOOKDOWN_EVAL=true")
message("and ensure basilisk environments are cached.")
