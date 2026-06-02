## Prepare a small SpatialExperiment for spatial chapter demos.
## Source: STexampleData::Visium_humanDLPFC() (Bioconductor)
## Output: ../data/visium-dlpfc-sub.rds
##
## The human DLPFC Visium sample (sample 151673) has ~3500 spots.
## We subset to ~500 spots + 500 genes for fast CI execution.

library(STexampleData)
library(SpatialExperiment)
library(scuttle)

set.seed(20241001)

spe <- Visium_humanDLPFC()

keep_spots <- sample(seq_len(ncol(spe)), 500)
spe <- spe[, keep_spots]

keep_genes <- sample(seq_len(nrow(spe)), 500)
spe <- spe[keep_genes, ]

spe <- logNormCounts(spe)

dir.create("../data", showWarnings = FALSE)
saveRDS(spe, "../data/visium-dlpfc-sub.rds")
