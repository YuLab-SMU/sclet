## This script converts a Seurat object into a SingleCellExperiment (SCE) object.
## The raw data is from https://github.com/zhanghao-njmu/SCP.
##
## Input:
## a Seurat object (pancreas_sub) containing the raw count matrix, log-transformed count matrix,
##   cell metadata, and gene metadata.
##
## Output:
## a SingleCellExperiment object (sce) saved as an RDS file ("pancreas_sub_sce.rds").
##
## Import from SingleCellExperiment, Seurat

# load data
library(SingleCellExperiment)
library(Seurat)
rm(list = ls())
load("pancreas_sub.rda")
srt <- pancreas_sub

# extract counts_matrix and logcounts_matrix
counts_matrix <- as.matrix(srt@assays$RNA@counts)
logcounts_matrix <- as.matrix(srt@assays$RNA@data)

# extract cell_metadata
cell_metadata <- srt@meta.data

# extract gene_metadata
gene_metadata <- data.frame(gene_short_name = rownames(counts_matrix))
rownames(gene_metadata) <- rownames(counts_matrix)

# create sce object
sce <- SingleCellExperiment(
  assays = list(counts = counts_matrix, logcounts = logcounts_matrix),
  colData = cell_metadata,
  rowData = gene_metadata
)

# add reduction result
reducedDims(sce)$PCA <- as.matrix(Embeddings(srt, "PCA"))
reducedDims(sce)$UMAP <- as.matrix(Embeddings(srt, "UMAP"))

# add clustering result
colLabels(sce) <- colData(sce)$SubCellType

saveRDS(sce, file = "pancreas_sub_sce.rds")