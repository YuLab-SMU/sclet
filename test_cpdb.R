devtools::load_all()
library(SingleCellExperiment)
library(SeuratObject)

# Let's create a small dummy SCE object
mat <- matrix(rpois(1000, lambda = 1), nrow = 50)
rownames(mat) <- paste0("GENE", 1:50)
# To make CPDB work, we might need some real gene names, let's use some known pairs if possible.
# Actually, CellPhoneDB maps by gene symbol.
# Known pairs: IL2 - IL2RA, CD40 - CD40LG
rownames(mat)[1:4] <- c("IL2", "IL2RA", "CD40", "CD40LG")
colnames(mat) <- paste0("CELL", 1:20)

sce <- SingleCellExperiment(assays = list(counts = mat))
# CPDB expects normalized counts usually
logcounts(sce) <- log1p(mat)

# cell type
colData(sce)$cell_type <- rep(c("T_cell", "B_cell"), each = 10)

# Make some fake expression for IL2 in T cells and IL2RA in B cells to ensure a p-value < 0.05
# Wait, for a small matrix it might not find significant pairs, but let's try.
mat["IL2", 1:10] <- 10
mat["IL2RA", 11:20] <- 10
mat["CD40", 11:20] <- 10
mat["CD40LG", 1:10] <- 10
assay(sce, "counts") <- mat
logcounts(sce) <- log1p(mat)

# Run CPDB
message("Running CellPhoneDB...")
sce <- RunCCI(sce, method = "CellPhoneDB", group = "cell_type", species = "human", threads = 2, iterations = 10)

message("Interactions:")
cci_res <- get_cci(sce)
print(head(cci_res))
