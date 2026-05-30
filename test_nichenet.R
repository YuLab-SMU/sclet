library(SingleCellExperiment)
devtools::load_all()

# Create mock SCE object
counts <- matrix(0, nrow = 10, ncol = 100)
rownames(counts) <- paste0("Gene", 1:10)
colnames(counts) <- paste0("Cell", 1:100)

# Make some cells express specific genes
# Sender cells (Cell1-50): express Gene1 (Ligand)
counts["Gene1", 1:50] <- 10
# Receiver cells (Cell51-100): express Gene2 (Receptor), Gene3 (Target)
counts["Gene2", 51:100] <- 10
counts["Gene3", 51:100] <- 10
counts["Gene4", 51:100] <- 10
counts["Gene5", 51:100] <- 10
counts["Gene6", 51:100] <- 10
counts["Gene7", 51:100] <- 10

sce <- SingleCellExperiment(assays = list(counts = counts))
sce <- NormalizeData(sce)

colData(sce)$CellType <- rep(c("Sender", "Receiver"), each = 50)

# Mock NicheNet prior models
# 1. lr_network: Gene1 -> Gene2
lr_network <- data.frame(
    from = "Gene1",
    to = "Gene2",
    source = "mock",
    database = "mock",
    stringsAsFactors = FALSE
)

# 2. weighted_networks: weights for lr_network
weighted_networks <- list(
    lr_sig = data.frame(
        from = "Gene1",
        to = "Gene2",
        weight = 1.0,
        stringsAsFactors = FALSE
    )
)

# 3. ligand_target_matrix: columns are ligands, rows are targets
ligand_target_matrix <- matrix(
    c(1.0, 0.1, 0.5, 0.2, 0.0),
    nrow = 5,
    ncol = 1,
    dimnames = list(c("Gene3", "Gene4", "Gene5", "Gene6", "Gene7"), "Gene1")
)

# Run NicheNet
sce <- RunCCI(
    sce,
    method = "NicheNet",
    group = "CellType",
    receiver = "Receiver",
    sender = "Sender",
    genes_oi = "Gene3",
    ligand_target_matrix = ligand_target_matrix,
    lr_network = lr_network,
    weighted_networks = weighted_networks,
    name = "mock_nichenet",
    top_n_ligands = 1,
    exprs_pct = 0.05
)

# Extract CCI
cci_df <- get_cci(sce)
print("NicheNet CCI Interactions:")
print(cci_df)

if (!is.null(cci_df) && nrow(cci_df) > 0) {
    print("NicheNet test passed!")
} else {
    print("NicheNet test failed!")
}
