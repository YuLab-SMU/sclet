# sclet 0.0.9

+ `DimPlot` wrapper for `ggsc::sc_dim` (2025-06-02, Mon)
+ `RunSingleR` wrapper for automatic cell type annotation (2025-06-02, Mon)

# sclet 0.0.8

+ Fix `Idents` and `FindClusters` consistency (2025-05-30, Fri)
+ Fix `subset_cell` bug (2025-05-30, Fri)
+ Fix `gene_summary_table` bugs (2025-05-30, Fri)
+ Add Suggests package checks (2025-05-30, Fri)
+ Improve `NormalizeData` and `FindVariableFeatures` (2025-05-30, Fri)

# sclet 0.0.7

+ `RunSuperCell()` to merge similar cells to metacell (2025-05-28, Wed)

# sclet 0.0.6

+ `sc_merge()` to merge a named list of SingleCellExperiment objects to a single one (2025-04-17, Thu)
+ `BatchRemover()` for batch correction (2025-03-20, Thu)
+ `FindVariableFeatures()` and `VariableFeatures()` supports `method = 'scran'` (2025-03-17, Mon)

# sclet 0.0.5

+ Slingshot supports (2025-03-12, Wed)
    - genecurve_plot
    - RunSlingshot
    - lineage_plot
    - pseudo_heatmap
    - pseudo_plot

# sclet 0.0.4

+ re-export `scuttle::logNormCounts()` and `scuttle::perCellQCMetrics(()` (2024-10-08, Tue)
+ Import 'Matrix' for calculate stats of SparseMatrix (2024-10-07, Mon)
+ re-export `aplot::plot_list()` (2024-10-06, Sun)

# sclet 0.0.3

+ Adding cell-cell communication module. (2024-10-06, Sun)
    - `runCellChat()` wrapper function for running CellChat.
+ `gene_summary_table()` to add gene summary information to marker gene table (2024-10-02, Wed)

# sclet 0.0.2

+ `FindMarkers()` and `FindAllMarkers()` (2024-10-02, Wed)
+ remove vignettes and host online book (2024-10-01, Tue)
    - <https://yulab-smu.top/sclet>
+ re-export `scater::plotColData()` 
+ `QCMetrics()` function to add seurat-like QC metrics ('nFeature_RNA' and 'nCount_RNA' to 'SingleCellExperiment' object)

# sclet 0.0.1

+ Seurat-like functions for SingleCellExperiment (2024-09-20, Fri)
    - `Read10X()`
    - `PercentageFeatureSet()`
    - `VlnPlot()`
    - `FeatureScatter()`
    - `subset()`
    - `NormalizeData()`
    - `FindVariableFeatures()`
    - `VariableFeatures()`
    - `VariableFeaturePlot()`
    - `ScaleData()`
    - `ElbowPlot()`
    - `FindNeighbors()`
    - `FindClusters()`
    - `Idents()`
    - `RenameIdents()`
    - `RunUMAP()`
+ re-export:
    - `SummarizedExperiment::colData()`
    - `SummarizedExperiment::rowData()`
    - `SummarizedExperiment::assay()`
    - `scater::runPCA()`
