
# sclet 0.0.2.002

* Adding cell-cell communication module.
  * runCellChat( ) wrapper function for running CellChat.
* add demo data for CellChat.
  * hmskin_sce: singlecellexperiment object, derived from cellchat demo data.
  * hmskin.r show the process of preparation.

# sclet 0.0.2.001

+ `gene_summary_table()` to add gene summary information to marker gene table (2024-10-02, Wed)

# sclet 0.0.2

+ `FindMarkers()` and `FindAllMarkers()` (2024-10-02, Wed)
+ remove vignettes and host online book (2024-10-01, Tue)
  - [https://yulab-smu.top/sclet](https://yulab-smu.top/sclet)
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
