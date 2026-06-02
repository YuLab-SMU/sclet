# sclet 1.0.0



# sclet 0.99.7

+ **Data Cleaning**: Introduced `RunDoubletFinder()` to detect doublets using `scDblFinder`, and `RunDecontX()` to remove ambient RNA contamination using `celda::decontX`. Both functions natively support the `SingleCellExperiment` object and integrate seamlessly with the `analysis-state` architecture, storing clean matrices safely as new assays (e.g., `decontXcounts`) without overwriting raw data.
+ **Imputation**: Added `RunImputation()` to perform zero-preserving imputation using the `ALRA` method. The function handles sparse matrices and integrates with the analysis-state machine by storing the imputed results in a separate assay (`alra_imputed`).
+ **Diffusion Map**: Added `RunDiffusionMap()` as a robust dimensional reduction alternative for complex trajectory inference. It uses the `destiny` package to calculate diffusion components, storing them in `reducedDim(sce, "DM")`.
+ **Integration**: Refined `RunIntegration(method = "scVI")` to avoid dense matrix materialization, converting sparse matrix inputs directly to `scipy.sparse.csr_matrix` via `reticulate` for improved performance and reduced memory usage on large datasets.
+ **Memory Safety**: Audited other large-matrix workflows and now keep sparse representations as long as possible in `RunSpatialDeconvolution()`, `RunSCENIC()`, `RunCellRank()`, `RunInSilicoPerturbation()`, `RunKNNPredict()`, while also avoiding full-assay transpose in `pseudo_heatmap()` and `genecurve_plot()` unless only the requested features are needed.
+ **Phenotype Association**: Added `RunPhenotypeAssociation()` to directly link single-cell transcriptomics with bulk clinical data (e.g., survival analysis) using the `scPAS` package. Supports Cox regression, binomial, and gaussian models, tracking the association results and parameters directly in the `analysis-state`.
+ **Fate Mapping Workflow**: Refined the `velocity -> fate` chain around `SingleCellExperiment` state records. `RunVelocity()` now registers velocity outputs through the standard `analysis-state` contract, `RunCellRank()` stores terminal states, named fate probabilities, and lineage driver statistics in a structured way, and `RunCellFate()` provides a thin user-facing entry point that currently dispatches to `CellRank`.
+ **Fate Visualization**: Added native fate plotting helpers following the `plot_<module>_<type>` convention: `plot_fate_probability()`, `plot_fate_terminal_states()`, and `plot_fate_driver_trends()`. These functions consume the stored CellRank artifacts directly and expose fate probabilities, terminal state assignments, and lineage driver expression trends on top of the tracked reductions/layers.
+ **Mainline Workflow Shells**: Introduced three high-level `Run<Mainline>Workflow()` entry points as semantic shells over existing atomic interfaces: `RunTrajectoryWorkflow()` (diffusion map + Slingshot + velocity + CellRank/CellFate), `RunProgramWorkflow()` (gene-set scoring + SCENIC regulon inference), and `RunReferenceWorkflow()` (SingleR/KNN reference mapping). These shells are designed to make the recommended analysis path easier to follow and document without hiding the underlying methods.
+ **Reference Mapping Layer**: Added `RunReferenceMapping()` as a unified entry point for reference-based annotation, dispatching to `SingleR` or lightweight `KNN` label transfer with a consistent state contract and accessor layer. Added `plot_reference_label_transfer_heatmap()` and `plot_reference_label_confidence()` for downstream inspection of label transfer quality across query groups.
+ **Program & Regulon Activity Visualization**: Added `plot_program_activity()` and `plot_program_heatmap()` to unify the visualization layer across `RunGeneSetScoring()` and `RunSCENIC()` outputs. Both consumers resolve program-level activity from the typed state records, making the source (geneset scoring vs. regulon) transparent to the plotting call. Added `plot_scenic_activity()` for regulon-level embedding visualization, and `plot_trajectory_overview()` as the first native overview plot for the trajectory mainline.
+ **State Contract Refinements**: `RunGeneSetScoring()` now writes a typed `geneset_scoring` state record; `RunSCENIC()`, `RunVelocity()`, `RunSingleR()`, and `RunCellRank()` all accept a `name` parameter and register results through the unified analysis-state layer instead of ad hoc metadata slots. `get_geneset_scoring()` now reads both new typed records and legacy formats transparently.
+ **State Priority & Perturbation Sensitivity Mainline**: Introduced `RunPerturbationPriority()` wrapping the Augur framework for ranking cell types by perturbation responsiveness, `RunRareCellDetection()` as a reserved interface for rare-cell identification, and `RunStatePriorityWorkflow()` as the semantic shell over the state-priority analysis mainline. Added `get_perturbation_priority()` / `has_perturbation_priority()` and `get_rare_cells()` / `has_rare_cells()` accessors.

# sclet 0.99.6

+ **Pathway Scoring & Marker Detection**: Enhanced `RunGeneSetScoring()` with `as_altExp = TRUE` to store pathway scores as an alternative experiment. Added `FindMarkerPathways()` to detect cell-type specific pathways robustly using `FindAllMarkers`.
+ **KEGG Integration**: Added `RunKEGG()` as a one-stop wrapper to download KEGG pathways via `clusterProfiler`, score cells, and perform differential marker pathway analysis seamlessly.
+ **Cell-Cell Communication**: Introduced `RunCCI()` as a unified interface for cell-cell communication analysis, seamlessly supporting `CellChat`, `CellPhoneDB`, and `NicheNet` backends. This eliminates the need for manual data formatting and environment setup (using `basilisk` for Python backends), storing results in a standardized format within the analysis state.
+ **Cell-Cell Communication Visualization**: Added two native CCI visualization methods: `plot_cci_sigmoid` (S-curve connection plot) and `plot_cci_arrow` (bidirectional arrow plot with dynamic expression mapping). These native functions eliminate the need for heavy third-party visualization dependencies.
+ **Unified Plotting Names**: Established a `plot_<module>_<type>` naming convention for visualizations (e.g., `plot_cci_bubble`, `plot_velocity`, `plot_lineage`), while retaining classic aliases like `DimPlot` for Seurat users.

# sclet 0.99.5

+ **AI Copilot & Evidence Governance**: Introduced `sclet_copilot()` as an intelligent diagnostic agent. Powered by the `aisdk` framework, it reads the complete `analysis-state` provenance of the `SingleCellExperiment` object to provide context-aware diagnostics.
+ **Cross-Chain Auditing**: Added `AuditAnalysisChain()` to perform rigorous cross-chain error control. It calculates a "State-Dependency Confidence Score" to help distinguish true biological signals from computational artifacts (e.g., over-alignment during integration), actively combating "Transcriptomic Overload".
+ **Spatial Deconvolution**: Added `RunSpatialDeconvolution()` to leverage the `cell2location` deep learning model via `basilisk`. It maps single-cell references onto spatial transcriptomics spots for high-resolution abundance inference.
+ **In-Silico Perturbation**: Added `RunInSilicoPerturbation()` to wrap the `celloracle` framework. This allows for predictive simulation of cell fate trajectory shifts following virtual gene knockouts.
+ **Advanced Fate Mapping**: Added `RunCellRank()` to combine directed RNA velocity with transcriptomic similarity, enabling robust inference of terminal states and fate probabilities for complex developmental processes.

# sclet 0.99.4

+ **Python Interoperability**: Introduced `basilisk` infrastructure (`sclet_env`) to seamlessly manage Python environments and dependencies, enabling robust integration with the Python single-cell ecosystem without polluting user's local environments.
+ **RNA Velocity**: Added `RunVelocity()` and `VelocityPlot()` to support RNA velocity analysis via `velociraptor` (`scVelo`), with results seamlessly integrated into the `velocity` analysis state.
+ **Gene Regulatory Networks**: Added `RunSCENIC()` to run the complete `pySCENIC` pipeline (GRN inference, motif enrichment, and AUCell scoring) via `basilisk`. The resulting AUC matrix is stored as an `altExp` and tracked in the `scenic` analysis state.
+ **Gene Set Scoring**: Added `RunGeneSetScoring()` to calculate gene set or pathway activity scores using `UCell`, `AUCell`, or `GSVA`. Scores are added to `colData` and registered in the `geneset_scoring` analysis state.
+ **State Management**: Updated `Status()` and accessors (`get_velocity()`, `get_scenic()`, `get_geneset_scoring()`) to fully support the newly added high-level analysis modules.

# sclet 0.99.3

+ `analysis-state contract`: Introduced a unified, lightweight analysis-state registry to manage downstream analysis states (annotation, integration, mapping, trajectory, communication, detest, enrichment) and their dependencies.
+ **Integration Workflow**: Added `RunIntegration()` as a unified entry point for multi-backend data integration, currently supporting `fastMNN` and `Harmony`. Downstream functions like `FindNeighbors` and `RunUMAP` can now automatically perceive and consume the integration states.
+ **Reference Mapping**: Added `RunKNNPredict()` and `RunReferenceMapping()` for lightweight label transfer based on `BiocNeighbors::queryKNN`. Added `ProjectionPlot()` for query-reference visualization.
+ **DE & Enrichment**: Added `RunDEtest()` (wrapping `FindMarkers`/`FindAllMarkers`) and `RunEnrichment()` (wrapping `clusterProfiler`). Included corresponding visualization tools: `DEtestPlot()` (volcano plot) and `EnrichmentPlot()` (dotplot).
+ **Visualization Ecosystem**: Enhanced visualization capabilities by leveraging the `ggsc` package. Added `CellDimPlot()`, `FeatureDimPlot()`, `GroupHeatmap()`, and `CellStatPlot()`.
+ **Technical Debt & Bug Fixes**:
  - Migrated to `scrapper::aggregateAcrossCells.se` for pseudobulk aggregation and `scrapper` for HVG selection to eliminate deprecated warnings from `scuttle` and `scran`.
  - Fixed missing dependencies (`enrichplot`, `scran`) and non-ASCII character warnings, achieving a clean `R CMD check` (0 errors, 0 warnings, 0 notes).

# sclet 0.99.2

+ Performance improvement: `FindAllMarkers()` uses `presto` for faster calculation (2026-02-06, Fri)
+ Bug fix: `BatchRemover()` robustness improvements (2026-02-06, Fri)
+ Internal: Renamed source files from `.r` to `.R` (2026-02-06, Fri)
+ `RunMilo` and `refit_milo` for Milo differential abundance analysis (2026-01-30, Fri)

# sclet 0.99.1

+ `RunSuperCell` for metacell construction (2026-01-16, Fri)
+ `RunSlingshot` for trajectory inference
+ `RunEnrichment` for functional enrichment analysis
+ `RunExplorer` for interactive data exploration using iSEE

# sclet 0.99.0

+ `ReadH5AD` and `WriteH5AD` for AnnData interoperability (2025-06-04, Wed)
+ `as.Seurat` and `as.SCE` for object conversion (2025-06-04, Wed)
+ Performance improvements: support DelayedMatrix in normalization (2025-06-04, Wed)

# sclet 0.1.0

+ `AggregateExpression` for pseudobulk analysis (2025-06-03, Tue)
+ `FindMarkers` supports pseudobulk DE via `DESeq2` (2025-06-03, Tue)
+ `BatchRemover` now records parameters in metadata (2025-06-03, Tue)

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
    - `RunCellChat()` wrapper function for running CellChat.
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
    - `RunPCA()`
+ re-export:
    - `SummarizedExperiment::colData()`
    - `SummarizedExperiment::rowData()`
    - `SummarizedExperiment::assay()`
