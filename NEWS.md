# sclet 1.0.1

+ **BiocNeighbors Compatibility Fix**: `RunKNNPredict()` (and thus `RunReferenceMapping(method = "KNN")`) no longer fails with `unable to find an inherited method for function 'queryKNN' for signature 'X = "dgCMatrix"'` on BiocNeighbors 2.0.x-2.4.x (Bioconductor 3.20-3.22, R 4.4/4.5), as reported in [#25](https://github.com/YuLab-SMU/sclet/issues/25). Matrices reaching `BiocNeighbors::queryKNN()`/`findKNN()` now pass through a dispatch-aware helper that keeps sparse inputs sparse on BiocNeighbors >= 2.6 and transparently falls back to ordinary matrices where sparse dispatch is unavailable, so no user-side version management is needed. Also fixed the `RunSymphonyMapping()` install hint: symphony has been removed from CRAN, so the message now points to `remotes::install_github('immunogenomics/symphony')`.
+ **fastMNN / `BatchRemover` semantic alignment**: `BatchRemover()` (and thus `RunIntegration(method = "fastMNN")`) now graft the correction artifacts back onto the *original* object instead of replacing it with the stripped batchelor output. This fixes three issues reported in [#24](https://github.com/YuLab-SMU/sclet/issues/24):
  + The `corrected` embedding (fastMNN's intended clustering/UMAP artifact) is kept in `reducedDim("corrected")`, registered as a `reduction` state, and made the active reduction, so `RunUMAP` / `FindNeighbors` consume it directly.
  + The original `counts`/`logcounts` and **all genes** are preserved, so marker detection and differential expression work again after integration.
  + The "corrected" layer is no longer backed by the low-rank `reconstructed` matrix (which batchelor marks *not for quantitative analysis* and which spans only the HVGs when `correct.all = FALSE`). It now points to the full gene-resolution `reconstructed` assay when `correct.all = TRUE`, otherwise to the original source expression; the integration record reports `corrected_expression` and `embedding_only` flags so downstream tools can distinguish an expression-corrected integration from an embedding-only one.
+ **In-Silico Perturbation State Contract**: `RunInSilicoPerturbation()` now registers results through the unified `analysis-state` contract: each run creates a typed `perturbation` state record (id `celloracle_<gene>`) with inputs/artifacts/params plus a command log entry, while the shift reduction stays in `reducedDim`. `get_perturbation()` reads typed records by id, keeps the most recent record active across multiple target genes, and falls back to legacy metadata for old objects.
+ **RegVelo GPU Backend**: Added `RunRegVelo()` as an RNA velocity backend that runs RegVelo through either the packaged Python stack or a user-managed `reticulate` Python, stores the inferred velocity as a `SingleCellExperiment` assay, and registers it in the standard velocity state layer.
+ **RegVelo to CellRank Workflow**: Extended `RunCellRank()` so CellRank can consume velocity assays produced by `RunRegVelo()`, including the GPU/SLURM-friendly `reticulate` path. Added CellRank diagnostics via `CellRankSummary()`, `VelocityFateCorrelation()`, and `plot_velocity_fate_correlation()`.
+ **RegVelo Smoke-Test Kit**: Added a self-contained `inst/regvelo-smoke-test` workflow for pyenv-based server setup, small h5ad preparation, SLURM RegVelo execution, and CellRank-after-RegVelo validation without conda/mamba.
+ **Basilisk Python Environment Reliability**:
  + Fixed `RunCellFate()`/`RunCellRank()` failing with an unsatisfiable conda dependency: `sclet_cellrank_env` pinned `python=3.12` together with `scvelo=0.3.2`, whose conda build requires `scikit-learn <1.2.0`, but no such scikit-learn build exists for Python 3.12. Bumped `scvelo=0.3.4`, which relaxed the `scikit-learn` upper bound.
  + **`GLIBCXX_... / CXXABI_... not found` — clarified host requirement.** A `GLIBCXX_3.4.31 not found` (or the `CXXABI_...` sibling, e.g. `CXXABI_1.3.15`) import failure is a **host/toolchain** problem, not a sclet bug: a compiled Python extension was built against a newer C++ runtime than the executing process provides. sclet cannot upgrade a system `libstdc++` (that needs root), so it makes every basilisk run (a) **fully activate** the environment (`full.activation = TRUE`) and (b) prepend the environment's own `lib`/`lib64` to `LD_LIBRARY_PATH` for the launched Python process. When the environment actually *bundles* a newer `libstdc++` (conda-backed basilisk environments ship `libstdc++-ng`), this makes the extension resolve the environment's library instead of an outdated system one — including when the import runs inside a basilisk worker (fork/PSOCK cluster). It is a harmless no-op for the default pip/venv backend, which bundles no `libstdc++`. Running these backends therefore still requires the host to satisfy **either** of: (a) a system C++ runtime providing the required symbol — e.g. Ubuntu ≥ 22.04 with GCC ≥ 12 (Ubuntu **20.04 and 22.04 do not** provide `GLIBCXX_3.4.31`/`CXXABI_1.3.15` by default), **or** (b) a conda-backed basilisk environment, which bundles its own newer `libstdc++`. A failed Python import now raises an actionable `sclet_glibcxx_not_found` error (matching both symbol families) pointing to these two options instead of a raw `ImportError`.
  + Fixed `RunSCENIC()` failing with `PackagesNotFoundError` for `pyscenic=0.12.1` and `loompy=3.0.7`, as reported in [#23](https://github.com/YuLab-SMU/sclet/issues/23). Neither is resolvable from conda-forge (pyscenic is a bioconda/PyPI package; `loompy=3.0.7` only exists on PyPI), so `sclet_scenic_env` now installs both via `pip`. This lets them resolve on both the conda and the pip/venv basilisk backends.
  + Added an internal `.env()` accessor so the documented troubleshooting command `basilisk::obtainEnvironmentPath(sclet:::.env("sclet_cellrank_env"))` actually works (the earlier guidance referenced a helper that did not exist).
  + Fixed `RunCellRank()`/`RunCellFate()` failing with `'Categorical' object has no attribute 'get_values'` on pandas >= 2.2, which removed the legacy `get_values()` API from `Series`/`Categorical`/`Index` (reported in [#23](https://github.com/YuLab-SMU/sclet/issues/23)). Every basilisk run now restores `get_values()` as an alias for `to_numpy()` in the launched Python process — a no-op on older pandas — so any dependency in the sclet Python stack that still calls it keeps working on current pandas.
+ **Sparse-assay gene curve / heatmap fix**: `plot_genecurve()` (`genecurve_plot`) and `pseudo_heatmap()` no longer fail with `t.default(...): argument is not a matrix` (or `cannot coerce class 'dgCMatrix' to a data.frame`) when the expression assay is a sparse `dgCMatrix` (the default `logcounts` sclet produces on large datasets), as reported in [#26](https://github.com/YuLab-SMU/sclet/issues/26). The requested feature block is now materialized with `as.matrix()` before the transposition, which keeps zero-expression cells in the data (so the smoothed/GAM gene curves are not biased upward) while preserving the fast positional `cbind` + `melt` reshape. The same latent bug in `pseudo_heatmap()` is fixed alongside it.

# sclet 1.0.0

+ **State Priority & Perturbation Sensitivity Mainline**: Introduced `RunPerturbationPriority()` wrapping the Augur framework for ranking cell types by perturbation responsiveness, `RunRareCellDetection()` as a reserved interface for rare-cell identification, and `RunStatePriorityWorkflow()` as the semantic shell over the state-priority analysis mainline. Added `get_perturbation_priority()` / `has_perturbation_priority()` and `get_rare_cells()` / `has_rare_cells()` accessors.
+ **Spatial Context & Niche Mainline**: Modernized `RunSpatialDeconvolution()` with `name` parameter, typed state record, and command logging. Introduced `RunSpatialColocalization()` and `RunSpatialNiche()` as reserved interfaces for colocalization and niche detection, and `RunSpatialWorkflow()` as the semantic shell over the spatial-context mainline. Updated `get_spatial()` to support id-based retrieval across deconvolution, colocalization, and niche records. Added `get_colocalization()` / `has_colocalization()` and `get_spatial_niche()` / `has_spatial_niche()` accessors.
+ **Multimodal Expansion Mainline**: Introduced `has_multimodal_adt()` / `has_multimodal_atac()` as modality detection helpers, `register_multimodal_adt()` / `register_multimodal_atac()` to record altExp ADT/ATAC metadata in the typed state layer, and `RunMultimodalWorkflow()` as the semantic shell that auto-detects and registers available modalities. This mainline is intentionally reserved — the object contract is established now so that modality-specific analysis backends can be plugged in later without breaking the navigation model.
+ **Visualization Layer Additions**: Added `plot_perturbation_ranking()` to visualize Augur perturbation priority AUCs as a ranked bar chart. Added `plot_program_dotplot()` for program activity summarized across groups (dot size = percent expressed, color = mean activity). Added `get_program()` / `has_program()` as a unified program-level accessor that transparently resolves signatures, pathways, and regulons from gene set scoring and SCENIC sources.
+ **Unified GRN Entry Point**: Added `RunGRN()` as a semantic frontend for gene regulatory network analysis. Currently dispatches to `RunSCENIC()`, keeping "GRN analysis" as the user-facing concept and leaving room for future backends (decoupleR, FigR, etc.) through the same interface.
+ **Spatial Visualization Layer**: Added `plot_spatial_deconvolution()` (ranked bar chart of mean cell type proportions from cell2location) and `plot_spatial_composition()` (heatmap of cell type proportions across spatial spots). Both consume the typed spatial state records produced by `RunSpatialDeconvolution()`.
+ **Velocity Latent Time Visualization**: Added `plot_velocity_latent_time()` to color cells on an embedding by velocity pseudotime (latent time), showing the predicted temporal ordering inferred by scVelo's dynamical mode.
+ **Symphony Reference Mapping**: Added `RunSymphonyMapping()` as a full Symphony backend for `RunReferenceMapping(method = "Symphony")`. The workflow builds a harmonized reference atlas from the reference SCE with `symphony::buildReference()`, maps query cells with `symphony::mapQuery()`, and predicts labels with `symphony::knnPredict()`. Supports batch variable integration (`vars` parameter) and confidence scoring. Predicted labels are stored in `colData$symphony_predicted` and confidence scores in `colData$symphony_confidence`, fully integrated with the analysis-state contract.
+ **Rare Cell Detection with Density**: Replaced the `RunRareCellDetection()` RareQ/Seurat backend with a fully SCE-native `"density"` method. The algorithm builds a kNN graph via `BiocNeighbors`, computes mutual connectivity Q-values, identifies high-density seeds, propagates labels, and flags small clusters as rare. No external package dependencies beyond `BiocNeighbors`. The RareQ backend has been removed.
+ **Spatial Colocalization with SVP**: Replaced the `RunSpatialColocalization()` stub with a full SVP backend. `SVP::runGLOBALBV()` computes global bivariate spatial autocorrelation between feature columns (e.g., deconvolution-derived cell type proportions), producing a pairwise correlation matrix. Features are auto-detected from prior deconvolution results if available, or specified explicitly. The correlation matrix is stored in the colocalization state record.
+ **Spatial Niche Detection with SVP LISA**: Replaced the `RunSpatialNiche()` stub with a full SVP backend. `SVP::runLISA()` computes Local Indicators of Spatial Association for each feature, categorizing spots as High-High (hotspot/niche), Low-Low, High-Low, Low-High, or not significant. LISA categories are stored in `colData` with a `lisa_` prefix. Features are auto-detected from prior deconvolution results.


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
