<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sclet`: A State-Aware Single-Cell & Spatial Toolkit for SingleCellExperiment

[![License:
Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
[![](https://img.shields.io/badge/devel%20version-1.1.0-blue.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/languages/code-size/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/last-commit/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet/commits/devel)
[![R-CMD-check](https://github.com/YuLab-SMU/sclet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YuLab-SMU/sclet/actions/workflows/R-CMD-check.yaml)

`sclet` is a state-aware single-cell and spatial transcriptomics
analysis toolkit built around `SingleCellExperiment`. It does not treat
`SCE` as a low-level backend hidden behind another abstraction. Instead,
it takes `SCE` as the primary analysis object and builds a unified
`Run*` API around it, with `Run<Mainline>Workflow()` shells that tie
independent modules into coherent, inspectable journeys.

Several verbs carry the same names as their Seurat counterparts
(`NormalizeData()`, `FindVariableFeatures()`, `FindClusters()`,
`FindMarkers()`, `DimPlot()`), so the syntax will look familiar. `sclet`
is not a compatibility layer over Seurat: the object of record is the
`SingleCellExperiment` throughout. One backend —
`RunPhenotypeAssociation()` — converts to Seurat internally because
scPAS requires it, and returns the results back into the `SCE`; nothing
else on the standard path does.

The package introduces an explicit **Analysis-State Contract**. In
practice, the object records what assay was used, which reduction was
active, what clustering was chosen, and what downstream results were
derived — as structured, queryable metadata rather than scattered ad hoc
slots. This provenance backbone makes workflows reproducible by default
and ready for downstream AI-assisted reasoning.

Key integration pattern: Python backends (scVI, scVelo, CellRank,
SCENIC, cell2location, CellOracle) are all isolated via `basilisk` and
exposed as plain R functions. No manual conda configuration. The first
call into a given backend builds its environment, which takes a few
minutes; later calls reuse it. The first argument is always the
`SingleCellExperiment`; the return value is always an `SCE` with the
analysis state updated.

## Quickstart

A complete analysis from raw counts, with the state recorded at every
step:

``` r
library(sclet)

## 1. Read data (Read10X also computes nFeature_RNA / nCount_RNA) and add mito fraction
sce <- Read10X("filtered_feature_bc_matrix/")
sce[["percent.mt"]] <- PercentageFeatureSet(sce, "^MT-")

## 2. Filter, then normalize -> HVG -> scale -> PCA -> neighbors -> clusters -> UMAP
sce <- subset_feature(sce, mincell = 3, peek = FALSE)
sce <- RunBasicWorkflow(sce, n_features = 2000, n_pcs = 30, cluster_resolution = 0.8)

## 3. What is on the object now
Status(sce)

## 4. Visualize and test
DimPlot(sce, reduction = "UMAP")
sce <- RunDEtest(sce, ident.1 = "1", ident.2 = "2")
get_detest(sce)$artifacts$result

## 5. Functional interpretation of those markers
sce <- RunEnrichment(sce, id = "detest", keyType = "SYMBOL", ont = "BP")
EnrichmentPlot(sce, showCategory = 5)
```

Anything that needs a class the object does not yet record tells you
what is missing rather than guessing — `check_de_readiness()`,
`check_enrichment_readiness()`, `check_fate_readiness()` and eleven more
cover every analysis stage.

## Analysis Mainlines

The package is organized around 13 user-facing mainlines, each answering
a concrete question in a single-cell or spatial project.

1.  **Data ingress & interoperability** — `Read10X()`, `ReadH5AD()`,
    `WriteH5AD()`, `as.SCE()`, `as.Seurat()`. Move between `SCE`,
    `Seurat`, and `AnnData` without losing the object story.
2.  **Core single-cell workflow** — `QCMetrics()`, `NormalizeData()`,
    `FindVariableFeatures()`, `RunPCA()`, `FindNeighbors()`,
    `FindClusters()`, `RunUMAP()`. QC through clustering in one coherent
    pipeline. `RunBasicWorkflow()` chains the whole route;
    `RunStandardPipeline()` is the shorter default-parameter form.
3.  **Quality control & data cleaning** — `RunDoubletFinder()`,
    `RunDecontX()`, `RunImputation()`. Remove doublets (scDblFinder),
    strip ambient RNA (DecontX), and fill dropouts with zero-preserving
    imputation (ALRA) before the results travel downstream.
4.  **Integration & batch correction** —
    `RunIntegration(method = c("fastMNN", "Harmony", "scVI"))`, plus
    `RunIntegrationRoutes()` to run several strategies and compare them
    on a uniform metric contract instead of eyeballing two UMAPs.
5.  **Cell identity & reference mapping** — `RunReferenceMapping()`,
    `RunSingleR()`, `RunKNNPredict()`, `RunSymphonyMapping()`. Annotate
    cells through multiple backends with confidence inspection.
6.  **Cell fate & dynamic processes** — `RunDiffusionMap()`,
    `RunSlingshot()`, `RunVelocity()`, `RunCellRank()`, `RunRegVelo()`,
    `RunCellFate()`. Trajectory, RNA velocity, and fate probability in
    one dynamics route.
7.  **Program, regulon & mechanistic interpretation** —
    `RunGeneSetScoring(method = c("UCell", "AUCell", "GSVA"))`,
    `RunGRN()`, `RunSCENIC()`, `RunKEGG()`. From gene-set scoring to
    full regulon activity networks.
8.  **Differential analysis & functional interpretation** —
    `RunDEtest()`, `FindMarkers_pseudobulk()`, `RunEnrichment()`. Marker
    discovery, pseudobulk testing, and enrichment in one chain.
9.  **Cell-cell communication & microenvironment** —
    `RunCCI(method = c("CellChat", "CellPhoneDB", "NicheNet"))`.
    Ligand-receptor signaling with multiple backend adapters.
10. **State priority & perturbation sensitivity** —
    `RunPerturbationPriority()`,
    `RunRareCellDetection(method = "density")`. Prioritize
    perturbation-sensitive (Augur) or rare cell states with SCE-native
    algorithms.
11. **Malignant & CNV inference** — `RunCNVInference()`,
    `ClassifyMalignantEpithelium()`, `CallMalignantCells()`,
    `FindMalignantSeedGenes()`. Infer copy-number profiles from
    expression with inferCNV, then call malignant cells through an
    explicit four-state decision rather than a single threshold.
12. **Spatial context & niche analysis** — `RunSpatialDeconvolution()`,
    `RunSpatialColocalization()`, `RunSpatialNiche()`. Spatial
    deconvolution (cell2location), colocalization (SVP GLOBALBV), and
    niche detection (SVP LISA).
13. **Multimodal expansion** — `RunMultimodalWorkflow()`. Reserved
    documentation path for future ADT / ATAC / multiome support.

The `Run<Mainline>Workflow()` shells sequence the modules of one
mainline and record the result as a single journey:
`RunBasicWorkflow()`, `RunReferenceWorkflow()`,
`RunTrajectoryWorkflow()`, `RunProgramWorkflow()`,
`RunSpatialWorkflow()`, `RunStatePriorityWorkflow()`,
`RunMultimodalWorkflow()`.

## The State Contract, and how to read it

`Status()` gives the one-screen summary. Everything below it is a paired
accessor: `has_*()` answers whether a result exists, `get_*()` returns
it.

``` r
Status(sce)                    # active assay/layer/reduction/graph/ident, analyses present, command count
get_analysis_context(sce)      # the full provenance DAG
has_integration(sce)           # TRUE/FALSE guard, one per analysis type
get_integration(sce)           # the recorded run, with its parameters and artifacts

ActiveIdent(sce)               # current clustering identity
PipelineSummary(sce)           # human-readable log of every command, in order
CommandLog(sce)                # the same log as a data.frame
GetAnalysisLedger(sce)         # every citable evidence node, with ids
```

Results live in the record’s `artifacts`, not in ad hoc slots —
`get_detest(sce)$artifacts$result` is the marker table,
`get_enrichment(sce)$artifacts$result` the enrichment result.

## AI layer

The recorded state is what makes an analysis plan verifiable, and
`sclet` ships an optional layer on top of it:

- `AIPlanAnalysis()` / `RunAIPlan()` — turn a question plus the object’s
  state into a step-wise plan.
- `ValidateAIPlan()` — check every step against the outputs that earlier
  steps declare, not only against the object as it stands.
- `ExecuteAIPlan()` — run the plan against the execution registry,
  recording each result as an evidence node.
- `RenderAIReport()` — render a report from the ledger and re-check it;
  a claim that the evidence cannot support is reported as such rather
  than printed as a finding.
- `check_claim_ceiling()` — ask how strong a claim the recorded evidence
  permits, without recording anything.
- `AIStatus()`, `GetAIProfile()`, `AISkills()` — inspect coverage, the
  object profile the planner sees, and the domain knowledge injected
  into planning.

The layer is inert unless used: no model is called during normal
analysis, and `aisdk` is a `Suggests` dependency, not a `Depends`.

To use it, install the optional package once and set a provider key:

``` r
install.packages("aisdk")
# then in ~/.Renviron (restart R afterwards):
#   OPENAI_API_KEY=...
```

Then describe what you want in one call:

``` r
scletAI(sce, "cluster the cells and find markers")
```

`scletAI()` prints the proposed plan and asks for approval in an interactive
session. The same engine is available with its switches exposed as
`RunAIAnalysis()`, as a read-only conversation with `AskAI()`, and a finished
round can be read back as a graph with `AnalysisGraph()`.

## Documentation

The [online book](https://yulab-smu.top/sclet) covers these mainlines as
chapters, organized into three parts:

- **Part I — How sclet Thinks** (2 chapters): the state contract and
  provenance-aware automation.
- **Part II — From Data to Discovery** (11 chapters): the analysis
  mainlines in project order. QC, doublet/ambient cleanup and imputation
  share the core workflow chapter rather than having one of their own.
- **Part III — Specialized Techniques** (7 chapters, plus a references
  appendix): visualization, metacells (SuperCell), differential
  abundance (Milo), phenotype association (scPAS), interactive
  exploration (iSEE), spatial variance (SVP), and AI-assisted analysis.

Malignant/CNV inference is documented in the reference manual
(`?RunCNVInference`, `?ClassifyMalignantEpithelium`) and does not yet
have a book chapter.

All mainline code chunks in the book run inside GitHub Actions, which
builds the book on demand at a requested `devel` ref: the book acts as
both a tutorial and a lightweight regression test suite. Pre-prepared
demo data (PBMC, Hermann spermatogenesis, DLPFC Visium, SCENIC
references) is generated from Bioconductor data packages and cached
during the CI build.

## :writing_hand: Authors

Guangchuang YU, Kaijia HUANG, Nan HE, Ming LI, Junrui LI

School of Basic Medical Sciences, Southern Medical University

<https://yulab-smu.top>

## :arrow_double_down: Installation

Get the development version from github:

``` r
## install.packages("remotes")
remotes::install_github("YuLab-SMU/sclet")
```

`sclet` declares 54 packages in `Suggests` and 4 in `Enhances`, so the
analysis verbs that wrap them are opt-in: `RunGeneSetScoring()` needs
`UCell`/`AUCell`/`GSVA`, `RunPerturbationPriority()` needs `Augur`,
`RunPhenotypeAssociation()` needs `scPAS`, `RunMilo()` needs `miloR`,
and so on. Each one stops with a `Please install ...` message naming the
package it needs, so nothing fails silently. Install the ones your
project uses:

``` r
BiocManager::install(c("scDblFinder", "miloR", "SingleR", "slingshot", "infercnv"))
```

## :book: Vignette

For more details, please refer to the online documents:

- [vignette](https://yulab-smu.top/sclet)
