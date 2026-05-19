<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sclet`: A State-Aware Single-Cell Analysis Toolkit for SingleCellExperiment

[![License:
Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
[![](https://img.shields.io/badge/devel%20version-0.99.5-blue.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/languages/code-size/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/last-commit/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet/commits/devel)

`sclet` is a state-aware analysis toolkit built around
`SingleCellExperiment`, the core Bioconductor container for single-cell
data. Rather than treating `SingleCellExperiment` as a low-level backend
hidden behind another framework, `sclet` takes it as the primary
abstraction for modern single-cell and spatial transcriptomics
workflows.

The package introduces an explicit **Analysis-State Contract (Provenance
DAG)**. In practice, this means the object records active layers,
reductions, graphs, identities, and higher-level analysis outputs as
part of a coherent workflow state. Upstream assumptions such as
batch-corrected layers, dimensional reductions, and derived analysis
results are therefore inspectable and reusable instead of being
scattered across ad hoc slots or implicit conventions.

In addition to core steps (preprocessing, dimensionality reduction,
clustering, annotation, and visualization), `sclet` provides a broad
workflow surface for both single-cell and spatial transcriptomics, with
selected Python tools managed through the `basilisk` sandbox:

- **Core Analysis**: State-aware preprocessing, dimensionality
  reduction, clustering, annotation, and provenance-aware workflow
  management.
- **Integration**: Multi-backend support (`fastMNN`, `Harmony`, `scVI`).
- **Python Interoperability**: Safe and isolated execution of `scVelo`
  (RNA Velocity), `pySCENIC` (Gene Regulatory Networks), `cellrank`
  (Advanced Fate Mapping), `cell2location` (Spatial Deconvolution), and
  `celloracle` (In-Silico Perturbation).
- **Downstream Extensions**: Gene set scoring (`UCell`, `AUCell`,
  `GSVA`), pseudobulk DE, enrichment analysis, cell-cell communication,
  and Milo differential abundance.
- **AI Copilot**: An intelligent diagnostic agent (`sclet_copilot`)
  powered by the `aisdk` framework that reads the object’s provenance to
  perform rigorous cross-chain error control.

For user-facing analysis verbs, `sclet` uses the `Run*` naming style
throughout, such as `RunIntegration()`, `RunSCENIC()`, and
`RunVelocity()`. This is now the public analysis API of the package.

## :writing_hand: Authors

Guangchuang YU

School of Basic Medical Sciences, Southern Medical University

<https://yulab-smu.top>

## :arrow_double_down: Installation

Get the development version from github:

``` r
## install.packages("remotes")
remotes::install_github("YuLab-SMU/sclet")
```

## :book: Vignette

For more details, please refer to the online documents:

- [vignette](https://yulab-smu.top/sclet)
