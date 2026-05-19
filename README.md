<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sclet`: A State-Aware Single-Cell Analysis Toolkit for SingleCellExperiment

[![License:
Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
[![](https://img.shields.io/badge/devel%20version-0.99.3-blue.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/languages/code-size/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/last-commit/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet/commits/devel)

Seurat is a mainstream tool for single-cell analysis, but its data
structure can change between versions and may be less explicit for
teaching and extension. In contrast, `SingleCellExperiment` is a
well-defined Bioconductor class with a rich ecosystem in R and a
corresponding Python counterpart, `AnnData`.

The goal of `sclet` is to provide a lightweight yet powerful set of
Seurat-like helpers for `SingleCellExperiment`. More importantly,
`sclet` introduces an innovative **Analysis-State Contract (Provenance
DAG)**. This means the object automatically tracks and manages upstream
assumptions (like batch-corrected layers or dimensional reductions),
eliminating “parameter hell” and allowing downstream functions (and even
AI Copilots) to seamlessly consume the correct contextual data.

In addition to core steps (preprocessing, dimensionality reduction,
clustering, and visualization), `sclet` offers a robust ecosystem of
wrappers for both R and Python top-tier tools via the `basilisk`
sandbox:

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

For user-facing analysis verbs, `sclet` uses `Run*` as the canonical API
naming style, such as `RunIntegration()`, `RunSCENIC()`, and
`RunVelocity()`. Legacy `run*` names are still available as
compatibility aliases.

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
