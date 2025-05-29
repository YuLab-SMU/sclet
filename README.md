<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sclet`: A Lightweight Toolkit for Single-Cell Data Analysis

[![License:
Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
[![](https://img.shields.io/badge/devel%20version-0.0.7-blue.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/languages/code-size/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/last-commit/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet/commits/devel)

Seurat is a mainstream tool for single-cell analysis, but its data
structure lacks a clear definition and often changes between versions.
This makes it easy to get started but difficult to master. In contrast,
`SingleCellExperiment` has a well-defined structure, a rich ecosystem in
R, and a corresponding Python counterpart, `AnnData`. Clearly, building
upon such a data structure is highly beneficial for in-depth learning
and even development. This is the motivation behind my development of
`sclet`—to help students in my team transition to
`SingleCellExperiment`.

`SingleCellExperiment` itself boasts a robust ecosystem, covering nearly
all types of analyses. We will also integrate existing single-cell R
packages or develop new functionalities to make them easier to use
within the `SingleCellExperiment` ecosystem.

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

For more details, please refer to the [online
vignette](https://yulab-smu.top/sclet).
