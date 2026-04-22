<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sclet`: A Lightweight Toolkit for Single-Cell Data Analysis

[![License:
Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
[![](https://img.shields.io/badge/devel%20version-0.0.7-blue.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/languages/code-size/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet)
[![](https://img.shields.io/github/last-commit/YuLab-SMU/sclet.svg)](https://github.com/YuLab-SMU/sclet/commits/devel)

Seurat is a mainstream tool for single-cell analysis, but its data
structure can change between versions and may be less explicit for
teaching and extension. In contrast, `SingleCellExperiment` is a
well-defined Bioconductor class with a rich ecosystem in R and a
corresponding Python counterpart, `AnnData`. The goal of `sclet` is to
provide a lightweight set of Seurat-like helpers for
`SingleCellExperiment`, making common workflows easier to learn and
apply.

In addition to core steps (preprocessing, dimensionality reduction,
clustering, visualization, batch correction, and pseudobulk differential
expression), `sclet` also offers optional wrappers for popular
downstream tools, including trajectory inference, enrichment analysis,
cell type annotation, cell-cell communication, and Milo differential
abundance.

For user-facing analysis verbs, `sclet` now uses `Run*` as the
canonical API naming style, such as `RunPCA()`, `RunMilo()`, and
`RunCellChat()`. Legacy `run*` names are still available as
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

- [vignette](https://yulab-smu.top/sclet).
- [中文文档](https://mp.weixin.qq.com/mp/appmsgalbum?__biz=MzI5NjUyNzkxMg==&action=getalbum&album_id=3906595076079616009&subscene=126&scenenote=https%3A%2F%2Fmp.weixin.qq.com%2Fs%3F__biz%3DMzI5NjUyNzkxMg%3D%3D%26mid%3D2247497179%26idx%3D1%26sn%3Df7f7bba43b50ad3ee0178c5433c342a5%26chksm%3Ded7973e4430f91e4f714add65636a6fcd36caff85da9272f90b9fa3aca5f5efc67ffbf24681e%26scene%3D126%26sessionid%3D1748485834%26subscene%3D91%26clicktime%3D1748485845%26enterid%3D1748485845%26key%3Ddaf9bdc5abc4e8d003cb227738599c77c9fce87c8f671e25c0864bc93f93b02f6418ab7aafd990f3a4e872f4b9d1574e2f3d8c2ef5f007b08aa017a0f7707dc524269249302384add00cfd21694b501cbab0ea9e94c3ad1a8b6ae1f70c8266956e0e58c7846fba753c485a19a4f1e22a06cb5e6caad2e85916e43e39973f64a9%26ascene%3D0%26uin%3DMTMxNjc4OTY2Mg%253D%253D%26devicetype%3DWindows%2B11%2Bx64%26version%3D63090c33%26lang%3Den%26countrycode%3DCN%26exportkey%3Dn_ChQIAhIQ8W%252B8eoJgciVToRoGlBZ9jBLjAQIE97dBBAEAAAAAABWFLmbKJU0AAAAOpnltbLcz9gKNyK89dVj0uzkEcwpT0bZ2zuDWiUxHxXUO0ahGgovowdMQSbAG3qDXBwZ5BxX5PP3WdcK7G0hSbdguCbCXV63%252BTtFDdm0crf1Rrzt3THBC26RPjsb0rZkeN3o4A9pQmaXM0skzNbVArNcC0vaVS53HuA5iuW5TlZL3FsRs2MFhdWlCeJsmIx2JQtfX%252Br5f54Fbtq%252FDSD3EJM4ySt9eIOQvD2rNS1o%252BwftM64AbYzqkp7fMmHhG%252Bs5j8OgTn9nBZDGZsqFa%26acctmode%3D0%26pass_ticket%3DI2kYj%252BzyOCQXjbBjcqegQiT9Zp3pOx2fQzQX1foRUOrxoZSBb3AlNjginGhf8WUe%26wx_header%3D1%26fasttmpl_type%3D0%26fasttmpl_fullversion%3D7752341-en_US-zip%26fasttmpl_flag%3D1&nolastread=1&sessionid=434166369#wechat_redirect)
