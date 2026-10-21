# G2 — Rd examples for the core analysis surface

Group: `.dev/tools/groups/G2.tsv` (48 Rd).
Sources touched: `R/dimred.R`, `R/clustering.R`, `R/preprocessing.R`, `R/markers.R`,
`R/Enrichment.R`, `R/pathways.R`, `R/geneset_scoring.R`, `R/annotation.R`,
`R/batch.R`, `R/integration.R`, `R/data_cleaning.R`, `R/imputation.R`,
`R/pseudobulk.R`, `R/refit_milo.R`.

Verification (all with `export R_LIBS_USER=~/Rlib`):

```
Rscript .dev/tools/check-examples.R --src-group G2       -> 48 PASS, 0 SKIP(dontrun), 0 FAIL, 0 NOEXAMPLE
Rscript .dev/tools/lint-example-symbols.R --src-group G2 -> no undefined symbols found in 48 examples
Rscript .dev/tools/lint-example-placement.R                -> none misplaced
```

**What "PASS" means with the current checker.** The checker evaluates the part of an
example outside `\dontrun{}`; SKIP(dontrun) is only reported when the whole example is
inside it. So:

* **runnable** — no `\dontrun{}`; the checker executed the whole example (18 files).
* **setup verified** — the setup outside `\dontrun{}` was executed by the checker; the
  heavy call is inside `\dontrun{}` (30 files). I additionally **executed every
  `\dontrun{}` body by hand** here (all Suggests packages are installed in this
  library) — 24 of the 30 run clean, 6 do not and are listed under "Not runnable".

## Per-file

| Rd | status | example shows | `\dontrun{}` body executed here |
|---|---|---|---|
| ActiveIdent.Rd | PASS (runnable) | `colLabels` -> `"colLabels"`; setter -> `"Cell_Cycle"` | n/a |
| AggregateExpression.Rd | PASS (setup verified) | pseudobulk by `sample`; `dim(agg)` | **runs OK** (needs `scrapper`) |
| BatchRemover.Rd | PASS (setup verified) | fastMNN over a `batch` column; corrected embedding | **runs OK** (needs `batchelor`) |
| CommandLog.Rd | PASS (runnable) | `log$command` -> "NormalizeData" | n/a |
| DefaultAssay.Rd | PASS (runnable) | after `NormalizeData()` -> `"logcounts"`; setter -> `"counts"` | n/a |
| DefaultGraph.Rd | PASS (runnable) | after `FindNeighbors()` -> `"knn_graph"` | n/a |
| DefaultReduction.Rd | PASS (runnable) | after `RunPCA()` -> `"PCA"`; setter round trip | n/a |
| ElbowPlot.Rd | PASS (runnable) | elbow curve over the 5 PCs | n/a |
| FindAllMarkers.Rd | PASS (setup verified) | markers per cluster | **runs OK** (needs `presto`) |
| FindClusters.Rd | PASS (runnable) | Louvain on the KNN graph -> 2 clusters (12 + 8 cells) | n/a |
| FindMarkerPathways.Rd | PASS (setup verified) | gene set scores as an altExp -> pathway markers | **runs OK** (needs `UCell`, `presto`) |
| FindMarkers.Rd | PASS (setup verified) | markers for one identity vs the rest | **runs OK** (needs `presto`) |
| FindMarkers_pseudobulk.Rd | PASS (setup verified) | DESeq2 on 4 pseudobulk samples, 2 per condition | **runs OK** (needs `DESeq2`) |
| FindNeighbors.Rd | PASS (runnable) | `dims = 1:5` -> `"knn_graph"` | n/a |
| FindVariableFeatures.Rd | PASS (setup verified) | `method = "scrapper"` -> HVG table | **runs OK** (needs `scrapper`) |
| Idents.Rd | PASS (runnable) | named membership vector per cell | n/a |
| NormalizeData.Rd | PASS (setup verified) | `logcounts` assay created | **runs OK** |
| PercentageFeatureSet.Rd | PASS (runnable) | per-cell % of genes matching a pattern | n/a |
| QCMetrics.Rd | PASS (setup verified) | `nCount_RNA` / `nFeature_RNA` in colData | **runs OK** (needs `scrapper`) |
| RenameIdents.Rd | PASS (runnable) | `c1`/`c2` -> `"TypeA"`/`"TypeB"` | n/a |
| RunDEtest.Rd | PASS (setup verified) | marker table recorded as a detest state | **runs OK** (needs `presto`) |
| RunDecontX.Rd | PASS (setup verified) | `decontXcounts` assay added | **runs OK** (needs `celda`) |
| RunDiffusionMap.Rd | PASS (setup verified) | `"DM"` reduction added | **runs OK** (needs `destiny`) |
| RunDoubletFinder.Rd | PASS (setup verified) | `scDblFinder.class` per cell | **runs OK** (needs `scDblFinder`) |
| RunEnrichment.Rd | PASS (setup verified) | `RunDEtest()` -> `RunEnrichment()` -> `n_enriched` | **fails** — see below |
| RunGeneSetScoring.Rd | PASS (setup verified) | UCell scores added to colData | **runs OK** (needs `UCell`) |
| RunImputation.Rd | PASS (setup verified) | `alra_imputed` assay added | **runs OK** (needs `ALRA`, and `k = 10`) |
| RunIntegration.Rd | PASS (setup verified) | Harmony on a `batch` column | **runs OK** (needs `harmony`) |
| RunKEGG.Rd | PASS (setup verified) | KEGG download -> scoring -> markers | **fails** — see below |
| RunKNNPredict.Rd | PASS (setup verified) | label transfer by KNN; labels in colData | **runs OK** (needs `BiocNeighbors`) |
| RunPCA.Rd | PASS (runnable) | `ncol(reducedDim(sce, "PCA"))` = 5 | n/a |
| RunReferenceMapping.Rd | PASS (setup verified) | `method = "KNN"` backend of the unified entry point | **runs OK** (needs `BiocNeighbors`) |
| RunSingleR.Rd | PASS (setup verified) | SingleR labels added to colData | **runs OK** (needs `SingleR`) |
| RunSymphonyMapping.Rd | PASS (setup verified) | Symphony label transfer | **fails** — see below |
| RunUMAP.Rd | PASS (runnable) | `ncol(reducedDim(sce, "UMAP"))` = 2 | n/a |
| ScaleData.Rd | PASS (runnable) | `scaled` assay added | n/a |
| VariableFeatures.Rd | PASS (runnable) | 20 HVGs named | n/a |
| get_detest.Rd | PASS (setup verified) | `rec$summary$n_markers` from the recorded detest state | **runs OK** (needs `presto`) |
| get_enrichment.Rd | PASS (setup verified) | `rec$summary$n_enriched` | **fails** — see below |
| has_detest.Rd | PASS (setup verified) | `TRUE` once a detest run is recorded | **runs OK** (needs `presto`) |
| has_enrichment.Rd | PASS (setup verified) | `TRUE` once an enrichment run is recorded | **fails** — see below |
| plot_program_activity.Rd | PASS (setup verified) | program score on the embedding | **runs OK** (needs `UCell`) |
| plot_program_dotplot.Rd | PASS (setup verified) | program activity per group | **runs OK** (needs `UCell`) |
| plot_program_heatmap.Rd | PASS (setup verified) | program activity per group as a heatmap | **runs OK** (needs `UCell`) |
| plot_reference_label_confidence.Rd | PASS (runnable) | confidence per transferred label, 20 cells | n/a |
| plot_reference_label_transfer_heatmap.Rd | PASS (runnable) | query group x predicted label table | n/a |
| refit_milo.Rd | PASS (setup verified) | Milo refit after dropping a neighbourhood | **fails** — see below |
| sce_merge.Rd | PASS (runnable) | two mock objects merged -> 40 cells | n/a |

## Not runnable, and why (honest list)

* **RunEnrichment / get_enrichment / has_enrichment** — `clusterProfiler::compareCluster(enrichGO)`
  rejects `mockSCE()` gene names with *"None of the keys entered are valid keys for
  'SYMBOL'"*. `mockSCE()` genes are `Gene_0001...Gene_0050` and are not `org.Hs.eg.db`
  symbols, so **no** enrichment table can be produced from a mock object at any
  threshold. The examples show the call chain (`RunDEtest()` -> `RunEnrichment()` ->
  accessor`) and will work on a real object with real gene symbols.
* **RunKEGG** — needs the `clusterProfiler` Suggests package *and* network access to
  KEGG; with a 50-gene mock universe UCell then fails with *"One or more signatures
  contain more genes than maxRank parameter"*. Shape demonstration.
* **RunSymphonyMapping** — Symphony's internal decomposition fails on a 20-cell object:
  *"max(nu, nv) must be strictly less than min(nrow(A), ncol(A))"*. It also failed at
  100 cells x 200 genes, so this is not only a size problem — I did not chase it
  further. Shape demonstration.
* **refit_milo** — `RunMilo()` on the mock object aborts with *"subscript out of
  bounds"* (it also fails at 100 cells x 200 genes), so the example cannot be taken
  past the `milo` construction step. Shape demonstration.

Everything else that is `\dontrun{}`-protected was executed end to end and works.

## Exported functions in G2 with no `man/<name>.Rd`

**None.** All 48 Rd files listed in `G2.tsv` exist under `man/`.
The replacement functions `DefaultReduction<-`, `DefaultAssay<-`, `DefaultGraph<-` and
`ActiveIdent<-` are exported but have no Rd of their own; they use `@rdname`, so they
share the getter's topic and inherit its examples.

## Notes for the Lead

* **The `suggests_needed` column under-reports for this group.** These functions
  need a Suggests package that the column records as `-`, so I wrapped them in
  `\dontrun{}` anyway and said why in the example:
  `FindAllMarkers`, `FindMarkers`, `RunDEtest`, `get_detest`, `has_detest`
  (**`presto`**), `FindMarkerPathways` (**`presto`** + a `RunGeneSetScoring()` altExp),
  `RunDiffusionMap` (**`destiny`**), `RunReferenceMapping` (**`BiocNeighbors`**),
  `RunSymphonyMapping` (**`symphony`**), `get_enrichment`/`has_enrichment`
  (**`clusterProfiler`**), `plot_program_activity`/`plot_program_dotplot`/
  `plot_program_heatmap` (**`UCell`**). `presto` deserves a note: it is a Suggests
  package, but `FindAllMarkers()`/`FindMarkers()` call
  `yulab.utils::install_zip_gh("immunogenomics/presto")` when it is missing — a
  **network install** inside an example, so `\dontrun{}` is the only safe shape.
* **`symphony` is still not declared in DESCRIPTION.** `RunSymphonyMapping()` calls
  `requireNamespace("symphony")` and the package is not on CRAN, so
  `Suggests: symphony` (with a GitHub install note) would be honest. `destiny` and
  `Augur` were added to `Suggests` by someone else while I was working, so my
  `RunDiffusionMap` example comment now says "needs the 'destiny' Suggests package".
* **`NormalizeData` is `\dontrun{}`-protected per the work list** (`suggests_needed =
  DelayedArray`). Reading the body, the `DelayedArray::sweep()` branch is only taken
  for DelayedArray-backed assays, and `mockSCE()`'s counts are a plain matrix — so
  the example would in fact run on any machine. I followed the work list rather than
  my own judgement; flipping that one to a runnable example is a one-line change if
  you prefer it.
* **Two examples use `mockSCE()`'s `Cell_Cycle` column** (`CellStatPlot`,
  `GroupHeatmap`) instead of `Treatment`, to avoid using `treat1`/`treat2` as cluster
  labels as instructed.
* **`as.Seurat()`/`as.SCE()` in G1b need `SingleCellExperiment::altExps(sce) <- list()`**
  — `mockSCE()` ships a "Spikes" altExp and Seurat's SCE converter swaps to every
  altExp it finds, after which the object has no `logcounts` and Seurat aborts.
* No `man/*.Rd`, `NAMESPACE` or `DESCRIPTION` file was edited by me.
