# G1b — Rd examples for the plotting and IO surface

Group: `.dev/tools/groups/G1b.tsv` (23 Rd).
Sources touched: `R/visualization.R`, `R/io.R`, `R/interop.R`, `R/interactive.R`,
`R/plot-cci.R`, `R/pseudo_heatmap.R`.

Verification (all with `export R_LIBS_USER=~/Rlib`):

```
Rscript .dev/tools/check-examples.R --src-group G1b     -> 21 PASS, 2 SKIP(dontrun), 0 FAIL, 0 NOEXAMPLE
Rscript .dev/tools/lint-example-symbols.R --src-group G1b -> no undefined symbols found in 23 examples
Rscript .dev/tools/lint-example-placement.R                -> none misplaced
```

**What "PASS" means with the current checker.** The checker evaluates the part of an
example that sits *outside* `\dontrun{}` and reports SKIP(dontrun) only when the whole
example is inside it. So the table below distinguishes two kinds of PASS:

* **runnable** — no `\dontrun{}` at all; the checker executed the entire example.
* **setup verified** — the setup lines outside `\dontrun{}` were executed by the
  checker; the heavy call sits inside `\dontrun{}`. I additionally **executed every
  `\dontrun{}` body by hand** in this environment (all the Suggests packages are
  installed here) and report the outcome in the last column.

## Per-file

| Rd | status | example shows | `\dontrun{}` body executed here |
|---|---|---|---|
| CellDimPlot.Rd | PASS (runnable) | PCA + `colLabels` -> `CellDimPlot()`; `nrow(p$data)` = 20 cells | n/a |
| CellStatPlot.Rd | PASS (runnable) | `group.by = "Cell_Cycle"` -> stacked proportion bar, 20 cells | n/a |
| DEtestPlot.Rd | PASS (setup verified) | `RunDEtest()` -> volcano plot; `class(p)[1]` | **runs OK** (needs `presto`) |
| DimPlot.Rd | PASS (runnable) | `NormalizeData()` -> `RunPCA(5)` -> `DimPlot(reduction = "PCA")`, 20 cells | n/a |
| EnrichmentPlot.Rd | PASS (setup verified) | `RunDEtest()` -> `RunEnrichment()` -> `EnrichmentPlot()` | **fails** — see "Not runnable" below |
| FeatureDimPlot.Rd | PASS (runnable) | three features painted on the PCA; `class(p)[1]` = "ggsc" | n/a |
| FeatureScatter.Rd | PASS (runnable) | two synthetic `colData` columns -> 20 points | n/a |
| GroupHeatmap.Rd | PASS (runnable) | 5 features grouped by `Cell_Cycle`; `class(p)[1]` = "pheatmap" | n/a |
| ProjectionPlot.Rd | PASS (runnable) | query + reference, both PCA'd, shared `label.main` -> one ggplot | n/a |
| Read10X.Rd | SKIP(dontrun) | whole example inside `\dontrun{}` (needs a real 10x folder) | placeholder path, as intended |
| ReadH5AD.Rd | SKIP(dontrun) | whole example inside `\dontrun{}` (needs a real `.h5ad` file) | placeholder path, as intended |
| RunExplorer.Rd | PASS (setup verified) | `RunExplorer()` launches the iSEE app | **runs OK** (needs `iSEE`) |
| VariableFeaturePlot.Rd | PASS (runnable) | after `FindVariableFeatures(method = "scran")`; 50 genes scored | n/a |
| VlnPlot.Rd | PASS (runnable) | one numeric `colData` feature -> violin; `class(p)[1]` = "patchwork" | n/a |
| WriteH5AD.Rd | PASS (setup verified) | `WriteH5AD()` into `tempdir()` | **runs OK** (needs `zellkonverter`) |
| as.SCE.Rd | PASS (setup verified) | SCE -> Seurat -> SCE round trip; `dim()` = 50 x 20 | **runs OK** (needs `Seurat`) |
| as.Seurat.Rd | PASS (setup verified) | SCE -> Seurat; `class(seurat_obj)[1]` | **runs OK** (needs `Seurat`) |
| plot_cci_arrow.Rd | PASS (setup verified) | `plot_cci_arrow()` between two cell types | **fails** — needs a recorded CCI run |
| plot_cci_bubble.Rd | PASS (setup verified) | `plot_cci_bubble()` | **fails** — needs a recorded CCI run |
| plot_cci_chord.Rd | PASS (setup verified) | `plot_cci_chord()` | **fails** — needs a recorded CCI run |
| plot_cci_network.Rd | PASS (setup verified) | `plot_cci_network()` | **fails** — needs a recorded CCI run |
| plot_cci_sigmoid.Rd | PASS (setup verified) | `plot_cci_sigmoid()` | **fails** — needs a recorded CCI run |
| pseudo_heatmap.Rd | PASS (setup verified) | `pseudo_heatmap()` over a pseudotime trajectory | **fails** — needs slingshot pseudotime |

## Not runnable, and why

* **EnrichmentPlot** — `RunEnrichment()` -> `clusterProfiler::compareCluster(enrichGO)`
  rejects `mockSCE()` gene names: *"None of the keys entered are valid keys for
  'SYMBOL'"*. `mockSCE()` calls its genes `Gene_0001...Gene_0050`, which are not
  `org.Hs.eg.db` symbols, so **no** enrichment can be produced from a mock object at
  any threshold. The example therefore demonstrates the call sequence, not a populated
  result. A real object with real gene symbols is required.
* **plot_cci_arrow / bubble / chord / network / sigmoid** — all five read the
  standardised CCI interaction table through `get_cci()`, which returns `NULL` unless a
  `communication` state record exists. The only ways to create one are `RunCCI()` ->
  `RunCellChat()` / `RunCellPhoneDB()` / `RunNicheNet()`; the first needs the `CellChat`
  Suggests package (and a Python env), the other two need Python backends, and — more
  decisively — `mockSCE()` has no ligand/receptor signal at all, so the recorded table
  would be empty even on a successful run. These five are shape demonstrations.
* **pseudo_heatmap** — requires `colData(sce)$slingPseudotime`, i.e. a real
  `RunSlingshot()` result. Shape demonstration.
* **Read10X / ReadH5AD** — the arguments are necessarily placeholders (a 10x folder, an
  `.h5ad` file); the whole example is inside `\dontrun{}`.

## Exported functions in G1b with no `man/<name>.Rd`

**None.** All 23 Rd files listed in `G1b.tsv` exist under `man/`.

## Notes for the Lead

* `as.Seurat()` / `as.SCE()` need one extra line in the setup —
  `SingleCellExperiment::altExps(sce) <- list()`. Reason: `mockSCE()` ships a
  **"Spikes" altExp**, and `Seurat:::as.Seurat.SingleCellExperiment()` swaps to every
  altExp it finds; after the swap the object has no `logcounts` assay and Seurat
  aborts with *"No data in provided assay - logcounts"*. This is a `mockSCE()` +
  Seurat interaction, not an sclet bug, but it will bite any user who converts a
  mock or spike-containing SCE, so it is worth a line in the interop docs.
* The `plot_cci_*` functions and `pseudo_heatmap()` have no way to be exercised
  without a real upstream analysis. If you want runnable examples for them, the
  prerequisite state has to become constructible from public API (e.g. an exported
  helper that records a CCI table); that is a package change, not a docs change.
* The `suggests_needed` column in `G1b.tsv` was accurate for this group. All eleven
  rows with a Suggests package are `\dontrun{}`-protected with a one-line reason, and
  the two rows whose column says `-` (`plot_cci_arrow`, `plot_cci_bubble`) are also
  `\dontrun{}`-protected for the reason given above.
* No `man/*.Rd`, `NAMESPACE` or `DESCRIPTION` file was edited by me.
