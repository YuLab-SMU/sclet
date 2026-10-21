# G3 — Rd examples for the domain layer

**Scope:** `SuperCell.R, cci.R, celloracle.R, cellrank.R, malignancy-cnv.R,
malignancy-combine.R, malignancy.R, multimodal.R, phenotype-evidence.R,
phenotype.R, regvelo.R, runMilo.R, runSlingshot.R, scenic.R, sling_plot.R,
spatial.R, velocity.R` — 49 exported functions, one `@examples` block each.

**Verification (all three tools re-run after the last edit):**

```
export R_LIBS_USER=~/Rlib
cd /home/wang/data/source/omics/sclet-wt
Rscript .dev/tools/check-examples.R --src-group G3        -> 49 PASS, 0 FAIL, 0 NOEXAMPLE
Rscript .dev/tools/lint-example-symbols.R --src-group G3  -> no undefined symbols
Rscript .dev/tools/lint-example-placement.R --src-group G3-> none misplaced
```

**Counts:** 49/49 rows PASS. 25 of the 49 examples are **fully executable**; the
other 24 are shape demonstrations whose heavy call sits inside `\dontrun{}`
with a one-line comment naming the prerequisite. Every example keeps at least
one real expression outside the `\dontrun{}`, so the checker actually executes
the setup and a broken setup line would be reported as FAIL.

Nothing outside the 17 assigned `R/` files was touched: no `man/*.Rd`, no
`NAMESPACE`, no `DESCRIPTION`, no `roxygenise()` run.

---

## Per file

| File | Functions | Result | What the examples show |
|---|---|---|---|
| `R/SuperCell.R` | `RunSuperCell` | 1 dontrun | metacell aggregation call shape (`gamma`, `k.knn`) |
| `R/cci.R` | `RunCCI`, `RunCellChat`, `RunCellPhoneDB`, `RunNicheNet` | 4 dontrun | the three CCI backends and the unified `RunCCI(method=)` dispatch |
| `R/celloracle.R` | `RunInSilicoPerturbation` | 1 dontrun | in-silico KO call shape (`target_gene`, `base_grn_path`) |
| `R/cellrank.R` | `RunCellRank`, `RunCellFate` | 2 dontrun | trajectory/fate run call shapes |
| | `CellRankSummary`, `VelocityFateCorrelation`, `plot_velocity_fate_correlation`, `plot_fate_probability`, `plot_fate_terminal_states`, `plot_fate_driver_trends` | 6 **runnable** | real CellRank state record → terminal-state counts, fate ranges, magnitude↔fate correlation, and the four CellRank plots |
| `R/malignancy-cnv.R` | `check_cnv_readiness` | **runnable** | the readiness report on a bare object: `ready` FALSE and the reason ("only 0 genes got a chromosome position…") |
| | `RunCNVInference` | 1 dontrun | inferCNV call shape with `reference_cells` |
| `R/malignancy-combine.R` | `CallMalignantCells` | **runnable** | a completed CNV record + signature → the four-state call (`malignant` / `likely_malignant` / `undetermined` / `non_malignant`), not a NULL |
| `R/malignancy.R` | `FindMalignantSeedGenes`, `ClassifyMalignantEpithelium` | 2 dontrun | bulk matrix + `sample_info` → seed; seed → iterative classification |
| `R/multimodal.R` | `has_multimodal_adt`, `has_multimodal_atac`, `register_multimodal_adt`, `register_multimodal_atac`, `RunMultimodalWorkflow` | 5 **runnable** | altExp presence, registration into the multimodal state, and the workflow record (`"ADT"`) |
| `R/phenotype-evidence.R` | `PhenotypeEvidenceProviders`, `check_phenotype_provider`, `CallPhenotypeCells` | 3 **runnable** | the provider registry (`list`, `get`), provider validation, and a real per-cell call (`probable` for signature-positive cells, `undetermined` elsewhere — correct, because only one line of evidence exists) |
| `R/phenotype.R` | `RunPhenotypeAssociation` | 1 dontrun | scPAS call shape with bulk matrix + survival phenotype |
| `R/regvelo.R` | `RunRegVelo` | 1 dontrun | GRN-matrix → velocity call shape |
| `R/runMilo.R` | `RunMilo` | 1 dontrun | `sample_col` / `condition_col` / `contrasts` shape (contrast written as condition levels, `"treat2 - treat1"`) |
| `R/runSlingshot.R` | `RunSlingshot` | 1 dontrun | `group` + `reduction` call shape |
| `R/scenic.R` | `RunSCENIC`, `RunGRN` | 2 dontrun | pySCENIC / unified GRN entry point shapes |
| | `plot_scenic_activity` | **runnable** | SCENIC altExp + state record → regulon activity on an embedding |
| `R/sling_plot.R` | `genecurve_plot` | 1 dontrun | pseudotime + features call shape |
| | `pseudo_plot` | **runnable** | `slingPseudotime` colData + a reduction → the pseudotime facets |
| `R/spatial.R` | `RunSpatialDeconvolution`, `RunSpatialColocalization`, `RunSpatialNiche`, `RunSpatialWorkflow` | 4 dontrun | cell2location / SVP LISA / SVP GLOBALBV / workflow shapes |
| | `plot_spatial_deconvolution`, `plot_spatial_composition` | 2 **runnable** | a deconvolution state record → composition bar chart and spot heatmap |
| `R/velocity.R` | `RunVelocity`, `VelocityPlot` | 2 dontrun | scVelo run and stream-plot shapes |
| | `plot_velocity_latent_time`, `VelocityMagnitude`, `plot_velocity_magnitude`, `TopVelocityGenes`, `plot_top_velocity_genes` | 5 **runnable** | a velocity state record → latent-time embedding; a velocity assay → per-cell magnitude, top genes and their plots |

---

## The 24 `\dontrun{}` examples and exactly what each needs

| Function | Needs |
|---|---|
| `RunSuperCell` | the `SuperCell` package |
| `RunCCI` | the `CellChat` package (the `method = "CellChat"` branch) |
| `RunCellChat` | the `CellChat` package |
| `RunCellPhoneDB` | the `reticulate` package **and** the CellPhoneDB Python environment |
| `RunNicheNet` | the `nichenetr` package **and** the three pre-built NicheNet model matrices |
| `RunInSilicoPerturbation` | the `reticulate` package **and** the CellOracle Python environment |
| `RunCellRank` | a velocity result **and** the `reticulate` package with a CellRank Python environment |
| `RunCellFate` | a velocity result **and** the `reticulate` package with a CellRank Python environment |
| `RunCNVInference` | the `infercnv` package **and** the JAGS system library |
| `FindMalignantSeedGenes` | the `limma` package |
| `ClassifyMalignantEpithelium` | the `limma` package |
| `RunPhenotypeAssociation` | the `scPAS` package |
| `RunRegVelo` | the `reticulate` package **and** a RegVelo Python environment |
| `RunMilo` | the `miloR` package |
| `RunSlingshot` | the `slingshot` package |
| `RunSCENIC` | the `reticulate` package **and** the pySCENIC Python environment |
| `RunGRN` | the `reticulate` package **and** the pySCENIC Python environment (it delegates to `RunSCENIC`) |
| `genecurve_plot` | the `reshape2` package (the body calls `reshape2::melt` directly) |
| `RunSpatialDeconvolution` | the `reticulate` package **and** the cell2location Python environment |
| `RunSpatialColocalization` | the `SVP` package |
| `RunSpatialNiche` | the `SVP` package |
| `RunSpatialWorkflow` | the `reticulate` package **and** the cell2location Python environment (it dispatches to the spatial backends) |
| `RunVelocity` | the `velociraptor` package **and** a scVelo Python environment |
| `VelocityPlot` | the `velociraptor` package **and** a scVelo Python environment |

### Where I departed from the TSV `suggests_needed` column

* **`check_cnv_readiness` — made runnable on purpose.** The TSV lists
  `infercnv,rjags`, but the function only *probes* for them
  (`requireNamespace(..., quietly = TRUE)`) and reports what is missing; it
  never needs them to run. Rule 10 prefers a real example for a
  readiness/probe function, so this one now shows a user exactly what their
  object is missing (`readiness$ready` is `FALSE`, `readiness$issues` names
  the missing gene order).
* **`RunCCI`, `RunGRN`, `RunCellFate`, `RunSpatialWorkflow`,
  `RunPhenotypeAssociation` — wrapped although the TSV says `-`.** Their own
  bodies reach a Suggests package indirectly (`RunCCI` → `RunCellChat` →
  `CellChat`; `RunGRN` → `RunSCENIC` → `reticulate`; `RunCellFate` →
  `RunCellRank` → `reticulate`; `RunSpatialWorkflow` → the spatial backends;
  `RunPhenotypeAssociation` → `requireNamespace("scPAS")`). They cannot run
  where Suggests is absent.

### What could not be made runnable, and why

Nothing in the group is blocked by something I could have fixed cheaply; the 24
`\dontrun{}` cases are exactly the ones whose *backend* is a Suggests package
or a Python environment. Two of them need an input that cannot be synthesised:

* `RunNicheNet` needs the three NicheNet model matrices
  (`ligand_target_matrix`, `lr_network`, `weighted_networks`). The example binds
  all three with `readRDS("…rds")` and a comment pointing at the NicheNet repo,
  so no name is read without being bound, but a reader still has to download
  them. That is a property of NicheNet, not of the example.
* `RunSCENIC` / `RunGRN` / `RunInSilicoPerturbation` / `RunCellPhoneDB` need
  local resource files (TF list, motif annotations, cisTarget `.feather`
  databases, base GRN CSV) as well as the Python environment; the example shows
  the argument names and the file names the function expects.

---

## Exported functions in G3 with no `man/<name>.Rd`

**None.** All 49 names in `G3.tsv` have a `man/<name>.Rd`, and every one of
them carries the `% Generated by roxygen2` header (i.e. roxygen owns them, not
a hand-written file).

---

## Notes for the Lead before the roxygenise pass

1. **One whitespace change outside `@examples`.** `R/runSlingshot.R` had a
   blank line between the `RunSlingshot` roxygen block and its definition.
   Roxygen binds the block either way, but the checker's source-side extractor
   only walks back over *contiguous* `#'` lines, so the example was invisible
   to it (`NOEXAMPLE`). I removed the blank line; the block is now contiguous
   and the example is picked up. Nothing else about the block changed.

2. **`lint-example-symbols.R` cannot see inside `\dontrun{}`.** It calls
   `parse(text = code)` on the raw example text, and `\dontrun{` is not valid
   R, so the parse fails and both the "bound names" and the "read names"
   walkers return empty — the linter then reports *"no undefined symbols
   found"* for **every** `\dontrun` example, silently. It only actually checks
   examples that contain no `\dontrun` at all. Reproduced:
   `parse(text = "…\n\\dontrun{\nseed <- 1\n}")` → *"unexpected symbol"*.
   Suggested fix: normalise the text before parsing, e.g.
   `code <- gsub("\\\\dontrun\\s*\\{", "if (FALSE) {", code)` (and the same for
   `\donttest`), then parse. With that, the linter would have found the `seed`
   defect the Lead spotted.

3. **I fixed five unbound-symbol defects by hand, inside `\dontrun{}` bodies**,
   which the linter as it stands cannot see:
   * `ClassifyMalignantEpithelium` — `seed` was read but never created (the
     defect the Lead flagged). The example now builds
     `seed <- list(tumor = rownames(sce)[1:5], normal = rownames(sce)[6:10])`
     inside the `\dontrun{}` and adds the `sample` column the function needs.
   * `RunRegVelo` — `grn` is now built as a 50×2 matrix with dimnames.
   * `RunSpatialDeconvolution` / `RunSpatialWorkflow` — `ref` is now a second
     `mockSCE` with a `cell_type` column.
   * `RunPhenotypeAssociation` — `bulk_matrix` and `phenotype_data` are now
     built.
   * `RunNicheNet` — the three model matrices are now bound with `readRDS`.

   I re-ran the same bound/read analysis with the `\dontrun{}` markers
   normalised to `if (FALSE) {`: **0 unbound names across all 49 examples.**

4. **How the state-reading functions get their prerequisite.** There is no
   exported state builder (`sclet_set_state_record()` is internal and `:::`
   is out of the question), so the examples that need a recorded result write
   the record straight into the object's own state:

   ```r
   S4Vectors::metadata(sce)$sclet$states$records$trajectory$cellrank <- record
   S4Vectors::metadata(sce)$sclet$states$active$trajectory <- "cellrank"
   ```

   This is the same convention the G1a group used in `R/analysis-accessors.R`,
   it preserves the rest of the state (unlike replacing `metadata(sce)$sclet`
   wholesale), and it keeps the examples honest: the record is built, so the
   call returns a real answer rather than `NULL`. The `record` list carries the
   same `method` / `inputs` / `artifacts` shape that `RunCellRank()`,
   `RunSCENIC()` and `RunVelocity()` write.

5. **Two judgement calls worth a second pair of eyes:**
   * `CallPhenotypeCells` and `CallMalignantCells` emit the package's own
     "evidence provider is unavailable" warnings (the CNV line is not there).
     That is the correct behaviour and the warning is informative, but it does
     mean those two examples print warnings. `R CMD check` treats example
     warnings as non-fatal.
   * `CallMalignantCells`'s example fabricates a *completed* CNV record
     (`status = "completed"`) so that the two-line call can actually be
     exercised. It is clearly commented as "a completed CNV run, as
     `RunCNVInference()` would record it", but a reader should understand the
     burden column is synthetic.

6. **`RunCCI` is a unified entry point over three backends**, so its example
   names the `CellChat` branch specifically (`method = "CellChat"`); the
   `CellPhoneDB` and `NicheNet` branches have their own examples.
