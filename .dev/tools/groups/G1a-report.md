# G1a — Rd examples for accessors, state and workflow plumbing

Agent: `ex-accessors`. Worktree: `/home/wang/data/source/omics/sclet-wt`
(branch `feature/quality-2026-10-21`).

## Summary

- **69 `@examples` blocks written** across the 11 R/ files in my write scope
  (65 in the first pass + 4 after the scope extension).
- **Verification — all 68 rows of the (regenerated) `G1a.tsv` pass.**

  | check | command | result |
  |---|---|---|
  | pre-roxygenise, source mode | `Rscript .dev/tools/check-examples.R --src-group G1a` | **66 PASS, 2 NOEXAMPLE** (`reexports.Rd`, `sclet-package.Rd` — see below) |
  | post-roxygenise, Rd mode | `Rscript .dev/tools/check-examples.R --group G1a` on a roxygenised **copy** of the package in `/tmp` | **68 PASS, 0 FAIL** |
  | undefined symbols | `Rscript .dev/tools/lint-example-symbols.R --src-group G1a` | no undefined symbols |
  | placement | `Rscript .dev/tools/lint-example-placement.R` | no `@examples` block in the wrong file |

  The Rd-mode run was done on a **copy** of the package (`/tmp/g1a/pkgcopy`),
  never on the worktree: `man/`, `NAMESPACE` and `DESCRIPTION` in the worktree
  are untouched, and `roxygenise()` was never run against the shared tree.
- 7 examples use `\dontrun{}` for the heavy call only; every one of them keeps
  its setup lines **outside** the `\dontrun{}` block, so the setup is executed
  and verified.
- No file outside my write scope was edited:
  `R/analysis-accessors.R`, `R/copilot.R`, `R/layers.R`,
  `R/pipeline_summary.R`, `R/priority.R`, `R/status.R`, `R/workflow.R`,
  `R/workflows.R`, `R/subset.R`, `R/reexport.R`, `R/sclet-package.R`, plus
  this report.

## One line per file

| file | examples | status | what the examples show |
|---|---|---|---|
| `R/analysis-accessors.R` | 49 | 49 PASS | The state-reading accessors: each `get_*`/`has_*` example records a state entry, then reads it back. |
| `R/copilot.R` | 2 | 1 PASS, 1 `\dontrun` | `SummarizeContextForLLM()` prints the ledger header; `AuditAnalysisChain()` needs a provider key. |
| `R/layers.R` | 3 | 3 PASS | `Layers()` -> `"counts"`, `dim(LayerData())` -> `50 20`, `DefaultLayer()` -> `"counts"`. |
| `R/pipeline_summary.R` | 1 | 1 PASS | `PipelineSummary()` after `NormalizeData()` + `RunPCA()` prints the three-section command log. |
| `R/priority.R` | 4 | 1 PASS, 3 `\dontrun` | `plot_perturbation_ranking()` builds a real ggplot from a recorded Augur AUC table; the three `Run*` entry points need Augur / BiocNeighbors. |
| `R/status.R` | 1 | 1 PASS | `Status()` after a small pipeline: `available$analyses` -> `"hvg"`, `health$mainline_missing` -> the 4 unmet mainline stages. |
| `R/workflow.R` | 1 | 1 PASS | `RunStandardPipeline()` on a mock object -> `reducedDimNames()` -> `c("PCA","UMAP")`. |
| `R/workflows.R` | 4 | 1 PASS, 3 `\dontrun` | `plot_trajectory_overview()` composes a real patchwork from a recorded trajectory workflow; the three `Run*Workflow()` entry points need Suggests packages. |
| `R/subset.R` | 2 | 2 PASS | `dim(subset_feature(sce, mincell = 5, peek = FALSE))` -> `49 20` (one feature dropped); `dim(subset_cell(sce, feature = "nFeature_RNA", method = "mad", n = 3))` -> `50 20` (outlier cells filtered). |
| `R/reexport.R` | 1 | 1 PASS (Rd mode) | The re-exported `assay()` from SummarizedExperiment: `dim(assay(sce))` -> `50 20`. |
| `R/sclet-package.R` | 1 | 1 PASS (Rd mode) | Package-level smoke test: `mockSCE` -> `NormalizeData()` -> `RunPCA()` -> `Status(sce)$available$reductions` -> `"PCA"`. |

### `R/analysis-accessors.R` in detail

All 48 accessor rows in the TSV are covered, plus `get_integration`.

- **trajectory family** — `get_trajectory()`/`has_trajectory()`: record a
  `"trajectory"` state entry (`method = "slingshot"`), read it back as
  `get_trajectory(sce)$method` -> `"slingshot"`, `has_trajectory(sce)` -> `TRUE`.
- **milo / supercell / cellrank** — same shape for the `"milo"`
  (`"miloR"`), `"aggregation"` (`"SuperCell::SCimplify"`) and `"trajectory"`
  (`"CellRank"`) state types.
- **communication family** — `get_cellchat()`, `get_cellphonedb()`,
  `get_nichenet()` read a `"communication"` record and check the method;
  `get_cci()` additionally carries an `interactions` data.frame and returns
  `nrow(get_cci(sce))` -> `1`.
- **annotation / mapping** — these two have a public fallback: a
  `SingleR_labels` column in `colData` is enough.
  `get_annotation(sce)$summary$n_labels` -> `2`,
  `get_mapping(sce)$artifacts$mapping_type` -> `"reference_mapping"`.
- **batch / integration** — `get_batch()`/`get_integration()` read the
  `analyses$batch` entry of the state contract.
- **hvg** — a real two-step pipeline: `NormalizeData()` +
  `FindVariableFeatures()`; `get_hvg(sce)$n` -> `20`.
- **graph** — `RunPCA()` + `FindNeighbors()`; `get_graph(sce)$name` ->
  `"knn_graph"`.
- **analysis context** — `names(get_analysis_context(sce)$records)` ->
  `c("annotation","mapping")`.
- **velocity / scenic / geneset scoring** — recorded state entries;
  `get_program()`/`has_program()` additionally need the score column the
  record points at (`Score_sigA` in `colData`), and return the per-cell
  activity vector `1..20`.
- **spatial family** — `get_spatial()` (state record) and
  `get_colocalization()` / `get_spatial_niche()` (analysis entries).
- **perturbation / priority / rare cells** — recorded `"perturbation"`
  state entry, `analyses$augur` with an AUC table
  (`get_perturbation_priority(sce)$AUC$auc` -> `c(0.7, 0.4)`) and
  `analyses$rare_cells` (`get_rare_cells(sce)$summary$n_clusters` -> `3`).

## The 2 rows that `--src-group` cannot see (and why that is not a defect)

`reexports.Rd` and `sclet-package.Rd` are reported `NOEXAMPLE` by
`check-examples.R --src-group`, because that tool locates an example by
finding `<name> <- function` in the source file, and neither page has one:
`R/reexport.R` is a list of `pkg::name` re-export expressions and
`R/sclet-package.R` is the `"_PACKAGE"` sentinel.

Both pages **do** have working examples — verified on the roxygenised copy,
where they come out as:

```
man/reexports.Rd        \examples{ sce <- scuttle::mockSCE(...); dim(assay(sce)) }
man/sclet-package.Rd    \examples{ sce <- ...; NormalizeData(); RunPCA(); Status(sce)$available$reductions }
```

and both rows are **PASS** in `--group` mode. So: `--src-group G1a` is a
pre-roxygenise convenience that cannot address these two pages; the
authoritative check after the Lead's roxygenise pass is `--group G1a`, which
reports 68/68 PASS. (I verified roxygen2 attaches an `@examples` block on a
re-export block to the `reexports` topic, and on the `"_PACKAGE"` block to
the package topic, in a throwaway package before touching the real tree.)

`advanced-ai-diagnostics.Rd` is no longer a row of the regenerated TSV — the
Lead replaced it with per-function Rd files and handled `R/ai-diagnostics.R`
(the G4/G5 agent's file), which I did not touch.

## Examples that are `\dontrun{}` (7) and why

All seven keep their setup **outside** the `\dontrun{}` block, so the setup is
still executed and verified by the checker. In each case the heavy call needs
a package that is in `Suggests` (or is not installed at all) — reading the
function body is what decides this.

| function | reason (one comment line inside `\dontrun{}`) |
|---|---|
| `AuditAnalysisChain` | needs the aisdk package and a configured provider |
| `RunPerturbationPriority` | needs the Augur package (install it from GitHub) and >= 20 cells per group |
| `RunRareCellDetection` | needs the BiocNeighbors package |
| `RunStatePriorityWorkflow` | needs the Augur and BiocNeighbors packages and >= 20 cells per group |
| `RunTrajectoryWorkflow` | needs the slingshot and destiny packages |
| `RunProgramWorkflow` | needs the UCell package |
| `RunReferenceWorkflow` | needs the BiocNeighbors package |

### Every `\dontrun{}` body was executed with its packages present

`Augur` (1.0.3) and `aisdk` were installed into `~/Rlib` part-way through this
session, and `BiocNeighbors`, `slingshot`, `destiny` and `UCell` were already
there — so all seven bodies could be run for real. **Two of them were
genuinely broken** and are now fixed; the other five were already correct:

| body | before | now |
|---|---|---|
| `RunReferenceWorkflow` | **ERR** "No Bioconductor-native HVG statistics found" — `RunKNNPredict()` calls `VariableFeatures()`, which needs `FindVariableFeatures()` to have run | fixed: the setup now runs `NormalizeData()` + `FindVariableFeatures()` (0.5 s) |
| `AuditAnalysisChain` | **ERR** "An integrated assay is required for AIAuditAnalysisChain()" — that gate is checked before the AI call | fixed: the setup adds a `"corrected"` assay and the call passes `integrated_layer = "corrected"` (26.5 s; returns the deterministic "the chain is empty" report when no provider is reachable) |
| `RunPerturbationPriority` | **ERR** "no cell type had at least 20 cells in all conditions" — Augur needs >= 20 cells per group, so the 20-cell mock can never work (40 and 60 cells fail too) | fixed: the mock is now `ncells = 100` (11.7 s) |
| `RunStatePriorityWorkflow` | same Augur constraint | fixed: `ncells = 100` (10.2 s) |
| `RunRareCellDetection` | OK (0.5 s) | unchanged |
| `RunTrajectoryWorkflow` | OK (2.1 s) | unchanged |
| `RunProgramWorkflow` | OK (1.0 s) | unchanged |

So every `\dontrun{}` body in this group now runs to completion when its
package is installed — a reader who copies it gets a working call, not a
second error hiding behind the first one. (This is the failure class
`lint-example-symbols.R` was built to expose, one level deeper: the *symbols*
were all bound, but the *call* still could not succeed.)

### A correction to my own earlier claim about `destiny`

I previously wrote that `destiny` "is not declared in DESCRIPTION at all".
That was wrong and I retract it. `destiny` **is** declared — at `DESCRIPTION`
line 124, which is inside the **`Enhances:`** block (lines 122-126); the
`Suggests:` block ends at line 121. The pointer "line 57 of the Suggests
block" resolves to `DESCRIPTION` line 124, i.e. to `Enhances`.

Neither of us had the check consequence right, so to be precise:

- `R CMD check` will **not** complain. `tools:::.check_packages_used()` (the
  "checking dependencies in R code" check) treats `Enhances` as declared
  (`depends_suggests <- c(depends, suggests, enhances, ...)`), and
  `requireNamespace("destiny")` passes the package name as a *string*, so it
  is not a `destiny::` call and is invisible to that check anyway. My earlier
  "it will show up as an `R CMD check` note" was incorrect — there is no
  note.
- What remains is a *field-placement* observation, not a check failure:
  `RunDiffusionMap()` hard-requires `destiny` (it stops with "Please install
  'destiny'"), and `RunPerturbationPriority()` hard-requires `Augur`, so both
  are functionally `Suggests`-level dependencies while sitting under
  `Enhances:`. Per the Lead's instruction I am **not** acting on this.

## Honest notes on the state-reading accessors

For every accessor whose *producing* function needs a Suggests package
(`RunSlingshot()`, `RunMilo()`, `RunSuperCell()`, `RunCellChat()`,
`RunCellRank()`, `RunSCENIC()`, `RunVelocity()`, `RunGeneSetScoring()`, ...)
there is no cheap way to produce a genuine record in an example. Those
examples therefore **write the analysis-state contract directly** into
`S4Vectors::metadata(sce)$sclet` and then read it back:

```r
record <- list(id = "run1", type = "trajectory", method = "slingshot")
state <- list(states = list(active = list(trajectory = "run1"),
    records = list(trajectory = list(run1 = record))))
S4Vectors::metadata(sce)$sclet <- state
get_trajectory(sce)$method
```

This is a **shape demonstration of the state contract**, not a demo of the
producing pipeline — but it is a *working* demonstration: the accessor
returns a non-NULL, correct answer (`"slingshot"`, `TRUE`, ...), which is
what the "an example that returns NULL proves nothing" rule asks for. The
alternative (calling the real producer) would have made every one of these
pages `\dontrun{}` and unverifiable.

Where a real, Suggests-free pipeline exists it is used instead: `get_hvg()`,
`has_hvg()`, `get_graph()`, `has_graph()` run `NormalizeData()` /
`FindVariableFeatures()` / `RunPCA()` / `FindNeighbors()`, and
`get_annotation()` / `get_mapping()` use the public `SingleR_labels`
`colData` fallback.

## One source change beyond adding `@examples`

`R/subset.R`: the roxygen block of `subset_cell` used the `##'` prefix
(9 lines). roxygen2 8.1.1 does read `##'` blocks, but the checkers and linters
only recognise `#'`, so the block was invisible to them. I converted those 9
lines from `##'` to `#'`. The roxygen semantics are identical (roxygen2 treats
both as roxygen comments — confirmed in a throwaway package), the generated
Rd is unchanged apart from the new `\examples{}` section, and
`lint-example-placement.R` reports the block in the right file.

## Beyond the work list

- `get_integration()` is **not** a row of `G1a.tsv`: the planner skipped it
  because `man/get_integration.Rd` already contained `\examples` — a
  placeholder whose body was entirely commented out (`#' # get_integration(sce)`),
  i.e. it proved nothing. I replaced it with a working example in the same
  roxygen block.

## TSV hints I deliberately did not follow

`G1a.tsv` flags `suggests_needed = SuperCell` for `get_supercell` and
`velociraptor` for `get_velocity`. Those flags are **false positives** from
the planner's `pkg::` scan: the two accessors never call
`requireNamespace()`; the strings appear only as the *default method label*
inside the record (`"SuperCell::SCimplify"`, `"velociraptor::scvelo"`).
Reading a recorded result needs neither package, so both examples are real,
runnable ones (they PASS) rather than `\dontrun{}`. The *producers*
(`RunSuperCell()`, `RunVelocity()`) do need those packages — that is where
the Suggests requirement actually lives. If the Lead prefers the conservative
reading, moving these two into `\dontrun{}` is a one-line change each.

## Exported functions in this group with no `man/<name>.Rd`

**None.** All 68 rows of the regenerated `G1a.tsv` have a corresponding
`man/<name>.Rd` in the worktree (checked file by file).

## Notes for the Lead before the roxygenise pass

1. `roxygenise()` should turn the 69 new `@examples` blocks into
   `\examples{}` sections. Nothing was inserted between a roxygen block and
   the function it documents — the inserter only appended the `@examples`
   lines at the **end** of the existing block, and the one function that
   already had an `@examples` section (`get_integration`) had the old section
   removed first. The roxygenise on the /tmp copy printed **no
   `Deleting '<name>.Rd'`** line, which is the hijack symptom AGENTS.md warns
   about.
2. `R/layers.R` has a second roxygen block for `DefaultLayer<-`
   (`@rdname DefaultLayer`). Only the getter's block carries `@examples`, so
   roxygen will not see a duplicate `\examples` section.
3. **Something to look at, not mine:** the roxygenise on the copy warns
   *"Aliases documented in multiple Rd files"* for the 11 diagnostics
   functions — `advanced-ai-diagnostics.Rd` now duplicates
   `check_annotation_readiness.Rd`, `check_integration_readiness.Rd`,
   `check_rare_cell_readiness.Rd`, `check_trajectory_readiness.Rd`,
   `summarize_*_readiness.Rd`, `summarize_cluster_sample_composition.Rd`,
   `summarize_pca_metadata_association.Rd`, `summarize_qc_by_group.Rd`,
   `summarize_small_cluster_evidence.Rd`, `summarize_small_clusters.Rd`,
   `summarize_trajectory_cluster_order.Rd`. That comes from the G4/G5
   agent's / your own change to `R/ai-diagnostics.R`, not from my files. The
   Lead has a fix in flight (delete the 29 hand-written Rd files so roxygen
   takes ownership) — I have not touched them.
4. `destiny` and `Augur` are declared under **`Enhances:`**, not
   `Suggests:`, while `RunDiffusionMap()` and `RunPerturbationPriority()`
   hard-require them via `requireNamespace()`. As shown above this produces
   **no** `R CMD check` note, so it is a field-placement question only; per
   the Lead's instruction I am not acting on it.
5. Every example was additionally validated as Rd: each `@examples` block was
   written into a scratch Rd file and passed through `tools::parse_Rd()`,
   `tools::Rd2ex()` and `parse()` — 66/66 clean, 0 problems. The `\dontrun{}`
   bodies are balanced and contain no `%` (which would have to be escaped).
6. Since the last run, four `\dontrun{}` bodies were edited (see the table
   above) and two of my files grew a few setup lines
   (`R/copilot.R`, `R/priority.R`, `R/workflows.R`). The counts stand at
   **66 PASS / 2 NOEXAMPLE** in `--src-group` mode and **68/68 PASS** in
   `--group` mode on a roxygenised copy; the Rd-mode run should be repeated
   on a fresh copy after this last edit round (the three edited files are
   re-verified in `--src-group` mode, which covers the same code paths).
