# G5 report — Rd examples for the AI domain adapters

Work list: `.dev/tools/groups/G5.tsv` (27 Rd — 26 plus the `advanced-ai-diagnostics.Rd`
topic page that was added to the group when the TSVs were regenerated).
Sources edited: `R/ai-communication.R`, `R/ai-comparison.R`, `R/ai-context.R`,
`R/ai-diagnostics.R`, `R/ai-enrichment.R`, `R/ai-fate-drivers.R`, `R/ai-fate.R`,
`R/ai-integration-routes.R`, `R/ai-interpretation-gaps.R`, `R/ai-investigate.R`,
`R/ai-multimodal.R`, `R/ai-profile.R`, `R/ai-spatial.R`,
`R/ai-trajectory-root-comparison.R`, `R/ai-trajectory-velocity-crosscheck.R` — and nothing else.
`man/*.Rd`, `NAMESPACE` and `DESCRIPTION` were **not** touched. No `roxygenise()` was run in
this worktree (a throw-away copy under `/tmp` was used to check what roxygen would do).

## Verification

```
export R_LIBS_USER=~/Rlib
cd /home/wang/data/source/omics/sclet-wt
Rscript .dev/tools/check-examples.R --src-group G5   # 27 / 27 PASS  (0 FAIL, 0 NOEXAMPLE, 0 SKIP)
Rscript .dev/tools/lint-example-symbols.R --src-group G5   # no undefined symbols found in 27 examples
Rscript .dev/tools/lint-example-placement.R          # none misplaced
```

Whole group runs in ~19 s wall clock including package loading.

## One line per file

| Rd | example | status | what the example shows |
|---|---|---|---|
| `AIInvestigate.Rd` | `R/ai-investigate.R` | PASS | `mockSCE` setup outside, the model call inside `\dontrun{}` (provider key). |
| `CompareAIAnalyses.Rd` | `R/ai-comparison.R` | PASS | **Working demo**: records a real integration route first, then compares what is on record (`status = "available"`, `baseline = "raw"`). |
| `GetAIProfile.Rd` | `R/ai-profile.R` | PASS | Builds the bounded, value-free dataset profile and reads its cell count and schema version. |
| `GetAnalysisLedger.Rd` | `R/ai-context.R` | PASS | Builds the bounded ledger and reads `dataset$n_cells` and the object fingerprint. |
| `RunIntegrationRoutes.Rd` | `R/ai-integration-routes.R` | PASS | **Working demo**: runs the raw baseline route under a confirmed design and reports the executed route. |
| `advanced-ai-diagnostics.Rd` | `R/ai-diagnostics.R` | PASS | The topic page's example lives on `summarize_qc_by_group()` (the page's `\name`): adds a numeric QC column, groups by `Treatment`, and reads `status` and the group count. See the caveat below. |
| `check_communication_readiness.Rd` | `R/ai-communication.R` | PASS | Readiness report: `status` plus the question naming the missing group column and the available `colData` columns. |
| `check_de_readiness.Rd` | `R/ai-interpretation-gaps.R` | PASS | `mockSCE` setup outside, the call inside `\dontrun{}` (limma Suggests). |
| `check_enrichment_readiness.Rd` | `R/ai-enrichment.R` | PASS | `mockSCE` setup outside, the call inside `\dontrun{}` (clusterProfiler Suggests). |
| `check_fate_drivers_readiness.Rd` | `R/ai-fate-drivers.R` | PASS | Readiness report: `status` plus the question naming the missing CellRank trajectory record. |
| `check_fate_readiness.Rd` | `R/ai-fate.R` | PASS | Readiness report: `status` plus the question naming the missing velocity record. |
| `check_milo_readiness.Rd` | `R/ai-interpretation-gaps.R` | PASS | `mockSCE` setup outside, the call inside `\dontrun{}` (miloR Suggests). |
| `check_multimodal_readiness.Rd` | `R/ai-multimodal.R` | PASS | Readiness report: `status` plus the question naming the unattached ADT modality. |
| `check_spatial_readiness.Rd` | `R/ai-spatial.R` | PASS | Readiness report: `status` plus the question naming the missing spatial analysis. |
| `check_velocity_readiness.Rd` | `R/ai-diagnostics.R` | PASS | Readiness report: `status` plus the question naming the missing `spliced`/`unspliced` assays. |
| `compare_fate_lineages.Rd` | `R/ai-fate-drivers.R` | PASS | Typed contract: `status = "not_available"`, `reason = "at_least_two_lineages_required"`. |
| `compare_multimodal_evidence.Rd` | `R/ai-multimodal.R` | PASS | Typed contract: `status = "not_available"`, `reason = "at_least_two_modalities_required"`. |
| `compare_rare_cell_evidence.Rd` | `R/ai-comparison.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_rare_cell_evidence_records"`. |
| `compare_spatial_analyses.Rd` | `R/ai-spatial.R` | PASS | Typed contract: `status = "not_available"`, `reason = "at_least_two_analyses_required"`. |
| `compare_trajectory_roots.Rd` | `R/ai-trajectory-root-comparison.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_trajectory_runs"`. |
| `compare_trajectory_velocity_evidence.Rd` | `R/ai-trajectory-velocity-crosscheck.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_trajectory_or_velocity_evidence_records"`. |
| `summarize_communication_evidence.Rd` | `R/ai-communication.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_communication_record"`. |
| `summarize_enrichment_evidence.Rd` | `R/ai-enrichment.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_enrichment_record"`. |
| `summarize_fate_lineage_drivers.Rd` | `R/ai-fate-drivers.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_cellrank_record"`; no lineage is chosen for the user. |
| `summarize_milo_evidence.Rd` | `R/ai-interpretation-gaps.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_milo_record"`. |
| `summarize_multimodal_evidence.Rd` | `R/ai-multimodal.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_multimodal_record"`. |
| `summarize_spatial_evidence.Rd` | `R/ai-spatial.R` | PASS | Typed contract: `status = "not_available"`, `reason = "no_spatial_record"`. |

## Examples that need a provider API key (wrapped in `\dontrun{}`)

| function | why |
|---|---|
| `AIInvestigate()` | Model call through the AI adapter. |

`mockSCE` setup sits outside the `\dontrun{}`; the reason line inside is
`# needs a provider API key`.

## Examples that need a Suggests package (wrapped in `\dontrun{}`)

Per the work list's `suggests_needed` column, these three bodies call `requireNamespace()` on
a Suggests package, so the example is not runnable where that package is absent:

| function | Suggests package | reason line inside the `\dontrun{}` |
|---|---|---|
| `check_de_readiness()` | `limma` | `# needs the limma Suggests package` |
| `check_enrichment_readiness()` | `clusterProfiler` | `# needs the clusterProfiler Suggests package` |
| `check_milo_readiness()` | `miloR` | `# needs the miloR Suggests package` |

**Honest caveat:** all three actually *degrade gracefully* — they call
`requireNamespace(..., quietly = TRUE)` and report the package's availability as a check
(`de_package_available`, `enrichment_packages_available`, `milo_available`) rather than erroring.
They would run on a machine without the Suggests package. I wrapped them anyway because the
work list flags them, and a reader who copies the example on a machine without the package
gets a different `status`/checks than the one shown. If the Lead prefers runnable pages, moving
the call out of the `\dontrun{}` is a one-line change per function.

## Examples that are a real, runnable demo

* `CompareAIAnalyses()` — records a real integration route (`RunIntegrationRoutes(...,
  routes = "raw", confirm = "always")`) and then compares what is on record: `status =
  "available"`, `baseline = "raw"`. The prerequisite state is genuinely built.
* `RunIntegrationRoutes()` — executes the raw baseline route under a confirmed design and
  reports `status = "completed"` and the executed route. The raw route needs no integration
  package, so this runs anywhere; the other routes (`fastMNN`, `Harmony`, `scVI`) are
  previewed instead of executed (the example says so in a comment).
* `GetAnalysisLedger()`, `GetAIProfile()` — read-only, no state needed.
* `summarize_qc_by_group()` — adds one numeric `colData` column to the `mockSCE` and gets a
  real `"available"` summary with two anonymized groups.

## Examples that are a **shape demonstration**, not a working demo

These are honest about it: the function reads a *recorded analysis result*, and the record
can only be produced by an analysis that needs a Suggests package (slingshot, velociraptor,
CellRank, CellChat, clusterProfiler, miloR, fastMNN, …) or a long pipeline, so a portable
example cannot build the state. The examples show the functions' documented **typed**
`not_available` contract — a status and a reason, never `NULL`, never a crash — and each
comment says what is missing. They are not NULL-returning no-ops.

`compare_rare_cell_evidence()`, `compare_trajectory_roots()`,
`compare_trajectory_velocity_evidence()`, `compare_fate_lineages()`,
`compare_multimodal_evidence()`, `compare_spatial_analyses()`,
`summarize_communication_evidence()`, `summarize_enrichment_evidence()`,
`summarize_fate_lineage_drivers()`, `summarize_milo_evidence()`,
`summarize_multimodal_evidence()`, `summarize_spatial_evidence()`,
`summarize_pca_metadata_association()`, `summarize_cluster_sample_composition()`,
`summarize_small_clusters()`, `summarize_small_cluster_evidence()`,
`summarize_trajectory_cluster_order()`.

The `check_*_readiness()` examples are a third category: on an object with nothing recorded
they answer `"not_ready"` / `"clarification_required"` **and name the missing prerequisite in
`questions`** — that question list is the function's product, so the example is meaningful,
not a placeholder.

## Caveats the Lead must know

### 1. 25 of the 27 G5 Rd files are HAND-WRITTEN — roxygen will skip them

`man/AIInvestigate.Rd`, `man/CompareAIAnalyses.Rd`, `man/RunIntegrationRoutes.Rd`,
`man/advanced-ai-diagnostics.Rd`, `man/check_communication_readiness.Rd`,
`man/check_de_readiness.Rd`, `man/check_enrichment_readiness.Rd`,
`man/check_fate_drivers_readiness.Rd`, `man/check_fate_readiness.Rd`,
`man/check_milo_readiness.Rd`, `man/check_multimodal_readiness.Rd`,
`man/check_spatial_readiness.Rd`, `man/check_velocity_readiness.Rd`,
`man/compare_fate_lineages.Rd`, `man/compare_multimodal_evidence.Rd`,
`man/compare_rare_cell_evidence.Rd`, `man/compare_spatial_analyses.Rd`,
`man/compare_trajectory_roots.Rd`, `man/compare_trajectory_velocity_evidence.Rd`,
`man/summarize_communication_evidence.Rd`, `man/summarize_enrichment_evidence.Rd`,
`man/summarize_fate_lineage_drivers.Rd`, `man/summarize_milo_evidence.Rd`,
`man/summarize_multimodal_evidence.Rd`, `man/summarize_spatial_evidence.Rd`.

They have no `% Generated by roxygen2` header, so `roxygen2::roxygenise()` reports *"already
exists and was not generated by roxygen2"* and **leaves them alone**. My `@examples` are in
the source and verified, but they will not reach the shipped pages until those files are
deleted once so roxygen takes ownership — which is what the project's own `AGENTS.md` says to
do with a hand-written Rd.

Verified in a scratch copy: after deleting the 25 hand-written files and re-running
`roxygenise()`, all 25 pages are regenerated from the roxygen blocks **with** the examples,
and `Rscript .dev/tools/check-examples.R --group G5` reports **26 PASS** (the 26th row is
`advanced-ai-diagnostics.Rd`, see next point). Nothing else in the group regressed.

### 2. `man/advanced-ai-diagnostics.Rd` is a hand-written page for 11 exported functions

`R/ai-diagnostics.R` documents 11 exported functions (`summarize_qc_by_group()`,
`summarize_pca_metadata_association()`, `summarize_cluster_sample_composition()`,
`summarize_small_clusters()`, `check_integration_readiness()`,
`check_annotation_readiness()`, `check_rare_cell_readiness()`,
`summarize_small_cluster_evidence()`, `check_trajectory_readiness()`,
`summarize_trajectory_cluster_order()`, `check_velocity_readiness()`), each with its own
roxygen block and `@export`. The shipped page is a hand-written `man/advanced-ai-diagnostics.Rd`
that gives all 11 the same `\name`/`\alias` set.

Two consequences the Lead has to decide about:

* `roxygenise()` writes **11 per-function pages** for these (it already does — that is a
  pre-existing condition, not caused by my change: `HEAD:R/ai-diagnostics.R` already carried
  all 11 `@export` tags), while the hand-written page stays. The result is
  `Rd files with duplicated alias` in `R CMD check` — I reproduced the exact warning on a
  minimal package, so it is a **WARNING, not an ERROR** (the build does not break), but it is
  a real defect.
* Because the page is hand-written, it will never receive an `\examples{}` section.

I added `@examples` to **all 11** roxygen blocks in `R/ai-diagnostics.R` (not just the one the
TSV points at), so that whichever way the page is resolved — one combined page or 11
per-function pages — every exported function on it carries a runnable example. The TSV row
`advanced-ai-diagnostics.Rd` (name `summarize_qc_by_group`) is satisfied by the
`summarize_qc_by_group()` block.

Recommended resolution for the Lead: delete `man/advanced-ai-diagnostics.Rd` (hand-written, so
roxygen never regenerates it) and let the 11 roxygen blocks own their pages. Then
`man/advanced-ai-diagnostics.Rd` no longer exists and the TSV row for it has no file — the
group's row would need to be re-pointed at the 11 pages (or at `summarize_qc_by_group.Rd`).
The `--group G5` run then reports that one row `MISSING` and 26 `PASS`, which is the honest
state.

### 3. A literal `%` breaks the Rd example text — avoid `%in%` in examples

I originally wrote `"normalize_data" %in% names(registry)`. Roxygen correctly escapes it to
`\%in\%` in the Rd, and `R CMD check` un-escapes it, so the shipped page is fine — but the
local checker's `--group` mode reads the Rd as plain text and then fails to parse
`"normalize_data" \%in\% ...`. I replaced both occurrences with
`any(names(registry) == "normalize_data")` and `any(ids == "ai_copilot_example")`.
**Worth telling the other three agents**: avoid `%` in examples, or the `--group` checker
fails on a page that is actually fine.

### 4. Exported functions in G5 with no `man/<name>.Rd` at all

**None.** All 27 rows in `G5.tsv` have a tracked `man/<name>.Rd`.
