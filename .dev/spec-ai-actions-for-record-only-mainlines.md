# AIP Spec: close the "has a record, has no action" gaps

> Status: **rev 2, post-review — for implementation** (2026-10-23)
> Rev 1 was reviewed by a fact-check (`/tmp/gap-spec/facts.md`) and a design
> review (`/tmp/gap-spec/design.md`). Both ran code; rev 2 adopts their blocking
> corrections. Revision note at the end.
> Baseline: `devel`, package 1.1.0, `include = "all"` → **40** actions.
> Recorded baseline: AI suite 74 files / 6670 assertions green.

## 1. Problem

Eight exported mainlines write a `sclet` state record and no action reaches them,
so the AI can *see* the result in the ledger and never *produce* it:

| function | record | why it matters |
|---|---|---|
| `RunImputation` (`R/imputation.R:22`) | `preprocess::imputation` (id literal `"imputation"`, `:95`) | the package's own imputation step |
| `RunRegVelo` (`R/regvelo.R:42`) | `velocity::<name>` | `run_velocity` wraps `RunVelocity` (scVelo) only |
| `ClassifyMalignantEpithelium` (`R/malignancy.R:378`) | `malignancy::<name>` | only the older `CallMalignantCells` is reachable |
| `RunKNNPredict` (`R/annotation.R:363`) | `mapping::<name>` | reference mapping covers SingleR only |
| `RunSymphonyMapping` (`R/annotation.R:524`) | `mapping::<name>` | same |
| `RunSpatialDeconvolution` (`R/spatial.R:26`) | `spatial::<name>` | the `spatial` group has **no executing action** |
| `RunSpatialColocalization` (`R/spatial.R:193`) | `spatial::<name>` | same |
| `RunSpatialNiche` (`R/spatial.R:484`) | `spatial::<name>` | same |

Two of these are **not** a registration gap, they are callee defects that make the
action impossible to write honestly. Rev 1 proposed wrapping them; the reviews
proved the wrapper cannot succeed (verified here):

- **`RunSpatialColocalization` / `RunSpatialNiche` reject gene features.** Both
  require `features %in% colnames(colData(sce))` (`R/spatial.R:223-231`,
  `:510-517`) yet hand that vector to `SVP::runGLOBALBV(features1=)`, which
  indexes `assay` rows. Measured: a real `SpatialExperiment` with gene features is
  refused with "Feature column(s) not found in colData"; with colData columns it
  dies inside SVP with `subscript out of bounds`. **No input shape works.**
  Fix the callee's feature contract first, then register the action.
- **`ai-data-support.R:508`** suggests `run_spatial_deconvolution` while checking
  `has_analysis_type("spatial_deconvolution")`, but deconvolution records are
  typed `"spatial"` — so the suggestion fires forever. Same class of defect.

`RunDiffusionMap` and `AggregateExpression` write no state record and are outside
this contract.

## 2. Scope, restated after review

Three tiers, so the milestone ships what is true and states the rest:

**Tier A — register now (5 actions).** The callee works and the input is
expressible: `run_imputation`, `classify_malignant_epithelium`,
`run_spatial_deconvolution` (driver), `map_reference_knn` (driver),
`map_reference_symphony` (driver).

**Tier B — register, gated on a solver input (2 actions).** `run_regvelo` needs a
`grn` object. The privacy layer denies `*_path`/`*_file` field names
(`R/ai-privacy.R:388-396`), so a path parameter cannot reach the planner either;
the honest answer is a driver-only action plus a documented way to produce the
prior (`RunSCENIC` + `RegVelo`'s own prior builder).

**Tier C — fix the callee, then register (2 actions).** The two SVP spatial
analyses. **Not in this milestone's acceptance**; §6 records the exact repair so
it is not lost.

## 3. Design — Tier A and B

### 3.1 Where each action lives

| group | actions | file |
|---|---|---|
| `preprocess` | `run_imputation` | `R/ai-registry-preprocess.R` |
| `malignant` | `classify_malignant_epithelium` | `R/ai-registry-malignant.R` |
| `annotation` | `map_reference_knn`, `map_reference_symphony` | `R/ai-registry-annotation.R` |
| `spatial` | `run_spatial_deconvolution` | new `R/ai-registry-spatial-actions.R` |
| `velocity` | `run_regvelo` | `R/ai-registry-trajectory-velocity-fate.R` |

`preprocess` currently has **4** actions (`normalize_data`,
`find_variable_features`, `scale_data`, `run_decontx`); it becomes 5. The new
spatial file must be wired into the composition
(`R/ai-execution.R:618-625`), and it is named `-spatial-actions.R` to distinguish
it from the read-only `R/ai-registry-spatial-multimodal.R` (whose read-only
property is per-action `mutates_object = FALSE`, not the file header).

### 3.2 Actions

| id | required params | optional | record | planner_visible |
|---|---|---|---|---|
| `run_imputation` | — | `method` (ALRA), `assay` | `preprocess::imputation` | TRUE |
| `classify_malignant_epithelium` | `seed_tumor`, `seed_normal` (`character_vector`) | `method`, `sample_col`, `n_top`, `padj`, `name` | `malignancy::<name>` | TRUE |
| `map_reference_knn` | `ref` (`object`) | `labels`, `features`, `k`, `name` | `mapping::<name>` | FALSE |
| `map_reference_symphony` | `ref` (`object`) | `labels`, `vars`, `k`, `name` | `mapping::<name>` | FALSE |
| `run_spatial_deconvolution` | `sce_ref` (`object`), `ref_group_key` | `ref_batch_key`, `max_epochs`, `name` | `spatial::<name>` | FALSE |
| `run_regvelo` | `grn` (`object`) | `spliced_assay`, `unspliced_assay`, `reduction`, `max_epochs`, `lam`, `lam2`, `name` | `velocity::<name>` | FALSE |

**The malignant seed is two vectors, not one.** `ClassifyMalignantEpithelium`'s
`seed` accepts either `FindMalignantSeedGenes()` output or a named list with a
`tumor`/`malignant` and a `normal` element (`R/malignancy.R:433-470`); a flat
character vector matches **neither** and would be rejected. The action takes
`seed_tumor` and `seed_normal` as two required parameters and assembles
`list(tumor = , normal = )` — explicit, and impossible to get half-right.

### 3.3 Prerequisites

| action | gates |
|---|---|
| `run_imputation` | `assay_missing` (named `assay` absent) |
| `classify_malignant_epithelium` | `seed_missing` (either vector empty), `seed_not_found` (no seed gene in `rownames`), `sample_column_missing`, `too_few_cells_for_split` (`ncol < 2 * min_cells_per_group`), `optional_package_missing: limma` |
| `map_reference_knn` | `reference_missing`, `reference_counts_missing`, `annotation_levels_missing`, `optional_package_missing: BiocNeighbors` |
| `map_reference_symphony` | `reference_missing`, `annotation_levels_missing`, `optional_package_missing: symphony` |
| `run_spatial_deconvolution` | `optional_package_missing: basilisk`, `counts_assay_missing` (spatial and reference), `reference_group_missing`, `no_common_genes` |
| `run_regvelo` | `spliced_unspliced_missing`, `grn_invalid` (not a list/matrix carrying `regulators` and `W`), `optional_package_missing: basilisk`, `reduction_missing` (honouring `planned$reductions`) |

No `spatial_coordinates_missing`: that condition is already reported by
`check_spatial_readiness` and `sclet_ai_spatial_coordinates()`
(`R/ai-spatial.R:78-92`), and the actions that would need it are Tier C.

### 3.4 Evidence (one node, `ev:` ids, numbers/logicals only)

| action | values | claim_level | note |
|---|---|---|---|
| `run_imputation` | `n_cells`, `n_genes`, `assay_imputed = TRUE` | measured | the record's `summary` is empty; sizes come from the object |
| `classify_malignant_epithelium` | `n_malignant`, `n_non_malignant`, `malignant_fraction`, `iterations`, `converged` | associated | id `ev:malignant_iterative_<name>` — the existing `call_malignant_cells` node is `ev:malignant_<name>`, and two calls must not share one id |
| `map_reference_knn` / `_symphony` | `n_labels`, `n_cells`, `n_reference_labels` | measured | `n_reference_labels` from the record's `n_reference_labels` when present, else `length(unique(labels))` |
| `run_spatial_deconvolution` | `n_common_genes`, `n_columns`, `n_cells` | measured | the record has **no** `summary`; `n_common_genes` from `inputs`, columns from `artifacts$columns` |
| `run_regvelo` | `n_targets`, `n_regulators`, `n_edges` | measured | the GRN counts live in `inputs$grn`, not `summary` |

### 3.5 Planner visibility: four hidden, two visible

Visible: `run_imputation` (no external input), `classify_malignant_epithelium`
(the two seed gene vectors are user choices the planner may carry). Hidden:
`map_reference_knn`, `map_reference_symphony` (a whole second SCE), `run_regvelo`
(a GRN object), `run_spatial_deconvolution` (a second SCE plus a group column).
Rev 1 had KNN visible; that was wrong — it needs the same `ref` object as
Symphony.

### 3.6 Planning tables

`sclet_ai_goal_action_map()` gains **one** entry, `run_imputation = "impute|ALRA"`.
`classify_malignant_epithelium` gets **no** entry: the existing
`call_malignant_cells = "malignant|malignancy|tumou?r cells?"` already matches, and
a second pattern would make the coverage judge demand both actions
(`R/ai-planning.R:1650-1672`). `map_reference_*` and `run_regvelo` are hidden, so
they get no entry either.

`sclet_ai_planned_state_records()` gains `run_imputation` (`preprocess::
"imputation"`) and `classify_malignant_epithelium` (`malignancy::<name>`), so
`run_imputation → find_variable_features …` validates in one plan.

### 3.7 Existing text that becomes wrong (must be updated in the same change)

- `R/ai-data-support.R:605-618` asserts "There is still no imputation action, and
  this table does not pretend otherwise", and the remedy table lacks it. Add
  `ambient`/`imputation` remedies once `run_imputation` exists.
- `R/ai-data-support.R:508` suggests `run_spatial_deconvolution` against
  `has_analysis_type("spatial_deconvolution")`, a type never written. Change the
  lookup to `"spatial"` (with the deconvolution `analysis_key`), so the
  suggestion means what it says.
- `inst/skills/sclet-analysis/SKILL.md` "Five mainline actions … are absent from
  available_actions" becomes "Six …" and names the new driver-only ids.

## 4. Non-goals

- **No Tier C in this milestone.** The SVP callee's feature contract is its own
  fix (§6), because a wrapper over a callee that rejects both input shapes cannot
  be tested honestly.
- No multimodal registration actions (`RunMultimodalWorkflow`,
  `register_multimodal_*)`: they move altExps the object already carries, so an
  evidence node would report a bookkeeping count with no analysis behind it. The
  gap is recorded, not hidden; if the read-only multimodal actions need a creator,
  that is a separate, cheaper design (an evidence-free registration action).
- No new readiness algorithms; gates read what exists.
- No Python execution in tests (`sclet_basilisk_run` mocked for RegVelo and
  cell2location).
- No change to `RunAIAnalysis()`'s default `include`.

## 5. Acceptance criteria

- **A1** `include = "all"` yields **46** actions (40 + 6); each new action is in
  exactly its own group and no other; `include = "preprocess"` yields 5.
- **A2** All 46 are discoverable through `sclet_ai_action_groups()`.
- **A3** The planner catalog contains `run_imputation` and
  `classify_malignant_epithelium` and none of the four driver-only new actions.
- **A4** Each action refuses a fixture lacking its input, naming the gate; no
  handler runs, no record, no evidence.
- **A5** Each action, on a fixture that satisfies its gates (real for imputation,
  `sclet_basilisk_run`-mocked for RegVelo/cell2location, mocked callee for the
  two mappings and malignant), writes its record and exactly one evidence node
  with numeric/logical values, `source` = record id, claim level per §3.4.
- **A6** A driver plan per action validates, executes and lands in
  `AnalysisGraph()` with a `produced` edge to the record and to the evidence node.
- **A7** The seed is two required vectors: a plan supplying one, or a flat
  vector, is refused with `seed_missing`/`seed_invalid`.
- **A8** Skill text updated and pinned: six driver-only ids, seed/reference
  inputs listed as user-owned.
- **A9** No existing test changes expectation.
  **This is false for two tests and they must be updated with the change:**
  `test-ai-registry-mainlines.R:107` pins `length(all_actions) == 40L`;
  `test-ai-execution.R:527-547` asserts that **no** `run_*` action in the
  spatial/multimodal/regvelo domains exists (its own comment says those domains
  were paused) — adding `run_regvelo` and `run_spatial_deconvolution` makes it
  fail by design, so the assertion becomes "no `run_*` action in those domains
  **without** `planner_visible`/`mutates_object` declarations", i.e. it keeps
  testing that the paused domains were not silently shipped an *unreviewed*
  mutation action, and names the two that this milestone reviewed.
  `test-ai-available-actions.R`'s hand-written `hidden` list gains the four ids.

## 6. Tier C repair, recorded for the next change

`RunSpatialColocalization` / `RunSpatialNiche` must accept **assay gene
features** (the SVP contract: `features1` are rows of `logcounts`), on a
`SpatialExperiment` or an object with a `spatial` reducedDim. Concretely:

1. Resolve features from `rownames` (genes), not `colData`; keep the deconvolution
   columns path only for abundance-based co-localisation, and pass those through
   `SVP`'s `gsvaexp`/`group.by`, not `features1`.
2. Gate on `sclet_ai_spatial_coordinates()` (already exists) instead of a
   colData-name check.
3. Then register `run_spatial_colocalization` (claim level `associated`) and
   `run_spatial_niche` (`measured`), both visible, and give
   `sclet_ai_planned_state_records()` a `run_spatial_deconvolution` entry so
   `deconvolution → colocalization` validates in one plan.

## 7. Tests

New `tests/testthat/test-ai-registry-record-only-actions.R`:

| case | asserts |
|---|---|
| registry shape | A1, A2 |
| planner catalog | A3 |
| gates | A4, one case per action |
| imputation real | A5 (ALRA on a normalized fixture) |
| regvelo + cell2location mocked | A5 via `sclet_basilisk_run` mock |
| mappings + malignant mocked | A5 with a callee stub that writes the real record shape |
| seed contract | A7: one vector, a flat vector, and both vectors |
| plan round trips | A6, one per action |
| skills | A8, in `test-ai-skills.R` |
| updated expectations | A9, in the two files named there |

The malignant-calibration logic (limma, iterative k-means) is **not** exercised:
the mock writes the record the real function writes. The test comment says so, so
a green run is not read as "the classifier works".

## 8. Risks

- **`grn` is an object**: a driver must build it (`RunSCENIC` output or
  `RegVelo`'s prior builder). The action's description will say which.
- **Two mappings, two `ref` objects**: both hidden, so the planner cannot
  propose an action it cannot fill.
- **`ev:malignant_*` collision** handled by the `_iterative_` infix.
- **Group counts** change (`preprocess` 4→5, `all` 40→46); the two pins that
  hard-code them are updated deliberately, not weakened.

## 10. Revision note (rev 3 — post-implementation review)

- **RegVelo's prior contract was backwards.** `RunRegVelo` accepts a
  matrix/Matrix (targets × regulators, with dimnames) or a data.frame with
  `target`/`regulator(/weight)`; a named list is rejected outright
  (`R/regvelo.R:260-313`). The gate now accepts exactly those two shapes and the
  test asserts it accepts what the callee accepts and rejects a named list.
- **`RunImputation` writes `alra_imputed`, not `ALRA`** (`R/imputation.R:60`).
  Schema, description and the test now use the real name, and the test runs
  ALRA for real (mocking it is what hid the wrong name).
- `run_imputation` forwards `k`: ALRA's automatic rank selection can pick a `k`
  that does not fit the matrix and `choose_k` fails on any shape without one
  (measured on both a tall and a wide fixture).
- `run_imputation` refreshes its evidence node (`replace = TRUE`) instead of
  raising `evidence_id_conflict` on a re-run with changed values.
- `map_reference_symphony` evidences `n_query_cells`/`n_ref_cells`, which is what
  the real record carries (`R/annotation.R:660-663`); `n_labels` was a stub
  artifact.
- The stale "there is still no imputation action" comment in
  `R/ai-data-support.R` is corrected.
- Baseline: `include = "all"` was 34 at HEAD and is 46 now (12 added across the
  two milestones: 6 mainline + 6 record-only).

## 9. Open questions for review

1. `run_imputation` in `preprocess` vs its own `imputation` group? (Recommend
   `preprocess`: it writes `type = "preprocess"` and is a one-step assay
   transform like `run_decontx`.)
2. Should `classify_malignant_epithelium` be visible at all, given the seed
   vectors? (Recommend yes: the genes are user choices the planner can carry,
   unlike a whole reference object; and the alternative is a visible gap next to
   `call_malignant_cells`.)
3. RegVelo driver-only, or also accept `grn` as a file handle through a
   non-path-named parameter (e.g. `grn_source`)? (Recommend driver-only for now —
   a renamed path field is a privacy workaround, not a design.)
