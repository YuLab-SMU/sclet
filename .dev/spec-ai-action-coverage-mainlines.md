# AIP Spec: AI actions for the six uncovered mainlines

> Status: **rev 2, post-review — for implementation** (2026-10-23)
> Rev 1 was reviewed by a fact-check (`/tmp/coverage-spec/facts.md`) and a design
> review (`/tmp/coverage-spec/design.md`). Rev 2 adopts their blocking
> corrections; the revision note at the end lists what changed.
> Baseline: `devel`, package 1.1.0. `AIDefaultExecutionRegistry(sce, include =
> "all")` has **34** actions (measured); these six mainlines have **0**.

## 1. Problem

Six exported mainlines have no AI action, so a plan naming them fails with
"action is not registered": `RunGeneSetScoring` (`R/geneset_scoring.R:30`),
`RunSuperCell` (`R/SuperCell.R:26`), `RunPerturbationPriority`
(`R/priority.R:54`), `RunPhenotypeAssociation` (`R/phenotype.R:23`),
`RunSCENIC` (`R/scenic.R:41`), `RunInSilicoPerturbation`
(`R/celloracle.R:26`). (Two of the six — gene-set scoring and phenotype
association — have no test either; that gap is not this spec's subject but it is
recorded here so nobody reads "exported" as "exercised".)

All six state types already exist (`R/state.R:686-721`: `aggregation`,
`scenic`, `geneset_scoring`, `perturbation`, `priority`, `phenotype`), so no
state-schema change is needed. The group vocabulary, however, is hardcoded in
**five** places (§4.1), none of which has a slot for them.

## 2. Feasibility, one row per mainline

| mainline | function | input the object cannot supply | locally runnable? |
|---|---|---|---|
| gene-set scoring | `RunGeneSetScoring` | `gene_sets` (named list) — a user modelling choice | yes (UCell/AUCell/GSVA) |
| metacells | `RunSuperCell` | nothing beyond a **`logcounts`** assay | yes (SuperCell) |
| perturbation priority | `RunPerturbationPriority` | a `condition_col` naming a colData column with ≥2 levels | yes (Augur) |
| phenotype association | `RunPhenotypeAssociation` | `bulk_matrix` + `phenotype_data` (external matrices) | **no** (`scPAS` absent) |
| regulon inference | `RunSCENIC` | three file inputs (TFs, motif annotations, databases) | no (Python `sclet_scenic_env`) |
| in-silico perturbation | `RunInSilicoPerturbation` | `target_gene` + `base_grn_path` (no default) | no (Python `sclet_celloracle_env`) |

Two consequences the design respects:

1. **The planner cannot invent any of these inputs.** A gene set, a target
   gene, a base GRN, a phenotype matrix or a database path is the user's
   scientific choice. Every such input is a **required parameter with no
   default**, and the gate or the schema names it.
2. **Five of the six are driver-mode only** (driver mode is documented in
   `inst/skills/sclet-driver/SKILL.md`). The in-package planner sees a catalog
   built from the registry, and `AIPlanAnalysis()` defaults to `include =
   "all"` (`R/ai-planning.R:1992`) — so registering them without a visibility
   rule would offer the model actions it can never fill and grow
   `available_actions` from 34 to 40 (see §4.2).

## 3. What "coverage" means here

The six become **registered, gated, recorded and evidence-bearing** actions,
executable by the driver path and by an explicit `ValidateAIPlan()`/
`ExecuteAIPlan()` call with a full registry. Only
`run_perturbation_priority` is also offered to the in-package planner, because
its one external input is a colData column name the planner can read from
`capabilities`. The other five stay out of the planner catalog until a
clarification mechanism can collect their inputs; pretending otherwise would
recreate the failure class this milestone exists to remove.

## 4. Design

### 4.1 The five tables that enumerate groups

| table | location |
|---|---|
| `allowed_groups` validator | `R/ai-execution.R:587` |
| `sclet_ai_action_groups()` hint vector | `R/ai-execution.R:627-631` |
| roxygen `@param include` list → `man/AIDefaultExecutionRegistry.Rd` | `R/ai-execution.R:571-574` |
| registry composition | `R/ai-execution.R:598-610` |
| `RunAIAnalysis()` default `include` | `R/ai-functions.R:401` (must **not** grow) |

New groups: `geneset`, `metacell`, `priority`, `phenotype`, `regulon`,
`perturbation`. The new registry block `R/ai-registry-mainlines.R` is appended
last, so existing action order is unchanged.

Two more name-keyed tables need entries, both in `R/ai-planning.R`:
`sclet_ai_goal_action_map()` (`:1646`) and `sclet_ai_planned_state_records()`
(`:316`).

### 4.2 Planner visibility (new, small)

`AIAction()` gains `planner_visible = TRUE`; `sclet_ai_describe_registry()`
(`R/ai-context.R:814`) skips descriptors with `planner_visible = FALSE`.
Validation and execution are untouched: a driver plan can use every action.

Marked `planner_visible = FALSE`: `run_geneset_scoring` (user-owned sets),
`run_supercell` (replaces the object, see §4.5), `run_phenotype_association`
(matrices), `run_scenic` and `run_in_silico_perturbation` (paths — the privacy
classifier denies `*_path`/`*_file` field names, `R/ai-privacy.R:357`, so their
required parameters cannot reach the model at all; keeping them out of the
catalog also keeps the zero-redaction pin in `test-ai-available-actions.R`
true).

`run_perturbation_priority` stays visible.

### 4.3 Actions

All six: `returns = "sce"`, `mutates_object = TRUE`,
`requires_confirmation = TRUE`, `idempotent = FALSE` (each writes an evidence
node under a stable id; the recorder refuses changed values without
`replace = TRUE`, so a retry must not be treated as free).

| registry id | group | required params | optional params | record type / id |
|---|---|---|---|---|
| `run_geneset_scoring` | geneset | `gene_sets` (list) | `method` ∈ UCell/AUCell/GSVA, `assay_use`, `ncores`, `name`, `as_altExp` | `geneset_scoring` / `name` (default `"Score"`) |
| `run_supercell` | metacell | — | `assay`, `layer`, `nHVG`, `hvg_method`, `gamma`, `k.knn`, `name` | `aggregation` / `name` (default `"supercell"`) |
| `run_perturbation_priority` | priority | `condition_col` | `group_col`, `assay`, `layer`, `name` | `priority` / `name` (default `"augur"`) |
| `run_phenotype_association` | phenotype | `bulk_matrix`, `phenotype_data` (both `object`) | `family` | `phenotype` / `"scPAS"` |
| `run_scenic` | regulon | `tfs_path`, `motif_annotations_path`, `database_paths` | `assay_use`, `num_workers`, `seed`, `name`, `output_dir` | `scenic` / `name` (default `"scenic"`) |
| `run_in_silico_perturbation` | perturbation | `target_gene`, `base_grn_path` | `perturbation_value`, `cluster_key`, `reduction` | `perturbation` / `celloracle_<gene>` |

Cost: geneset/supercell `low`; priority/phenotype `medium`; scenic/oracle
`high`. `planner_visible` as in §4.2.

### 4.4 Prerequisites

Each returns `"<gate_name>: <what is missing and what to do>"`.

| action | gates |
|---|---|
| geneset | `gene_sets_empty`, `gene_sets_not_found` (no member in `rownames`), `geneset_package_missing` (the selected method's package), `counts_assay_missing` |
| supercell | `logcounts_assay_missing` (or the named `assay` absent), `optional_package_missing: SuperCell`, `too_few_metacells` (`ncol(object) < 2 * gamma`; `gamma` is cells per metacell, `R/SuperCell.R:14`) |
| priority | `condition_col_missing`, `condition_levels_insufficient` (<2 distinct values, `R/priority.R:96-100`), `group_column_missing` (`group_col` absent and no `colLabels`), `optional_package_missing: Augur` |
| phenotype | `optional_package_missing: scPAS`, `optional_package_missing: Seurat` (`as.Seurat()`, `R/phenotype.R:35`), `bulk_matrix_missing`, `phenotype_data_missing` |
| scenic | `optional_package_missing: basilisk`, `counts_assay_missing`, `scenic_files_missing` (each of the three), `scenic_database_missing` (empty `database_paths`) |
| perturbation | `optional_package_missing: basilisk`, `counts_assay_missing`, `cluster_identity_missing` (no explicit `cluster_key`, no `ActiveIdent`), `reduction_missing`, `base_grn_missing`, `target_gene_not_found` |

### 4.5 `run_supercell` replaces the object

Measured: 60 cells → 6 metacells, and the state records went from
`preprocess` to `aggregation` **only** — `RunSuperCell` rebuilds the internal
state (`sclet_rebuild_internal_state()`, `R/state.R:150-165`). In a plan it
therefore discards the ledger and every earlier record. Hence
`planner_visible = FALSE`, a description that says so in the first sentence, and
a note in the skill. The action is still worth having for the driver path (a
metacell object is what the next analysis needs); it is not a mid-plan step.
`AIAction(replaces_object = TRUE)` marks it, and `ValidateAIPlan()` refuses a
plan in which such a step shares the plan with any other step
(`object_replacement_not_standalone`), because the round's record is written on
the replacement object and the other steps' records would be gone.

### 4.6 Evidence

One node per action, `kind = "deterministic_summary"`, `source` = the record id,
ids following the `ev:` convention (`R/ai-execution.R:762,811`); a re-run with
a different name gets a different node, which is correct.

| action | evidence id | values (numeric/logical only) | claim_level |
|---|---|---|---|
| geneset | `ev:geneset_<name>` | `n_sets`, `n_cells`, `n_genes_used` | measured |
| supercell | `ev:supercell_<name>` | `n_metacells`, `n_cells`, `gamma`, `k_knn` | measured |
| priority | `ev:augur_<name>` | `n_cell_types`, `auc_min`, `auc_median`, `auc_max` | associated |
| phenotype | `ev:phenotype_association` | `n_columns_added`, `n_cells`, `external_inputs = TRUE` | observed |
| scenic | `ev:scenic_<name>` | `n_regulons`, `n_cells`, `python_backend = TRUE` | measured |
| perturbation | `ev:perturbation_<gene>` | `n_cells`, `perturbation_value`, `simulated = TRUE` | associated |

`n_genes_used` is computed in the handler from the resolved `gene_sets`
(the record stores only set names). `auc_*` are computed from
`record$artifacts$AUC` (columns `cell_type`, `auc`). `n_prioritized` is dropped:
no such field exists and a threshold would be an invented definition.

Gamma names, paths and column names stay in the record's `inputs`/`params`; no
free text enters an evidence node.

### 4.7 Skills

Extend the existing "Modelling choices that belong to the user" section
(`inst/skills/sclet-analysis/SKILL.md:171-196`) with the six inputs, state that
five of the actions are driver-only, and say what a plan should do instead (ask
or omit). The whole skill is injected, so no selection change is needed; pin via
the existing `system_prompt` capture in `test-ai-skills.R`.

## 5. Non-goals

- No change to `RunAIAnalysis()`'s default `include`.
- No clarification/collect-input mechanism (this spec hides what it cannot
  collect).
- No Python backend execution in tests; `sclet_basilisk_run` and the exported
  `Run*` functions are mocked for the four non-local actions.
- No `...` passthrough into the actions; the schema lists explicit parameters.
- No claim that a `returns = "value"` action can contribute record edges: it
  cannot, and (verified 2026-10-23) none of the ten value-returning actions
  writes a state or analysis record, so nothing is lost by that rule.
- No claim that a simulated perturbation is measured (`associated`, with
  `simulated = TRUE`).

## 6. Acceptance criteria

- **A1** `include = "all"` yields 40 actions; each new action belongs to its own
  group (single-group actions only — `confirm_design`'s triple membership is
  pre-existing).
- **A2** Every one of the 40 is discoverable through
  `sclet_ai_action_groups()`.
- **A3** Each action refuses a fixture that lacks its input, naming the gate;
  no handler runs.
- **A4** With the backend mocked (or real where local), each action writes its
  record, registers exactly one evidence node with `source` = that record id,
  all-numeric/logical values, and the claim level from §4.6.
- **A5** A driver-mode plan per action validates, executes, and appears in
  `AnalysisGraph()` with a `produced` edge to its record.
- **A6** Every user-owned input is required: omitting it refuses with the
  parameter/gate named.
- **A7** The planner catalog (`sclet_ai_describe_registry()` over `include =
  "all"`) contains `run_perturbation_priority` and none of the five
  driver-only actions; `available_actions` sanitizes with no redaction.
- **A8** The skill text reaches the prompt, pinned.
- **A9** No existing test changes expectation.

## 7. Tests

New `tests/testthat/test-ai-registry-mainlines.R`:

| case | asserts |
|---|---|
| registry shape | A1: ids, groups, cost, `planner_visible`, `allowed_state_types` |
| group invariant | A2 across all 40 |
| planner catalog | A7: five hidden, one visible; zero redactions |
| gates | A3, one case per action, fixture lacking exactly that input |
| handlers + evidence | A4: geneset and supercell real; the other four mocked |
| plan round trip | A5: one mocked plan per non-local action |
| user-owned inputs | A6: missing `gene_sets` / `target_gene` / `condition_col` / matrices named in the refusal |
| default include | A9 helper: `RunAIAnalysis`'s default still five groups |

`run_supercell`'s real test uses a normalized 60-cell fixture (`logcounts`
present) and asserts the metacell count and the evidence values; it does not run
after another plan (the object is replaced).

## 8. Risks

- **Mocking exported functions** (`with_mocked_bindings(RunSCENIC = ...)`,
  `.package = "sclet"`) hides argument-order mistakes in the real call. Mitigated
  by a signature assertion in the handler test namespace: the mock records
  `names(list(...))` and the test compares it with `names(formals(RunSCENIC))`.
- **scPAS absent** means the phenotype real path is untested; the comment in the
  test says so rather than implying coverage.
- **Evidence conflict on re-run**: ids include the record `name`/gene, and the
  recorder refuses silently changed values; `idempotent = FALSE` makes retries
  explicit.
- **Supercell destroys state**; documented in description + skill, hidden from
  the planner.

## 10. Revision note (rev 2)

### Rev 4 — final coverage and trace audit (2026-10-23)

- Identity and reduction gates honour the plan's own earlier steps
  (`planned$active_ident`, `planned$reductions`/`active_reduction`), so
  `find_clusters → run_perturbation_priority` and `run_pca → … → run_in_silico_perturbation`
  validate in one round instead of demanding state the object does not have yet.
- The four state-record writers that omitted `created_at` were reduced to one
  real omission (`R/malignancy-cnv.R`'s failed-run cnv record; `celloracle`
  fixed here, `cellrank`/`regvelo` already stamped at runtime). Without a stamp,
  the trace diff cannot see a same-id rewrite and the graph keeps attributing
  the changed record to the old producer.
- The coverage audit found **eight** mainlines still without actions beyond the
  six this spec closed: `RunImputation`, `RunRegVelo`,
  `ClassifyMalignantEpithelium`, `RunKNNPredict` / `RunSymphonyMapping`
  (`RunReferenceMapping`'s non-SingleR branches), `RunSpatialDeconvolution`,
  `RunSpatialColocalization`, `RunSpatialNiche`, and the multimodal registration
  shell (`register_multimodal_adt/_atac`). `RunDiffusionMap` and
  `AggregateExpression` write no state record and are outside the AI action
  contract. These are the next milestone, not part of this one.

### Rev 3 — post-implementation review

- `replaces_object` + the standalone-plan guard (§4.5).
- The identity sentinel: `ActiveIdent()` returns `"colLabels"`, which is not a
  colData column on the objects this package produces (the labels live in
  `label`). The action resolves the real column and passes `cluster_key`
  explicitly; the gates mirror the callees (an explicit key must exist, an
  absent one falls back to `ActiveIdent()`).
- The goal→action map carries only planner-visible actions: an entry for a
  hidden action made the repair message demand a step the planner cannot propose,
  and paid for a repair round.
- `run_perturbation_priority` refuses a record with no usable AUC table
  (`priority_no_result`) instead of registering an empty evidence node.
- Gates validate named `assay_use` / `cluster_key` / `group_col`.
- Evidence `source` resolution prefers the analysis store: a state record and an
  analysis record can share a bare id (`supercell`), which made the edge
  ambiguous and dropped it.
- Tests gained: group exclusivity, positive gate cases, the replacement refusal,
  AUC median, and the evidence→analysis `source_of` edge.

## 9. Revision note (rev 2)

Adopted from the reviews: `planner_visible` (§4.2) instead of "driver-mode by
convention"; path parameters cannot reach the planner (privacy deny list); five
group tables, not three; `RunSuperCell` needs `logcounts` and replaces the
object; `gamma` is a ratio, so the cell gate is `2 * gamma`; evidence ids are
`ev:` and `supercell`/`scenic` were already taken; `n_genes_used` computed in the
handler and `n_prioritized` dropped; phenotype counts are `observed`;
`planner_visible`/catalog pins; `cluster_key` default can be `NULL`; perturbation
needs a counts assay; phenotype needs Seurat.
