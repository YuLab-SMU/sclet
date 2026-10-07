# P5 Spec (Slice 1): Velocity Readiness + Action Adapter + Evidence

> Status: draft for implementation (2026-10-07)
> Parent roadmap: `.dev/ai-product-roadmap.md`, P5 ("受控扩展动态、空间和多模态"), first item in
> the recommended order: "velocity readiness -> velocity execution -> trajectory / velocity
> interpretation -> CellRank / fate -> spatial -> multimodal".
> Precondition satisfied: P0 (`92332ca`), P0.1 (`607e48d`), P1 (`9bf701d`), P2 (`1688260`), P3
> (`ce1f3b4`), P4 (`43f9f2f`) are all committed on `devel`.
> Baseline before this change: `devtools::test(filter = "ai-")` ->
> `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1202 ]`. `make check` -> `0 errors | 0 warnings | 0 notes`.

## 1. Problem

The roadmap's P5 explicitly forbids starting with "just a `run_*` action": every new domain
must ship readiness diagnostics, an action adapter, prerequisites, success criteria, an
evidence summary, interpretation-report wiring, and failure/uncertainty tests together. This
spec covers the first slice only -- **velocity readiness and a minimal, evidence-producing
`run_velocity` action** -- not the full velocity interpretation report (deferred to the next
slice, see section 2) and not CellRank/fate/spatial/multimodal (deferred per the roadmap's own
order).

**Verified while writing this spec** (every claim below was checked by reading the actual code,
not assumed):

1. **`RunVelocity()` already exists as a real, deterministic R workflow function**
   (`R/velocity.R:19-133`), basilisk/`velociraptor`-backed. It requires `"spliced"` and
   `"unspliced"` assays (hard `stop()` at line 35-37) and a dimensionality reduction (`stop()`
   at line 31-33, defaulting to `DefaultReduction(sce)`). On success it writes
   `velocity_<name>`-prefixed reductions and up to four `colData` columns
   (`velocity_pseudotime`, `velocity_confidence`, `root_cells`, `end_points`, intersected with
   what `velociraptor` actually returned), and records a `velocity` state record via
   `sclet_set_state_record()` (lines 86-114) with `method`/`inputs`/`artifacts`/`params`/
   `summary`/`created_at` fields -- the same state-record shape every other domain already
   uses. **Confirmed `velociraptor` is installed in this environment**
   (`requireNamespace("velociraptor", quietly = TRUE)` returns `TRUE`), so this spec's tests can
   exercise the real handler, not a mock.

2. **No `check_velocity_readiness()` exists.** Confirmed with
   `grep -n 'check_.*_readiness' R/ai-diagnostics.R`: only `check_rare_cell_readiness`,
   `check_trajectory_readiness`, `check_annotation_readiness`, `check_integration_readiness`
   exist. This spec adds a fifth, modeled on `check_annotation_readiness()`
   (`R/ai-diagnostics.R:834-906`) rather than `check_trajectory_readiness()`, because velocity
   has no analogous human-judgment "root/semantic" parameter to gate on: `RunVelocity()`'s
   `mode` argument (`"deterministic"`/`"stochastic"`/`"dynamical"`) is a numerical-method
   choice, not a biological assumption the user must declare the way `run_trajectory`'s
   `start_cluster` or `run_integration`'s `batch` column is. **Confirmed by reading
   `RunCellRank()` and `RunRegVelo()` signatures too**: neither has a design-semantic parameter
   either (`RunCellRank` takes `reduction`/`cluster_key`/`velocity_id`/`backend`;
   `RunRegVelo` takes `grn`/`regulators`/assay names/`reduction`/hyperparameters) -- so this
   spec does not add a `ConfirmAIDesignSemantics()` gate for velocity, matching what the real
   handler actually needs, not inventing a confirmation step the domain does not require.

3. **No `run_velocity` (or any velocity-related) action exists in the registry.** Confirmed
   with `grep -n '"velocity"' R/ai-execution.R`: the string does not appear at all.
   `AIDefaultExecutionRegistry()`'s `allowed_groups` (`R/ai-execution.R:179`) is
   `c("read", "preprocess", "dimred", "graph", "cluster", "integration", "annotation",
   "rare_cell", "trajectory", "all")` -- `"velocity"` is not among them and must be added.

4. **No evidence-recording pattern exists for velocity.** The nearest analogs are
   `sclet_ai_record_trajectory_evidence()` (`R/ai-execution.R:1038-1093`, confirmed in the P4
   spec to record `claim_level = "consistent_with"` with `pseudotime_is_absolute_time = FALSE`)
   and `run_rare_cell_detection`'s evidence recording. Velocity evidence must follow the same
   discipline: bounded aggregate values only, an honest `claim_level`, and an explicit
   limitation note (velocity direction/pseudotime is a model-fitted estimate, not a measured
   biological fact -- this is the same "do not let a derived quantity be read as ground truth"
   principle `run_trajectory`'s evidence note already states for pseudotime).

5. **This spec's scope is a vertical slice through all seven of the roadmap's required
   deliverables, scaled down to what one new domain actually needs**:
   - readiness diagnostics: `check_velocity_readiness()` (new);
   - action adapter: `run_velocity` registered under a new `"velocity"` include-group (new);
   - prerequisites: reuse `RunVelocity()`'s own two hard requirements
     (spliced/unspliced assays, a reduction) inside the action's `prerequisites` function, so a
     bad plan is rejected by `ValidateAIPlan()` before `RunVelocity()` ever raises its own
     `stop()`;
   - design confirmation: **deliberately not needed** for this domain (justified in item 2
     above) -- state this explicitly rather than silently omitting it;
   - success criteria: not a new mechanism -- a plan author can already attach
     `success_criteria` (P0.1) referencing `run_velocity`'s own step id/output, so this spec
     only needs to confirm that path works for the new action, not build a new one;
   - evidence summary: a new `sclet_ai_record_velocity_evidence()`, modeled on the trajectory
     equivalent;
   - interpretation report: covered implicitly, since P3's `RunAIAnalysis(interpret = TRUE)`
     and `AIExplainAnalysis()` are already domain-agnostic (confirmed in P4: `interpret = TRUE`
     already worked for `run_de_test` with zero velocity-specific code) -- this spec adds one
     test proving the same holds for `run_velocity`, not new interpretation code;
   - failure/uncertainty tests: readiness-not-ready and prerequisite-failure scenarios routed
     through the full `RunAIAnalysis()` orchestrator, matching P4's own pattern.

## 2. Non-goals (do not implement these here)

- No `RunCellRank`/`RunRegVelo` action adapters, no CellRank/fate action, no spatial action
  catalog, no multimodal action catalog. The roadmap's own recommended order places all of
  these strictly after velocity readiness and execution; do not get ahead of it.
- No trajectory/velocity *combined* interpretation (the roadmap's second-listed item,
  "trajectory / velocity interpretation", is a distinct, later slice -- e.g. correlating
  velocity-derived pseudotime with `run_trajectory`'s slingshot pseudotime, or a joint report).
  This spec's "interpretation" scope is limited to confirming the existing, domain-agnostic
  `interpret = TRUE` path works for `run_velocity`, nothing more.
- No `ConfirmAIDesignSemantics()` gate for velocity (justified in section 1 item 2). Do not add
  one "for consistency" with trajectory/integration -- that would invent a confirmation
  requirement the real `RunVelocity()` handler does not have.
- No change to `R/velocity.R`, `R/cellrank.R`, or `R/regvelo.R` themselves. If this spec's
  tests reveal a genuine defect in `RunVelocity()`'s own behavior, fix it narrowly and document
  it with the same precision as every prior round's bug fixes -- but do not refactor the
  deterministic workflow layer as part of adding its AI-facing wrapper.
- No change to `R/ai-planning.R`'s success-criteria mechanism, `R/ai-privacy.R`'s allowlist
  logic itself (only its allowlist *regex contents* may need extending if this spec introduces
  a new field name -- check against the existing allowlist before assuming an extension is
  needed), or any other cross-cutting AI subsystem file beyond what section 3 names.
- Do not attempt multiple velocity `mode` values' numerical correctness (that is
  `velociraptor`'s own concern, already covered by its own upstream tests). This spec's tests
  check the AI-facing contract (readiness, prerequisites, evidence shape, orchestration), not
  velocity's scientific correctness.
- Do not make the new readiness check or action handle `RunCellRank`'s `velocity_id` lookup or
  any CellRank-specific state. That is explicitly out of scope for this slice.

## 3. What to build

### 3.1 `check_velocity_readiness(object, reduction = NULL)` (new, in `R/ai-diagnostics.R`)

Model directly on `check_annotation_readiness()` (`R/ai-diagnostics.R:834-906`), reading it in
full before writing this function. Required checks, each mapped to `RunVelocity()`'s own two
hard requirements (section 1 item 1):

```r
check_velocity_readiness <- function(object, reduction = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    assay_names <- SummarizedExperiment::assayNames(object)
    has_spliced <- "spliced" %in% assay_names
    has_unspliced <- "unspliced" %in% assay_names
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    resolved <- sclet_ai_diag_embedding_name(object, reduction)
    has_reduction <- !is.null(resolved)
    questions <- character()
    if (!has_spliced || !has_unspliced) {
        questions <- c(questions, paste0(
            "RNA velocity requires 'spliced' and 'unspliced' assays (e.g. from a splicing-aware",
            " quantification pipeline such as velocyto or alevin-fry); the current assays are: ",
            paste(assay_names, collapse = ", ")
        ))
    }
    if (!has_reduction) {
        questions <- c(questions, paste0(
            "No usable embedding is available; run the dimensionality reduction first. ",
            "Available reductions: ", paste(reductions, collapse = ", ")
        ))
    }
    status <- if (!has_spliced || !has_unspliced || !has_reduction) "not_ready" else "ready_for_diagnostic"
    list(
        status = status,
        checks = list(
            spliced_assay_available = has_spliced,
            unspliced_assay_available = has_unspliced,
            reduction_resolved = resolved,
            available_reductions = reductions,
            raw_values_included = FALSE
        ),
        blocked_actions = if (identical(status, "not_ready")) "velocity" else character(),
        questions = questions,
        notes = paste(
            "No velocity mode (deterministic/stochastic/dynamical) is recommended here: that is a",
            "numerical-method choice for run_velocity's own params, not a readiness gate. No design",
            "confirmation is required for this domain -- velocity has no sample/batch/root semantic",
            "parameter the way trajectory or integration do."
        )
    )
}
```

Use `sclet_ai_diag_embedding_name()` (already used by `check_trajectory_readiness()`, confirmed
at `R/ai-diagnostics.R:684`) rather than inventing a different reduction-resolution helper.

### 3.2 `run_velocity` action (new, in `R/ai-execution.R`)

Add a new `"velocity"` include-group (extend `allowed_groups` at `R/ai-execution.R:179` from
`c("read", ..., "trajectory", "all")` to also include `"velocity"`), and register the action
inside an `if ("velocity" %in% include) { ... }` block, following the exact structural pattern
of the `"trajectory"` block (`R/ai-execution.R:931-1030`):

```r
run_velocity = AIAction(
    name = "run_velocity",
    description = "Estimate RNA velocity with velociraptor::scvelo and record a bounded aggregate evidence node.",
    handler = function(object, params) {
        name <- params$name %||% "velocity"
        result <- RunVelocity(
            object,
            mode = params$mode %||% "deterministic",
            use.dimred = params$use.dimred %||% NULL,
            subset.row = if (is.null(params$subset.row)) NULL else sclet_ai_vector_param(params$subset.row, "character"),
            name = name
        )
        sclet_ai_record_velocity_evidence(result, name = name)
    },
    input_schema = list(
        mode = "character",
        use.dimred = "character",
        subset.row = "character_vector",
        name = "character"
    ),
    prerequisites = function(object, params, planned = NULL) {
        available <- planned$assays %||% SummarizedExperiment::assayNames(object)
        if (!all(c("spliced", "unspliced") %in% available)) {
            return(paste0(
                "spliced_unspliced_missing: run_velocity requires 'spliced' and 'unspliced' assays; available: ",
                paste(available, collapse = ", ")
            ))
        }
        reduction <- params$use.dimred %||% planned$active_reduction %||%
            tryCatch(DefaultReduction(object), error = function(e) NULL)
        available_reductions <- planned$reductions %||% tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
        if (is.null(reduction) || !nzchar(reduction) || !reduction %in% available_reductions) {
            return(paste0(
                "reduction_missing: a dimensionality reduction is required for run_velocity; available: ",
                paste(available_reductions, collapse = ", ")
            ))
        }
        if (!base::requireNamespace("velociraptor", quietly = TRUE)) {
            return("optional_package_missing: velociraptor is required for run_velocity; install via BiocManager::install('velociraptor')")
        }
        TRUE
    },
    returns = "sce",
    requires_confirmation = TRUE,
    output_schema = list(
        required_states = list(velocity = "any"),
        note = paste(
            "Velocity direction and velocity_pseudotime are model-fitted estimates from scVelo,",
            "not measured biological facts; they depend on the chosen mode and the spliced/unspliced",
            "ratio and must not be read as ground-truth cell-state transitions. No cell is removed",
            "or filtered by this action."
        )
    ),
    mutates_object = TRUE,
    allowed_state_types = c("velocity", "ai_evidence"),
    estimated_cost = "high",
    idempotent = FALSE
)
```

Use `estimated_cost = "high"` (not `"medium"` like trajectory/rare-cell) because
`velociraptor::scvelo()` runs a basilisk-backed Python model fit, which is substantially more
expensive than the other domains' handlers -- state this honestly rather than copying
trajectory's cost label by default. Verify `planned$assays`/`planned$active_reduction`/
`planned$reductions` are the correct field names by reading how `run_pca`'s and
`run_trajectory`'s own `prerequisites` functions already read `planned` (confirmed at
`R/ai-execution.R:343-345` for the assay pattern); do not invent different field names.

### 3.3 `sclet_ai_record_velocity_evidence()` (new, in `R/ai-execution.R`)

Model directly on `sclet_ai_record_trajectory_evidence()` (`R/ai-execution.R:1038-1093`), read
it in full first. Read back the just-written state record via
`sclet_get_state_record(object, "velocity", name)` (matching the trajectory recorder's own
pattern of reading back its own state record rather than trusting handler-local variables), and
build a bounded aggregate summary from whichever of the four possible `colData` columns
(`velocity_pseudotime`, `velocity_confidence`, `root_cells`, `end_points`) are actually present
-- `RunVelocity()`'s own code intersects against what `velociraptor` returned
(`R/velocity.R:70-73`), so not all four are guaranteed. **Resolved while writing this spec by
actually running `RunVelocity()` on a real synthetic fixture (`mode = "deterministic"`,
40 cells, 30 genes, ~13.5s wall time)**: all four columns were present and all four are
`numeric` (confirmed via `class()`), including `root_cells` and `end_points` -- both are
continuous values from scVelo's Markov diffusion terminal-state estimation (observed sample
values like `0, 0, 0, 0, 0` for `root_cells` and `0.34, 0.63, 0.49, 0.17, 0.43` for
`end_points` in that run), **not** boolean flags or counts. Do not implement the
"bounded presence/count summary" this spec's own earlier draft speculated for them. Instead,
summarize **all four** present columns identically, with the same
`n_observed`/`mean`/`quantile_25`/`median`/`quantile_75`/`max` shape the trajectory recorder
already uses for pseudotime columns (`R/ai-execution.R:1056-1064`) -- reuse that exact
summarizing pattern uniformly across `velocity_pseudotime`, `velocity_confidence`, `root_cells`,
and `end_points`, whichever subset is actually present, rather than special-casing any of them.
Evidence shape:

```r
evidence <- list(
    id = paste0("ev:velocity_", name),
    kind = "deterministic_summary",
    values = values,  # the bounded aggregate list described above, plus:
                       # mode, use.dimred, reductions_written, fields_written,
                       # velocity_is_model_estimate = TRUE, raw_values_included = FALSE
    claim_level = "consistent_with"
)
tryCatch(
    RecordAIEvidence(object, evidence, source = name, parents = character(), scope = NULL),
    error = function(e) object
)
```

Use `claim_level = "consistent_with"` (matching trajectory's own evidence claim level,
confirmed at `R/ai-execution.R:1087`) rather than a stronger level: velocity is itself a
model-fitted estimate one step further removed from direct observation than clustering or even
pseudotime ordering, so it should not claim more confidence than trajectory's own pseudotime
evidence does.

### 3.4 Privacy: one field rename, no allowlist regex change needed

**This section was substantially corrected while writing this spec by actually running the
privacy scanner against realistic nested payloads, not by guessing from field-kind lookups in
isolation.** An initial draft of this spec assumed every new field name needed an allowlist
regex extension (the same shape as P1/P2's fixes). Direct testing against the real nesting
depth this evidence actually lives at (`state_records.ai_evidence.<id>.summary.values.*`,
confirmed by reading `R/ai-context.R:405-408` and `R/ai-privacy.R:112-128`) showed this
assumption was wrong for almost every field:

- `R/ai-privacy.R:122-128` has an escape hatch: once nested at `depth >= 2` under an
  already-`"aggregate"`-or-`"container"`-classified parent, an `"unknown"`-kind field is
  automatically promoted to `"aggregate"` as long as its value is numeric/logical/`NULL`/a
  list/or a short (`<= 12`-element) character vector. `state_records` and `ai_evidence` are
  both already allowlisted as `"aggregate"` (confirmed interactively). This means
  `velocity_is_model_estimate` (logical), `use.dimred` (short character), `mode`,
  `reductions_written`, and the per-column quantile summaries (`n_observed`/`mean`/
  `quantile_25`/`median`/`quantile_75`/`max`, matching the trajectory recorder's own field
  names) **all pass through with zero redaction and zero warnings already**, confirmed by
  constructing the exact realistic nested shape this spec's evidence will produce and running
  it through `sclet_ai_payload_scan()` directly -- not by reasoning about field-kind lookups on
  bare names in isolation, which is misleading at depth 0 (bare-name lookups report `"unknown"`
  for all of these, but that is irrelevant once the depth-escape-hatch applies).
- The one genuine exception: `coldata_fields_written` is hard-redacted (as `"metadata"`, not
  merely `"unknown"`) **regardless of nesting depth**, because the `"metadata"` classification
  in `sclet_ai_payload_field_kind()` is checked before the depth-escape-hatch logic runs, and
  it matches on the substring `coldata` anywhere in the field name. Confirmed directly: both
  `coldata_fields_written` and `coldata_columns_written` resolve to `"metadata"` even at depth
  3; renaming to **`fields_written`** (no `coldata` substring) resolves to plain `"unknown"`,
  which the depth-escape-hatch then correctly promotes to `"aggregate"`.
- **Fix: name the field `fields_written` in section 3.3's evidence `values`, not
  `coldata_fields_written`.** This is a one-word rename in the spec itself (already reflected
  in section 3.3 above); no change to `R/ai-privacy.R` is needed for this slice.
- The one case that does require caution, confirmed by direct testing: a bare top-level field
  literally named `velocity` (e.g. `context$velocity <- ...`) is **not** safe -- it sits at
  depth 0/1 with no aggregate parent yet, so the depth-escape-hatch does not apply, and it
  hard-errors as `"Unknown outbound payload field: velocity"`. **This spec does not add any
  top-level `velocity` context field** (readiness and evidence both live nested under existing
  allowlisted containers -- `check_velocity_readiness()`'s own result is not attached to
  `GetAnalysisLedger()`'s context at all in this slice, matching how
  `check_rare_cell_readiness()`/`check_integration_readiness()` are also never attached
  verbatim, only consulted internally or passed into `AIInvestigate()`'s own `diagnostics` sub-list
  which is itself already nested under an allowed path). If a future round ever adds a
  top-level `velocity` field to the ledger, it would need its own allowlist entry at that point
  -- not a concern for this slice, but worth naming explicitly so a future implementer does not
  have to rediscover it.

## 4. Required test scenarios (`tests/testthat/test-ai-velocity.R`)

Write one new test file. Build a real spliced/unspliced fixture (synthetic counts are fine, as
long as both assays exist with matching dimensions) and run it through `NormalizeData()` +
`FindVariableFeatures()` + `RunPCA()` first, matching the construction style of the existing
domain fixtures in `tests/testthat/test-ai-execution.R` and
`tests/testthat/test-ai-domain-adapter-e2e.R`. `skip_if_not_installed("velociraptor")` on every
scenario that actually executes `run_velocity` (confirmed installed in this environment, so
these will genuinely run here, but the guard keeps the suite portable per the P4 round's own
convention).

1. **`check_velocity_readiness()` on an object missing spliced/unspliced**: `status ==
   "not_ready"`, `checks$spliced_assay_available == FALSE` and/or
   `checks$unspliced_assay_available == FALSE`, `"velocity" %in% blocked_actions`.
2. **`check_velocity_readiness()` on an object with spliced/unspliced but no reduction**:
   `status == "not_ready"`, `checks$reduction_resolved` is `NULL`.
3. **`check_velocity_readiness()` on a fully-ready object**: `status == "ready_for_diagnostic"`,
   `length(questions) == 0`.
4. **`run_velocity` prerequisites reject a plan missing spliced/unspliced** through
   `ValidateAIPlan()` directly (matching the existing per-action prerequisite-test pattern in
   `tests/testthat/test-ai-execution.R`, not through the full orchestrator for this one
   low-level check): error contains `"spliced_unspliced_missing"`.
5. **`run_velocity` prerequisites reject a plan with no available reduction**: error contains
   `"reduction_missing"`.
6. **`run_velocity` survives `ExecuteAIPlan()` directly on a real, ready fixture** (matching
   `test-ai-execution.R`'s own direct-execution-layer coverage style for every other domain):
   `status == "completed"`, the resulting evidence (`evidence[["ev:velocity_<name>"]]`) is
   non-null, `claim_level == "consistent_with"`, and no `values` entry is an atomic vector of
   length `>= ncol(object)` (the same bounded-value check pattern reused in every prior round).
7. **`run_velocity` end-to-end through `RunAIAnalysis()`** (matching P4's own
   `sclet.ai.call`-mock pattern from `tests/testthat/test-ai-beginner.R:71`, reused exactly as
   `tests/testthat/test-ai-domain-adapter-e2e.R` already does for the other four domains):
   `result$status == "completed"`, evidence present with the correct `claim_level`.
8. **A velocity readiness failure routed through the full orchestrator**: mock a plan whose
   `run_velocity` step targets an object with no spliced/unspliced assay. Assert
   `result$status == "invalid_plan"` (not a raised R error) via `RunAIAnalysis()`, matching
   P4's own prerequisite-failure-through-orchestrator pattern (scenario 5 in
   `.dev/spec-p4-domain-adapter-e2e.md`).
9. **`RunAIAnalysis(interpret = TRUE)` works for `run_velocity`**: reuse the single
   `sclet.ai.call` mock dispatching on `task` (`"analysis_plan"` vs `"analysis_explanation"`),
   exactly as P4's own scenario 6 does. Assert `result$report$interpretation` is a valid
   `sclet_ai_result`, `result$report$interpretation_error` is `NULL`, and the interpretation was
   actually recorded on the returned object's ledger (checking for an `ai_interpretation_*`
   analysis id) -- do not repeat the P3-round mistake of only checking that `interpretation` is
   non-null without confirming it was actually persisted.
10. **Evidence never contains a raw per-cell vector**: construct a fixture, run `run_velocity`
    directly, and assert every entry in the recorded evidence's `values` is either a scalar or a
    bounded-length vector (length `< ncol(object)`), confirming
    `sclet_ai_record_velocity_evidence()` only ever writes aggregate summaries.
11. **No design confirmation is required**: construct a ready fixture with no
    `ConfirmAIDesignSemantics()` call at all, and confirm `run_velocity` still validates and
    executes successfully through `ValidateAIPlan()`/`ExecuteAIPlan()` -- a regression guard
    proving this spec did not accidentally add a design-confirmation requirement velocity does
    not need (per section 1 item 2 / section 2's explicit non-goal).
12. **Full suite regression guard**: after this change, every existing test still passes
    unmodified.
13. **Privacy regression guard (verifies section 3.4's finding, not a new allowlist fix)**:
    record real `run_velocity` evidence on an object, then call `AIStatus()` or
    `AIPlanAnalysis()` on that object (a real `sclet_ai_call()`-reachable path, matching the
    exact mistake pattern from the P1 round -- the implementer's own test file back then only
    called the state accessor directly, never the real call path, which is why 4 warnings went
    undetected until a manual reproduction caught them) and assert **0 warnings** fire. This
    confirms section 3.4's finding holds for real recorded evidence, not just the synthetic
    nested payload this spec used to verify it. If this test unexpectedly fails with a
    redaction warning for any field other than the ones already known-and-handled in section
    3.4 (i.e. if the depth-escape-hatch does not behave the way this spec verified it would for
    the *real* evidence shape `RunVelocity()`/`sclet_ai_record_velocity_evidence()` actually
    produce, as opposed to the hand-constructed shape this spec tested), extend the allowlist
    regex at `R/ai-privacy.R:347` for that specific field and document exactly which field and
    why -- do not assume this will be needed, but do not skip verifying it either.

## 5. Documentation updates

- `man/AIDefaultExecutionRegistry.Rd`, `man/check_velocity_readiness.Rd` (new): hand-edit/create
  (do **not** run `devtools::document()` or `roxygen2::roxygenise()` -- this caused an
  unauthorized `DESCRIPTION` change in the P0.1 round and has been avoided by hand-editing in
  every round since). Document the new `"velocity"` include-group and the new readiness
  function's return shape.
- `NEWS.md`: one new top entry describing exactly what was added: `check_velocity_readiness()`,
  the `run_velocity` action and its `"velocity"` include-group, `sclet_ai_record_velocity_evidence()`,
  the `fields_written` naming choice (and why `coldata_fields_written` was avoided -- it trips
  the privacy scanner's metadata classification regardless of nesting depth), whether the
  fallback privacy allowlist extension from section 3.4 was actually needed or not (state this
  explicitly either way, with evidence), and the explicit statement that no
  `ConfirmAIDesignSemantics()` gate was added for this domain and why. State plainly what was
  *not* done: no CellRank/fate/spatial/multimodal, no combined trajectory/velocity
  interpretation, no change to `RunVelocity()`/`RunCellRank()`/`RunRegVelo()` themselves unless
  a genuine defect was found and fixed (document precisely if so).
- `.dev/ai-advanced-analysis-spec.md`: update the top status line to note velocity readiness +
  action + evidence are implemented, following the exact precedent of how P4's entry was
  phrased (state what shipped, do not overclaim P5 as a whole is done -- CellRank/fate/spatial/
  multimodal and the trajectory/velocity joint interpretation remain explicitly open).
- Do not edit `.dev/ai-product-roadmap.md`.

## 6. Acceptance checklist (self-verify before reporting done)

- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, **0 warnings**, and a pass count
      **greater than or equal to** 1202 (the baseline measured before this change).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`. If a transient, unrelated
      environment error occurs (this happened once in the P2 round with a missing-package
      error unrelated to any touched file, confirmed transient by a clean re-run), re-run once
      before reporting a result.
- [ ] `git diff --check` is clean.
- [ ] `LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-diagnostics.R R/ai-execution.R
      tests/testthat/test-ai-velocity.R` finds nothing (also check `R/ai-privacy.R` if scenario
      13 ends up requiring the fallback allowlist extension described in section 3.4).
- [ ] `git status --short` shows changes only to: `R/ai-diagnostics.R`, `R/ai-execution.R`,
      `man/AIDefaultExecutionRegistry.Rd`, a new `man/check_velocity_readiness.Rd` (or
      equivalent), `NEWS.md`, `.dev/ai-advanced-analysis-spec.md`, and the new test file.
      `R/ai-privacy.R` must NOT appear unless scenario 13 actually failed and required the
      fallback fix section 3.4 describes -- if it appears, the report must say precisely why.
      No `DESCRIPTION` side effect. `R/velocity.R`/`R/cellrank.R`/`R/regvelo.R` must NOT appear
      unless a genuine, documented defect was found and fixed in one of them.
- [ ] Scenario 13 (section 4) is run and shown to pass with **0 warnings** through a real
      `sclet_ai_call()`-reachable path (not a direct call to the state accessor, the exact
      P1-round mistake this scenario exists to prevent) -- confirming section 3.4's verified
      finding (no allowlist regex change needed, only the `fields_written` rename) holds for
      the real evidence shape, not just the spec's own synthetic test payload.
- [ ] Scenario 11 (section 4) confirms no `ConfirmAIDesignSemantics()` requirement was
      accidentally introduced for velocity.
- [ ] `grep -c 'AIAction(' R/ai-execution.R` is exactly 18 (the pre-change baseline of 17, plus
      exactly one new `run_velocity` action -- not two, not zero).
- [ ] If (and only if) scenario 13 required the fallback allowlist extension described in
      section 3.4, re-run the dangerous-field-name safety check (patient/donor/token/secret/etc.
      combinations) to confirm they all still resolve to `deny` after the extension.
- [ ] No file outside the authorized list in this section's `git status` item was touched.

## 7. Implementation order

1. Read this entire spec once before writing any code.
2. Read `check_annotation_readiness()` (`R/ai-diagnostics.R:834-906`) and
   `sclet_ai_diag_embedding_name()` in full before writing `check_velocity_readiness()`.
3. Read the `"trajectory"` action block (`R/ai-execution.R:931-1030`) and
   `sclet_ai_record_trajectory_evidence()` (`R/ai-execution.R:1038-1093`) in full before writing
   `run_velocity` and `sclet_ai_record_velocity_evidence()`.
4. Read `RunVelocity()` in full (`R/velocity.R:19-133`) to confirm every parameter name and the
   exact `colData`/reduction fields it can produce before writing the handler and evidence
   recorder.
5. Implement `check_velocity_readiness()` (section 3.1).
6. Extend `allowed_groups` and register `run_velocity` (section 3.2).
7. Implement `sclet_ai_record_velocity_evidence()` (section 3.3).
8. Use `fields_written` (not `coldata_fields_written`) in the evidence `values` per section
   3.3/3.4 -- no privacy allowlist change is expected to be needed for this slice, confirmed by
   this spec's own direct testing against the real nesting depth. Write scenario 13 (section 4)
   to confirm this holds for the real evidence shape; only extend `R/ai-privacy.R` if that
   scenario actually fails and surfaces a field this spec did not anticipate.
9. Write the test file (section 4), 13 scenarios.
10. Run `devtools::test(filter = "ai-", reporter = "progress")` in the background (allow extra
    time beyond the usual ~80-100s: `velociraptor::scvelo()` runs a real basilisk-backed Python
    model fit per scenario that executes it, which is slower than any prior round's handlers;
    a foreground call will certainly time out), iterate until green.
11. Run `make check` in the background (allow extra time for the same reason).
12. Update documentation per section 5.
13. Self-verify every item in section 6, listing a concrete yes/no and evidence for each.
14. Report: paste the actual terminal output of the test run and `make check`, list every file
    changed with `git status --short` and `git diff --stat`, and do not commit -- leave the
    changes in the worktree for review.


