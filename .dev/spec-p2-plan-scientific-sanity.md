# P2 Spec: Plan Stop Conditions, Data-Scale Sanity, and Human Confirmations

> Status: draft for implementation (2026-10-07)
> Parent roadmap: `.dev/ai-product-roadmap.md`, P2 ("补强 Plan / Scientific Sanity")
> Precondition satisfied: P0 (commit `92332ca`), P0.1 (commit `607e48d`), P1 (commit `9bf701d`)
> are all committed on `devel`.
> Baseline before this change: `devtools::test(filter = "ai-")` ->
> `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1114 ]`. `make check` -> `0 errors | 0 warnings | 0 notes`.

## 1. Problem

P0.1 (commit `607e48d`) already added `assumptions`/`candidate_routes`/`selected_route`/
`success_criteria`/`risks` to `sclet_ai_plan`, and `ValidateAIPlan()` already validates
criterion resolvability and route consistency. **Verified while writing this spec**:
`grep -rn "stop_condition" R/ tests/ .dev/` finds `stop_condition` only in
`.dev/ai-product-roadmap.md`'s prose -- nowhere in code or tests. The roadmap's P2 names three
concrete remaining gaps:

1. No `stop_conditions` field exists anywhere on the plan.
2. `ValidateAIPlan()` has no deterministic, generic "does this plan's scale make sense for this
   object" check. Per-action `prerequisites` functions do narrow, action-specific checks (e.g.
   `R/ai-execution.R:636-640` requires at least 2 cells per group for one DE-style action), but
   nothing checks a parameter like `ncomponents` against the object's actual `ncol()`/`nrow()`
   before execution. **Verified**: `run_pca`'s own `prerequisites` (`R/ai-execution.R:342-351`)
   only checks that the source assay/layer exists; it does not compare `params$ncomponents` to
   `ncol(object)` or `nrow(object)`. If a plan requests `ncomponents = 50` on a 10-cell object,
   `ValidateAIPlan()` currently returns `valid = TRUE`, and the failure only surfaces as an
   opaque error from inside `RunPCA()`/`prcomp()`/`irlba()` at actual execution time.
3. Clarification answers already become `user_decision` evidence (via `ResolveAIClarifications()`,
   `R/ai-functions.R:526`) and already appear in `analysis_story$user_decisions` (P1,
   `R/ai-context.R:404-422`), but nothing on the plan itself lists which human confirmations a
   *specific* plan depends on before it can run. The roadmap asks for "将 clarification 结果直接
   转为 plan 的 human_confirmations" -- a plan-level view of what the user already decided and
   what the plan still needs decided, not a new confirmation mechanism.

## 2. Non-goals (do not implement these here)

- No new action handlers. No velocity/spatial/multimodal/CellRank work.
- No per-action scale metadata field on `AIAction()`/action descriptors. **Verified**: the
  current `AIAction()` signature (`R/ai-execution.R:21-34`) has no `min_cells`/`min_features`/
  `required_modality` slot, and no action anywhere declares one. Do not invent one. The
  data-scale check added by this spec (section 3.2) is deliberately generic and self-contained:
  it reads `ncol(object)`/`nrow(object)` directly inside `ValidateAIPlan()` and checks one
  specific, already-named parameter (`run_pca`'s `ncomponents`) against them. It is not a
  general "every action declares its scale requirements" framework -- that would require
  touching every `AIAction()` call site in `R/ai-execution.R`, which is out of scope.
- No change to `prerequisites` functions on any existing `AIAction()`. The new check lives
  entirely inside `ValidateAIPlan()` in `R/ai-planning.R`.
- No change to the clarification catalog (`sclet_ai_clarification_catalog()`,
  `R/ai-functions.R`) or to `ResolveAIClarifications()`'s own logic. `human_confirmations` reads
  from a plan's own content (which design roles / which `check_*_readiness()` clarification
  paths it would trigger) and from already-recorded `user_decision` evidence; it does not add a
  new clarification type or a new confirmation mechanism.
- No change to `execution$status`'s existing value set (`"completed"`, `"failed"`,
  `"completed_with_errors"`, `"dry_run"`). A stop condition that is met during dry-run is
  reported as data on the dry-run result, not as a new status value.
- Do not touch `R/ai-comparison.R`, `R/ai-evidence.R`, `R/ai-privacy.R`, `R/ai-diagnostics.R`,
  `R/ai-design-confirmation.R`, `R/ai-integration-routes.R`, `R/ai-context.R` (P1's
  `analysis_story` is already shipped and out of scope for this change), or any `run_*`/action
  registration code. This spec is scoped to `R/ai-planning.R` and `R/ai-execution.R` only (the
  latter only for the `ExecuteAIPlan()` dry-run branch, not any action handler).
- Do not add any new privacy-sensitive field. Every new field introduced by this spec (section
  3) must consist only of: step ids already present in the plan, criterion/condition ids already
  defined by the plan author, parameter names/values already present in `plan$actions`, and
  booleans/counts. No raw colData values, no per-cell data, no secrets.

## 3. What to build

### 3.1 `stop_conditions` on the plan

Add a sixth optional field to `new_sclet_ai_plan()` (`R/ai-planning.R:152-?`, alongside the five
P0.1 fields), defaulting to `list()`, fully backward compatible:

```r
stop_conditions = list(
    list(
        id = "excessive_failure",
        description = "stop if more than one step fails",
        source = "action_output",    # or "evidence", same two sources as success_criteria
        step_id = "<id>",            # required when source == "action_output"
        evidence_id = "<id>",        # required when source == "evidence"
        check = "exists",            # one of: exists, equals, gte, lte, in
        field = NULL,                # required when check != "exists"
        value = NULL                 # required when check != "exists"
    )
)
```

This is intentionally the same shape as `success_criteria` (added in P0.1, see
`sclet_ai_normalize_success_criteria()` at `R/ai-planning.R:142-189` -- read it first and copy
its exact normalization pattern, do not invent a different shape for symmetry's sake). Implement
`sclet_ai_normalize_stop_conditions()` as a near-identical sibling function; the only semantic
difference is in how `ExecuteAIPlan()` later interprets a met condition (section 3.4), not in
how the field is shaped or normalized.

`ValidateAIPlan()` must validate `stop_conditions` with the same two checks `success_criteria`
already gets (reuse the identical logic, do not write a second divergent implementation):
- `stop_condition_unresolvable`: `source == "action_output"` referencing a `step_id` not present
  in `plan$actions`, or `source == "evidence"` with an empty/missing `evidence_id`.
- `stop_condition_incomplete`: `check != "exists"` with a missing `field` or `value`.

### 3.2 Generic data-scale sanity check in `ValidateAIPlan()`

Add one new, narrowly-scoped deterministic check inside `ValidateAIPlan()`'s existing
per-step prerequisite loop (`R/ai-planning.R:420-439`, the `for (step in normalized)` loop that
already calls `sclet_ai_call_prerequisites()`). Immediately after that existing prerequisite
call for each step, add:

```r
if (!is.null(step$action) && identical(step$action, "run_pca")) {
    requested_ncomp <- step$params$ncomponents
    if (!is.null(requested_ncomp)) {
        requested_ncomp <- suppressWarnings(as.integer(requested_ncomp))
        max_possible <- min(ncol(object), nrow(object))
        if (!is.na(requested_ncomp) && requested_ncomp > max_possible) {
            errors <- c(errors, paste0(
                "data_scale_incompatible: step ", step$id,
                " requests ncomponents=", requested_ncomp,
                " but the object only has min(ncol, nrow)=", max_possible,
                " (ncol=", ncol(object), ", nrow=", nrow(object), ")"
            ))
        }
    }
}
```

This is deliberately a single, narrow, named check (`run_pca` vs `ncomponents`) rather than a
generic framework, per the non-goals in section 2. Use the literal error prefix
`data_scale_incompatible:` so a future, broader version of this check can be found and extended
by searching for that prefix. Do not extend this to any other action in this pass. Position this
check so it still runs even when `object` is supplied but the step's `descriptor` is `NULL` (an
unregistered action name) -- it must not depend on `descriptor` being resolvable, only on
`object` and `step$params` being available, since the existing loop already does
`next` when `descriptor` is `NULL` or has no `prerequisites` function at
`R/ai-planning.R:423-425`; the new check must run before that `next`, not after it.

### 3.3 `human_confirmations` on the validation result

Add a `human_confirmations` field to the `sclet_ai_plan_validation` object returned by
`ValidateAIPlan()` (`R/ai-planning.R:461-472`, the `validation <- list(...)` block). This is a
read-only, derived summary -- it does not call `ResolveAIClarifications()` or
`sclet_ai_format_clarification()`, and it does not add a new confirmation mechanism. It only
inspects two things already available inside `ValidateAIPlan()`:

1. **Design confirmations already on record.** If `object` is supplied, call
   `GetAnalysisLedger(object)$analysis_story$design_confirmations` (added in P1, confirmed
   present at `R/ai-context.R`; read it first to confirm the exact field names
   `id`/`created_at`/`design` before using them -- do not assume). For each entry, add:
   ```r
   list(role = "<role>", column = "<column name>", status = "confirmed",
        confirmed_at = "<created_at or NULL>")
   ```
   one list entry per role/column pair inside that record's `design` list.

2. **Design confirmations the plan's own errors indicate are still missing.** If
   `errors` (the vector already being built earlier in `ValidateAIPlan()`) contains any string
   matching the literal substring `"design_semantics_not_confirmed"` (this is the existing error
   produced by `run_integration`'s prerequisites when design is unconfirmed -- verify the exact
   substring by reading `R/ai-execution.R`'s `run_integration` action prerequisites before
   relying on it; do not invent a different substring if the real one differs), add:
   ```r
   list(role = NA_character_, column = NA_character_, status = "required_not_confirmed",
        confirmed_at = NULL)
   ```
   one entry per such error (do not deduplicate across multiple steps requesting the same thing
   in this pass; that refinement is not required for this spec).

If neither case applies, `human_confirmations` is `list()`. This field must never contain a
colData value, only role names, column names (schema metadata, not data), and timestamps already
present in the records.

### 3.4 Dry-run shows success criteria, stop conditions, and confirmations

Extend `ExecuteAIPlan()`'s `dry_run` branch (`R/ai-execution.R:1592-1615`) to surface the plan's
own `success_criteria`, `stop_conditions`, and the validation's `human_confirmations`, without
evaluating any of them against real execution output (there is no real output in dry-run; do
not call `sclet_ai_assess_success_criteria()` here, since that function assumes a completed
execution and real `outputs`/`object_after` -- calling it on dry-run data would silently produce
meaningless `not_available` results that look like real assessments). Add these as new,
purely-descriptive fields alongside the existing dry-run return list:

```r
list(
    object = object,
    plan = plan,
    validation = validation,
    results = results,
    status = "dry_run",
    dry_run = TRUE,
    recorded = FALSE,
    execution_id = NULL,
    success_assessment = list(),
    success_criteria_preview = plan$success_criteria %||% list(),
    stop_conditions_preview = plan$stop_conditions %||% list(),
    human_confirmations = validation$human_confirmations %||% list()
)
```

Keep `success_assessment = list()` exactly as it already is (do not populate it in dry-run --
that would misrepresent an unevaluated criterion as assessed). The three new `*_preview`/
`human_confirmations` fields are additive; no existing dry-run consumer reads a field by that
exact name today (verify this with `grep -rn "success_criteria_preview\|stop_conditions_preview"
R/ tests/` before relying on it -- it should return nothing, confirming these are genuinely new
names with no collision).

## 4. Required test scenarios (`tests/testthat/test-ai-plan-stop-conditions.R`)

Write one new test file. Each scenario must be a real call through the standard pipeline
(`new_sclet_ai_plan()` -> `ValidateAIPlan()` -> `ExecuteAIPlan()` where applicable), never a
direct call to an internal helper that bypasses `ValidateAIPlan()`'s own dispatch:

1. `new_sclet_ai_plan()` with no `stop_conditions` argument: `plan$stop_conditions` is
   `list()`, and every other P0.1 field default is unchanged (regression guard: construct a
   plan and assert `assumptions`/`candidate_routes`/`selected_route`/`success_criteria`/`risks`
   still default exactly as before this change).
2. A stop condition with `source = "action_output"` referencing a `step_id` not present in
   `plan$actions`: `ValidateAIPlan()` returns `valid = FALSE` with an error containing
   `"stop_condition_unresolvable"`.
3. A stop condition with `source = "evidence"` and an empty `evidence_id`: same
   `stop_condition_unresolvable` error.
4. A stop condition with `check = "gte"` and no `field`/`value`: error containing
   `"stop_condition_incomplete"`.
5. A well-formed `stop_conditions` list on an otherwise-valid plan: `ValidateAIPlan()` returns
   `valid = TRUE` (confirms the new validation does not reject correct input).
6. A plan with one `run_pca` step requesting `ncomponents` larger than
   `min(ncol(object), nrow(object))` on a small real `SingleCellExperiment`: `ValidateAIPlan()`
   returns `valid = FALSE` with an error containing `"data_scale_incompatible"`.
7. The same plan with `ncomponents` within range: `ValidateAIPlan()` returns `valid = TRUE`
   (confirms the check only fires when it should).
8. A plan with no `run_pca` step at all: confirms the data-scale check never fires for
   unrelated actions (no `data_scale_incompatible` error appears).
9. An object with a real `ai_design_confirmation` record (via `ConfirmAIDesignSemantics()`,
   same pattern as `test-ai-analysis-story.R`'s design-confirmation test): a plan validated
   against that object has `validation$human_confirmations` containing one entry with
   `status == "confirmed"` and the correct `role`/`column`.
10. A plan whose steps would trigger an unconfirmed-design error from `run_integration`'s real
    prerequisites (construct this the same way `test-ai-execution.R` already does for its
    design-confirmation tests -- read that file first to copy the exact unconfirmed-design
    setup, do not invent a different one): `validation$human_confirmations` contains one entry
    with `status == "required_not_confirmed"`.
11. `ExecuteAIPlan(..., dry_run = TRUE)` on a plan with non-empty `success_criteria` and
    `stop_conditions`: the dry-run result's `success_criteria_preview` and
    `stop_conditions_preview` equal the plan's own fields verbatim, and `success_assessment` is
    still `list()` (not populated during dry-run).
12. `ExecuteAIPlan(..., dry_run = TRUE)` result's `human_confirmations` field equals
    `validation$human_confirmations` (confirms it is threaded through, not recomputed
    differently).
13. Full suite regression guard: after this change, every existing test in
    `tests/testthat/test-ai-plan-success-criteria.R` and `tests/testthat/test-ai-execution.R`
    still passes unmodified (do not edit either file as part of this spec unless a scenario in
    this list requires a genuinely new fixture that does not already exist there).

## 5. Documentation updates

- `man/ValidateAIPlan.Rd` and `man/ExecuteAIPlan.Rd`: regenerate or hand-edit (consistent with
  how P0.1 did it, by hand editing rather than running `roxygen2::roxygenise()` for the whole
  package) to document `stop_conditions`, `human_confirmations`, and the three new dry-run
  fields. Do not run `make rd`.
- `NEWS.md`: one new top entry describing exactly what was added (stop_conditions field,
  data-scale check scoped to `run_pca`/`ncomponents` only, human_confirmations derived from
  existing design-confirmation records and existing unconfirmed-design errors, three new
  dry-run preview fields). State plainly what was *not* done: no per-action scale metadata
  framework, no new clarification type, no change to `execution$status` values.
- `.dev/ai-advanced-analysis-spec.md`: update the top status line and the plan-schema section
  (section 8, already updated once by the P0 commit) to reflect `stop_conditions` and
  `human_confirmations` as implemented. Do not claim P2 is fully done if any item in section 6
  below does not actually hold -- state precisely what shipped and what of the roadmap's P2
  wish-list remains (if anything).

## 6. Acceptance checklist (self-verify before reporting done)

- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, **0 warnings**, and a pass count
      **greater than or equal to** 1114 (the baseline measured before this change).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`.
- [ ] `git diff --check` is clean (no whitespace/line-ending errors).
- [ ] `LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-planning.R R/ai-execution.R
      tests/testthat/test-ai-plan-stop-conditions.R` finds nothing.
- [ ] `git status --short` shows changes only to: `R/ai-planning.R`, `R/ai-execution.R`,
      `man/ValidateAIPlan.Rd`, `man/ExecuteAIPlan.Rd`, `NEWS.md`,
      `.dev/ai-advanced-analysis-spec.md`, and the new test file. No other file changed; no
      `DESCRIPTION` side effect from accidentally running full roxygen regeneration (this
      happened in the P0.1 round -- avoid it the same way P0.1's follow-up fix did, by hand
      editing the two `.Rd` files instead of running `devtools::document()`).
- [ ] `new_sclet_ai_plan()` with no `stop_conditions` argument still defaults every field
      exactly as it did before this change (regression, scenario 1).
- [ ] The `run_pca`/`ncomponents` data-scale check never fires for any other action name
      (scenario 8).
- [ ] `human_confirmations` never contains a raw colData value, only role/column names
      (schema metadata) and timestamps.
- [ ] Dry-run's `success_assessment` is still always `list()` -- the three new preview fields
      are additive, not a replacement for the existing (correct) "nothing was assessed yet"
      signal.
- [ ] No `AIAction()` registry entries were added, removed, or had their `prerequisites`
      function body changed. Confirm with `grep -c 'AIAction(' R/ai-execution.R` before and
      after -- the count must be identical to the count on `HEAD` before this change.
- [ ] No file outside the authorized list in this section's `git status` item was touched.

## 7. Implementation order

1. Read this entire spec once before writing any code.
2. Read `R/ai-planning.R` in full (548 lines) and `R/ai-execution.R` lines 1548-1823 (the
   `ExecuteAIPlan` signature through its end) to confirm current behavior before changing it.
   Also read `sclet_ai_normalize_success_criteria()` (`R/ai-planning.R:142-189`) in full before
   writing `sclet_ai_normalize_stop_conditions()`, since section 3.1 requires copying its exact
   pattern rather than inventing a divergent shape.
3. Read `R/ai-context.R`'s `analysis_story$design_confirmations` builder (added in P1) to
   confirm the exact field names (`id`/`created_at`/`design`) before using them in section 3.3.
   Read `R/ai-execution.R`'s `run_integration` action prerequisites to confirm the exact
   `design_semantics_not_confirmed` error substring before relying on it in section 3.3.
4. Read `test-ai-execution.R`'s existing design-confirmation test setup before writing test
   scenario 10, so the new test reuses the same real unconfirmed-design construction rather
   than inventing a different one.
5. Implement `sclet_ai_normalize_stop_conditions()`, wire into `new_sclet_ai_plan()` (section 3.1).
6. Extend `ValidateAIPlan()`: stop-condition validation (3.1), `run_pca` ncomponents check (3.2),
   `human_confirmations` derivation (3.3).
7. Extend `ExecuteAIPlan()`'s dry-run branch (3.4).
8. Write the test file (section 4), 13 scenarios.
9. Run `devtools::test(filter = "ai-", reporter = "progress")` in the background (it takes
   roughly 70 seconds; a foreground call will likely time out), iterate until green.
10. Run `make check` in the background (it takes several minutes).
11. Update documentation per section 5.
12. Self-verify every item in section 6, listing a concrete yes/no and evidence for each.
13. Report: paste the actual terminal output of the test run and `make check`, list every file
    changed with `git status --short` and `git diff --stat`, and do not commit -- leave the
    changes in the worktree for review.









