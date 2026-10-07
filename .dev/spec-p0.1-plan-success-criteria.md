# P0.1 Spec: Plan Success Criteria / Candidate Routes Contract

> Status: draft for implementation (2026-10-06)
> Parent roadmap: `.dev/ai-product-roadmap.md`, confirmed by independent review (subagent-74, 2026-10-06)
> Scope: this spec covers P0.1 only. P0.2 (interpretation report) and P0.3 (claim ceiling on
> read-only entry points) are separate specs that depend on this one landing first.

## 1. Problem

`sclet_ai_plan` currently only carries `plan_id`, `task`, `rationale`, `context_fingerprint`,
`actions`, `requires_confirmation`, `ai_result`, `metadata` (`R/ai-planning.R:152-176`). There is
no field anywhere in the codebase for *why* a route was chosen or *how success is measured*:
`grep -rn "success_criteri" R/` returns nothing. `ValidateAIPlan()` checks action registration,
parameter schemas, dependency ordering, and prerequisites, but nothing about whether the plan can
ever be judged to have succeeded or failed.

This spec adds that missing layer without touching execution semantics, the action registry, or
any existing passing test's expectations.

## 2. Non-goals (do not implement these here)

- No new action handlers. No velocity/spatial/multimodal/CellRank work.
- No automatic route selection by the AI. `selected_route` is always a value the plan declares;
  nothing in this spec picks it algorithmically.
- No change to `ExecuteAIPlan()`'s action-execution loop itself (the per-step `for` loop). Only
  what happens with the already-computed `results`/`outputs` after that loop, to assess criteria.
- No rendering/report/narrative layer (`R/ai-tools.R`'s read-only entry points are untouched).
  `AIExplainAnalysis()` consuming `success_assessment` is P0.2, not this spec.
- No change to `sclet_ai_claim_ceiling()` or `sclet_ai_audit_result_claims()` semantics.
- Do not require `success_criteria` on every plan. A plan with zero criteria stays valid; this is
  additive, not a hard gate that breaks the current 488-pass AI test baseline.

## 3. Schema additions

### 3.1 `new_sclet_ai_plan()` (`R/ai-planning.R`)

Add five new parameters, all optional, all defaulting to an empty list so every existing caller
(including every existing test) keeps working unchanged:

```r
new_sclet_ai_plan <- function(
    task = "analysis_plan",
    actions = list(),
    context_fingerprint = NULL,
    rationale = NULL,
    requires_confirmation = TRUE,
    ai_result = NULL,
    metadata = list(),
    assumptions = list(),
    candidate_routes = list(),
    selected_route = NULL,
    success_criteria = list(),
    risks = list()
)
```

Field shapes (all are plain lists checked structurally, not S4/R6):

- `assumptions`: character vector or list of character strings. Each entry is a short
  human-readable statement the plan depends on (e.g. `"batch column confirmed by user"`).
- `candidate_routes`: list of named lists, each with at minimum `id` (character) and
  `description` (character). Example: `list(list(id = "fastmnn", description = "..."), list(id =
  "harmony", description = "..."))`. May be empty when there is only one obvious route.
- `selected_route`: `NULL`, or a single character matching one `candidate_routes[[i]]$id` when
  `candidate_routes` is non-empty. If `candidate_routes` is empty, `selected_route` must also be
  `NULL` or a free-text character (no cross-check possible).
- `success_criteria`: list of named lists. Each entry has:
  - `id`: character, unique within the plan.
  - `description`: character, human-readable.
  - `source`: one of `"action_output"` or `"evidence"`.
  - `step_id`: required when `source == "action_output"`; must match a `plan$actions[[i]]$id`.
  - `evidence_id`: required when `source == "evidence"`; a string identifier the plan author
    expects to exist in the ledger/evidence graph after execution (checked at validation time
    only insofar as the referenced step is expected to produce it — see 3.3).
  - `check`: one of `"exists"`, `"equals"`, `"gte"`, `"lte"`, `"in"`. Defaults to `"exists"`.
  - `field`: optional character, a dotted path into the step's output summary or evidence values
    (e.g. `"n_rare_clusters"`). Required when `check` is not `"exists"`.
  - `value`: required when `check` is not `"exists"`; the comparison target.
- `risks`: character vector or list of character strings describing known risks or trade-offs of
  the selected route (e.g. `"fastMNN may over-correct a true biological batch-confounded signal"`).

No new field is permitted to carry executable code, R expressions, or arbitrary closures — every
value must be a plain atomic/character/list value so it survives `jsonlite`/structured-output
round-trips untouched.

### 3.2 Normalizer

Add `sclet_ai_normalize_success_criteria(success_criteria)` in `R/ai-planning.R` (near
`sclet_ai_normalize_plan_actions()`), following the exact normalization pattern already used for
plan actions: coerce to a list, fill missing optional fields with their defaults, leave invalid
shapes for `ValidateAIPlan()` to reject rather than silently dropping them. Mirror the existing
style of `sclet_ai_normalize_plan_actions()` (read it first; this spec's reviewer should check the
new function reuses that pattern, not a divergent one).

Call this normalizer from inside `new_sclet_ai_plan()` so `plan$success_criteria` is always in
normalized shape by the time `ValidateAIPlan()` sees it, same as `actions` already is.

### 3.3 `ValidateAIPlan()` additions (`R/ai-planning.R:187-375`)

Add checks, appended to the existing `errors`/`warnings` accumulation — do not restructure the
existing control flow, only insert new checks at the appropriate points:

1. **Uniqueness**: `success_criteria[[i]]$id` must be unique within the plan (same pattern as the
   existing `anyDuplicated(ids)` check for actions, line ~220).
2. **Resolvability** (the core new rule): for each criterion with `source == "action_output"`,
   `step_id` must match a known `ids` value (the same `ids` vector already computed from
   `plan$actions`). If it does not, add error:
   `"success_criterion_unresolvable: <id> references unknown step <step_id>"`.
   For `source == "evidence"`, this spec does **not** require the evidence to already exist (it
   won't, before execution) — only that `evidence_id` is a non-empty character. Add error
   `"success_criterion_unresolvable: <id> has empty evidence_id"` if blank.
3. **Field requirement**: when `check != "exists"`, `field` and `value` must both be present and
   non-NULL. Error: `"success_criterion_incomplete: <id> requires field and value for check '<check>'"`.
4. **`selected_route` consistency**: if `candidate_routes` is non-empty, `selected_route` must be
   `NULL` or match one `candidate_routes[[i]]$id`. Error:
   `"selected_route_unknown: '<selected_route>' is not in candidate_routes"`.
5. These are errors (not warnings) — an unresolvable criterion is as invalid as an unregistered
   action, by the same logic already applied to actions in this function.

Keep `valid <- !length(errors)` as the final gate (existing line ~341) — do not change how
`valid` is computed, only what feeds into `errors`.

### 3.4 `ExecuteAIPlan()` additions (`R/ai-execution.R:1354-1559`)

After the existing per-step `for` loop finishes (after line ~1531, before the `record` block),
compute a `success_assessment` list, one entry per `plan$success_criteria[[i]]`:

```r
list(
    id = criterion$id,
    description = criterion$description,
    status = "met" | "not_met" | "not_available",
    reason = <character, always present>,
    observed_value = <value actually found, or NULL>
)
```

Resolution rules:

- `source == "action_output"`: look up `outputs[[criterion$step_id]]` (the existing `outputs`
  list already populated by the loop, keyed by step id — reuse it, do not build a parallel
  structure). If that step's `results` entry has `status == "failed"` or never ran (e.g. loop
  `break`-ed before reaching it), status is `"not_available"`, reason explains which (`"step
  <step_id> did not complete"`). Otherwise extract `field` from the output summary (when `check !=
  "exists"`) or just check presence (when `check == "exists"`), and compare per `check`.
- `source == "evidence"`: look up the referenced id via the same evidence accessor used elsewhere
  (`sclet_ai_evidence_get_all()` / `sclet_ai_evidence_get()`, confirm exact name in `R/ai-evidence.R`
  before writing) restricted to nodes written by this execution (i.e. only evidence recorded
  during this `ExecuteAIPlan()` call, not any pre-existing evidence on the object — read
  `sclet_ai_record_execution()` to see what, if anything, is already tracked as "evidence written
  this run"; if nothing tracks that yet, scope the check to `source` matching one of this plan's
  step ids, same restriction as rare-cell comparison's `sclet_ai_rare_cell_run_ids()` pattern in
  `R/ai-comparison.R`). If not found: `"not_available"`, reason `"no evidence node found for
  <evidence_id>"`.
- All comparisons (`equals`/`gte`/`lte`/`in`) must be `NA`-safe: a missing or non-comparable value
  resolves to `"not_available"`, never a cryptic R error, and never silently to `"met"`.
- This computation must be pure (no mutation, no new evidence, no new state record) — it only
  reads `outputs`/evidence and produces a plain list.

Attach the result:

```r
execution$success_assessment <- success_assessment  # possibly empty list, always present
```

on both the `dry_run` early-return branch (as an empty list, since nothing executed) and the real
execution return branch (as the computed list). This keeps the `sclet_ai_execution` shape
consistent regardless of `dry_run`.

Do not let a `"not_met"` or `"not_available"` criterion change `execution$status` — `status`
(`"completed"` / `"failed"` / `"completed_with_errors"`) continues to reflect only action
execution outcomes, exactly as today. Success-criteria evaluation is a separate, additional field,
not a new way to fail the run. This matters: conflating the two would silently change behavior
for every existing caller of `ExecuteAIPlan()` that branches on `status`.

## 4. Tests (must add, must pass)

New test file `tests/testthat/test-ai-plan-success-criteria.R`. At minimum:

1. `new_sclet_ai_plan()` with no new args produces a plan with
   `success_criteria == list()`, `candidate_routes == list()`, `selected_route == NULL`,
   `assumptions == list()`, `risks == list()` — i.e. fully backward compatible.
2. A plan whose criterion references an unknown `step_id` fails `ValidateAIPlan()` with an error
   containing `"success_criterion_unresolvable"`.
3. A plan whose criterion has `check = "gte"` but no `value` fails validation with
   `"success_criterion_incomplete"`.
4. A plan with `selected_route = "x"` and `candidate_routes = list(list(id = "y", ...))` fails
   validation with `"selected_route_unknown"`.
5. A plan with `selected_route` matching a candidate route, and criteria that reference real
   step ids, passes validation.
6. `ExecuteAIPlan()` dry run returns `execution$success_assessment` as `list()`.
7. `ExecuteAIPlan()` real execution with a criterion on a step that completes successfully (e.g.
   `check = "exists"` against an existing action's known output field) returns status `"met"`.
8. `ExecuteAIPlan()` real execution where the referenced step fails (plan has 2 steps, second
   references the first's output, first's handler is forced to error) returns `"not_available"`
   for any criterion on the un-run/failed step, and `execution$status` is still `"failed"`
   (unaffected by the criterion outcome).
9. `ExecuteAIPlan()` with an `evidence`-sourced criterion whose evidence_id was never written
   returns `"not_available"`, not an error.
10. Confirm existing `test-ai-execution.R` and `test-ai-diagnostics.R`/`test-ai-claims.R` tests
    are unaffected (same pass count or higher, zero new failures).

Use the project's existing SCE test fixtures / helpers (check `tests/testthat/helper-*.R` if
present) rather than inventing a new fixture pattern.

## 5. Documentation updates required

- `man/new_sclet_ai_plan.Rd` (if hand-written; check with `glob man/*.Rd` whether one exists —
  if `new_sclet_ai_plan` is currently undocumented/internal, confirm via `@noRd`/absence of
  `@export` before deciding whether a new `.Rd` is needed at all) — update parameter list only if
  the function is exported and documented today.
- `man/ValidateAIPlan.Rd` and `man/ExecuteAIPlan.Rd`: add `@return` detail noting
  `success_assessment` field and the new validation error codes, matching existing doc style.
- `NEWS.md`: one bullet, in the existing style used for the rare-cell-comparison entry
  (`NEWS.md` top section), describing the new plan fields and that they are additive/optional.
- `.dev/ai-advanced-analysis-spec.md`: update the plan contract section to say
  `success_criteria`/`candidate_routes`/`risks`/`assumptions` are implemented, matching the
  existing convention of marking spec items done only when code backs the claim.
- Do **not** touch `.dev/ai-product-roadmap.md` itself in this change — only `update_goal`/future
  conversation turns should re-mark roadmap phase status, not this implementation commit.

## 6. Acceptance criteria (binary, checked at review)

- [ ] `git diff --check` clean (no non-ASCII, no trailing whitespace issues) in `R/ai-*.R`.
- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, **0 warnings**, and a pass count
      **greater than or equal to** the pre-change baseline, which is confirmed as
      `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1027 ]` (measured 2026-10-06, before this change).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`.
- [ ] No new exported function changes its existing default-argument behavior for any caller that
      does not pass the five new plan fields.
- [ ] No execution action handler, registry entry, or `R/ai-execution.R` action list was added or
      removed. `grep -c 'AIAction(' R/ai-execution.R` count before and after must be identical.
- [ ] No file outside `R/ai-planning.R`, `R/ai-execution.R`, their man pages, `NEWS.md`,
      `.dev/ai-advanced-analysis-spec.md`, and the new test file was modified.
- [ ] All temporary scratch/debug files removed before final report; `git status --short` clean
      except for the intended changed files.

## 7. Suggested implementation order (for the delegated subagent)

1. Read `R/ai-planning.R` in full and `R/ai-execution.R:1354-1611` (already read for this spec;
   re-read before coding to avoid stale assumptions).
2. Read `R/ai-evidence.R` to confirm the exact evidence-accessor function name/signature used by
   `source == "evidence"` criteria.
3. Implement `sclet_ai_normalize_success_criteria()` and wire it into `new_sclet_ai_plan()`.
4. Extend `ValidateAIPlan()` with the five checks in order given in section 3.3.
5. Extend `ExecuteAIPlan()` with the `success_assessment` computation.
6. Write the test file from section 4.
7. Run `devtools::test(filter = "ai-")`, iterate until green.
8. Run `make check`, fix any new notes/warnings.
9. Update documentation per section 5.
10. Report: paste the actual test and check output, list every file changed, and self-verify every
    item in section 6 before declaring done.
