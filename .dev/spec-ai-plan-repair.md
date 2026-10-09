# Spec: Plan Repair Loop for `AIPlanAnalysis()`

Status: ready for implementation
Predecessor: `de4e1c1` (registry vocabulary in context + eager plan validation)
Baseline: all 21 `tests/testthat/test-ai-*.R` files pass

## 1. Problem (grounded in a real API run, not hypothesis)

A live DeepSeek call (`model = "deepseek:deepseek-chat"`, `DEEPSEEK_API_KEY` set,
`SCLET_RUN_ONLINE_TESTS=true`) on a 40-cell x 60-gene synthetic `SingleCellExperiment`
produced a `sclet_ai_plan` whose action names were **all real registry keys**
(`normalize_data`, `find_variable_features`, `scale_data`, `run_pca`, `find_neighbors`,
`find_clusters`, `run_de_test`, `run_umap`, `run_doublet_detection`, plus one unregistered
step), yet `plan$validation$valid` was `FALSE`, with errors of two distinct kinds:

- **Parameter type mismatches** (verified from `plan$validation$errors`):
  `parameter scale.factor has invalid type`, `parameter nfeatures has invalid type`,
  `parameter ncomponents has invalid type`, `parameter dims has invalid type`,
  `parameter k has invalid type`, `parameter resolution has invalid type`,
  `parameter min.pct has invalid type`, `parameter logfc.threshold has invalid type`.
- **Unmet prerequisites**: `prerequisites not met for action step_3: assay is not available:
  character`, `... step_4: PCA source assay/layer is not available: character`,
  `... step_5: reduction is not available: character`,
  `... step_7: cluster identity is not set`.
- Plus one unregistered action name: `action is not registered: QC threshold
  inspection/setting`.

Reading those errors closely shows the underlying cause is almost certainly a **parameter
encoding mismatch, not model incompetence**: several errors say the value's *type* is wrong
(e.g. `step_3: assay is not available: character` means the model passed the literal string
`"character"` where an assay name was expected, and `parameter ncomponents has invalid type`
where an integer was expected). The model was given `input_schema` as bare type labels
(`ncomponents = "integer"`, `assay = "character"`) in `context$available_actions`, and appears
to have echoed those type labels back as values. That makes repair information-rich and
therefore viable: the fix is to tell the model what it got wrong, in the same vocabulary it
was given.

**The loop currently stops at that point.** `AIPlanAnalysis()` returns the invalid plan and
nothing consumes `plan$validation$errors`. The model is never told why its plan was rejected
and never gets a chance to fix it, even though the errors are already computed, precise, and
actionable. This spec closes that loop.

## 2. Design

### 2.1 Public API change

Add one parameter to `AIPlanAnalysis()` (`R/ai-planning.R`, current signature ends at the
`registry = NULL, ...` parameters around line 653-661):

```r
max_repair_rounds = 1L
```

- `0L` preserves today's exact behaviour (one proposal call, no repair).
- `1L` (the default) allows one repair attempt when the first proposal fails validation.
- Values are clamped to `0:3`; a non-integer or out-of-range value is an error, not silently
  clamped — a caller asking for 99 rounds should be told no, not quietly given 3.

Default `1L` is a deliberate choice: unlike P3's `interpret` (an optional extra narrative),
repair is what makes the function's primary output usable, and the evidence above shows the
first proposal is frequently invalid. The spec author flags this as the single most
challengeable decision here — if review concludes the default should be `0L` for cost
reasons, that is a one-line change, but it should be argued explicitly rather than assumed.

### 2.2 Internal helper

Extract the current body of `AIPlanAnalysis()` (from the `registry` default through the
`plan$validation <- ValidateAIPlan(...)` assignment) into an internal helper, so the loop and
the single call cannot drift apart:

```r
sclet_ai_plan_once <- function(object, goal, model, structured_output,
                               fallback_on_structure_error, registry,
                               repair_errors = NULL, ...)
```

Behaviour:
- Identical to today's `AIPlanAnalysis()` body when `repair_errors` is `NULL`.
- When `repair_errors` is a non-empty character vector, append a repair block to the
  `system_prompt` (see 2.3). Everything else — context construction, `sclet_ai_call` task
  name `"analysis_plan"`, `new_sclet_ai_plan()` construction, attaching `plan$validation` — is
  unchanged, so a repair call produces the same object shape as a first call.
- Returns the plan (with `$validation` attached), never throws on an invalid plan
  (`strict = FALSE` is retained).

### 2.3 Repair prompt content

Build the repair block with a small **pure helper**, not inline string assembly:

```r
sclet_ai_repair_prompt_block <- function(errors)  # character vector -> character vector
```

Rationale: verified by reading `R/ai-adapter.R` that the `options(sclet.ai.call = ...)` mock
branch calls `mock(task = task, context = context, tools = tools, schema = schema,
model = model, ...)` — it does **not** receive `system_prompt`, and the mock branch returns
before the real path ever assembles `prompt`. So no mock-based test can assert anything about
the repair prompt's contents. A pure helper makes that content directly unit-testable, which
is the only way to actually verify the most safety-critical part of this change (that the
model is told not to invent actions, and is shown the real errors verbatim).

`AIPlanAnalysis()` (via `sclet_ai_plan_once()`) appends the helper's output to the existing
`system_prompt`, which already points the model at `context$available_actions`, only when
`repair_errors` is non-empty. The block must contain:

- A clear statement that the previous plan was **rejected by deterministic validation** and
  was not executed.
- The exact `plan$validation$errors` strings, verbatim, one per line. Do not paraphrase,
  summarise, or truncate them — they are the actionable signal and the model must see the
  actual text.
- Explicit instruction to fix **only** what the errors identify, preserve the rest of the
  plan, and — critically — to **omit an unfixable step rather than invent a replacement
  action name**, echoing the first-call instruction. This matters because the observed failure
  mode includes an invented action (`QC threshold inspection/setting`); a repair prompt that
  merely says "fix these errors" invites a second invented name.
- A note that `input_schema` values describe the *type* a parameter expects and must not be
  passed as parameter values themselves. This directly targets the type-label-echo failure
  mode diagnosed in section 1.

### 2.4 Loop, selection, and metadata

```r
plan <- sclet_ai_plan_once(..., repair_errors = NULL)
rounds_used <- 0L
best <- plan
while (!isTRUE(best$validation$valid) && rounds_used < max_repair_rounds) {
    rounds_used <- rounds_used + 1L
    candidate <- sclet_ai_plan_once(..., repair_errors = best$validation$errors)
    if (length(candidate$validation$errors) < length(best$validation$errors)) best <- candidate
    else break   # no improvement: stop early, do not burn further rounds
}
```

Rules:
- Validate after every round; stop as soon as `valid` is `TRUE`.
- **Keep the plan with the fewest validation errors seen** (ties go to the earliest attempt).
  Rationale: returning the *last* attempt would sometimes hand the caller a plan worse than an
  earlier one. Deterministic and directly testable.
- **Stop early when a round does not reduce the error count**, even if rounds remain.
  A model that is not improving on round 1 is unlikely to improve on round 2, and each round
  is a real, billed pair of provider calls (see 2.5). This bound must be tested.
- Record on the returned plan, in `plan$metadata`:
  - `repair_rounds_used` (integer),
  - `validation_error_counts` (integer vector, one entry per attempt, first attempt first),
    so a caller or benchmark can see the trajectory without re-running.

### 2.5 Cost note (must appear in the roxygen docs)

Confirmed by reading `sclet_ai_call()` (`R/ai-adapter.R`): because `AIPlanAnalysis()` passes
`tools = sclet_ai_tool_registry(object)`, every call takes the tool-loop path, which is an
agent round-trip **plus** a second `generate_object()` call to convert the text to the
structured result. A repair round therefore roughly doubles the provider calls of an
`AIPlanAnalysis()` invocation. The `max_repair_rounds` roxygen block must state this plainly
so a caller can make an informed cost decision, and `plan$metadata$repair_rounds_used` makes
the actual number observable after the fact.

## 3. Non-goals

- No change to `ValidateAIPlan()`'s checks or strictness. Repair consumes its output; it does
  not relax or extend it.
- No auto-execution. Repair never calls `ExecuteAIPlan()`/`RunAIPlan()`. The
  "never execute an action" contract in `AIPlanAnalysis()`'s docs stays true.
- No repair for `AIRecommendNextStep()`. It returns an `sclet_ai_result` with recommendations,
  not a validated `sclet_ai_plan`, so there is nothing to validate against. Do not add a
  speculative repair path there.
- No change to `R/ai-privacy.R`, `R/ai-context.R`, or `R/ai-tools.R`. If repair needs a new
  context field, stop and report rather than adding one — the previous round showed that new
  context fields interact with the privacy classifier in non-obvious ways, and this design
  deliberately carries the repair signal in the `system_prompt` (assembled from internal
  validation strings, containing no cell-level data) instead.
- No retry/backoff logic, no streaming, no change to `sclet_ai_call()` itself.
- Do not "fix" the model's type-label echo by rewriting `input_schema` rendering in
  `sclet_ai_describe_registry()` in this change. That is a plausible follow-up (the schema
  vocabulary shipped in `de4e1c1` may itself be part of the cause), but it is a separate
  change with its own verification burden, and mixing it in would make it impossible to tell
  whether repair or the schema rewrite did the work.

## 4. Required tests (`tests/testthat/test-ai-plan-repair.R`)

All tests use the existing `options(sclet.ai.call = function(task, context, ...) ...)` mock
pattern (see `tests/testthat/test-ai-phase3.R`'s `AIPlanAnalysis` test) so they are
deterministic and free of provider cost.

**Verified by reading `R/ai-adapter.R`'s mock branch**: the mock is invoked as
`mock(task = task, context = context, tools = tools, schema = schema, model = model, ...)` and
the branch returns before the real path assembles `prompt`. The mock therefore **cannot** see
`system_prompt`, and no mock-based test can assert the repair prompt's contents. Consequences:

- Scenarios 1-9 distinguish first-call from repair-call by **invocation count** (invalid plan
  on call 1, valid on call 2). State this in a comment in the test file.
- Prompt *content* is covered separately by scenario 10, which tests the pure
  `sclet_ai_repair_prompt_block()` helper directly. Do not attempt to assert prompt content
  through the mock; it is not observable there.

1. **First proposal valid → no repair.** Mock always returns a valid plan. Assert `valid` is
   `TRUE`, `repair_rounds_used == 0L`, mock called exactly once.
2. **Invalid then valid → exactly one repair.** Invalid on call 1, valid on call 2. Assert
   final `valid` `TRUE`, `repair_rounds_used == 1L`, mock called exactly twice, and
   `validation_error_counts` has two entries with the first `> 0`.
3. **`max_repair_rounds = 0L` → no repair attempted.** Invalid-only mock. Assert
   `repair_rounds_used == 0L` and exactly one call (the backward-compatibility guarantee).
4. **Never improves → stops early, no infinite loop.** Mock always returns the same invalid
   plan. With `max_repair_rounds = 3L`, assert exactly 2 calls (the "stop when a round does not
   reduce the error count" rule), not 4.
5. **Returns the best attempt, not the worst.** Mock returns invalid-with-3-errors, then
   invalid-with-1-error, then invalid-with-5-errors, with `max_repair_rounds = 3L`. Assert the
   returned plan is the 1-error one and `validation_error_counts` records all three in order.
6. **Repair executes nothing.** On the invalid-then-valid path, assert no `ai_execution_*`
   record appears in `GetAnalysisLedger()` and no `ExecuteAIPlan` command-log entry exists.
7. **Out-of-range `max_repair_rounds` errors.** `-1L`, `99L`, `NA`, and `"two"` each raise an
   error rather than being silently clamped or ignored.
8. **No-repair path keeps the plan shape.** With `max_repair_rounds = 0L`, assert the plan
   still carries every field the pre-change version had (`plan_id`, `task`, `rationale`,
   `actions`, `validation`, `metadata`, `assumptions`, `candidate_routes`, `selected_route`,
   `success_criteria`, `risks`, `stop_conditions`), class `sclet_ai_plan`, plus the new
   metadata entries.
9. **Full-suite regression guard** — the whole `test-ai-*` suite still passes.
10. **Repair prompt block content** (tests `sclet_ai_repair_prompt_block()` directly, since the
    mock cannot observe `system_prompt` — see the preamble):
    - Every error string passed in appears **verbatim** in the returned block (assert with
      `expect_true(any(grepl(error, block, fixed = TRUE)))` per error, using `fixed = TRUE` so
      regex metacharacters in real error text cannot cause a false pass).
    - The block states the plan was rejected and not executed.
    - The block instructs against inventing action names / omitting unfixable steps.
    - The block warns that `input_schema` values are types, not values.
    - An empty or `NULL` `errors` input yields an empty block (so callers can append
      unconditionally without special-casing), or `NULL` — pick one and assert it, but do not
      return a non-empty block for empty input.

## 5. Acceptance checklist

- [ ] `R/ai-planning.R` is the only production file changed (plus `.Rd`/`NEWS`/spec docs). If a
      change to another `R/` file turns out to be necessary, stop and report rather than making
      it — the previous round demonstrated that adjacent subsystems carry their own failure
      modes.
- [ ] All 21 existing `tests/testthat/test-ai-*.R` files still pass.
- [ ] The new test file's 9 scenarios all pass.
- [ ] `git status --short` shows only `R/ai-planning.R`, the new test file,
      `man/AIPlanAnalysis.Rd`, `NEWS.md`, and this spec. No `DESCRIPTION` side effect (do not
      run `devtools::document()`; hand-edit the `.Rd`).
- [ ] `sclet_ai_plan_once()` is not exported; `NAMESPACE` is unchanged.
- [ ] The roxygen block for `max_repair_rounds` states the cost implication from section 2.5.

## 6. Implementation order

1. Read `R/ai-planning.R` lines ~604-712 (`AIPlanAnalysis()` and `print.sclet_ai_plan()` in
   their post-`de4e1c1` state) in full before editing.
2. Read `R/ai-adapter.R`'s mock branch in `sclet_ai_call()` to determine what the
   `options(sclet.ai.call = ...)` mock can and cannot observe (see section 4's note).
3. Extract `sclet_ai_plan_once()`; keep `AIPlanAnalysis()`'s default path behaviour identical.
4. Add the loop, best-attempt selection, early-stop rule, and metadata.
5. Add `max_repair_rounds` validation (reject non-integer / out of `0:3`).
6. Update roxygen (new param, cost note, `$metadata` fields) and hand-edit
   `man/AIPlanAnalysis.Rd`.
7. Write the 9 test scenarios.
8. Run the new test file first, then the full `ai-` suite; iterate until green.
9. Report: paste the actual per-file test output, list files changed with `git diff --stat`,
   and do **not** commit — leave the changes in the worktree for review.

## 7. Post-implementation (not part of the delegated task)

After review, the loop will be verified against the **real** DeepSeek API from the
`sclet-benchmark-ai` repo, recording first-attempt validity, rounds used, the error-count
trajectory, and wall-clock cost. That measurement — not the mock tests — is what establishes
whether repair actually rescues real model plans, and it is the first honest answer to "can a
real model drive a validated plan end to end".
