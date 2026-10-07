# P5 Spec (Slice 2): Velocity Execution Robustness

> Status: draft for implementation (2026-10-07)
> Parent roadmap: `.dev/ai-product-roadmap.md`, P5, second item in the recommended order:
> "velocity readiness -> **velocity execution** -> trajectory / velocity interpretation ->
> CellRank / fate -> spatial -> multimodal".
> Precondition satisfied: P0 (`92332ca`) through P5 slice 1 (`e4c5e75`) are committed on `devel`.
> Baseline before this change: `devtools::test(filter = "ai-")` ->
> `[ FAIL 0 | WARN 24 (benign velociraptor/scuttle deprecations) | SKIP 1 | PASS 1249 ]`.
> `make check` -> `0 errors | 0 warnings | 0 notes`.

## 1. Problem

Slice 1 (`e4c5e75`) delivered velocity readiness, the `run_velocity` action, and bounded
evidence. The roadmap's second item, "velocity execution", is distinct from slice 1: slice 1
proved the action runs once with `mode = "deterministic"` and produces evidence; this slice
proves the already-generic execution machinery (`success_criteria`, `stop_conditions`,
dry-run preview, failure classification) genuinely works for `run_velocity`, and that
`run_velocity` is robust across its own parameter space, not just happy-path.

**Verified while writing this spec** (every claim below was checked by reading the actual
code, not assumed):

1. **`sclet_ai_assess_success_criteria()` (`R/ai-execution.R:1450-1608`, confirmed by reading
   it in full) is fully generic**: it keys off `crit$step_id` and `outputs[[crit$step_id]]`, with zero
   action-name-specific branching. The same holds for the stop-conditions check in
   `ValidateAIPlan()` (P2) and the dry-run preview fields (P2's
   `success_criteria_preview`/`stop_conditions_preview`/`human_confirmations`). **None of this
   needs new code for `run_velocity`** -- it already works by construction, the same way it
   already works for `run_pca`/`run_de_test`/`run_trajectory`. What is missing is a **test**
   proving this, since `tests/testthat/test-ai-velocity.R` (slice 1) has zero `success_criteria`
   or `stop_conditions` scenarios -- confirmed with
   `grep -n "success_criteria\|stop_conditions" tests/testthat/test-ai-velocity.R` returning
   nothing.

2. **`run_velocity` has not been exercised with `mode = "stochastic"` or `mode = "dynamical"`
   anywhere in the test suite.** Confirmed with
   `grep -rn "mode.*stochastic\|mode.*dynamical" tests/testthat/test-ai-velocity.R` returning
   nothing -- slice 1's tests only ever used the default `"deterministic"` mode. `velociraptor`
   is confirmed installed in this environment (section 1 of the slice 1 spec already verified
   this), so these can be exercised for real, not mocked.

3. **The dry-run branch (`ExecuteAIPlan(..., dry_run = TRUE)`) has never been exercised for
   `run_velocity`.** Confirmed with `grep -n "dry_run" tests/testthat/test-ai-velocity.R`
   returning nothing. Slice 1's scenario 5 (prerequisite failure) goes through
   `RunAIAnalysis()`'s validation-rejection path, which is a different code path from a
   validated plan's dry-run preview (P2's `success_criteria_preview`/
   `stop_conditions_preview`/`human_confirmations` fields, confirmed present on every dry-run
   result since P2, `R/ai-execution.R:1728-1729`).

4. **`run_velocity`'s own `output_schema$note`** (`R/ai-execution.R:1084-1089`, verified by
   direct read, correcting an earlier draft of this item that misquoted trajectory's note
   instead) actually states: "Velocity direction and velocity_pseudotime are model-fitted
   estimates from scVelo, not measured biological facts; they depend on the chosen mode and the
   spliced/unspliced ratio and must not be read as ground-truth cell-state transitions. No cell
   is removed or filtered by this action." -- this is prose guidance for an LLM reading the
   context, not a machine-checkable invariant. This slice adds exactly one machine-checkable
   assertion matching the one falsifiable claim in that prose ("No cell is removed or filtered
   by this action"): that executing `run_velocity` leaves `ncol()`/`colnames()` unchanged
   (mirroring the identical regression guard P4's rare-cell/doublet test already uses,
   confirmed at `tests/testthat/test-ai-domain-adapter-e2e.R`'s
   `expect_equal(ncol(result$object), ncol(sce))` pattern), since this has not actually been
   tested for `run_velocity` specifically.

5. **Invalid-`mode` handling**: direct inspection of `RunVelocity()` (`R/velocity.R:19-25`)
   shows it already calls `mode <- match.arg(mode)` immediately, so a typo'd mode string
   already fails loudly with R's standard `match.arg` error -- it does **not** silently fall
   back to `"deterministic"`. This is correct, existing behavior with no gap to fix. This slice
   adds a test asserting this holds through the full orchestrator layer (a typo'd `mode`
   surfaces as a structured action failure via `RunAIAnalysis()`, not that the underlying
   behavior itself needs changing), since no existing test exercises this path above
   `RunVelocity()` itself.

## 2. Non-goals (do not implement these here)

- No new mechanism in `R/ai-planning.R` or `R/ai-execution.R`'s success-criteria/stop-conditions
  machinery. Section 1 item 1 confirmed it is already fully generic. This slice is tests only,
  unless a genuine defect is found while writing them (same discipline as every prior round:
  fix narrowly, document precisely what and why).
- No trajectory/velocity combined interpretation (still the roadmap's third item, still
  deferred). No CellRank/fate/spatial/multimodal (still deferred).
- No change to `R/velocity.R`'s `RunVelocity()` itself unless a genuine, reproduced defect is
  found. Section 1 item 5 already confirmed its `match.arg(mode)` behavior is correct as-is.
- No change to `check_velocity_readiness()`, `sclet_ai_record_velocity_evidence()`, or
  `run_velocity`'s own descriptor fields (`prerequisites`, `output_schema`, `estimated_cost`,
  etc.) in `R/ai-diagnostics.R`/`R/ai-execution.R` unless a genuine defect surfaces while
  writing the tests in section 3. If no defect is found, this slice's entire diff is one new
  test file plus `NEWS.md`/`.dev/ai-advanced-analysis-spec.md`.
- Do not add tests for `mode = "dynamical"`'s numerical correctness or scientific validity --
  that is `velociraptor`'s own concern. This slice only confirms the AI-facing orchestration
  (success criteria, stop conditions, dry-run, failure classification) behaves correctly
  regardless of which valid `mode` was requested.
- Do not touch `R/ai-privacy.R`. Slice 1 already confirmed and tested that velocity evidence
  passes the privacy scanner cleanly; this slice adds no new field shape.

## 3. Required test scenarios (`tests/testthat/test-ai-velocity-execution.R`)

Write one new test file, separate from slice 1's `tests/testthat/test-ai-velocity.R` (do not
edit that file unless a scenario below requires a genuinely new fixture that does not already
exist there -- prefer adding a new helper in the new file). Reuse slice 1's fixture-construction
style exactly (`NormalizeData()` + `FindVariableFeatures()` + `RunPCA()` on a real
spliced/unspliced object), and `skip_if_not_installed("velociraptor")` on every scenario that
executes `run_velocity` for real.

1. **`success_criteria` end-to-end for `run_velocity`**: build a plan with one `run_velocity`
   step and a `success_criteria` entry referencing that step's output (e.g.
   `field = "required_states"`, `check = "exists"`, matching the pattern already used for other
   actions in `tests/testthat/test-ai-plan-success-criteria.R` -- read that file's existing
   action_output-source scenario first and reuse its exact shape). Run through the real
   `RunAIAnalysis()` orchestrator (`confirm = "yes"`), assert
   `result$execution$success_assessment` contains one entry with `status == "met"`.
2. **`stop_conditions` end-to-end for `run_velocity`**: same pattern, using a `stop_conditions`
   entry instead (reuse `tests/testthat/test-ai-plan-stop-conditions.R`'s existing shape).
   Assert `ValidateAIPlan()` accepts the well-formed condition (`valid == TRUE`) and that it
   round-trips through dry-run's `stop_conditions_preview` (scenario 3 below covers the
   dry-run assertion itself; this scenario covers the plan-construction/validation side).
3. **Dry-run preview for `run_velocity`**: `ExecuteAIPlan(..., dry_run = TRUE)` on a validated
   plan containing a `run_velocity` step with real `success_criteria`/`stop_conditions`. Assert
   `success_criteria_preview`/`stop_conditions_preview` equal the plan's own fields verbatim
   (matching the exact assertion pattern already used in
   `tests/testthat/test-ai-plan-stop-conditions.R`'s own dry-run scenario), and that
   `success_assessment` is still `list()` (unevaluated in dry-run, per P2's own design).
4. **`mode = "stochastic"` executes successfully end-to-end**: a real `run_velocity` step with
   `params = list(mode = "stochastic", ...)`, run through `RunAIAnalysis(confirm = "yes")`.
   Assert `status == "completed"` and that evidence is recorded (reuse slice 1's evidence
   assertions).
5. **`mode = "dynamical"` executes successfully end-to-end**: same as scenario 4 with
   `mode = "dynamical"`. If this mode is substantially slower or requires additional
   `velociraptor` arguments to converge on a tiny synthetic fixture, use a slightly larger
   fixture (e.g. 60-80 cells, matching slice 1's trajectory-style fixture size) rather than
   skip the scenario -- the point is proving the orchestration path works for this mode, which
   requires it to actually complete.
6. **An invalid `mode` string fails as a structured action failure, not a silent fallback**:
   a plan with `params = list(mode = "not_a_real_mode")`, run through
   `RunAIAnalysis(confirm = "yes")`. Assert the step's result status is `"failed"` (not
   `"completed"`) and that the failure reason is traceable to the `match.arg()` error (contains
   something like "should be one of" or references the invalid value -- assert on whatever the
   real error text actually is, read it first rather than guessing the exact string).
7. **`run_velocity` never changes `ncol()`/`colnames()`**: run a real `run_velocity` step
   through `RunAIAnalysis(confirm = "yes")`, assert `ncol(result$object) == ncol(sce)` and
   `identical(colnames(result$object), colnames(sce))`, matching the P4-round pattern cited in
   section 1 item 4.
8. **Full suite regression guard**: after this change, every existing test still passes
   unmodified (do not edit `tests/testthat/test-ai-velocity.R` or any other existing file unless
   a scenario above genuinely requires it).

## 4. Documentation updates

- `NEWS.md`: one new top entry stating plainly this slice is test-only (no new mechanism, no
  new field, no change to `run_velocity`'s descriptor) unless section 3 surfaced a genuine
  defect, in which case describe that fix with the same precision as every prior round's entry.
- `.dev/ai-advanced-analysis-spec.md`: update the top status line to note that velocity
  execution (success criteria, stop conditions, dry-run, multi-mode, invalid-mode failure
  classification) now has end-to-end coverage, following the roadmap's second P5 item. Do not
  claim trajectory/velocity interpretation (the roadmap's third item) is covered by this --
  it is explicitly deferred per section 2.
- Do not edit `.dev/ai-product-roadmap.md`.

## 5. Acceptance checklist (self-verify before reporting done)

- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, and the only warnings present
      are the same benign `velociraptor`/`scuttle` deprecation notices already present in the
      baseline (no new privacy/redaction/consent warnings), with a pass count greater than the
      1249 baseline by at least 8 (one per scenario in section 3, possibly more if a scenario
      naturally splits into multiple assertions-based `test_that` blocks).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`.
- [ ] `git status --short` shows changes only to: the new test file, `NEWS.md`,
      `.dev/ai-advanced-analysis-spec.md`, and (only if a genuine defect was found and fixed)
      `R/ai-execution.R`, `R/ai-diagnostics.R`, or `R/velocity.R` -- documented precisely if so.
      `R/ai-privacy.R` and `DESCRIPTION` must NOT appear.
- [ ] `grep -c 'AIAction(' R/ai-execution.R` is still exactly 18 (unchanged from slice 1) unless
      a genuine defect required changing `run_velocity`'s own descriptor (not adding a new one).
- [ ] Scenario 6 (invalid mode) asserts on the actual real error text observed, not a guessed
      string.
- [ ] Scenario 3 (dry-run) asserts `success_assessment` is still `list()` in the dry-run result
      -- the preview fields are additive, not a replacement for the "nothing was assessed yet"
      signal (same discipline as the P2/P3 rounds' own dry-run tests).
- [ ] No file outside the authorized list above was touched.

## 6. Implementation order

1. Read this entire spec once before writing any code.
2. Read `tests/testthat/test-ai-plan-success-criteria.R` and
   `tests/testthat/test-ai-plan-stop-conditions.R` in full to copy their exact plan-construction
   and assertion idioms before writing scenarios 1-3.
3. Read `tests/testthat/test-ai-velocity.R` (slice 1) in full to reuse its fixture-construction
   helpers rather than inventing new ones for the same spliced/unspliced setup.
4. Write the test file (section 3), 8 scenarios.
5. Run `devtools::test(filter = "ai-", reporter = "progress")` in the background (expect
   120-180 seconds given multiple real `velociraptor::scvelo()` executions across different
   modes; a foreground call will certainly time out), iterate until green.
6. Run `make check` in the background (allow extra time for the same reason; if it reports an
   error unrelated to any file this spec touches, re-run once to check for transience).
7. Update documentation per section 4.
8. Self-verify every item in section 5, listing a concrete yes/no and evidence for each.
9. Report: paste the actual terminal output of the test run and `make check`, list every file
   changed with `git status --short` and `git diff --stat`, and do not commit -- leave the
   changes in the worktree for review.




