# P4 Spec: End-to-End RunAIAnalysis Coverage for Existing Domain Adapters

> Status: draft for implementation (2026-10-07)
> Parent roadmap: `.dev/ai-product-roadmap.md`, P4 ("稳定已有领域 adapter")
> Precondition satisfied: P0 (`92332ca`), P0.1 (`607e48d`), P1 (`9bf701d`), P2 (`1688260`), P3
> (`ce1f3b4`) are all committed on `devel`.
> Baseline before this change: `devtools::test(filter = "ai-")` ->
> `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1172 ]`. `make check` -> `0 errors | 0 warnings | 0 notes`.

## 1. Problem

The roadmap's P4 asks to prove, before any new domain is added, that the four existing domain
adapters (integration route comparison, marker/DE/annotation, rare-cell/doublet, trajectory)
complete a real `context -> plan -> execute -> report` loop. **Verified while writing this
spec** (every claim below was checked by reading the actual code, not assumed):

1. **`RunAIAnalysis()` -- the actual orchestrator that chains plan -> validate -> dry-run ->
   confirm -> execute -> record -> report -- has only ever been driven end-to-end for two
   action names.** Confirmed with `grep -rln "RunAIAnalysis(" tests/testthat/*.R`: only
   `tests/testthat/test-ai-beginner.R` and `tests/testthat/test-ai-evidence-linked-report.R`
   call it, and the only action names appearing in their mocked plans are `inspect_status` and
   `run_integration` (confirmed with `grep -n 'action\s*=\s*"' tests/testthat/test-ai-beginner.R`).
   `run_de_test`, `run_annotation`, `run_rare_cell_detection`, `run_doublet_detection`, and
   `run_trajectory` have **never** been driven through `RunAIAnalysis()`.

2. **Every one of those five actions already has real, substantive coverage -- but only at the
   `ValidateAIPlan()`/`ExecuteAIPlan()` layer, bypassing `AIPlanAnalysis()` and `RunAIAnalysis()`
   entirely.** For example, `tests/testthat/test-ai-execution.R` (around line 595-639) builds a
   real clustered+normalized object, calls `run_annotation` through `ValidateAIPlan()` +
   `ExecuteAIPlan()` directly, and correctly asserts the resulting evidence's `claim_level` is
   bounded to `c("associated", "consistent_with")` and never contains a per-cell vector. This
   proves the *execute* layer is solid for these actions. What has never been exercised is the
   layer above it: whether a plan shaped the way `AIPlanAnalysis()` actually produces (via
   `result$proposed_actions` from `sclet_ai_call()`, mapped verbatim into `plan$actions` --
   confirmed at `R/ai-planning.R:648`) survives `RunAIAnalysis()`'s full control flow, including
   its `confirm` logic, its `RecordAIResult()` call for the plan's own `ai_result`, and now (per
   P3) its optional `interpret` step and `report$success_assessment`.

3. **The right mocking mechanism already exists and is simple.**
   `tests/testthat/test-ai-beginner.R:71` uses `options(sclet.ai.call = function(task, context,
   ...) list(answer = ..., proposed_actions = list(list(id = ..., action = ..., params =
   list(...)))))` -- a single-option mock that intercepts `AIPlanAnalysis()`'s call to
   `sclet_ai_call()` directly, in-process, with no JSON/structured-output round-trip. This means
   a mocked plan's `params` can carry real R objects (a matrix for `run_annotation`'s `ref`, a
   character vector for `labels`) with no serialization concerns -- unlike the
   `sclet.ai.create_agent`/`sclet.ai.generate_object` pair used in
   `tests/testthat/test-ai-integration.R` and `tests/testthat/test-ai-evidence-linked-report.R`,
   which is for a different call path (`AIStatus()`/`AIExplainAnalysis()`'s tool-loop-plus-
   structured-output flow, not `AIPlanAnalysis()`'s direct `sclet_ai_call()` usage). This spec
   uses the `sclet.ai.call` mock for all four new end-to-end tests, matching
   `test-ai-beginner.R`'s existing pattern, not the more complex mock.

4. **All domain packages this spec needs are installed in this environment.** Verified
   directly: `requireNamespace("slingshot", quietly = TRUE)`,
   `requireNamespace("scDblFinder", quietly = TRUE)`, and
   `requireNamespace("BiocNeighbors", quietly = TRUE)` all return `TRUE`. `harmony` is not
   installed, but is irrelevant here: the integration test fixture already in
   `tests/testthat/test-ai-execution.R` uses `method = "fastMNN"`, which has no optional
   dependency.

5. **The roadmap's acceptance bullets translate to four concrete, testable claims per domain**,
   none of which require new code in `R/`:
   - a mocked `AIPlanAnalysis()` plan naming the domain's real action, with real parameters,
     survives `RunAIAnalysis()` end to end (`confirm = "yes"`, `dry_run` preview, then real
     execution);
   - the resulting `execution$object`'s evidence (if the action records any) has the
     `claim_level`/bounded-value properties the domain's own execute-layer tests already prove,
     confirming nothing about the orchestration layer silently weakens or loses that guarantee;
   - a prerequisite failure (missing design confirmation, missing root cluster, missing
     reduction, etc. -- whichever each domain's own `prerequisites` already enforces) surfaces
     through `RunAIAnalysis()`'s `invalid_plan` status and `sclet_ai_format_clarification()`,
     not as a raw R error;
   - `RunAIAnalysis(..., interpret = TRUE)` (from P3) works for at least one of these domains,
     confirming P3's new parameter is not accidentally coupled to the two action names it was
     tested against.

6. **This spec adds no new action, no new field, no new validation rule.** It is purely a test
   file: `tests/testthat/test-ai-domain-adapter-e2e.R`. If, while writing these tests, a genuine
   defect in existing `R/` code is found (matching this session's P0.1/P1/P2/P3 pattern of the
   implementer finding and fixing real bugs along the way), fix it narrowly and document it the
   same way those rounds did -- but do not go looking for refactoring opportunities absent a
   concrete, reproduced failure.

## 2. Non-goals (do not implement these here)

- No new action handlers. No velocity/spatial/multimodal/CellRank work -- this is explicitly
  what the roadmap's P4 forbids before P4 itself is proven ("这一阶段原则上不新增领域 action").
- No change to any `AIAction()` descriptor, `prerequisites` function, or evidence-recording
  logic in `R/ai-execution.R`. If a test in section 4 reveals that an existing prerequisite
  check or evidence record is actually wrong, treat that as a real, narrowly-scoped bug fix
  (same discipline as every prior round), not a license to redesign the action.
- No change to `R/ai-planning.R`, `R/ai-functions.R`, `R/ai-privacy.R`, `R/ai-claims.R`,
  `R/ai-evidence.R`, `R/ai-diagnostics.R`, `R/ai-design-confirmation.R`,
  `R/ai-integration-routes.R`, `R/ai-comparison.R`, or `R/ai-context.R` unless a concrete,
  reproduced defect is found while writing the tests in section 4. If no defect is found, this
  spec's entire diff is one new test file plus `NEWS.md`/`.dev/ai-advanced-analysis-spec.md`.
- Do not use the `sclet.ai.create_agent`/`sclet.ai.generate_object` mock pair from
  `tests/testthat/test-ai-integration.R`/`test-ai-evidence-linked-report.R` for driving
  `AIPlanAnalysis()`. That pair mocks a different call path. Use the simpler
  `options(sclet.ai.call = function(task, context, ...) ...)` mock from
  `tests/testthat/test-ai-beginner.R:71`, which intercepts `AIPlanAnalysis()`'s own
  `sclet_ai_call()` call directly.
- Do not attempt to make every domain's end-to-end test independent of optional Bioconductor
  packages. `skip_if_not_installed("slingshot")` / `skip_if_not_installed("scDblFinder")` /
  `skip_if_not_installed("BiocNeighbors")` guards are correct and expected (matching the
  pattern already used throughout `tests/testthat/test-ai-execution.R` for these same
  packages) -- do not try to stub out the underlying algorithms themselves. All three are
  confirmed installed in this environment (section 1 item 4), so the tests will actually run
  here, but the guards keep the suite portable.
- Do not test every domain's full prerequisite-failure matrix (that already exists at the
  `ValidateAIPlan()` layer, confirmed by the existing coverage in `test-ai-execution.R` and
  `test-ai-diagnostics.R`). Section 4 asks for exactly one prerequisite-failure scenario,
  routed through the full `RunAIAnalysis()` orchestrator, per domain that has one meaningful to
  show (confirming the error surfaces as a structured `invalid_plan` report, not that every
  possible prerequisite is independently re-verified).

## 3. Fixture pattern (reuse across all four domains)

Every scenario in section 4 follows this shape. Read `tests/testthat/test-ai-beginner.R:67-102`
first and copy its exact mocking idiom; do not invent a different one.

```r
sce <- <build a small real SingleCellExperiment, run whatever deterministic R workflow
         steps the target action needs as prerequisites -- e.g. NormalizeData() +
         FindVariableFeatures() + RunPCA() + FindClusters() for annotation/DE/trajectory,
         or just a counts assay for doublet detection>

old <- options(sclet.ai.call = function(task, context, ...) {
    list(
        answer = "<a short plan rationale string>",
        proposed_actions = list(list(
            id = "<step id>",
            action = "<the real registered action name>",
            params = list(<real, correctly-typed parameters the action's own
                            prerequisites/handler expects>)
        ))
    )
})
on.exit(options(old), add = TRUE)

registry <- AIDefaultExecutionRegistry(sce, include = c("read", "<the domain's own
    include-group name>"))
result <- RunAIAnalysis(sce, goal = "<a short natural-language goal>", confirm = "yes",
    registry = registry)
```

Then assert on `result$status`, `result$execution$status`,
`GetAnalysisLedger(result$object)`'s evidence/analyses, and (for the prerequisite-failure
scenario) `result$report$clarification`.

**Known include-group names** (confirmed in `R/ai-execution.R:179`):
`"annotation"` (covers both `run_de_test` and `run_annotation`), `"rare_cell"` (covers both
`run_doublet_detection` and `run_rare_cell_detection`), `"trajectory"` (covers `run_trajectory`),
`"integration"` (covers `run_integration`, already has one `inspect_status`/`run_integration`
end-to-end test in `test-ai-beginner.R`, so section 4 does not require a new one for
integration -- see section 4 item 4).

## 4. Required test scenarios (`tests/testthat/test-ai-domain-adapter-e2e.R`)

Write one new test file. Every scenario must go through the real `RunAIAnalysis()` public
function using the fixture pattern in section 3 -- never call `ValidateAIPlan()`/
`ExecuteAIPlan()` directly (that layer is already covered; the point of this spec is the layer
above it):

1. **Marker/DE end-to-end**: build a real clustered object (`NormalizeData()` +
   `FindVariableFeatures()` + `RunPCA()` + `FindNeighbors()` + `FindClusters()`, matching the
   exact fixture already used in `tests/testthat/test-ai-execution.R` around line 595-604 --
   read it and reuse the same construction, do not invent a different one). Mock a plan whose
   single action is `run_de_test` with no `ident.1`/`ident.2` (triggers `FindAllMarkers`
   behavior per `R/ai-execution.R:558`). Run through `RunAIAnalysis(confirm = "yes")`. Assert
   `result$status == "completed"`, and that `GetAnalysisLedger(result$object)`'s evidence
   contains an `"ev:findallmarkers"` node (or the id the plan's `name` param would produce) with
   `claim_level == "associated"` (matching `R/ai-execution.R:608`'s hardcoded claim level for
   this action) and no per-cell-length atomic vector among its `values` (reuse the exact
   bounded-value check pattern from `tests/testthat/test-ai-execution.R:636-639`).

2. **Annotation end-to-end**: using the same clustered fixture, mock a plan whose action is
   `run_annotation` with real `ref`/`labels` parameters shaped exactly like
   `tests/testthat/test-ai-execution.R`'s existing annotation fixture (read lines 606-620 and
   reuse the same reference matrix construction, do not invent a different one -- the mock's
   `proposed_actions[[1]]$params$ref` can be the real matrix object directly, per section 1 item
   3). Run through `RunAIAnalysis(confirm = "yes")`. Assert `result$status == "completed"`,
   that the resulting evidence's `claim_level %in% c("associated", "consistent_with")` (per the
   existing execute-layer test's own assertion), and that `ActiveIdent()` is unchanged
   (annotation must not silently overwrite cluster identity, matching
   `tests/testthat/test-ai-execution.R:626`'s existing check).

3. **Rare-cell/doublet end-to-end**: `skip_if_not_installed("BiocNeighbors")` and
   `skip_if_not_installed("scDblFinder")`. Build a small real object with a `counts` assay, run
   `NormalizeData()` + `FindVariableFeatures()` + `RunPCA()` (a PCA reduction is a prerequisite
   for `run_rare_cell_detection`, confirmed at `R/ai-execution.R:897-904`). Mock a plan with
   **two** steps: `run_doublet_detection` (no required params) followed by
   `run_rare_cell_detection` with `reduction = "PCA"`. Run through
   `RunAIAnalysis(confirm = "yes")`. Assert `result$status == "completed"`, that
   `colData(result$object)` gained `scDblFinder.class`, and that no cell was removed
   (`ncol(result$object) == ncol(sce)` -- matching the "label only, never filter" contract
   documented in both actions' `output_schema$note`, confirmed at `R/ai-execution.R:857` and
   `921`).

4. **Trajectory end-to-end**: `skip_if_not_installed("slingshot")`. Using the same clustered
   fixture as scenario 1/2 (needs a real cluster column and a reduction), determine one real
   existing cluster label from `Idents()` (do not hardcode a guessed label; read it from the
   actual fixture at test-construction time) and mock a plan with action `run_trajectory`,
   `group` set to the active identity's column name, `start_cluster` set to that real cluster
   label, `reduction = "PCA"` (or `"UMAP"` if one was computed -- use whichever reduction the
   fixture actually has; read `R/ai-execution.R:941` to confirm the action's own default is
   `"UMAP"` and either compute one or pass `reduction` explicitly, do not assume). Run through
   `RunAIAnalysis(confirm = "yes")`. Assert `result$status == "completed"` and that
   `GetAnalysisLedger(result$object)`'s evidence contains an `"ev:trajectory_<name>"` node
   (matching the id format at `R/ai-execution.R:1084`, where `<name>` is the plan's `name` param
   or the action's default `"slingshot"`) with `claim_level == "consistent_with"` (the
   hardcoded value at `R/ai-execution.R:1087`) and `values$pseudotime_is_absolute_time ==
   FALSE` / `values$relative_ordering_only == TRUE` (confirmed fields at
   `R/ai-execution.R:1072-1073`, guarding against pseudotime being misread as real time).

5. **One prerequisite-failure routed through the full orchestrator**: pick exactly one domain
   from scenarios 1-4 (trajectory is the clearest, since `start_cluster_missing` produces a
   human-readable, catalog-mapped clarification per `sclet_ai_clarification_catalog()`'s
   existing `"trajectory_root"` entry -- confirmed at `R/ai-functions.R:362/367`). Mock a plan
   whose `run_trajectory` step omits `start_cluster`. Run through
   `RunAIAnalysis(confirm = "yes")`. Assert `result$status == "invalid_plan"` (not a raised R
   error) and that `result$report$clarification$status == "clarification_required"` with at
   least one question whose `id == "trajectory_root"` (confirming the full chain from a real
   prerequisite failure to a structured, catalog-mapped clarification works end to end, not
   just at the `ValidateAIPlan()` layer where `test-ai-beginner.R` already partially covers this
   pattern for a different error type).

6. **`interpret = TRUE` works for a domain other than the two already tested in P3**:
   **Resolved while writing this spec**: `sclet_ai_call()` checks `getOption("sclet.ai.call")`
   first (`R/ai-adapter.R:225,260`) before falling through to the real/structured-output path
   that the `sclet.ai.create_agent`/`generate_object` pair governs -- and *every* AI-facing
   entry point (`AIPlanAnalysis()`, `AIExplainAnalysis()`, `AIStatus()`, `AIInvestigate()`, ...)
   ultimately calls `sclet_ai_call()`. So one `sclet.ai.call` mock function covers both the plan
   call and the interpretation call; dispatch on the mock's own `task` argument (`"analysis_plan"`
   for the planning call, `"analysis_explanation"` for `AIExplainAnalysis()`'s call -- confirmed
   literal task strings at `R/ai-planning.R` and `R/ai-tools.R:154`) to return a different canned
   response for each. Reuse scenario 1's or 2's fixture, extend its single mock to branch on
   `task`, call `RunAIAnalysis(..., interpret = TRUE)`, and assert `result$report$interpretation`
   is a valid `sclet_ai_result` (via `validate_sclet_ai_result(result$report$interpretation,
   error = FALSE)`) and (per the P3 round's own lesson) that `result$report$interpretation_error`
   is `NULL` and the interpretation was actually recorded on `result$object`'s ledger -- do not
   repeat the P3 round's own silent-failure mistake of only checking that `interpretation` is
   non-null.

7. Full suite regression guard: after this change, every existing test still passes unmodified
   (do not edit any pre-existing test file as part of this spec unless a scenario above
   requires a genuinely new fixture that does not already exist, and even then prefer adding a
   new helper in the new file over editing an existing one).

## 5. Documentation updates

- `NEWS.md`: one new top entry stating plainly that this is test-only (no new action, no new
  field, no new validation rule) unless section 4 surfaced a genuine defect, in which case
  describe that fix with the same precision as every prior round's entry (exact file, exact
  bug, exact fix, exact regression test that would have caught it).
- `.dev/ai-advanced-analysis-spec.md`: update the top status line to note that P4's four named
  adapters now have `RunAIAnalysis()` end-to-end coverage. Do not claim more than what the five
  acceptance bullets in section 1 item 5 actually establish -- in particular, do not claim
  every domain's full prerequisite/failure matrix is covered through the orchestrator (only one
  example scenario is required per section 4 item 5), and do not claim this proves readiness
  for P5's new domains (P4's own stated purpose is narrower: proving the *existing* adapters
  work, not clearing P5 to start).
- Do not edit `.dev/ai-product-roadmap.md`.

## 6. Acceptance checklist (self-verify before reporting done)

- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, **0 warnings**, and a pass count
      **greater than or equal to** 1172 (the baseline measured before this change).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`. If a transient, unrelated
      environment error occurs (this happened once in the P2 round with a missing-package
      error unrelated to any touched file, and was confirmed transient by a clean re-run),
      re-run once before reporting a result.
- [ ] `git diff --check` is clean (no whitespace/line-ending errors).
- [ ] `LC_ALL=C grep -nP '[^\x00-\x7F]' tests/testthat/test-ai-domain-adapter-e2e.R` finds
      nothing, and the same check against any `R/` file touched (if section 4 surfaced a real
      defect) finds nothing either.
- [ ] `git status --short` shows changes only to: the new test file, `NEWS.md`,
      `.dev/ai-advanced-analysis-spec.md`, and -- only if a genuine, reproduced defect was
      found and fixed -- the specific `R/` file(s) and `man/*.Rd` file(s) that defect required.
      No `DESCRIPTION` side effect (do not run `devtools::document()` or
      `roxygen2::roxygenise()`; hand-edit any `.Rd` file if one genuinely needs updating).
- [ ] Every one of the six scenarios in section 4 uses `RunAIAnalysis()` as its entry point,
      not `ValidateAIPlan()`/`ExecuteAIPlan()` directly (grep the new test file for
      `RunAIAnalysis(` and confirm it appears at least six times, once per scenario, outside of
      comments).
- [ ] Scenario 3 confirms no cell was removed by either rare-cell action
      (`ncol(result$object) == ncol(sce)`), matching the "label only, never filter/remove"
      contract already documented on both actions.
- [ ] Scenario 5's failure is asserted as `result$status == "invalid_plan"`, never as a raised
      R error reaching the test via `expect_error()` -- the whole point is that
      `RunAIAnalysis()` already converts this into a structured, non-throwing report.
- [ ] No `AIAction()` registry entries were added, removed, or changed unless section 4
      surfaced a genuine defect requiring one. Confirm with `grep -c 'AIAction(' R/ai-execution.R`
      -- the count must be identical to the pre-change count of 17, or you must explicitly
      justify and document why it changed.
- [ ] No file outside the authorized list in this section's `git status` item was touched.

## 7. Implementation order

1. Read this entire spec once before writing any code.
2. Read `tests/testthat/test-ai-beginner.R:1-102` in full to internalize the exact
   `sclet.ai.call` mocking idiom this spec reuses throughout.
3. Read `tests/testthat/test-ai-execution.R`'s existing fixtures for annotation (lines
   ~595-639), rare-cell/doublet, and trajectory actions, to reuse their exact object
   construction rather than inventing new ones.
4. Section 4 scenario 4's evidence assertion is already resolved (points to
   `"ev:trajectory_<name>"`, `claim_level == "consistent_with"`); do not re-investigate.
5. Section 4 scenario 6's mocking mechanism is already resolved (one `sclet.ai.call` mock,
   dispatching on `task`); do not re-investigate, do not reach for the
   `sclet.ai.create_agent`/`generate_object` pair.
6. Write the test file (section 4), 6 scenarios (grouped as needed; it is fine to have more
   than 6 `test_that()` blocks if a scenario is naturally split into a setup check plus the
   main assertion, but every scenario's substance must be present).
7. Run `devtools::test(filter = "ai-", reporter = "progress")` in the background (it takes
   roughly 80-85 seconds; a foreground call will likely time out), iterate until green.
8. Run `make check` in the background (it takes several minutes).
9. Update documentation per section 5.
10. Self-verify every item in section 6, listing a concrete yes/no and evidence for each.
11. Report: paste the actual terminal output of the test run and `make check`, list every file
    changed with `git status --short` and `git diff --stat`, and do not commit -- leave the
    changes in the worktree for review.





