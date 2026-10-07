# P3 Spec: Evidence-Linked Report (Reference Validation, Success-State, Optional Post-Execution Interpretation)

> Status: draft for implementation (2026-10-07)
> Parent roadmap: `.dev/ai-product-roadmap.md`, P3 ("补强 Interpret / Evidence-linked Report")
> Precondition satisfied: P0 (`92332ca`), P0.1 (`607e48d`), P1 (`9bf701d`), P2 (`1688260`) are
> all committed on `devel`.
> Baseline before this change: `devtools::test(filter = "ai-")` ->
> `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1147 ]`. `make check` -> `0 errors | 0 warnings | 0 notes`.

## 1. Problem

The roadmap's P3 lists six work items and five acceptance criteria. **Verified while writing
this spec** (every claim below was checked by reading the actual code or running it, not
assumed):

1. **`RecordAIResult()`'s claim-ceiling audit already rejects over-claimed reports.**
   `R/ai-functions.R:93-102` calls `sclet_ai_audit_result_claims()` (default `audit_claims =
   TRUE`) and refuses to record when `audit$status == "overclaimed"`. This already satisfies
   the roadmap's acceptance item "claim ceiling 违规的 report 不能被记录". Nothing to build here.

2. **`ValidateAIEvidenceRefs()` exists, is exported, is fully implemented, but is called
   nowhere except its own test file.** Confirmed with
   `grep -rn "ValidateAIEvidenceRefs(" R/ tests/testthat/*.R` -- the only matches outside
   `R/ai-evidence.R` itself are in `tests/testthat/test-ai-evidence.R`. `AIExplainAnalysis()`,
   `AIInvestigate()`, `RecordAIResult()`, and `RunAIAnalysis()` never call it.

3. **This produces a real, empirically-confirmed gap**: a finding with `claim_level` outside
   `sclet_ai_findings_needing_support` (`c("consistent_with", "suggestive", "causal")`, defined
   at `R/ai-claims.R:58`) can cite a completely nonexistent `evidence_refs` id and
   `RecordAIResult()` will record it without any error. Verified directly:
   ```r
   sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, nrow = 3, ncol = 2)))
   result <- new_sclet_ai_result(task = "t", findings = list(list(
       statement = "x", severity = "info", claim_level = "associated",
       evidence_refs = c("ev:does_not_exist_at_all"))))
   RecordAIResult(sce, result)   # records without error -- the stale ref is never checked
   ```
   For `claim_level %in% c("consistent_with", "suggestive", "causal")` the stale ref *is*
   indirectly caught, but only as a side effect of `sclet_ai_claim_ceiling()` finding zero
   resolvable support (not through `ValidateAIEvidenceRefs()`, and not for any other claim
   level). This is the single clearest, most important gap this spec closes.

4. **`success_assessment` (from P0.1/P2) is never read by any report-facing function.**
   Confirmed with `grep -rn "success_assessment" R/ai-tools.R R/ai-investigate.R
   R/ai-functions.R R/ai-adapter.R` -- no match. The roadmap's "把 success criteria 的实际状态
   纳入 report" is unimplemented.

5. **`RunAIAnalysis()` already builds a thin status-echo `report` field after execution**
   (`R/ai-functions.R:325-345`: `status`/`goal`/`plan_id`/`actions`/`errors`/`warnings`/
   `execution_id`), and it already records `plan$ai_result` via `RecordAIResult()` on success
   (line 318-324) -- but that recorded result is the *planning-stage* `ai_result`, never a
   post-execution interpretation of what the execution actually produced. The roadmap's "让
   `RunAIAnalysis()` 在执行成功后可选生成 report" means something richer than the existing
   status echo: an actual evidence-linked interpretation of the execution's outputs. This spec
   adds that as an opt-in step (not forced on, to avoid a mandatory extra LLM call every time
   `RunAIAnalysis()` succeeds).

6. **"强制输出 uncertainty、alternative explanations 和 limitations" cannot be deterministically
   enforced in R** the way the other items can -- it is a system-prompt-level request to the
   LLM, not something an R function can verify the LLM actually did (the existing
   `system_prompt` text in `AIExplainAnalysis()`/`AIInvestigate()` already asks for this). This
   spec does not claim to solve this item; see section 2's non-goals.

7. **"把 report 的 context fingerprint、analysis refs 和 user decisions 写回 ledger" is
   partially done.** `RecordAIResult()` already writes `context_fingerprint`
   (`R/ai-functions.R:117-118`) and the full `findings` (which contain each finding's own
   `evidence_refs`, including any that happen to be `user_decision` nodes) via
   `artifacts$findings`. What is missing: a distinct, queryable list of which `user_decision`
   evidence ids a report's findings actually cited, rather than that information being buried
   inside the opaque `artifacts$findings` blob.

## 2. Non-goals (do not implement these here)

- No new action handlers. No velocity/spatial/multimodal/CellRank work.
- Do not attempt to force the LLM to output uncertainty/alternative explanations/limitations
  via any new validation gate. The existing `system_prompt` wording in `AIExplainAnalysis()` and
  `AIInvestigate()` already asks for this; this spec does not add a new enforcement mechanism
  for it, because there is no deterministic way to verify an LLM's prose actually discussed
  "alternative explanations" short of another LLM call, which is out of scope. If asked to
  self-check against the roadmap's acceptance list, state plainly that this one item is
  prompt-level guidance, not a verified contract, and was already present before this spec.
- Do not change `sclet_ai_findings_needing_support` (`R/ai-claims.R:58`) or any claim-ceiling
  logic in `R/ai-claims.R`. The fix in section 3.1 is a *reference-resolvability* check, wholly
  separate from and in addition to the existing claim-strength audit -- it does not change which
  claim levels require independent support.
- Do not change `ValidateAIEvidenceRefs()`'s own logic in `R/ai-evidence.R`. It is correct and
  complete already; this spec only wires its existing output into `RecordAIResult()`.
- Do not make the new post-execution interpretation step in `RunAIAnalysis()` (section 3.3)
  run by default. It must be strictly opt-in via a new parameter defaulting to `FALSE`, because
  it issues an additional LLM call that most callers (especially existing tests and the
  beginner-facing default path) do not expect and should not pay for automatically.
- Do not touch `R/ai-comparison.R`, `R/ai-diagnostics.R`, `R/ai-design-confirmation.R`,
  `R/ai-integration-routes.R`, `R/ai-planning.R`, `R/ai-execution.R`, `R/ai-privacy.R`, or any
  `run_*`/action registration code. This spec is scoped to `R/ai-evidence.R` (read-only,
  confirming the existing `ValidateAIEvidenceRefs()` signature, no edits expected there),
  `R/ai-functions.R` (`RecordAIResult()`, `RunAIAnalysis()`) only. `R/ai-tools.R` is read-only
  reference material (see section 3.2, resolved: no code change needed there) and must not be
  edited by this spec.
- Do not add a new privacy-sensitive field. Every new field this spec introduces (section 3)
  must consist only of: evidence/step/criterion ids already present in the result or plan,
  booleans, counts, and status strings already used elsewhere in this codebase (`"met"`,
  `"not_met"`, `"not_available"`, `"resolved"`, `"unresolved"`). No raw colData values, no
  per-cell data, no secrets.
- Do not change the five existing `sclet_ai_result` required fields
  (`task`/`findings`/`evidence`/`warnings`/`recommendations`/`proposed_actions`) or
  `validate_sclet_ai_result()`'s required-field list. New information is additive metadata, not
  a new required field, so that no existing caller or test breaks.

## 3. What to build

### 3.1 Wire `ValidateAIEvidenceRefs()` into `RecordAIResult()`

This is the primary fix (section 1, item 3). In `RecordAIResult()` (`R/ai-functions.R:93-138`),
after the existing claim-ceiling audit block and before building `record`, add a reference-
resolvability check that runs for **every** finding regardless of `claim_level` (unlike the
claim-ceiling audit, which only inspects `sclet_ai_findings_needing_support` levels):

```r
if (isTRUE(audit_claims)) {
    all_refs <- unique(unlist(lapply(result$findings, function(f) {
        if (is.list(f)) as.character(f$evidence_refs %||% character()) else character()
    }), use.names = FALSE))
    if (length(all_refs)) {
        ref_check <- ValidateAIEvidenceRefs(
            object, all_refs,
            fingerprint = result$context$fingerprint %||% NULL
        )
        if (length(ref_check$errors)) {
            stop("refusing to record an AI result citing unresolvable evidence: ",
                paste(ref_check$errors, collapse = "; "),
                call. = FALSE)
        }
    }
}
```

**Verified**: `ValidateAIEvidenceRefs()`'s actual return shape (`R/ai-evidence.R:259-269`) is
`list(valid, status, refs, resolved, errors, current_fingerprint, requesting_scope,
requesting_kind)` -- `errors` is confirmed to be the correct field name, matching the code
above exactly as written; no adjustment needed. This check is gated by the same `audit_claims`
parameter the claim-ceiling audit already uses -- it is not a separate on/off switch, since both are "does this result honestly
represent its evidence" checks and should turn off together for the same callers (e.g. tests
that intentionally construct unresolvable fixtures without wanting the claim-ceiling check
either).

### 3.2 `success_assessment` state is already visible to `AIExplainAnalysis()` -- no code change needed here

**Resolved while writing this spec, do not re-investigate**: `sclet_ai_normalize_record()`
(`R/ai-context.R:512-567`) already includes the record's `"summary"` field in its base `fields`
list (line 530) for every record at every detail level, via `sclet_ai_safe_value()`. Since
`AIExplainAnalysis()` (`R/ai-tools.R:150-166`) calls `sclet_ai_context(object, target = target,
detail = "full")`, and `RecordAIResult()` (`R/ai-functions.R:93-138`) writes
`n_findings`/`n_recommendations`/`n_proposed_actions` into that same `summary` field, any new
`summary` sub-field this spec adds in sections 3.3/3.4 (`cited_user_decisions`, and
`success_assessment` once it is recorded on an `ai_*` record) is **already** passed through to
the LLM context automatically, with zero code change required in `R/ai-tools.R` or
`R/ai-context.R`. Do not edit either file. Section 4's test scenario 6 and the man-page item in
section 5 reflect this: `man/AIExplainAnalysis.Rd` is **not** edited by this spec, and
`R/ai-tools.R` is **not** in this spec's authorized file list (removed from section 2's earlier
draft wording accordingly).

### 3.3 Record `success_assessment` and cited `user_decision` ids on `RunAIAnalysis()`'s report

In `RunAIAnalysis()` (`R/ai-functions.R:205-346`), the successful-execution branch
(`R/ai-functions.R:308-345`) already has `execution$success_assessment` available (from
`ExecuteAIPlan()`, confirmed present per P0.1/P2's existing contract) but never puts it into
`report`. Add it:

```r
report = list(
    status = execution$status,
    goal = goal,
    plan_id = plan$plan_id,
    actions = vapply(plan$actions, function(x) x$action, character(1)),
    errors = validation$errors,
    warnings = validation$warnings,
    execution_id = execution$execution_id,
    success_assessment = execution$success_assessment
)
```

Also add a new, strictly-opt-in parameter `interpret = FALSE` to `RunAIAnalysis()`'s signature.
When `interpret = TRUE` and the execution actually completed (`execution$status %in%
c("completed", "completed_with_errors")`), call `AIExplainAnalysis(updated, id =
execution$execution_id, model = model)` after the existing `RecordAIResult()` call, store its
result as `report$interpretation`, and if that interpretation itself is a valid
`sclet_ai_result` (it already will be, since `AIExplainAnalysis()`'s result passes through
`sclet_ai_call()`'s own validation), also attempt to record it via
`RecordAIResult(updated, interpretation, id = paste0("ai_interpretation_", plan$plan_id))`,
wrapped in the same `tryCatch` discipline the rest of this function already uses for optional
steps (do not let a failed interpretation call abort the whole `RunAIAnalysis()` return --
catch the error, put it in `report$interpretation_error`, and still return the already-completed
execution result). When `interpret = FALSE` (the default), `report$interpretation` is simply
absent -- do not add a placeholder `NULL` field that changes the shape callers already depend
on; use the same additive-field discipline as section 3.1/3.2.

### 3.4 Record cited `user_decision` evidence ids distinctly

In `RecordAIResult()` (`R/ai-functions.R:93-138`), after computing `all_refs` in section 3.1,
also compute which of those refs are `user_decision` nodes, and add a `cited_user_decisions`
field to `record$summary`:

```r
cited_user_decisions <- character()
if (length(all_refs)) {
    all_nodes <- tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())
    cited_user_decisions <- all_refs[vapply(all_refs, function(id) {
        node <- all_nodes[[id]]
        !is.null(node) && sclet_ai_evidence_kind_is_user_decision(node$kind %||% "")
    }, logical(1L))]
}
```

Add `cited_user_decisions = cited_user_decisions` to the existing `summary = list(...)` block
(`R/ai-functions.R:128-132`), alongside the existing `n_findings`/`n_recommendations`/
`n_proposed_actions`. This makes "which human decisions did this report lean on" a directly
queryable field on the record rather than requiring a caller to re-parse `artifacts$findings`
and cross-reference evidence kinds themselves.

## 4. Required test scenarios (`tests/testthat/test-ai-evidence-linked-report.R`)

Write one new test file. Each scenario must go through the real public functions
(`RecordAIResult()`, `RunAIAnalysis()`, `AIExplainAnalysis()`), never an internal helper that
bypasses them:

1. **Regression guard**: a `new_sclet_ai_result()` whose findings all have resolvable
   `evidence_refs` (construct real evidence via `RecordAIEvidence()` first) still records
   successfully via `RecordAIResult()` -- confirms section 3.1 does not reject valid input.
2. **The exact empirically-confirmed gap this spec closes**: a finding with
   `claim_level = "associated"` (chosen because section 1 proved this level currently bypasses
   all checking) and `evidence_refs = "ev:does_not_exist_at_all"` now causes `RecordAIResult()`
   to raise an error (any error -- assert with `expect_error()`, and assert the error message
   contains either `"unresolvable"` or `"evidence ref"`, matching whatever literal text the
   implementation in section 3.1 actually uses).
3. A finding with `claim_level = "hypothesis"` and **no** `evidence_refs` at all: still records
   successfully (confirms the new check only fires when `evidence_refs` is non-empty; a
   hypothesis citing nothing is unaffected, matching existing behavior proven in section 1).
4. A finding with `claim_level = "consistent_with"` and a stale ref: still produces the
   *original* claim-ceiling error (not a different, newly-introduced error) -- confirms section
   3.1's new check does not change behavior for the claim levels that were already covered.
5. `audit_claims = FALSE`: a finding with a stale ref now records without error (confirms the
   new check is gated by the same parameter as the existing claim-ceiling audit, not a separate
   always-on gate).
6. `RunAIAnalysis()` on a plan that completes successfully: `result$report$success_assessment`
   is identical to `result$execution$success_assessment` (confirms section 3.3's first change;
   use a minimal plan/action fixture, e.g. the same `inspect_status` action pattern already used
   throughout `tests/testthat/test-ai-plan-stop-conditions.R`).
7. `RunAIAnalysis(..., interpret = FALSE)` (the default): `result$report$interpretation` is
   `NULL`/absent -- confirms the opt-in default does not change existing behavior for any
   caller that does not ask for it.
8. `RunAIAnalysis(..., interpret = TRUE)` on a plan that completes successfully, with a mocked
   `sclet.ai.create_agent`/structured-output option (reuse the exact mocking pattern already
   used in `tests/testthat/test-ai-integration.R`'s `AIStatus()` tests -- read that file's setup
   before writing this scenario, do not invent a different mocking approach): confirms
   `result$report$interpretation` is populated and is a valid `sclet_ai_result` (check with
   `validate_sclet_ai_result(result$report$interpretation, error = FALSE)` returning `TRUE`).
9. `RunAIAnalysis(..., interpret = TRUE)` where the mocked interpretation call itself errors:
   confirms the error is caught, `result$report$interpretation_error` is populated with a
   character message, and `result$status`/`result$execution` are unaffected (the already-
   completed execution is still returned, not discarded because of the optional interpretation
   failure).
10. `RecordAIResult()` on a result whose findings cite a real `user_decision` evidence node
    (construct one via the same pattern as `tests/testthat/test-ai-analysis-story.R`'s
    `user_decisions` test, or via `RecordAIEvidence(..., kind = "user_decision")` directly):
    confirms the recorded record's `summary$cited_user_decisions` contains that evidence id.
11. `RecordAIResult()` on a result whose findings cite no `user_decision` nodes: confirms
    `summary$cited_user_decisions` is `character()` (empty, not absent -- the field must always
    be present on every record going forward, even when empty, so a caller can rely on its
    presence rather than checking for `NULL` first).
12. Full suite regression guard: after this change, every existing test in
    `tests/testthat/test-ai-claims.R`, `tests/testthat/test-ai-evidence.R`, and
    `tests/testthat/test-ai-execution.R` still passes unmodified (do not edit any of them as
    part of this spec unless a scenario above requires a genuinely new fixture that does not
    already exist there).

## 5. Documentation updates

- `man/RecordAIResult.Rd` and `man/RunAIAnalysis.Rd`: hand-edit (do **not** run
  `devtools::document()` or `roxygen2::roxygenise()` -- this caused an unauthorized
  `DESCRIPTION` change in the P0.1 round and was explicitly avoided in P1/P2 by hand-editing
  `.Rd` files directly) to document: `RecordAIResult()`'s new evidence-reference-resolvability
  check (gated by the existing `audit_claims` parameter) and the new `summary$cited_user_decisions`
  field; `RunAIAnalysis()`'s new `interpret` parameter, the new
  `report$success_assessment`/`report$interpretation`/`report$interpretation_error` fields.
- `man/AIExplainAnalysis.Rd`: **not edited**. Section 3.2 already resolved that no code change
  is needed in `R/ai-tools.R`, so there is nothing new to document there. State this plainly in
  the final report rather than padding the diff with an unnecessary documentation-only change.
- `NEWS.md`: one new top entry describing exactly what was added: the evidence-reference check
  in `RecordAIResult()` (and the specific empirically-confirmed gap it closes -- non-
  `consistent_with`/`suggestive`/`causal` claims previously citing nonexistent evidence
  unchecked), the `cited_user_decisions` field, `RunAIAnalysis()`'s `success_assessment` in its
  report and its new opt-in `interpret` parameter. State plainly what was *not* done: no
  enforcement of "uncertainty/alternative explanations/limitations" output (prompt-level only,
  pre-existing), no change to claim-ceiling logic, no change to `ValidateAIEvidenceRefs()`
  itself.
- `.dev/ai-advanced-analysis-spec.md`: update the top status line and whichever existing section
  documents `RecordAIResult()`/`RunAIAnalysis()`'s contract to reflect what shipped. State
  precisely which of the roadmap's P3 acceptance criteria now hold and which (if any) remain
  open -- in particular, be explicit that "强制输出 uncertainty..." was already prompt-level
  before this spec and remains so; do not claim it as newly solved.

## 6. Acceptance checklist (self-verify before reporting done)

- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, **0 warnings**, and a pass count
      **greater than or equal to** 1147 (the baseline measured before this change).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`. If a transient, unrelated
      environment error occurs (e.g. a missing-package error during the isolated install that
      is not reproducible and unrelated to any file this spec touches), verify it is transient
      by re-running `make check` once more before reporting a result -- this exact situation
      occurred once during the P2 round and was confirmed transient by a clean re-run.
- [ ] `git diff --check` is clean (no whitespace/line-ending errors).
- [ ] `LC_ALL=C grep -nP '[^\x00-\x7F]' R/ai-functions.R
      tests/testthat/test-ai-evidence-linked-report.R` finds nothing.
- [ ] `git status --short` shows changes only to: `R/ai-functions.R`,
      `man/RecordAIResult.Rd`, `man/RunAIAnalysis.Rd`, `NEWS.md`,
      `.dev/ai-advanced-analysis-spec.md`, and the new test file. `R/ai-tools.R` and
      `man/AIExplainAnalysis.Rd` must NOT appear (section 3.2 resolved that no change is needed
      there). No `DESCRIPTION` side effect.
- [ ] The new evidence-reference check in `RecordAIResult()` is gated by the same
      `audit_claims` parameter as the existing claim-ceiling audit, not a new separate flag
      (scenario 5).
- [ ] The exact empirically-confirmed gap from section 1 item 3 is closed and covered by a
      test that would have failed before this change (scenario 2) -- re-run that exact
      reproduction from section 1 against the fixed code to confirm it now raises an error.
- [ ] `RunAIAnalysis()`'s `interpret` parameter defaults to `FALSE` and changes no observable
      behavior for any existing caller that does not pass it (scenario 7, plus the full
      regression guard in scenario 12).
- [ ] `summary$cited_user_decisions` is present (possibly empty) on every record written by
      `RecordAIResult()` going forward, never silently absent (scenarios 10-11).
- [ ] No `AIAction()` registry entries were added, removed, or changed. Confirm with
      `grep -c 'AIAction(' R/ai-execution.R` -- the count must be identical to `HEAD`'s current
      count of 17 (this spec does not touch `R/ai-execution.R` at all, so this should be
      trivially true, but confirm it explicitly rather than assuming).
- [ ] No file outside the authorized list in this section's `git status` item was touched.

## 7. Implementation order

1. Read this entire spec once before writing any code.
2. Read `R/ai-functions.R` in full (693 lines) to confirm current `RecordAIResult()` and
   `RunAIAnalysis()` behavior before changing either.
3. Read `R/ai-evidence.R`'s `ValidateAIEvidenceRefs()` (`R/ai-evidence.R:215-269`) in full,
   including its exact return shape, before writing section 3.1's integration code. Correct the
   assumed field name in section 3.1's draft code if the real one differs.
4. Section 3.2 is already resolved (no code change needed in `R/ai-tools.R`); do not
   re-investigate it, do not edit that file.
5. Implement section 3.1 (`RecordAIResult()` evidence-reference check).
6. Implement section 3.4 (`cited_user_decisions`), in the same function, same change.
7. Implement section 3.3 (`RunAIAnalysis()`'s `success_assessment` in report, and the opt-in
   `interpret` parameter).
8. Write the test file (section 4), 12 scenarios. For scenario 8's mocking setup, read
   `tests/testthat/test-ai-integration.R`'s existing `sclet.ai.create_agent` mock pattern first
   and reuse it rather than inventing a new one.
9. Run `devtools::test(filter = "ai-", reporter = "progress")` in the background (it takes
   roughly 70-80 seconds; a foreground call will likely time out), iterate until green.
10. Run `make check` in the background (it takes several minutes). If it reports an error
    about a missing package unrelated to any file this spec touches, re-run once to check for
    transience before treating it as a real failure.
11. Update documentation per section 5.
12. Self-verify every item in section 6, listing a concrete yes/no and evidence for each.
13. Report: paste the actual terminal output of the test run and `make check`, list every file
    changed with `git status --short` and `git diff --stat`, and do not commit -- leave the
    changes in the worktree for review.






