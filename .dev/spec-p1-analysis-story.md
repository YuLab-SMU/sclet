# P1 Spec: Deterministic `analysis_story`

> Status: draft for implementation (2026-10-07)
> Parent roadmap: `.dev/ai-product-roadmap.md`, P1 ("补强 Understand / Analysis Story")
> Precondition satisfied: P0 (API responsibility map, `.dev/ai-advanced-analysis-spec.md`
> new "P0. API 职责映射" section) and P0.1 (`success_criteria`/`candidate_routes` plan
> contract, commit `607e48d`) are both committed on `devel`.
> Baseline before this change: `devtools::test(filter = "ai-")` ->
> `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 1073 ]`. `make check` -> `0 errors | 0 warnings | 0 notes`.

## 1. Problem

`GetAnalysisLedger()` (`R/ai-context.R`) already exposes `analyses`, `state_records`,
`workflows`, `lineage`, `health`, `capabilities`, and `blocked_actions`. An AI reading this
ledger can see *that* a `trajectory` state record exists and *that* `has_rare_cells` is `TRUE`,
but it cannot answer, in one bounded field, "what order did these happen in, what fed into
what, which steps are still missing, and which parts of this are a human decision versus a
computed result." That is the `analysis_story` gap named in the roadmap's P1.

Two concrete, verified defects compound this (found while reading the current context builder
for this spec, not assumed):

- `sclet_ai_blocked_actions()` (`R/ai-context.R:511-520`) unconditionally returns
  `execute_analysis = list(blocked = TRUE, reason = "Phase 0 context is read-only")` regardless
  of whether an execution registry exists or any action has ever run. This is stale text from
  before Phase C added real execution (`ExecuteAIPlan()`, `RunIntegrationRoutes()`, etc.) and
  nothing in the codebase actually reads this specific field (confirmed:
  `grep -rn "blocked_actions\$execute_analysis\|blocked_actions\\[\\[.execute_analysis" R/ tests/`
  returns nothing) — but any entry point that forwards the full ledger to an LLM
  (`AIStatus()`, `AskAI()`, `AIExplainAnalysis()`, all of which call `sclet_ai_context()`/
  `GetAnalysisLedger()`) is handing the model a factually wrong claim about its own
  capabilities. This spec removes/fixes that field as part of making the context honest.
- `sclet_ai_capabilities()$execution` (`R/ai-context.R:494-509`) is hardcoded `FALSE` with the
  same staleness and the same "nothing reads it" property
  (`grep -rn "capabilities\\$execution" R/ tests/` returns nothing). Fixed in the same pass,
  same reasoning.

## 2. Non-goals

- No new action handlers, no velocity/spatial/multimodal work.
- No change to `Status()`'s existing `health`/`mainline_missing` semantics (that mechanism is
  specific to the velocity -> fate -> perturbation mainline and stays as-is).
- No change to `sclet_ai_claim_ceiling()`, `RecordAIResult()`, or any evidence/claim-level rule.
- Do not infer biological meaning, suggest a next step, or rank analyses by importance.
  `analysis_story` is a chronological/structural summary, not an interpretation (that is P3's
  job). If a record order cannot be determined (missing `created_at`), say so; never guess.
- Do not change `sclet_ai_collect_records()`'s existing return shape
  (`analyses`/`state_records`/`workflows`/`lineage`/`warnings`) — `analysis_story` is a new,
  additional field built from that existing data, not a replacement for it.
- Do not touch `R/ai-tools.R`, `R/ai-functions.R`, `R/ai-planning.R`, `R/ai-execution.R`, or any
  action registry file. This spec is scoped to `R/ai-context.R` only, plus its tests and docs.

## 3. What `analysis_story` must contain

Add a new field to the object returned by `sclet_ai_build_context()`
(`R/ai-context.R:91-110`), computed by a new function `sclet_ai_build_analysis_story()`:

```r
result$analysis_story <- sclet_ai_build_analysis_story(object, records, status)
```

Shape:

```text
analysis_story
    schema_version        "1.0"
    timeline               list of entries, oldest first (see 3.1)
    n_steps                integer, length(timeline)
    unordered_steps        character vector of record keys with no usable created_at
    user_decisions         list, see 3.2
    design_confirmations   list, see 3.3
    evidence_gaps          list, see 3.4
    conflicts              list, see 3.5
    raw_values_included    FALSE (always)
```

### 3.1 `timeline`

One entry per record already present in `records$analyses` (the flattened, deduplicated view
`sclet_ai_collect_records()` already builds — reuse it, do not re-walk `state$states$records`
or `state$analyses` separately). Each entry:

```r
list(
    key = <character, the same key used in records$analyses>,
    type = <character>,
    id = <character or NULL>,
    created_at = <character, ISO-ish string from as.character(Sys.time())-style value, or NULL>,
    status = <character or NULL, e.g. "completed">,
    method = <character or NULL>,
    depends_on = <character vector, see below>
)
```

`created_at` ordering: use whatever `created_at` value the normalized record already carries
(`sclet_ai_normalize_record()` already copies this field through when present — confirmed at
`R/ai-context.R:347-350`). Sort `timeline` ascending by parsed `created_at`
(`as.POSIXct(x, tryFormats = ...)`, falling back to string comparison if parsing fails — do not
error on an unparseable timestamp, just keep it in original encounter order relative to other
unparseable ones and place all unparseable entries at the end). Records with no `created_at` at
all go into `unordered_steps` (their `key` only) and are excluded from `timeline`'s ordering
guarantee, but are still worth including in `timeline` at the end with `created_at = NULL` so
nothing silently disappears — the key point is `unordered_steps` tells the reader which entries'
position in `timeline` is not meaningful.

`depends_on`: derive from the existing `lineage[[key]]$parents` (`records$lineage` already
computed by `sclet_ai_collect_records()`, confirmed at `R/ai-context.R:305-319`) restricted to
parent keys that are themselves present in `records$analyses` (drop dangling references rather
than erroring — same defensive pattern already used elsewhere in this file, e.g.
`sclet_ai_target_matches()`'s graceful fallbacks).

### 3.2 `user_decisions`

Pull every record from `records$state_records[["ai_evidence"]]` (if present) whose `kind` is
`"user_decision"` (the convention already established by `sclet_ai_record_clarification_response()`
and `RecordAIEvidence()` — confirm the exact field name by reading `R/ai-evidence.R` before
coding). For each:

```r
list(
    id = <character>,
    created_at = <character or NULL>,
    summary = <character, a short bounded description already present on the node — do not
        invent one; if the node has no human-readable summary field, use "user decision recorded"
        as a fixed fallback string, never fabricate specifics>
)
```

### 3.3 `design_confirmations`

Pull every record from `records$state_records[["ai_design_confirmation"]]` (if present). Read
`R/ai-design-confirmation.R` first to confirm this exact shape (verified while writing this
spec): `ConfirmAIDesignSemantics()` writes `summary$design` as a role-to-column-name mapping,
e.g. `list(batch = "batch_id", condition = "condition")`. Column *names* are schema-level
identifiers already exposed elsewhere in the ledger (e.g. `dataset$columns`), not sensitive
per-cell data — only `summary$design_value_key` (the per-column value hash used for staleness
detection) and the actual `colData` values are sensitive, and neither of those is included here.
For each confirmed record:

```r
list(
    id = <character>,
    created_at = <character or NULL, from summary$confirmed_at>,
    design = <the summary$design list itself, i.e. role -> column name,
        e.g. list(batch = "batch_id", condition = "condition")>
)
```

Do not include `design_value_key` or `object_fingerprint` from the record's `summary` — those
are internal staleness-detection machinery for `ConfirmAIDesignSemantics()`'s own prerequisite
checks, not part of a human-facing story and not useful without the raw colData to compare
against.

### 3.4 `evidence_gaps`

This is the one genuinely new piece of reasoning, kept narrow and purely structural (no claim
about biology). **Verified while writing this spec**: of the four exported `check_*_readiness()`
functions in `R/ai-diagnostics.R`, only `check_rare_cell_readiness()` currently exposes a
per-signal "missing supporting evidence" flag (`checks$doublet_evidence_available`, confirmed at
`R/ai-diagnostics.R:275-283`). `check_integration_readiness()`, `check_annotation_readiness()`,
and `check_trajectory_readiness()` only expose `status`/`blocked_actions`/`questions` — they do
not have a comparable per-signal gap flag today. Do not invent one for them in this spec.

Scope `evidence_gaps` to exactly this one case for now (narrower than originally scoped, by
design, rather than fabricating flags that do not exist):

If `"rare_cells"` is present in `records$analyses` or `records$state_records`, call
`check_rare_cell_readiness(object)` and, if `checks$doublet_evidence_available` is `FALSE`, add:

```r
list(
    analysis_type = "rare_cells",
    gap = "doublet_evidence_available",
    note = <character, copy the existing note string check_rare_cell_readiness() already
        produces for this condition verbatim from its `notes` field — do not write new prose>
)
```

If `"rare_cells"` is not present, or `doublet_evidence_available` is `TRUE`, `evidence_gaps` is
`list()`. This must be read-only and must not call any `run_*`/action handler. (A future spec
can extend this to other readiness functions once they grow a comparable flag; this spec does
not add new flags to any `check_*_readiness()` function — that would be scope creep beyond
`R/ai-context.R`.)

### 3.5 `conflicts`

Narrow, deterministic duplicate-record detection only: if `records$analyses` contains more than
one entry whose `type` is identical (e.g. two `"trajectory"` records from different runs), report:

```r
list(
    type = <character>,
    keys = <character vector of the colliding record keys>,
    note = "multiple <type> records exist; only the active one (if any) is in active_states"
)
```

Do not attempt to compare their actual output values or decide which is "better" — that is
explicitly P3/Interpret territory, not this spec.

## 4. Function placement and signature

```r
sclet_ai_build_analysis_story <- function(object, records, status)
```

Add to `R/ai-context.R`, placed after `sclet_ai_collect_records()` and before
`sclet_ai_normalize_record()` (or any position in the file that keeps related helpers together
— match the file's existing ordering convention, do not scatter). Call it from
`sclet_ai_build_context()` (`R/ai-context.R:91-110`) after `records` is computed:

```r
result$analysis_story <- sclet_ai_build_analysis_story(object, records, status)
```

`result$fingerprint <- sclet_ai_fingerprint(result)` already runs after this (line 109) and
already serializes the whole `result` list — no change needed there; `analysis_story` becomes
part of what gets fingerprinted automatically, which is correct (a change in the story should
change the fingerprint).

### 4.1 Fix the two stale fields (same file, same commit)

In `sclet_ai_blocked_actions()` (`R/ai-context.R:511-520`): replace the hardcoded
`execute_analysis = list(blocked = TRUE, reason = "Phase 0 context is read-only")` with a value
that reflects reality without overclaiming. Since this ledger builder has no registry argument
and genuinely cannot know what the caller intends to execute, the honest answer is "unknown from
context alone, not a blanket block":

```r
execute_analysis = list(
    blocked = NA,
    reason = "execution capability depends on the AIExecutionRegistry supplied at call time, not on this read-only context"
)
```

In `sclet_ai_capabilities()` (`R/ai-context.R:494-509`): same reasoning, change
`execution = FALSE` to `execution = NA` with the same rationale (context alone cannot assert
whether execution is possible — that depends on the registry passed to `ValidateAIPlan()`/
`ExecuteAIPlan()`, which this function does not see). Do not delete the field; downstream code
that merely checks for its presence must not break, and no code currently branches on its truth
value (confirmed by the greps in section 1), so changing `FALSE` to `NA` is safe.

## 5. Tests (must add, must pass)

New test file `tests/testthat/test-ai-analysis-story.R`. At minimum:

1. An SCE with zero recorded analyses: `analysis_story$timeline == list()`,
   `analysis_story$n_steps == 0L`, `unordered_steps == character()`,
   `user_decisions == list()`, `design_confirmations == list()`, `evidence_gaps == list()`,
   `conflicts == list()`.
2. An SCE that goes through `NormalizeData() -> RunPCA() -> FindNeighbors() -> FindClusters()`
   (reuse whatever minimal fixture pattern `test-ai-execution.R`'s
   `sclet_ai_test_rare_object()`-style helpers already use for a small real object): confirm
   `analysis_story$timeline` entries are in non-decreasing `created_at` order, and
   `n_steps == length(timeline)`.
3. Two sequential `RunRareCellDetection()` calls with different `name=` producing two
   `"rare_cells"` state records: `conflicts` reports one entry with `type == "rare_cells"` and
   both keys present.
4. A record manually constructed via `sclet:::sclet_set_analysis_state()` with no `created_at`
   in its `summary`/top-level fields (i.e. a record whose normalized form lacks `created_at`):
   its key appears in `unordered_steps`.
5. After `RecordAIEvidence()` records a `kind = "user_decision"` node (construct this the same
   way the existing clarification tests do — read `test-ai-beginner.R`'s
   `sclet_ai_record_clarification_response()` tests for the exact pattern first), confirm it
   appears in `analysis_story$user_decisions` with no raw design-value leakage.
6. After `ConfirmAIDesignSemantics(object, design = list(batch = "batch_col"))` on an object
   whose `colData` has a `batch_col` column with real labels (e.g. `c("donor_1", "donor_2")`):
   confirm `design_confirmations[[1]]$design == list(batch = "batch_col")` (the column name is
   expected to appear — it is schema metadata, not secret), and confirm the actual data values
   (`"donor_1"`, `"donor_2"`) do **not** appear anywhere in
   `utils::capture.output(str(result$analysis_story))`. Also confirm `design_value_key` and
   `object_fingerprint` are absent from the reported `design_confirmations` entry.
7. An object with `scDblFinder.class` missing and a `rare_cells` record present: `evidence_gaps`
   contains one entry with `analysis_type == "rare_cells"` and `gap == "doublet_evidence_available"`
   (per section 3.4).
8. `sclet_ai_blocked_actions()$execute_analysis$blocked` is `NA`, not `TRUE` or `FALSE`.
9. `sclet_ai_capabilities()`'s returned `execution` field is `NA`.
10. `GetAnalysisLedger(object)$analysis_story$raw_values_included` is `FALSE` always, and no
    full assay/expression values appear anywhere in
    `utils::capture.output(str(GetAnalysisLedger(object)$analysis_story))` for an object with a
    real counts matrix.
11. Confirm `GetAnalysisLedger(object)$fingerprint` changes when a new analysis is recorded
    between two calls (sanity check that `analysis_story` participates in the existing
    fingerprint without needing special-casing).
12. Confirm existing `test-ai-diagnostics.R`, `test-ai-execution.R`, `test-ai-claims.R`,
    `test-ai-profile.R` test counts are unaffected (same pass count or higher, zero new
    failures) — these are the files most likely to incidentally depend on `GetAnalysisLedger()`'s
    exact shape.

## 6. Documentation updates required

- `man/GetAnalysisLedger.Rd` (if hand-written and exported, which it is — confirm with
  `glob man/*.Rd` first): add `analysis_story` to the `@return` description, matching the
  existing bullet-point style for `analyses`/`state_records`/etc.
- `NEWS.md`: one bullet, same style as the P0.1 entry, describing the new field and the two
  stale-field fixes (`blocked = NA` / `execution = NA`), explicit that this is additive and that
  nothing in the codebase previously branched on the old `FALSE`/`TRUE` values (so this is not a
  behavior-breaking change for any existing caller).
- `.dev/ai-advanced-analysis-spec.md`: update §6 ("AI Context 扩展") to describe
  `analysis_story` as implemented, and update the new P0 responsibility table's
  `GetAnalysisLedger()` row description to mention it now also exposes `analysis_story`.
- Do not touch `.dev/ai-product-roadmap.md` — leave roadmap phase bookkeeping to the
  conversation layer, not the implementation commit (same convention as P0.1).

## 7. Acceptance criteria (binary, checked at review)

- [ ] `git diff --check` clean, no non-ASCII in `R/ai-context.R` or the new test file.
- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, **0 warnings**, pass count
      **greater than or equal to** `1073` (confirmed baseline before this change, measured
      2026-10-07).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`.
- [ ] Only `R/ai-context.R`, its man page, `NEWS.md`, `.dev/ai-advanced-analysis-spec.md`, and
      the new test file are modified. No action registry, no `R/ai-execution.R`,
      no `R/ai-planning.R` changes. `git diff --stat` must show exactly these files.
- [ ] `analysis_story` never contains a full assay/expression matrix, a per-cell vector, or a raw
      metadata value from `colData`/`rowData` (only column *names*, record *keys*, and already-
      aggregated flags/notes already produced elsewhere).
- [ ] `evidence_gaps`/`conflicts` never call a mutating action, never call `RecordAIEvidence()`,
      never change `execution` state. The whole function is read-only, verified by confirming
      `identical(object, <same object after calling GetAnalysisLedger>)` in at least one test.
- [ ] All temporary scratch/debug files removed before final report; `git status --short` clean
      except for the intended changed files.

## 8. Suggested implementation order (for the delegated subagent)

1. Read `R/ai-context.R` in full (already read for this spec; re-read before coding).
2. Read `R/ai-evidence.R` to confirm the exact field name/shape of a `user_decision` node.
3. Read `R/ai-design-confirmation.R` to confirm the exact shape of an `ai_design_confirmation`
   state record (what key holds the confirmed design list).
4. Read `check_rare_cell_readiness()` in `R/ai-diagnostics.R` to confirm the exact
   `doublet_evidence_available` flag shape and its accompanying note string (per section 3.4,
   `evidence_gaps` is scoped to this one case only in this spec).
5. Implement `sclet_ai_build_analysis_story()` per section 3.
6. Wire it into `sclet_ai_build_context()` per section 4.
7. Apply the two stale-field fixes in section 4.1.
8. Write the test file from section 5.
9. Run `devtools::test(filter = "ai-")` in the background (it takes ~60s; a foreground call will
   time out), iterate until green.
10. Run `make check` in the background (it takes several minutes), fix any new notes/warnings.
11. Update documentation per section 6.
12. Report: paste the actual test and check output, list every file changed, and self-verify
    every item in section 7 before declaring done. Do not commit; leave changes in the worktree.
