# P5 Spec (Slice 3): Trajectory/Velocity Evidence Cross-Check

> Status: draft for implementation (2026-10-08)
> Parent roadmap: `.dev/ai-product-roadmap.md`, P5, third item in the recommended order:
> "velocity readiness -> velocity execution -> **trajectory / velocity interpretation** ->
> CellRank / fate -> spatial -> multimodal".
> Precondition satisfied: P5 slice 1 (`e4c5e75`) and slice 2 (`fd20cbb`) are committed on
> `devel`.
> Scope decision (direct human, this round): a **read-only cross-check function**, not a
> combined `interpret = TRUE` narrative. Narrower, lower-risk, no `report` layer changes.
> Baseline before this change: `devtools::test(filter = "ai-")` ->
> `[ FAIL 0 | WARN 30 (benign velociraptor/scuttle deprecations) | SKIP 1 | PASS 1274 ]`.
> `make check` -> `0 errors | 0 warnings | 0 notes`.

## 1. Problem

`sclet_ai_record_trajectory_evidence()` (`R/ai-execution.R:1100-1155`, confirmed by direct read)
and `sclet_ai_record_velocity_evidence()` (`R/ai-execution.R:1157-1206`, confirmed by direct
read) each record one bounded evidence node (`ev:trajectory_<name>` /
`ev:velocity_<name>`) when their respective action runs. Nothing today reads both nodes
together. Confirmed by `grep -n 'trajectory.*velocity\|velocity.*trajectory' R/ai-execution.R
R/ai-diagnostics.R R/ai-functions.R` returning only one unrelated match (the
`allowed_groups` registry string listing both names side by side).

**Verified while writing this spec** (every claim below checked by reading the actual code):

1. **Both evidence nodes carry a pseudotime-shaped quantile summary that can be compared
   without touching per-cell data.** Trajectory's node carries `first_lineage_median` /
   `first_lineage_iqr` / `first_lineage_max` (from its first `slingPseudotime_N` column,
   confirmed at `R/ai-execution.R:1139-1143`) plus `start_group_code` /
   `relative_ordering_only = TRUE`. Velocity's node carries a per-column summary keyed by
   column name (confirmed at `R/ai-execution.R:1190-1195`), and if `velocity_pseudotime` was
   written, `values[["velocity_pseudotime"]]$median` / `$quantile_25` / `$quantile_75` /
   `$max` are present with the identical shape trajectory uses
   (`n_observed`/`mean`/`quantile_25`/`median`/`quantile_75`/`max`, confirmed identical field
   names at `R/ai-execution.R:1119-1126` vs `1170-1177`).

2. **Both pseudotime values are explicitly documented as relative, not absolute, orderings**:
   trajectory's `relative_ordering_only = TRUE` / `pseudotime_is_absolute_time = FALSE`
   (`R/ai-execution.R:1134-1135`), velocity's `velocity_is_model_estimate = TRUE`
   (`R/ai-execution.R:1187`). A cross-check can therefore only compare **relative** shape
   (e.g. spread via IQR, relative magnitude of median vs max) between the two bounded
   summaries -- it must not claim either pseudotime is the "correct" one, and must not attempt
   per-cell rank correlation (that would require reading raw colData columns together, which
   is explicitly out of scope -- see non-goals).

3. **The existing read-only comparator `compare_rare_cell_evidence()`
   (`R/ai-comparison.R:320-436`, confirmed by direct read in full) is the right structural
   model**: returns a typed `status = "not_available"` result with a `reason` code when fewer
   than two usable inputs exist (here: trajectory evidence missing, velocity evidence missing,
   or both), otherwise `status = "available"` with anonymized/bounded fields, every nested list
   carrying `raw_values_included = FALSE`, and a `matching`/`independence`-style block
   documenting exactly what the comparison does and does not establish. The new function
   should match this convention, not invent a new one.

4. **Export/doc convention confirmed**: `compare_rare_cell_evidence()` uses a roxygen `@export`
   block (`R/ai-comparison.R:306-319`), a hand-written `man/compare_rare_cell_evidence.Rd`, and
   a manual one-line addition to `NAMESPACE` (`export(compare_rare_cell_evidence)`) --
   confirmed this package does not run `devtools::document()` to regenerate these (per the
   standing project constraint already followed in slice 1). The new function must follow the
   same three-file pattern by hand.

5. **No existing include-group or action registry entry is needed.** This is a free function
   like `compare_rare_cell_evidence()`, not an `AIAction`. It takes a `SingleCellExperiment`
   object directly and reads its ledger; it does not go through `ValidateAIPlan()`/
   `ExecuteAIPlan()`. Confirmed `compare_rare_cell_evidence()` itself is never registered as an
   `AIAction` (`grep -n 'compare_rare_cell_evidence' R/ai-execution.R` returns nothing).

## 2. Non-goals (do not implement these here)

- No `interpret = TRUE` / report-layer change. The direct human explicitly chose the narrower
  read-only comparator over a combined narrative this round.
- No new `AIAction`, no registry change, no `include` group change. This is a free function
  taking `object` directly, exactly like `compare_rare_cell_evidence()`.
- No per-cell rank correlation (e.g. Spearman between raw `slingPseudotime_1` and raw
  `velocity_pseudotime` columns). That would require reading two raw colData columns together
  rather than two bounded evidence summaries, which is a materially different privacy posture
  from every other comparator in `R/ai-comparison.R`. If a future round wants real per-cell
  correlation, it needs its own privacy review and is explicitly out of scope here.
- No claim about which pseudotime (trajectory's or velocity's) is "more correct". The output
  may only describe whether their bounded shape summaries are consistent or inconsistent, with
  a `claim_level` no stronger than `"consistent_with"` (matching both source evidence nodes'
  own `claim_level`, confirmed identical at `R/ai-execution.R:1149` and `1200`).
- No change to `sclet_ai_record_trajectory_evidence()` or `sclet_ai_record_velocity_evidence()`
  themselves unless a genuine defect is found while writing this function (same discipline as
  every prior round).
- Do not touch `R/ai-privacy.R`. The new function only ever reads already-recorded, already
  bounded evidence `values` lists (via `sclet_ai_evidence_get_all()`); it never reads raw
  colData directly, so it introduces no new privacy surface.

## 3. Required behavior: `compare_trajectory_velocity_evidence()`

New file `R/ai-trajectory-velocity-crosscheck.R` (new file rather than appending to
`R/ai-comparison.R`, to keep this slice's diff self-contained and easy to review/revert; model
its internal style on `compare_rare_cell_evidence()` but do not literally copy its rare-cell-
specific matching logic).

```r
compare_trajectory_velocity_evidence <- function(object, trajectory_id = NULL, velocity_id = NULL)
```

- `object`: a `SingleCellExperiment`. Validate with the same
  `if (!inherits(object, "SingleCellExperiment")) stop(..., call. = FALSE)` pattern
  `compare_rare_cell_evidence()` uses (`R/ai-comparison.R:321-323`).
- `trajectory_id` / `velocity_id`: optional analysis-run name strings (matching the `name`
  parameter both `run_trajectory` and `run_velocity` accept). When `NULL`, use whichever
  trajectory/velocity evidence node is present if there is exactly one of each; if more than
  one of either exists and no id was given, return `status = "ambiguous_run"` naming the
  available ids rather than silently picking one (do not guess).
- Read both evidence nodes via `sclet_ai_evidence_get_all(object)`, looking for
  `id == paste0("ev:trajectory_", trajectory_id)` and `id == paste0("ev:velocity_", velocity_id)`
  (or scanning for ids matching `^ev:trajectory_` / `^ev:velocity_` when no id was given).
- **Not-available cases** (mirror `compare_rare_cell_evidence()`'s typed-reason convention,
  `R/ai-comparison.R:334-349`):
  - Neither evidence node present: `status = "not_available"`,
    `reason = "no_trajectory_or_velocity_evidence_records"`.
  - Only trajectory present: `status = "not_available"`,
    `reason = "velocity_evidence_missing"`.
  - Only velocity present: `status = "not_available"`,
    `reason = "trajectory_evidence_missing"`.
  - Velocity evidence present but has no `velocity_pseudotime` column summarized (i.e.
    `values$velocity_pseudotime` is `NULL` -- this happens if `mode` never wrote that column,
    confirmed possible since `velocity_coldata` in `sclet_ai_record_velocity_evidence()` is an
    `intersect()` against whatever actually exists, `R/ai-execution.R:1160-1163`):
    `status = "not_available"`, `reason = "velocity_pseudotime_column_not_available"`.
  - Trajectory evidence present but has no lineage summarized (i.e.
    `values$first_lineage_median` is `NULL`, confirmed possible when `pseudotime_columns` was
    empty, `R/ai-execution.R:1138` guards this with `if (length(lineage_summaries))`):
    `status = "not_available"`, `reason = "trajectory_lineage_not_available"`.
- **Available case**: build a bounded comparison using only the two nodes' existing quantile
  summaries (`trajectory`'s `first_lineage_median`/`first_lineage_iqr`/`first_lineage_max`;
  `velocity`'s `values$velocity_pseudotime$median`/`$quantile_25`/`$quantile_75`/`$max`).
  Compute, on these bounded aggregates only (never on raw columns):
  - `trajectory_relative_spread = first_lineage_iqr / first_lineage_max` (guard divide-by-zero:
    if `first_lineage_max` is `0` or non-finite, set this field to `NA_real_` rather than
    erroring).
  - `velocity_relative_spread` computed the same way from velocity's quantiles.
  - `spread_ratio = velocity_relative_spread / trajectory_relative_spread` (same divide-by-zero
    guard), and a boolean `comparable_spread = is.finite(spread_ratio)`.
  - `spread_consistency_label`: one of three fixed strings depending on whether `spread_ratio` (when
    finite) falls inside `[0.5, 2.0]` ("relative spread of the two pseudotime estimates is of
    similar order of magnitude"), outside that band ("relative spread of the two pseudotime
    estimates differs by more than 2x; this may reflect genuinely different dynamics captured
    by each method, not necessarily an error"), or `comparable_spread == FALSE`
    ("relative spread could not be compared because one or both bounded summaries were
    degenerate"). Pick the `[0.5, 2.0]` band because it is a simple, symmetric, pre-registered
    threshold rather than one tuned after seeing data -- do not adjust this band based on any
    particular dataset's result.
  - `claim_level = "consistent_with"` (never anything stronger).
  - `raw_values_included = FALSE` on the top-level result and nowhere does the function return
    anything resembling per-cell values.
  - Returned list shape:
    ```r
    list(
        status = "available",
        trajectory_run = trajectory_id_actually_used,
        velocity_run = velocity_id_actually_used,
        trajectory_relative_spread = ...,
        velocity_relative_spread = ...,
        spread_ratio = ...,
        comparable_spread = ...,
        spread_consistency_label = ...,
        claim_level = "consistent_with",
        caveats = list(
            "Both pseudotime values are relative orderings, not absolute time, and are model-
             or algorithm-dependent (slingshot lineage ordering vs scVelo-estimated dynamics).",
            "This comparison only checks whether the bounded spread summaries are of similar
             magnitude; it does not establish per-cell agreement, directional (root-to-tip)
             agreement, or that either method is biologically correct."
        ),
        raw_values_included = FALSE
    )
    ```
    (Exact caveat wording may be adjusted for clarity while implementing, but must preserve
    both substantive points: relative-not-absolute, and spread-only-not-per-cell.)

## 4. Required test scenarios (`tests/testthat/test-ai-trajectory-velocity-crosscheck.R`)

Reuse fixture-construction helpers from `tests/testthat/test-ai-velocity.R` and the trajectory
tests (grep for the existing trajectory fixture helper, e.g. in
`tests/testthat/test-ai-execution.R` or a dedicated trajectory test file, and reuse rather than
reinvent). `skip_if_not_installed("velociraptor")` and `skip_if_not_installed("slingshot")` (or
whatever package `RunTrajectory()` actually depends on -- check `R/trajectory.R`'s own
`requireNamespace()` call first) on every scenario that executes real actions.

1. **Neither evidence present**: fresh object, no trajectory or velocity run yet. Assert
   `status == "not_available"`, `reason == "no_trajectory_or_velocity_evidence_records"`.
2. **Only trajectory evidence present**: run `run_trajectory` only (through the real action or
   `RunTrajectory()` + `sclet_ai_record_trajectory_evidence()` directly, whichever is simpler to
   set up correctly). Assert `status == "not_available"`,
   `reason == "velocity_evidence_missing"`.
3. **Only velocity evidence present**: symmetric case. Assert
   `reason == "velocity_evidence_missing"` would be wrong here -- assert
   `reason == "trajectory_evidence_missing"`.
4. **Both present, comparable**: run both `run_trajectory` and `run_velocity` for real on the
   same object (reuse/adapt whichever fixture already has both spliced/unspliced assays and a
   cluster column trajectory needs -- check whether `tests/testthat/test-ai-velocity.R`'s
   fixture already has a cluster column; if not, add one real `kmeans()`-style cluster
   assignment the way other trajectory tests do). Assert `status == "available"`,
   `claim_level == "consistent_with"`, `raw_values_included == FALSE`, and that
   `trajectory_relative_spread`/`velocity_relative_spread`/`spread_ratio` are all finite
   numeric scalars (do not assert a specific `spread_consistency_label` string value on real data --
   assert only that it is one of the three allowed fixed strings from section 3, since the
   actual spread ratio on a tiny synthetic fixture is not something to hand-predict).
5. **Ambiguous run selection**: run `run_velocity` twice under two different `name`s on the
   same object (with one `run_trajectory`). Call
   `compare_trajectory_velocity_evidence(object)` with no `velocity_id` given. Assert
   `status == "ambiguous_run"` and that both velocity run names are listed somewhere in the
   result (field name your choice, document it in the function's own roxygen).
6. **Explicit ids resolve the ambiguity**: same two-velocity-run object as scenario 5, but call
   with an explicit `velocity_id` matching one of the two runs. Assert `status == "available"`.
7. **Degenerate spread (divide-by-zero guard)**: construct a case where `first_lineage_max`
   or velocity's `quantile_75$max`-equivalent is `0` (e.g. a trivial fixture where all
   pseudotime values collapse to the same value -- check whether this is actually
   constructible without a contrived mock; if genuinely impractical to construct for real,
   it is acceptable to call the internal helper function directly with a hand-built evidence
   list for this one scenario only, clearly commented as doing so). Assert
   `comparable_spread == FALSE` and `spread_consistency_label` matches the degenerate-case string.
8. **No raw per-cell values anywhere in the output**: for the "available" scenario (4), assert
   `length(result) `is reasonably bounded (e.g. no field has `length() > 20`, mirroring the
   bounded-aggregate assertions already used in `tests/testthat/test-ai-velocity.R`'s own
   privacy-regression scenario) and that none of the result's field *names* match the
   `sclet_ai_evidence_bad_field()` pattern from `R/ai-evidence.R:1-7` (reuse that exact
   function via `sclet:::sclet_ai_evidence_bad_field` to check field names programmatically
   rather than eyeballing).
9. **Privacy regression guard**: with both evidence nodes present (reuse scenario 4's object),
   call `AIStatus(object)` and assert zero privacy/redaction/consent-related warnings fire,
   matching the pattern used in `tests/testthat/test-ai-velocity.R`'s own scenario 13.

## 5. Documentation updates

- New `man/compare_trajectory_velocity_evidence.Rd`, hand-written following the exact
  structure of `man/compare_rare_cell_evidence.Rd` (read that file first).
- One new line in `NAMESPACE`: `export(compare_trajectory_velocity_evidence)`, inserted
  alphabetically among the existing `export(compare_*)` lines if any, otherwise wherever the
  existing `export(...)` block's ordering convention places it (check the file's actual
  ordering before deciding; do not assume strict alphabetical if the existing file is not
  strictly alphabetical).
- `NEWS.md`: one new top entry describing this as a read-only comparator (not an `interpret`
  change), what it does and does not establish, following the precision level of every prior
  round's entry.
- `.dev/ai-advanced-analysis-spec.md`: update the top status line, following the existing
  running-sentence convention, to note the trajectory/velocity cross-check comparator is done
  as a read-only function (not combined interpretation), and that trajectory/velocity combined
  `interpret` narrative remains open (it was deliberately not chosen this round).
- Do not edit `.dev/ai-product-roadmap.md`.

## 6. Acceptance checklist (self-verify before reporting done)

- [ ] `devtools::test(filter = "ai-")` reports **0 failures**, only the same benign
      `velociraptor`/`scuttle` deprecation warnings as the current baseline (no new
      privacy/redaction/consent warnings), with a pass count greater than the 1274 baseline by
      at least 9 (one per scenario in section 4).
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`.
- [ ] `git status --short` shows changes only to: the new `R/ai-trajectory-velocity-
      crosscheck.R`, the new test file, `NAMESPACE` (one line), `man/compare_trajectory_
      velocity_evidence.Rd` (new), `NEWS.md`, `.dev/ai-advanced-analysis-spec.md`. No other
      `R/` file, no `DESCRIPTION` change, no `R/ai-privacy.R` change.
- [ ] `grep -c 'AIAction(' R/ai-execution.R` is unchanged (still 18) -- this function is not an
      action.
- [ ] Every field name in the function's return value passes
      `!sclet:::sclet_ai_evidence_bad_field(name)` (confirmed in section 4 scenario 8; also
      self-verify this directly against the real implementation's actual field names, not just
      the names proposed in section 3, in case implementation needs to deviate).
- [ ] The function never calls `RecordAIEvidence()` itself (it is read-only; it does not write
      a new evidence node) -- confirm by reading the implementation.
- [ ] No file outside the authorized list above was touched.

## 7. Implementation order

1. Read this entire spec once before writing any code.
2. Read `R/ai-comparison.R:306-436` (`compare_rare_cell_evidence()` plus its roxygen header) in
   full to copy its structural idioms exactly.
3. Read `tests/testthat/test-ai-velocity.R` and whichever existing file has the trajectory
   fixture helper (grep for `run_trajectory` across `tests/testthat/*.R` to find it) to reuse
   fixture-construction code rather than inventing new helpers.
4. Read `man/compare_rare_cell_evidence.Rd` to copy its `.Rd` structure exactly for the new man
   page.
5. Implement `compare_trajectory_velocity_evidence()` in the new
   `R/ai-trajectory-velocity-crosscheck.R` per section 3. If any detail in section 3 does not
   match what the real evidence node shapes actually contain once you re-read
   `sclet_ai_record_trajectory_evidence()`/`sclet_ai_record_velocity_evidence()` yourself,
   adjust the implementation to match reality and document the deviation precisely in your
   final report (same discipline as every prior round -- do not silently guess).
6. Add the roxygen `@export` block, hand-write `man/compare_trajectory_velocity_evidence.Rd`,
   add the one `NAMESPACE` line.
7. Write the test file (section 4), 9 scenarios.
8. Run `devtools::test(filter = "ai-", reporter = "progress")` in the background (expect
   60-120 seconds; scenarios 4-6 run real `run_trajectory`/`run_velocity` actions), iterate
   until green.
9. Run `make check` in the background.
10. Update documentation per section 5.
11. Self-verify every item in section 6, listing a concrete yes/no and evidence for each.
12. Report: paste the actual terminal output of the test run and `make check`, list every file
    changed with `git status --short` and `git diff --stat`, and do not commit -- leave the
    changes in the worktree for review.




