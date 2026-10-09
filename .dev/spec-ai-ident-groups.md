# Spec: Expose identity group codes to planning (and resolve them back)

Status: ready for implementation
Baseline: `86707db` (all 22 `test-ai-*.R` files pass; `make check` 0/0/0)

## 1. Problem (grounded in two real API runs)

The plan-repair benchmark (`sclet-benchmark-ai`, commit `0a8c365`) ran
`AIPlanAnalysis()` twice against DeepSeek on a 60-cell / 200-gene
`scuttle::mockSCE()` fixture whose active identity is a `cluster` column with
real labels `treat1` (28 cells) and `treat2` (32 cells). Neither run reached
`valid = TRUE`, and **both failed on the same thing**:

- run 1: `prerequisites not met for action step_7: ident.1 contains unknown
  group labels: 0` — the model guessed the label `"0"`.
- run 2: `... ident.1 contains unknown group labels: <cluster labels from
  step_6>` — the model emitted a literal placeholder.

Cause, verified by building the real planning context and inspecting it:
the context exposes the identity column **name** but never any label
**value**. Specifically (`sclet_ai_context(sce, target = "planning", detail =
"full")` on the fixture above):

- `dataset.columns` contains `"cluster"` (the column name, among other
  colData column names)
- `active_view$ident` is `"cluster"`
- the action schema lists `ident.1`/`ident.2` as type `character_vector`
- **no field anywhere contains `treat1`/`treat2`, or even the number of
  groups**

Confirmed separately that this is not a privacy-scanner stripping problem:
`sclet_ai_sanitize_context()` redacts **zero** fields matching
`ident|group|cluster|label` on that context, because there is no label value
in it to strip. `run_de_test`'s `ident.1` therefore cannot be filled
reliably, and the only currently-valid alternative the model has is to omit
`ident.1` entirely (the `FindAllMarkers` path, `R/markers.R:179`).

This is the same class of gap as the action-vocabulary gap fixed in
`de4e1c1`: the model is asked to produce a value from a vocabulary it was
never shown.

## 2. Design

### 2.1 One deterministic authority for group codes

`group_N` codes are **already an established convention** in this codebase,
but there is no shared helper and the three existing sites disagree on
ordering:

| site | ordering used |
|---|---|
| `R/ai-diagnostics.R:43-49` | `levels <- unique(values)` — order of first appearance |
| `R/ai-execution.R:710` (`run_annotation` handler) | `u <- unique(as.character(idents))` — order of first appearance |
| `R/ai-privacy.R:281-287` (k-anonymization) | `table(...)` — name-sorted |

Add one helper as the authority for the new surface:

```r
sclet_ai_ident_group_codes <- function(object)
```

Returns, for the object's active identity, a deterministic mapping computed
with **`sort(unique(as.character(labels)))`** — i.e. name-sorted, matching
the existing privacy-scanner convention rather than the `unique()`
appearance-order convention. Sorted order is chosen because it is stable
regardless of cell ordering in the object, which is the property that makes
a code reproducible across two independently-computed contexts (the context
builder and the resolver).

Return shape:

```r
list(
    col = "cluster",                    # identity column name, or NULL
    codes = c("group_1", "group_2"),    # group_N in sorted-label order
    labels = c("treat1", "treat2"),     # real values, NEVER placed in a context
    sizes = c(28L, 32L),
    n_groups = 2L
)
```

Returns a zero-row/`NULL`-valued result (not an error) when there is no
active identity, or when `Idents(object)` yields no values.

**Corrected during implementation — do not reintroduce the original
mistake.** An earlier draft of this spec also said the helper should return
an empty result when the active identity is `colLabels`. That is wrong:
`FindClusters()` sets `ActiveIdent` to exactly `"colLabels"`, and it is the
most common way an object acquires an identity. `Idents()` already resolves
the `colLabels` case (`R/clustering.R`), so a `colLabels`-identified object
has a perfectly usable identity and must produce a mapping like any other.
Bailing out on `colLabels` made group codes work only on objects with an
explicit colData ident column, and silently fail on the standard clustering
path — caught only because the test fixture for the execution scenario was
switched to a `FindClusters()`-produced object. Use `Idents(object)` and let
it handle the source; do not branch on the string `"colLabels"`.

Deliberately **not** done here: reconciling the three pre-existing sites to
use this helper. Their outputs may be asserted by existing tests and their
semantics differ subtly (`ai-diagnostics.R` encodes a per-cell vector; the
privacy scanner's mapping is internal to sanitization). That reconciliation
is listed in section 4 as a follow-up, not part of this change.

### 2.2 Expose level codes (never real labels), nested under `capabilities`

**The shape below was determined empirically, not by reading the
classifier's regex.** A first draft exposed a new top-level `ident_groups`
field. Running it through `sclet_ai_sanitize_context()` showed that shape is
broken: `ident_groups.groups.code` (the one value the model actually needs),
`ident_groups.col`, `ident_groups.n_groups` and
`ident_groups.all_groups_available` were **all redacted**, leaving only the
`n_cells` counts. A second attempt that renamed the top-level field to
`ident_levels` was redacted **entirely**.

Root cause, established by direct probes:
- The privacy classifier's "aggregate escape hatch"
  (`R/ai-privacy.R:122-128`) applies only at `depth >= 2L` **and** only when
  the *parent* was classified `aggregate` or `container`.
- A **new top-level** field (depth 1) therefore never gets the hatch: it must
  itself match the allowlist regex or it is dropped wholesale — which is why
  `ident_levels` vanished.
- A container field whose *name* matches the `group`/`cluster` rule is
  classified `"group"`, which is neither `aggregate` nor `container`, so the
  hatch stops applying to its children — which is why
  `ident_groups.groups.code` vanished.

The verified shape therefore **nests under the existing, already-allowlisted
`capabilities` container** and avoids `group`/`cluster` in the field names
that matter:

```r
# inside the list returned by sclet_ai_capabilities()
identity_levels = list(
    column = "cluster",
    n_levels = 2L,
    aggregation_threshold = 10L,
    codes = list(
        list(code = "group_1", n_cells = 28L),
        list(code = "group_2", n_cells = 32L)
    ),
    omitted_below_threshold = 0L,
    all_levels_available = TRUE,
    raw_values_included = FALSE
)
```

**Verified by running exactly this structure through
`sclet_ai_sanitize_context(privacy = "standard")`: zero fields redacted, the
entire subtree survived intact including every `code`, and no real label
appeared.** Re-run that check during implementation and paste the result. If
the shape is changed at all, re-verify — this behaviour is not predictable
from the regex.

Rules:
- **Only `group_N` codes are exposed. A real label value must never appear.**
- Levels with fewer cells than the codebase's established aggregation
  threshold of `10` are omitted from `codes` and counted in
  `omitted_below_threshold`; `all_levels_available` is `FALSE` in that case.
  This mirrors the scanner's own `counts < 10L` rule (`R/ai-privacy.R:282`)
  rather than inventing a threshold. Omitting (not redacting) is deliberate:
  a redacted entry would still tell the model that a level exists and invite
  it to reference a code it cannot use.
- `n_levels` reports the count of **exposed** levels, so the model cannot
  infer a hidden level's existence from a count mismatch.
- `sclet_ai_capabilities()` currently takes `(status, records)`. It will need
  the object (or the pre-computed mapping) to build this; thread it in the
  smallest way that works rather than restructuring the function. Do not
  change its existing `execution = NA` field.

### 2.3 Resolve codes back to real labels at use sites

Add:

```r
sclet_ai_resolve_ident_value <- function(object, value)
```

- If `value` is `NULL`/empty → return it unchanged.
- If every element matches `^group_[0-9]+$` → map each via
  `sclet_ai_ident_group_codes(object)` to the real label. An out-of-range
  code (e.g. `group_9` when only 2 groups exist) is **left unchanged**, so
  the existing validation reports its normal "unknown group labels" error
  rather than a confusing new one.
- Any other value (already a real label, or arbitrary text) is returned
  unchanged, preserving today's behaviour for hand-written plans and every
  existing test.

Apply it at both places where the value is consumed, because params flow to
validation and execution separately and both must agree:

1. `run_de_test`'s `prerequisites` (`R/ai-execution.R`) — resolve
   `ident.1`/`ident.2` before comparing against the real group table, so a
   plan using `group_1` validates.
2. `run_de_test`'s `handler` — resolve before calling `RunDEtest()`, since
   `RunDEtest()`/`FindMarkers()` match `ident.1` against the **real**
   identity values (`R/markers.R:178-183`). Skipping this would validate
   cleanly and then fail at execution, which is worse than the current state.
3. `run_trajectory`'s `prerequisites` **and** `handler` for `start_cluster`
   (the analogous "must name a group value" parameter). Note `run_trajectory`'s
   `group` parameter is a *column name*, not a value, and needs no
   resolution — do not touch it.

### 2.4 Prompt

Extend the `AIPlanAnalysis()` system prompt (and, for consistency,
`AIRecommendNextStep()`'s, since it shares the planning target) to state
that group-valued parameters (`ident.1`, `ident.2`, `start_cluster`) must be
taken from `context$capabilities$identity_levels$codes` (each entry's
`code` field), and that when
`all_levels_available` is `FALSE` or no code fits, the parameter should be
omitted where the action supports omission (`run_de_test` without `ident.1`
runs the all-cluster marker path) rather than guessed.

## 3. Non-goals

- **Never expose a real label value in a context.** Only `group_N` codes.
  This is the whole point of the design; a solution that puts `treat1` into
  the prompt is wrong even if it makes the benchmark pass.
- No change to `R/ai-privacy.R`. The existing `"group"` classification and
  its k=10 rule are reused as-is. If the empirical sanitize check (section
  2.2) shows the new field is damaged, fix it by reshaping the field
  (rename/renest), not by editing the privacy classifier.
- No change to `RunDEtest()`, `FindMarkers()`, `RunSlingshot()`, or any
  other deterministic analysis function (`R/markers.R`, `R/trajectory.R`).
  Resolution happens in the AI-facing layer only.
- No reconciliation of the three pre-existing `group_N` sites
  (`R/ai-diagnostics.R:43`, `R/ai-execution.R:710`, `R/ai-privacy.R:281`) to
  the new helper. Documented as a follow-up in section 4 with the reason.
- No exposure of group values for *any* colData column other than the active
  identity. Exposing arbitrary column value sets is a much larger privacy
  surface and is not needed for the observed failures.
- No change to `ValidateAIPlan()`'s checks or strictness. Resolution makes
  more plans *valid* by supplying a real value the model asked for by code;
  it does not relax any check.

## 4. Required tests (`tests/testthat/test-ai-ident-groups.R`)

Use `scuttle::mockSCE()` with a `cluster` column and `ActiveIdent(sce) <-
"cluster"` for the "has groups" cases (matching the benchmark fixture: two
groups, both ≥10 cells). Mock `options(sclet.ai.call = ...)` for anything
that would otherwise reach a provider; no test in this file may make a real
API call.

1. **Helper mapping is sorted and deterministic.** `sclet_ai_ident_group_codes()`
   returns `codes == c("group_1","group_2")` and `labels` in name-sorted
   order, with `sizes` matching `table(Idents(sce))`. Then **reorder the
   columns of the object** (`sce[, rev(seq_len(ncol(sce)))`) and assert the
   mapping is byte-identical — this is the property that makes the code
   usable across two independently-built contexts.
2. **No active identity is not an error.** On an object with no identity set,
   the helper returns a result with `n_groups == 0L` (or `NULL` codes) and
   does not throw; `identity_levels` in the context is correspondingly empty.
3. **Context exposes codes and counts, and no real label.** Assert
   `context$capabilities$identity_levels$codes` contains `group_1`/`group_2`
   with `n_cells` 28/32, and assert that **no serialized form of the context**
   (e.g. `paste(utils::capture.output(str(context)), collapse=" ")`) contains
   the strings `treat1` or `treat2`. The latter is the load-bearing
   assertion; do not omit it.
4. **Sanitized context is intact and still label-free.** Run the context
   through `sclet_ai_sanitize_context()` (`privacy = "standard"`), then
   assert (a) `payload$capabilities$identity_levels` has the same codes and
   counts as the unsanitized one, and (b) the serialized sanitized payload
   still contains neither `treat1` nor `treat2`. Paste the sanitized
   `identity_levels` subtree and the matching `policy$redacted_fields` in the
   report.
5. **Small levels are omitted, not redacted.** Construct a fixture with one
   level below the aggregation threshold of 10 (e.g. 32 + 4 cells). Assert
   the small level is absent from `codes`, `omitted_below_threshold == 1L`,
   `all_levels_available == FALSE`, and `n_levels` equals the number of
   *exposed* levels (not the total).
6. **Resolver maps codes, passes everything else through.**
   `sclet_ai_resolve_ident_value(sce, "group_1")` returns the real label for
   group 1; a vector `c("group_1","group_2")` maps element-wise; an
   already-real label is returned unchanged; `"group_9"` (out of range) is
   returned unchanged; `NULL` and `character(0)` are returned unchanged.
7. **A plan using a code passes validation.** Build a plan with a
   `run_de_test` step carrying `ident.1 = "group_1"` and assert
   `ValidateAIPlan(plan, object = sce, registry = ...)$valid` is `TRUE`,
   where the same plan with `ident.1 = "treat1"` also validates (backward
   compatibility for hand-written plans).
8. **Execution resolves too, not just validation.** Run that plan through
   `ExecuteAIPlan(..., dry_run = FALSE)` and assert the resulting evidence
   node id is the one built from the **real** label (`ev:detest_treat1`),
   proving the handler resolved the code before calling `RunDEtest()`. A
   handler that skipped resolution would validate fine and then fail (or
   produce a differently-named node), which is exactly the failure this
   scenario exists to catch.
9. **Out-of-range code still produces the normal error.** A plan with
   `ident.1 = "group_99"` must fail validation with the ordinary
   "unknown group labels" error, not a new/unclear one.
10. **`start_cluster` resolves the same way.** A plan with a `run_trajectory`
    step carrying `start_cluster = "group_1"` validates against the fixture
    (skip if `slingshot` is unavailable, matching the existing trajectory
    tests' guard), asserting the same resolver is wired there.
11. **Full-suite regression guard is external.** Do **not** nest a
    whole-suite run inside this file (see `.dev/spec-ai-plan-repair.md`
    section 4 scenario 9 for why: it was ~20x slower and downgraded
    file-level errors to skips). The suite is run separately.

## 5. Acceptance checklist

- [ ] `R/ai-context.R` (helper + context field) and `R/ai-execution.R`
      (`run_de_test` prerequisites + handler, `run_trajectory` prerequisites +
      handler for `start_cluster`) are the only production files changed, plus
      their `.Rd` files, `NEWS.md`, and this spec. `R/ai-privacy.R`,
      `R/markers.R`, `R/trajectory.R`, `DESCRIPTION` and `NAMESPACE` must be
      untouched.
- [ ] No new exported function (the helper and resolver are internal); if
      either is exported, say so explicitly and justify it.
- [ ] All existing `test-ai-*.R` files still pass.
- [ ] The new file's 11 scenarios pass.
- [ ] `make check` reports `0 errors | 0 warnings | 0 notes`. **Run it — do
      not commit on the test suite alone.** A previous commit in this
      session (`de4e1c1`) skipped `make check` and shipped an Rd
      code/documentation mismatch that only surfaced later; if a signature
      or roxygen block changes, hand-edit the matching `.Rd`.
- [ ] The report pastes: the sanitized `capabilities$identity_levels` subtree,
      the `redacted_fields` entries matching `identity_levels`, and a grep
      proving neither `treat1` nor `treat2` appears in the serialized
      sanitized payload.
- [ ] `git status --short` shows only the files listed in the first item.

## 6. Implementation order

1. Read `R/ai-context.R`'s `sclet_ai_build_context()` and `sclet_ai_capabilities()`
   in full to see how to thread the object (or the pre-computed mapping) into
   `sclet_ai_capabilities()` with the smallest change.
2. Read `R/ai-privacy.R`'s `sclet_ai_payload_atomic()` group branch
   (`counts < 10L` rule) so the exposure threshold matches it exactly rather
   than re-deriving a number.
3. Read `run_de_test`'s prerequisites and handler, and `run_trajectory`'s, in
   full before editing either.
4. Implement `sclet_ai_ident_group_codes()` and `sclet_ai_resolve_ident_value()`.
5. Add `identity_levels` under `capabilities` (section 2.2). The exact shape
   in section 2.2 was already verified against the privacy sanitizer; still
   re-run the check after implementing and paste the result, and re-verify if
   you deviate from that shape at all.
6. Wire resolution into the four points listed in section 2.3.
7. Update the two system prompts (section 2.4).
8. Write the 11 test scenarios; run the new file first, then the full `ai-`
   suite.
9. Run `make check`.
10. Report as required by section 5, and do **not** commit — leave changes in
    the worktree for review.

## 7. Follow-ups this change deliberately does not do

- Reconcile `R/ai-diagnostics.R:43`, `R/ai-execution.R:710` and
  `R/ai-privacy.R:281` onto the new helper. They use two different orderings
  (appearance vs sorted) and their outputs may be asserted by existing tests;
  unifying them is its own change with its own review.
- The "secondary error" recorded in `sclet-benchmark-ai`: a step whose params
  already failed schema validation still runs its prerequisite check, so an
  NA value produces a redundant "missing value where TRUE/FALSE needed"
  alongside the correct error.
- Re-running the benchmark with repeated samples to measure whether the
  repair loop now reaches `valid = TRUE` on a real model, and at what cost.
  That is the measurement that would close out the P5 planning line.
