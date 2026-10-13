# Spec: name `success_criteria` and `stop_conditions` in the structured-output schema

Status: ready for implementation
Baseline: `3cd9af6` (all `test-ai-*.R` pass; `make check` 0/0/0)

## 1. Problem, measured

**The premise of this spec is a correction.** `.dev/spec-ai-execution-benchmark.md`
section 6 records:

> The model never uses `success_criteria`. It declared zero success criteria in
> every executed run, so `success_assessment` was always empty and acceptance
> criterion 5 in section 5 could not be exercised at all.

That claim does not survive contact with the recorded data. Reading both logs
directly (`sclet-benchmark-ai/results/*.jsonl`, 2026-10-13):

| log | goal given to the model | records | `n_success_criteria` |
|---|---|---|---|
| `execute_plan.jsonl` | "Cluster the cells and then find marker genes for each cluster." | 16 runnable | **0 in all 16** |
| `plan_success_criteria.jsonl` | "Cluster the cells and **report which clusters are supported by marker genes**." | 7 | 4, 5, 5, 5, 5, 5, 8 |

All seven of the second log's plans are valid, all were executed, all were fully
assessed (`n_assessed == n_success_criteria`, `criteria_fully_assessed = TRUE`),
and most criteria came back `met`. The machinery works end to end on real
evidence; the model uses it **when the goal asks for something that can be
accepted or rejected**.

So the gap is not "the model never uses the feature" and not a broken
implementation. It is a vocabulary gap in the output contract, of exactly the
kind measured for `evidence_refs` and fixed in `3cd9af6`:

- `sclet_ai_native_result_schema()` (`R/ai-adapter.R`) declares
  `success_criteria = z_array(z_any())` and `stop_conditions =
  z_array(z_any())`. A model that decides mid-response to express "the clusters
  should be supported by markers" as a criterion has **no declared shape to
  write into**, so it falls back to prose in the rationale.
- The planning prompt (`R/ai-planning.R`, the `success_criteria` block) does
  carry the full vocabulary -- field names, the `source` rules, the
  `check` enum. That text reaches the tool-loop stage. It does not reach the
  **structured-summary stage**, which is where `proposed_actions` and every
  other top-level field is actually produced: the plan is built from
  `result$proposed_actions` / `result$success_criteria` after
  `sclet_ai_generate_object()`.
- The two failure modes differ in an important way, and this spec must not
  repeat the mistake of assuming they are the same. `evidence_refs` was absent
  from the context entirely; `success_criteria` is *described* in the prompt
  but unshaped in the schema. The measured 0/16 above is therefore a
  goal-sensitivity finding as much as a schema finding, and any fix must be
  measured against both a criterion-friendly goal and a bare pipeline goal.

## 2. Requirements

1. The structured-output schema must declare the real shape of both
   `success_criteria` and `stop_conditions`: `id`, `description`, `source`
   (enum: `action_output`, `evidence`), `step_id`, `evidence_id`, `check`
   (enum: `exists`, `equals`, `gte`, `lte`, `in`), `field`, `value`. All
   nullable, `.additional_properties = TRUE` retained so task-specific extras
   survive.
2. The two enums must share one helper with the validator, the way
   `sclet_ai_finding_claim_levels()` / `sclet_ai_finding_severities()` do, so
   the schema cannot accept a `source` or `check` that `ValidateAIPlan()`
   would reject.
3. The structured-summary prompt (and the tool-loop prompt, one line) must
   state when a criterion is worth proposing: when the plan produces something
   that can be verified afterwards, a criterion naming it should be included;
   a plan that genuinely produces nothing verifiable may omit both lists.
4. Nothing may **require** criteria. A plan with zero criteria stays valid --
   that is P0.1's rule and this spec does not change it. The instruction is
   permissive, never mandatory.
5. No change to `ValidateAIPlan()`, `new_sclet_ai_plan()`,
   `sclet_ai_assess_success_criteria()`, or any planning prompt wording beyond
   what requirement 3 adds.

## 3. Acceptance criteria

1. `sclet_ai_native_result_schema()` declares `success_criteria` and
   `stop_conditions` as arrays of objects whose `properties` include all eight
   field names above, with `source` and `check` as enums matching the values
   `sclet_ai_normalize_success_criteria()` / `ValidateAIPlan()` accept.
2. A round trip through the schema: a mock structured response carrying one
   well-formed criterion reaches `new_sclet_ai_plan()` with `id`, `source`,
   `step_id`, `check`, `field`, `value` intact, and survives
   `ValidateAIPlan()` (existing `test-ai-plan-criteria-wiring.R` covers the
   same shape through the mock path -- this criterion is about the schema
   actually naming the fields).
3. A criterion-shaped mock response is **rejected at the schema level** when it
   carries an out-of-vocabulary value, e.g. `check = "matches"` or
   `source = "analysis_output"`, rather than being silently accepted and then
   failing validation later. (This is what `z_any()` loses today.)
4. `make check` 0/0/0 and all existing `test-ai-*.R` files still pass,
   including `test-ai-plan-criteria-wiring.R` and `test-ai-privacy.R`.
5. Offline: a hand-built criterion object passed through
   `sclet_ai_normalize_success_criteria()` is accepted by `ValidateAIPlan()`
   for every `check` in the enum, and rejected for `check` outside it --
   proving the shared vocabulary is the same one the validator uses.

## 4. What this spec deliberately does not do

- It does not require the model to emit criteria (requirement 4), and it does
  not try to make a bare "cluster and find markers" goal produce them. The
  measured 0/16 is a legitimate model judgement about that goal, not a defect
  to engineer away.
- It does not re-run the benchmark. The measurement here is the correction of
  the recorded claim plus the schema-level tests; an online re-measurement
  (both goal shapes, via `sclet-benchmark-ai`) is a separate, billed step and
  should be decided on its own.
- It does not touch `stop_conditions` semantics beyond giving them the same
  declared shape they already have in the normalizer.

## 5. Risk

Narrowing `z_any()` to a declared object can only reject output the local
validator would have rejected anyway; the risk is a provider refusing a shape
it previously accepted. Mitigation is the same as `3cd9af6`: all fields
nullable, additional properties retained, and the full existing suite plus the
structured-summary mock tests run before landing.
