# Spec: end-to-end execution benchmark (does execution do what the plan says?)

Status: implemented (see section 6 for the outcome)
Baseline: `61f908b` (all `test-ai-*.R` pass; `make check` 0/0/0)

## 1. Problem

Every benchmark in `sclet-benchmark-ai` so far stops at `ValidateAIPlan()`.
That measures whether a real model can *propose* a plan the validator accepts.
It says nothing about the second half of the loop: whether executing an
accepted plan actually produces the artifacts the plan declared, and whether
the execution report can be trusted about that.

The gap is concrete, not theoretical. `ExecuteAIPlan()` returns a
self-assessment (`status`, per-step `status`, `success_assessment`,
`recorded`), and the standing lesson from this benchmark work is that a
self-report must be verified against the artifact rather than believed. So the
benchmark must check the *object*, not the report:

1. Does every step report `completed`?
2. For each completed step, is the output its own `output_schema` declares
   actually present on the final object? (assay / reduction / graph / HVG /
   identity)
3. Does the ledger record the execution, and is `recorded` truthful?
4. Does a dry run change nothing?
5. Does `success_assessment` agree with the observed object?

## 2. Fixture

`mockSCE()` does not work for this benchmark, and finding that out was part of
the work. Executing the standard pipeline on it fails at the last step:

```
step_1..step_6  completed
step_7 run_de_test  failed  "Must have at least 2 groups defined."
```

because `FindClusters()` finds **one** cluster in 60 cells of near-uniform
`mockSCE()` counts (`colLabels` has the single level `1`). That measures the
fixture, not the engine. A benchmark whose last step cannot succeed on any
input teaches nothing.

`R/harness.R` therefore gains `fixture_two_population()`: 60 cells / 100 genes
with two real populations (genes 1-20 up in cells 1-30, genes 21-40 up in cells
31-60, Poisson counts). Verified offline to yield exactly 2 clusters through
the real pipeline, with `RunDEtest()` succeeding afterwards.

Crucially the fixture is **raw counts only** -- no `ActiveIdent`, no
preprocessing. That is the realistic starting state and it exercises the full
plan from scratch.

## 3. sclet-side fix required first

The raw-counts fixture exposed a real limitation, found while probing:

`run_de_test`'s prerequisite read `ActiveIdent(object)` and ignored the
`planned` argument entirely. Since `find_clusters` declares
`output_schema$active_ident = "colLabels"`, a plan of *cluster -> find markers*
was rejected at validation time on any object that did not already have an
identity. The ordinary request "cluster the cells and find markers for each
cluster" was therefore **unplannable from a raw object** -- it failed with
"cluster identity is not set", before any model was involved.

The identity a later step will see is not knowable before the producing step
runs: neither its group sizes nor its label set. So the prerequisite now
defers when `planned$active_ident` is set and differs from the object's current
identity, returning `TRUE` and leaving the producing step's output contract
(`required_states`, `active_ident`) and the execution-time check as the
authority. Deferral is not a relaxation of verification, it is a move of the
check to the point where the fact exists.

Verified both directions: the deferral applies for cluster-then-markers on a
raw object and for re-clustering an already-clustered object (where the
object's *current* groups would be the wrong ones to check), while an object
with a single identity group and no planned change is still refused.

## 4. What the benchmark records

`R/actions/execute_plan.R`, one real provider call sequence:

1. `AIPlanAnalysis()` on the raw fixture, `max_repair_rounds = 1`.
2. If the plan is not valid, record `executed = FALSE, reason =
   "invalid_plan"` and stop. Executing an invalid plan is not a measurement.
3. **Dry run**: assert `status == "dry_run"`, `recorded == FALSE`,
   `execution_id` is `NULL`, and the returned object is unchanged
   (fingerprint equal). A dry run that mutates is a contract violation.
4. **Real run** with `confirmation = validation$confirmation_token` and
   `record = TRUE`.
5. Independently derive, from each action's own `output_schema` in the
   registry, the artifacts that step promised; then check the **final
   object** for each one.

Recorded fields: `status`, `dry_run_status`, `dry_run_mutated`,
`n_steps`, `step_statuses`, `n_completed`, `declared_artifacts`,
`missing_artifacts`, `recorded`, `execution_id_present`, `ledger_has_record`,
`success_assessment_keys`, `evidence_ids`, `n_clusters`, `elapsed_secs`, and
`failure_errors`.

## 5. Acceptance criteria

1. `make check` is 0 errors / 0 warnings / 0 notes.
2. The new deferral tests pass, and the existing suite does not regress.
3. The benchmark runs against the real provider and records a result; a
   `missing_artifacts` list that is non-empty is a *finding*, not a fixture
   problem, and must be reported as such.
4. `docs/action-coverage.md` gains the `execute_plan` row and the finding.

## 6. Outcome

PENDING -- to be filled in from `results/execute_plan.jsonl` after the real run,
and only with what that record shows. Written before the run, deliberately
left empty rather than guessed.

The limit that will apply regardless of the result: one plan on one fixture
across one provider establishes that the execution half of the loop works end
to end and is worth verifying, not that it is robust across action types.
Actions not in this plan (`run_integration`, `run_annotation`, `run_trajectory`,
`run_velocity`, the doublet/rare-cell pair) would still have no execution
benchmark.
