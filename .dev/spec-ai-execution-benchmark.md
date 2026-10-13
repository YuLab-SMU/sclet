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

**Superseded: the default is now real data.** The fixture below was the first
attempt; `mockSCE()` turned out to be unusable and the synthetic stand-in was
still a toy, so both the planning and execution benchmarks now default to
**pbmc3k** (2700 cells x 32738 genes, real clusters and real marker genes),
loaded from a local extracted 10x directory by `fixture_pbmc3k()` in
`R/harness.R`. `fixture =` is an argument on both benchmarks and is recorded on
every result. The reasoning below still explains why a toy fixture cannot
support this pipeline at all.

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

Recorded in `sclet-benchmark-ai`: `results/execute_plan.jsonl` and the
"execution benchmark" finding in `docs/action-coverage.md`. Sixteen trials in
four batches, each batch run against the fixes from the one before it.

What the benchmark established, with the numbers:

- **The execution report never disagreed with the object.** In all 14 executed
  runs, every artifact each step's `output_schema` declared was really present
  on the final object: `missing_artifacts` was empty 14/14. This is the check
  the benchmark exists for, and it passed every time it could run.
- **Dry runs are inert** 14/14 (`recorded = FALSE`, `execution_id = NULL`,
  fingerprint unchanged), and **the ledger carries the execution record**
  14/14 (`type = "ai_execution"`, retrievable by `execution_id`).
- **Executing a validated plan completed** in 9 of 14 executed runs, versus 5
  failures. But the *kind* of failure changed once the fixes landed: batches
  1-3 lost runs to failures discovered mid-execution, whereas batch 4 (after
  all four fixes) had **no mid-execution failure at all** -- 3/4 completed and
  the fourth was rejected at validation, before any work was done.

Three defects were found by running it, all now fixed in `sclet` (see section
3 and the commit list): the cluster-then-markers prerequisite that ignored
`planned`, the unresolvable gene selector that surfaced as an opaque
`subscript out of bounds`, and validation accepting references
(`${step_2.output.required_hvg}`) that execution could never resolve.

Two findings deliberately left open:

- **The model never uses `success_criteria`.** It declared zero success
  criteria in every executed run, so `success_assessment` was always empty and
  acceptance criterion 5 in section 5 could not be exercised at all. The
  feature is implemented (P0.1) and the model does not reach for it; that is a
  separate question from execution and is recorded rather than papered over.
- **No referenceable field exposes the variable-feature set.** Two batch-3
  failures were the model trying to pass "the HVGs from step_2" into
  `scale_data(features=)`. `find_variable_features`' `output_schema` has
  `required_hvg` (a flag, not a feature list), so there is nothing legitimate
  to point at, and the model invents prose instead. Omission works -- the
  handler defaults correctly -- but whether to expose the feature set as a
  referenceable output is a design decision left for later, on more evidence
  than two runs.

  **Resolved as the next step, with 4/4 evidence.** On pbmc3k all four planning
  runs tried to reference step 2's feature set and all four were rejected, in
  four different spellings (`${step_2.output.required_hvg}`,
  `${step_2.output.hvg}`, `${step_2.output}`, `${step_2}`). That is the single
  largest first-attempt blocker and it is a capability gap, not a bug:
  `VariableFeatures(object)` already returns the set, so the fix is a new
  referenceable output field. Because the value only exists after the producing
  step runs, it needs the same deferred-projection treatment as
  `planned$active_ident` -- validation must accept the reference without being
  able to check the value, and leave the execution-time contract as the
  authority. New API surface, so it gets its own spec.

  **Implemented.** `.dev/spec-ai-hvg-reference.md`, landed as
  `find_variable_features` declaring a `variable_features` output field plus
  `sclet_ai_is_deferred_gene_selector()` for the deferred acceptance. Two
  details were decided against the measurement rather than from it: the static
  `output_schema` value is the marker `"hvg"`, not the gene names, because
  `output_schema` reaches `available_actions` before any step runs and 2000
  symbols there would enter every planning context; the real set appears only
  in the execution summary. And `ScaleData()` with a subset keeps the full
  assay, writing only the selected rows (`R/preprocessing.R`), so the
  acceptance check is on the layer's recorded `features` and on which rows
  changed -- an earlier draft asserted the assay was narrowed, which the
  package does not do. All four rejected spellings above are still rejected,
  and their error now names `variable_features`.

The limit that applies regardless: one fixture, one provider. The plan is
always the same seven-action pipeline, so this does not test execution of
`run_integration`, `run_annotation`, `run_trajectory`, `run_velocity`, or the
doublet/rare-cell pair. Those still have no execution benchmark.
