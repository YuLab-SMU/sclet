---
name: sclet-driver
description: >
  How an external AI agent drives sclet WITHOUT any AI package configured:
  read the shipped analysis skills, write the plan as an R list, validate it
  with ValidateAIPlan(), execute with ExecuteAIPlan() or RunAIPlan(), then
  read the ledger and evidence back and report with calibrated claim
  levels. Covers the exact function signatures, the plan JSON shape, the
  confirmation-token requirement, how to recover from validation errors,
  how to advance an analysis in rounds (the engine executes one plan and
  never re-plans mid-run, so each round is a fresh plan the caller starts),
  and what must never be done (executing unvalidated, inventing action
  names, claiming results without ledger evidence). Use when the agent runs
  sclet from R on the user's behalf, when sclet.ai.call / aisdk is absent
  or deliberately unused, or when asked to analyse a SingleCellExperiment
  object through sclet's registered actions.
---

# Driving sclet as an external agent (no aisdk required)

sclet's planning model (`AIPlanAnalysis`) needs the `aisdk` package and an
API key, because it asks a provider to write the plan. **Everything else is
plain R.** If you are the agent, YOU are the planner: read these skills,
write the plan yourself, and let sclet validate and execute it under the
same guardrails a built-in model would face. The package must be installed
(`library(sclet)`); the plan functions do not touch the network.

## The contract, in one paragraph

You propose a plan as an R list of steps; `ValidateAIPlan()` checks it
against the registry of actions actually available on THIS object (names,
parameter schemas, prerequisites, dependency declarations) and, when valid,
returns a confirmation token; `ExecuteAIPlan()` runs it with that token,
records every step in the analysis ledger, registers evidence nodes, and
returns the mutated object plus a per-step report; you then report results
at the claim level the recorded evidence supports. The engine is the
authority on what is allowed; these skills are how you get it right on the
first attempt.

## The workflow

```r
library(sclet)

# 1. What can run on THIS object? The registry is object-aware: actions
#    whose prerequisites the object cannot meet are listed with reasons.
sce                                   # your SingleCellExperiment
reg <- AIDefaultExecutionRegistry(sce, include = "all")

# 2. Read the shipped skills (this file is one of them) and plan against
#    them -- stage order, producers before consumers, parameters that are
#    the user's modelling choice.
skills <- AISkills()                  # c("sclet-claims", "sclet-provenance", ...)
guidelines <- ReadAISkill("sclet-analysis")$body

# 3. Write the plan. Every action name must be a verbatim registry key;
#    every parameter a real value of the declared type; every reference
#    ${step_id.output.field} backed by depends_on.
plan <- new_sclet_ai_plan(
    task = "pbmc_markers",
    context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
    actions = list(
        list(id = "s1", action = "normalize_data"),
        list(id = "s2", action = "find_variable_features",
             params = list(nfeatures = 2000L), depends_on = "s1"),
        list(id = "s3", action = "scale_data", depends_on = "s2"),
        list(id = "s4", action = "run_pca",
             params = list(ncomponents = 50L), depends_on = "s3"),
        list(id = "s5", action = "find_neighbors",
             params = list(dims = 1:30), depends_on = "s4"),
        list(id = "s6", action = "find_clusters",
             params = list(resolution = 0.8), depends_on = "s5"),
        list(id = "s7", action = "run_de_test",
             params = list(name = "markers"), depends_on = "s6")
    )
)

# 4. Validate. Read the errors -- they name the received value and the fix.
validation <- ValidateAIPlan(plan, object = sce, registry = reg)
if (!isTRUE(validation$valid)) {
    # Fix the plan in R from validation$errors, then re-validate.
    # Do not execute without a valid validation object.
    stop(paste(validation$errors, collapse = "\n"))
}

# 5. Execute. The confirmation token comes from the validation object --
#    executing without it is refused by design (dry_run first is the safe
#    default and costs nothing).
executed <- ExecuteAIPlan(sce, plan, registry = reg,
    validation = validation, dry_run = FALSE,
    confirmation = validation$confirmation_token)

# 6. Read the object and the ledger, then report at the evidence's level.
executed$status                       # "completed" / "completed_with_errors"
executed$results                      # per-step status and output summaries
GetAnalysisLedger(executed$object)    # the full audit trail
```

`RunAIPlan(sce, plan, reg, dry_run = FALSE, confirm = TRUE)` collapses
steps 4-5 (it validates internally); use it for the common case and
`ExecuteAIPlan` when you need the validation object separately. A dry run
(`dry_run = TRUE`, the default) returns what WOULD run, changing nothing.

## What the engine does for you (do not reimplement, do not bypass)

- **Validation refuses what cannot work**: unregistered action names,
  parameters of the wrong shape, prerequisites the object cannot meet,
  dependencies that were not declared, references to fields no step
  exposes. The error text names the received value and, for malformed
  references, the corrected spelling.
- **The plan is fingerprint-bound**: the validation matches the exact plan
  content and object state, so a plan edited after validation is rejected
  (`sclet_ai_stale_plan`) rather than executed against stale assumptions.
- **Every executed step is recorded** in the ledger with its parameters,
  status and output summary; evidence-producing actions register evidence
  nodes with claim levels and parent lineage.
- **Empty plans are refused** (`plan_empty`): a plan with zero steps
  performs no analysis. If the goal genuinely needs no computation, answer
  directly instead of proposing an empty plan.

## Resuming across processes (no resident session required)

The ledger lives ON the object, so the whole analysis state survives
`saveRDS()`. An agent does not need a resident R session: each interaction
can be a fresh `Rscript` that reads the object, plans the next stage against
the fingerprint it carries, executes, and saves it back. Verified end to
end: phase 1 (pipeline + markers) in one process, `saveRDS`, phase 2
(enrichment) in a fresh process reading the object back -- fingerprint
identical, records readable, ceiling computable.

Two rules when resuming:

- **Reference records, not past steps.** A plan may only reference steps of
  ITS OWN plan; `${old_plan_step.output.detest_id}` across a save/load
  boundary is refused (`depends on unknown step(s)`). Read what is on the
  object (`GetAnalysisLedger(sce)`, or the state records) and name the
  record id literally: `detest_id = "markers"`.
- **Re-validate every time.** The fingerprint binds the plan to the object
  state it was validated against, which is exactly what makes the
  stateless pattern safe: an object edited between processes cannot accept
  a plan written against the old state.

A resident session (`callr::r_session`, or an R process the agent keeps
talking to) saves fixture load time and lets the agent inspect objects
interactively, but it is an optimization, not a requirement.

## Advancing in rounds: starting the next plan is yours

Nothing in the engine re-plans. `ExecuteAIPlan()` walks `plan$actions` in
order and stops at the first failure that is not `continue_on_error`; a
failed step is retried with the SAME action and the SAME parameters, and no
step rewrites the ones after it. A non-trivial analysis is therefore a
SEQUENCE of rounds, and starting each round is your decision:

```r
# one round: plan against this state, validate, execute
ledger <- GetAnalysisLedger(sce)   # the records to plan from, and the fingerprint to bind
plan <- ...                 # written for THIS state, or AIPlanAnalysis(sce, goal = goal)
validation <- ValidateAIPlan(plan, object = sce, registry = reg)
if (!isTRUE(validation$valid)) stop(paste(validation$errors, collapse = "\n"))
executed <- ExecuteAIPlan(sce, plan, registry = reg,
    validation = validation, dry_run = FALSE,
    confirmation = validation$confirmation_token)
sce <- executed$object      # the object the next round plans against
```

`RunAIAnalysis(sce, goal, confirm = "yes")$object` is one such round when
aisdk is configured: it plans, validates, previews, asks, executes and
optionally interprets in one call -- running it again is how you continue.
There is no loop option to look for.

Start a new round when:

1. **A stage finished and the next stage feeds on it.** Plan the next stage
   from the ledger you now hold rather than committing a whole pipeline as
   one plan; the later parameters stop being guesses that way.
2. **A step failed for a reason that changes the route** -- a prerequisite
   gate, a parameter the data cannot satisfy, a handler rejecting a value.
   Retrying the same plan replays the same failure. Read the error, plan the
   corrected route (different threshold, the missing upstream step, a
   smaller scope), or tell the user what the object lacks.
3. **The plan never covered the goal.** `plan$metadata$goal_coverage$missing`
   (and `RunAIAnalysis()$report$goal_coverage`) names goal items with no
   step behind them. `AIPlanAnalysis()` spends one round on that before
   returning; a hand-written plan gets no such help -- check it yourself.
4. **What came back decides whether the next step is worth doing** -- QC
   flags, an unexpected cluster count, zero DE genes past threshold, a batch
   effect you did not expect. That is what the round boundary is for: the
   next parameter is a function of the result. (A plan's `success_criteria`
   are evaluated into `success_assessment`, which is observational -- it
   never changes `status`, so judge the round from the results and the
   ledger, not from that list.)
5. **The user answered a question.** A recorded decision appears in
   `analysis_story$user_decisions`, so the next plan can act on it. Never
   invent the answer to an open choice (which reference, which condition,
   which contrast) just to avoid a round trip -- name the choice you need.

Two properties of the round boundary are worth using rather than fighting:

- **An old plan is not reusable after a round.** Validation binds a plan to
  the ledger fingerprint it was written against, and recording a round's
  results changes that fingerprint, so re-validating the previous plan fails
  with `plan context fingerprint does not match the current object`. That is
  the feature, not an obstacle: the plan you checked described a different
  object than the one you now hold.
- **Cross-round inputs are record ids, not `${step.output}`.** A reference
  reaches only the steps of its own plan; carry a value from an earlier
  round by naming the record that produced it (`detest_id = "markers"`),
  which the ledger in front of you supplies literally.

One thing not to attempt: a plan has no conditional steps. Its actions run
unconditionally in order, so "if this succeeds, also do X" is a second
round, not a longer plan.

### End every round by rendering the report (and read its self-check)

A round is not finished when execution returns. It is finished when the
report has been rendered and its self-check has agreed with the ledger:

```r
executed <- ExecuteAIPlan(sce, plan, registry = reg, validation = validation,
    dry_run = FALSE, confirmation = validation$confirmation_token)
sce <- executed$object
report <- RenderAIReport(sce, id = <the record this round produced>)
```

`RenderAIReport()` is not decoration. It re-derives every recorded claim
level from the evidence and re-resolves every citation, and it prints any
disagreement in its own verification section. **A round that is not
rendered loses that check** -- and the next round then plans from a raw
ledger, which grows with every round; a planner reading an ever-longer
ledger is exactly the shape in which citation coverage falls off.

The next round is planned from **this report plus
`analysis_story$data_support`**, not from a fresh reading of everything:

- the report says what this round established, at which claim level, with
  which limitations;
- `data_support$unexplored` says which available analyses the ledger still
  has no result for;
- `data_support$blocked` says what the data cannot answer and which gate
  blocks it.

When the report's self-check FAILS (unresolved citations, or a claim level
that does not re-derive), do not carry that round's findings into the next
plan as settled facts. Report the discrepancy, or fix the record with
`RecordAIEvidence(replace = TRUE)` if the numbers are known to be right.
`RunAIAnalysis()` now does this for you: it renders by default
(`render_report = TRUE`) and sets `report$verification_status` plus
`report$continue_blocked` when the check fails.

### Reading evidence back: two shapes

Evidence nodes are reachable in two forms, and they are not interchangeable:

- `GetAnalysisLedger(object, include_artifacts = TRUE)$state_records$ai_evidence[[id]]$summary`
  -- the ledger form, wrapped in `summary`;
- `sclet_ai_evidence_get_all(object)[[id]]` (internal) -- the bare node, with
  `values` at the top level.

Picking the wrong one yields `NULL` silently, which reads like "no evidence"
instead of "wrong accessor" (one report was generated with empty numbers that
way before the verification pass caught it). Prefer the ledger form, and treat
a `NULL` from the other one as a shape mismatch, not an absent node.

Also: the ledger drops `artifacts` unless `include_artifacts = TRUE`, so a
recorded result's `findings` live there.

## Recovering from validation errors

Validation errors are the repair loop, one round, done by you:

1. `value_reference_malformed: ... uses ${step_4.reduction}; write
   ${step_4.output.reduction} instead` -- fix the spelling to all three
   segments, and make sure the producing step is in `depends_on`.
2. `parameter X has invalid type: expected integer_vector, got the text
   "1:20"` -- pass `1:20` as an integer vector, not a string. Never pass a
   type label or placeholder as a value.
3. `prerequisites not met for action step_N: <gate>_missing ...` -- the
   gate is real. Either reorder the plan so an earlier step produces what
   is missing, drop the step, or tell the user what is missing. Do not
   plan around a gate with a substitute analysis.
4. `unknown parameter(s)` -- check the action's `input_schema`; do not
   pass parameters an action does not declare (e.g. `run_integration`
   takes no `dims`).
5. `plan_empty: ...` -- you returned zero actions. A zero-step plan is
   never the way to record "nothing to do". When the goal names N kinds of
   work, write at least one action per kind; when part of the goal is
   genuinely out of scope for the data, keep the steps you CAN plan and say
   in `success_criteria` (or the rationale) which part you are not doing and
   why. Measured 2026-10-19: an empty plan for a 43k-cell goal was refused,
   and re-proposing without that instruction produced the same shape again.
6. `action is not registered: <name>` -- `action` must be a verbatim
   registry key. Writing `"not_run"`, `"none"` or `"skip"` to drop a step is
   the wrong repair: remove the step from `actions`, and remove its id from
   any other step's `depends_on` and from every `${<id>.output.<field>}`
   reference and success criterion that named it.
7. `depends on unknown step(s): step_N` -- a step is depended on that is not
   in the plan. There are exactly two honest repairs: delete that id from
   `depends_on` (plus any reference or criterion that names it), or put the
   removed step back into `actions`. That error appears when a step was
   dropped without cleaning its dependents -- measured 2 of 9 trials on the
   2026-10-19 Kang repeats. Do not leave a dependency pointing at a step
   that does not exist.

Every message that ends in a `fix:` clause is telling you what to change;
apply exactly that change and re-validate, rather than rewriting the plan
from scratch.

## Reporting: claims are graded by the evidence, not by confidence

After execution, read the evidence nodes on the returned object and phrase
findings at the level their citations support -- see the `sclet-claims`
skill for the ladder and the single-cell fabrication traps. The
`sclet-provenance` skill covers what the ledger already records and what
your report must state (parameters, thresholds, design choices).

Two things you must never do:

- **Claim a result the ledger does not carry.** A step that failed, or a
  plan that was never executed, is not evidence of anything.
- **Execute without validation.** The confirmation token is not red tape;
  it is the check that what runs is what was checked.

## What still needs aisdk

Only `AIPlanAnalysis()` (provider-written plans), the interpretive
entry points that call a model (`AIStatus`, `AIReviewQC`,
`AIExplainAnalysis`, `AskAI`, ...), and online tests. If the user later
configures aisdk, the same plans you wrote by hand become the baseline the
provider-written plans are compared against.
