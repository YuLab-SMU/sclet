---
name: sclet-driver
description: >
  How an external AI agent drives sclet WITHOUT any AI package configured:
  read the shipped analysis skills, write the plan as an R list, validate it
  with ValidateAIPlan(), execute with ExecuteAIPlan() or RunAIPlan(), then
  read the ledger and evidence back and report with calibrated claim
  levels. Covers the exact function signatures, the plan JSON shape, the
  confirmation-token requirement, how to recover from validation errors,
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
