# Spec: expose the variable-feature set as a referenceable plan output

Status: implemented and verified (2026-10-13)
Baseline: `490836e` (all `test-ai-*.R` pass; `make check` 0/0/0)

## 1. Problem, measured

`.dev/spec-ai-execution-benchmark.md` section 6 records the finding, with the
follow-up measurement that settled it:

> On pbmc3k all four planning runs tried to reference step 2's feature set and
> all four were rejected, in four different spellings
> (`${step_2.output.required_hvg}`, `${step_2.output.hvg}`,
> `${step_2.output}`, `${step_2}`). That is the single largest first-attempt
> blocker and it is a capability gap, not a bug: `VariableFeatures(object)`
> already returns the set, so the fix is a new referenceable output field.

Four rejections in four spellings is a model reaching for the same absent
thing four times, which distinguishes it from a one-off mistake. The mechanism
is concrete: `find_variable_features`' `output_schema` is

```r
output_schema = list(
    required_commands = "FindVariableFeatures",
    required_hvg = TRUE
)
```

`required_hvg` is a boolean flag meaning "the object has an HVG set", not the
set itself. `scale_data(features =)` and `run_pca(subset_row =)` are declared
`character_vector`, so the model's prose description of what it wants passed
the schema check and then died in `sclet_ai_gene_selector()`, which has always
rejected anything that is not a row name (and now additionally reads R range
notation). The repair loop never saw the real cause either, because `scale_data`
rejects the prose in `prerequisites` -- that part is already fixed -- but a
model can only be repaired toward something that exists.

Omission works: dropping `features` makes `ScaleData()` scale every gene,
which is correct but expensive on a real object and silent about the fact that
the model wanted a subset.

## 2. What is being decided

Two options, and the choice is not obvious:

**(a) Expose the set as a value.** `output_schema` gains
`variable_features`, and the execution summary carries the real gene names once
they exist. Downstream steps then write
`features = "${step_2.output.variable_features}"` and the existing projection
machinery resolves it. This is the option the earlier spec anticipated, and it
is what the deferred-projection design exists for.

**(b) A dedicated promise.** A distinct literal such as `"$HVG"` that
`sclet_ai_gene_selector()` resolves against `VariableFeatures(object)` at
execution time, with no gene names ever entering the plan, the ledger, or the
context.

This spec implements **(a)**, with one constraint taken from option (b).

The reason is that (b) introduces a second value language beside the existing
`${step.output.field}` references, and the measured failure is precisely that
the model already tries four spellings of one thing. A second syntax would make
that five. Option (a) reuses the machinery the prompt already documents
("reference it with the plan's `${step_id.output.field}` syntax"), so the fix
teaches the model nothing new.

The constraint: **the set is carried in the execution summary, never in the
context**. `available_actions` exposes `output_schema` statically, before any
step has run, and the static declaration must not promise a value it cannot
have yet -- the same rule already enforced by `sclet_ai_referenceable_output_fields()`
refusing to offer SCE fields to a `returns = "value"` action. So
`variable_features` is registered as a referenceable field name whose static
form is the promise, and only the post-execution summary carries the genes.
A model can cite the field; it can never read 2000 gene names out of the
planning context.

Bounding, concretely: a pbmc3k run keeps the default 2000 features and the gene
symbols serialise to roughly 20 KB. That figure is acceptable once (in the
execution record), and is what (b) would also have to pay if it ever wanted to
report which features were scaled. It is not acceptable multiplied across every
downstream step's summary, so the summary must not copy the vector into steps
that did not reference it.

## 3. Requirements

1. `find_variable_features`' `output_schema` declares a `variable_features`
   field. Its static value must be a promise the projection can carry and that
   `sclet_ai_exposable_contract_keys()` will surface -- a character marker, not
   a list, so the sanitizer sees a plain scalar.
2. `sclet_ai_referenceable_output_fields()` must list `variable_features` for
   this action, so `${step_2.output.variable_features}` validates instead of
   being rejected as `value_reference_unresolvable`. The error message that
   lists what a step exposes must therefore include it, which is how the model
   learns the spelling.
3. Projection must resolve the reference to the real vector at execution time,
   through the existing `planned_outputs` path, with no new syntax.
4. `sclet_ai_gene_selector()` must accept an already-resolved character vector
   of real row names unchanged (it does today) and must keep rejecting prose
   (it does today). No change to the selector.
5. The static `output_schema` must not carry the gene set into the planning
   context; only the execution summary does. Assert this in a test that builds
   a context with a registry and greps for a gene name.
6. `run_pca(subset_row =)` must accept the same reference -- it already routes
   through the same selector, so this is a consequence of (2), not new code,
   but it must be tested because it is the second consumer.
7. No change to `ValidateAIPlan()`'s rules, `ScaleData()`, `RunPCA()`,
   `FindVariableFeatures()`, or the planning prompt's reference paragraph
   (which already documents the syntax). The one prompt change allowed: none.

## 4. Acceptance criteria

1. A plan whose `scale_data` step carries
   `features = "${step_2.output.variable_features}"` (with
   `find_variable_features` as `step_2` and a declared `depends_on`)
   validates cleanly, executes, and scales **exactly** the HVG rows: the
   `scaled` layer records `VariableFeatures(object_after)` in its `features`
   param, those rows differ from `logcounts`, and every other row still equals
   `logcounts`. (Note: `ScaleData()` with a subset keeps the full assay and
   writes only the selected rows -- measured, `R/preprocessing.R`. An earlier
   draft of this criterion asserted `nrow(assay) == length(hvg)`, which is not
   what the package does; the criterion above is the real contract.)
2. The same plan with the four spellings the model actually used
   (`required_hvg`, `hvg`, `output`, the bare step id) is still rejected, with
   an error that lists `variable_features` among the step's exposed fields.
   This is what makes the fix discoverable: the rejection teaches the spelling.
3. The planning context (`GetAnalysisLedger()` / `sclet_ai_context()` with a
   registry) contains the field name but no gene name.
4. A `run_pca` step using `subset_row = "${step_2.output.variable_features}"`
   validates and executes the same way.
5. Prose in `features` is still rejected at validation, with the actionable
   message the selector already produces -- the fix must not weaken the
   existing guard.
6. `make check` 0/0/0 and all existing `test-ai-*.R` files pass, including
   `test-ai-velocity-execution.R`, which asserts a stop condition naming
   `required_states`.

## 5. Risk

The static declaration in `output_schema` is also read by
`sclet_ai_assess_success_criteria()` (a criterion may name `variable_features`
as its `field`) and by the privacy classifier on the path
`available_actions.<action>.output_schema.variable_features`. Requirement 5 and
acceptance criterion 3 cover both by asserting on a real context. The
`value` shape is chosen as a scalar character marker rather than a list
precisely so the classifier sees nothing that looks like a label list.

## 6. Two corrections found while implementing

Both were found by running the code, not by reading it, and both changed the
implementation:

**(a) "Validation accepts, execution cannot" reappeared inside this fix.** The
first implementation accepted the deferred reference in `prerequisites` and
left the handler alone, on the theory that the marker was the handler's problem
later. It was not: the projection rewrites
`${step_2.output.variable_features}` into the static marker `"hvg"` *before*
`prerequisites` runs, so the value that reached validation was the literal, and
the handler then handed that literal to `sclet_ai_gene_selector()`, which
rejected it -- a plan that passed validation and died at execution, which is the
exact defect this codebase keeps re-encountering. Two changes fixed it:
`sclet_ai_project_plan_value()` no longer projects a field whose static value is
a deferred promise (the reference is passed through intact), and
`sclet_ai_resolve_deferred_selector()` resolves it against the real object in
the handler. As a consequence the deferred path is now driven purely by the
reference *shape*, so a model that types the literal `"hvg"` is rejected
consistently at both layers instead of being accepted and then failing.

**(b) The obvious end-to-end assertion was wrong about `ScaleData()`.** An
earlier draft checked `nrow(assay(object, "scaled")) == length(hvg)`.
`ScaleData()` with a subset keeps the full assay and writes only the selected
rows (`R/preprocessing.R`: `scaled_assay[features, ] <- scaled`), so that
assertion fails on correct behaviour. The check is now on the layer's recorded
`features` param and on which rows actually changed -- HVG rows differ from
`logcounts`, every other row still equals it.

## 7. Outcome

Implemented on `devel`. New tests: `tests/testthat/test-ai-hvg-reference.R`,
10 scenarios, 65 assertions, confirmed to have teeth by removing the
`variable_features` declaration and observing 9 failures across four scenarios
including both end-to-end ones. Full suite 2732 pass / 0 fail / 78 warn
(unchanged) / 3 skip. `make check` 0 errors / 0 warnings / 0 notes.

Verified end to end, not just at the schema level:

- `scale_data(features = "${step_2.output.variable_features}")` validates
  `valid = TRUE`, executes `completed`, and the `scaled` layer records exactly
  `VariableFeatures(object_after)` while every non-HVG row still equals
  `logcounts`.
- `run_pca(subset_row = "${step_2.output.variable_features}")` validates and
  executes, registering the `PCA` reduction.
- A plan-level criterion with `field = "variable_features"` on the
  `find_variable_features` step assesses as `met` with the real 20 genes
  observed.
- The planning context carries the field *name* and none of the real gene
  names (22 assertions over the serialised registry view).
- The four spellings the model actually used are all still rejected, and the
  `value_reference_unresolvable` message now lists `variable_features` among
  the step's exposed fields, which is how the model learns the right spelling.

## 8. Online measurement (2026-10-13, real provider) -- partial confirmation

Two fresh `plan_success_criteria` and `execute_plan` records on pbmc3k
(`deepseek-chat`), read back directly:

- **The field name is in the model's vocabulary now.** First attempt in the
  `plan_success_criteria` run (2026-10-13 22:33): `step_3 uses
  ${step_2.variable_features}` -- the model reaches for `variable_features`
  itself, missing only the `.output` segment. The historical four spellings
  (`required_hvg`, `hvg`, bare `${step_2.output}`, `${step_2}`) never named the
  field. The `value_reference_malformed` message told it the exact syntax
  ("write `${step_2.output.variable_features}` instead") -- acceptance
  criterion 2's teaching behaviour is working as designed -- but one repair
  round did not fix it and the plan was rejected (`plan_valid = FALSE`,
  `executed = FALSE`).
- **No reference has completed end to end on a real model yet.** The
  `execute_plan` success run (22:30) omitted `features` entirely
  (`scale_used_hvg_subset = FALSE`, ScaleData densified all 32,738 genes), so
  the resolved-promise path (acceptance criterion 1) is proven offline only.

So the fix removed the capability gap the 4/4 measurement showed -- the model
is no longer reaching for something that does not exist -- and the remaining
barrier is syntax precision (`${step.field}` vs `${step.output.field}`), which
the malformed-reference error now spells out verbatim. Whether that error plus
repair converges on a valid referencing plan is a measurement for a larger
repetition sample (this is 1 attempt of it), and the plan-repair loop's
behaviour on malformed references is now a known open thread for the
`spec-ai-plan-repair.md` owners.

## 9. Online update (2026-10-14, real provider)

`execute_plan` re-run x5 under a degraded environment (3 recorded, 2 lost to
provider-side hangs -- see sclet-benchmark-ai `docs/action-coverage.md`,
2026-10-14 finding). The decision-relevant data from the three records:

- The model **writes the correct reference syntax after one repair round**:
  lines 23/24 show `plan_steps` with
  `features=${step_2.output.variable_features}` and
  `subset_row=${step_2.output.variable_features}` (`n_reference_steps = 2`),
  after first attempts used `${step_2.variable_features}`. The teach-on-rejection
  error text is demonstrably read and applied.
- But both repaired plans stayed `invalid_plan` with 1 residual error that no
  field recorded -- a benchmark gap now closed (`remaining_errors`,
  sclet-benchmark-ai `31b48fa`). What the residual was remains unknown.
- The one executed plan (line 22) omitted `features` entirely
  (`scale_used_hvg_subset = FALSE`), so **the end-to-end resolved-promise
  execution on live model output is STILL unmeasured**, after 6 executed
  attempts across two batches. The offline end-to-end (test-ai-hvg-reference)
  remains the only proof of criterion 1.

Net: the capability-gap half of the claim is confirmed on live models (the
model can and does name the field correctly); the usage-completion half is
still pending, now blocked by an unrecorded residual error, whose recorder
exists as of `31b48fa` for the next batch.
