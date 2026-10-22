# AI-layer audit — sclet, `devel`

Independent audit of the optional AI layer of `sclet`, run on the live `devel`
checkout on **2026-10-22**. Nothing in `R/`, `man/`, `tests/`, `NAMESPACE`,
`DESCRIPTION` or `.dev/` was modified; this file is the only artefact written.

**Method.** Every claim below was produced by running code against the working
tree (`devtools::load_all()`), not by reading it. The scripts are listed in
Appendix A. The design contract read as the yardstick is `.dev/ai-product-roadmap.md`,
`.dev/spec-ai-evidence-authorship-v2.md`, `.dev/spec-ai-derived-claim-levels.md`
(there is no `spec-ai-claim-levels.md`; the file is `spec-ai-derived-claim-levels.md`)
and `.dev/spec-m12-plan-quality-and-reliability.md`, plus the project's own
`AGENTS.md`, which records a number of defects this audit does **not** repeat.

**Deliberately not done.** No billed provider call was made by this audit.
`aisdk::get_model()` resolves to `deepseek:deepseek-chat` on this machine, so
the provider path is *configured*; every number quoted for provider behaviour
below is taken from the recorded benchmarks in `../sclet-benchmark-ai/results/*.jsonl`
and is labelled with the date and commit it was measured at.

---

## 0. Headline

| Question | Answer |
|---|---|
| How many places must you touch to add one action? | **21** (8 mechanical, 13 judgement) — §1 |
| Is the validation token flow airtight? | **No.** Two breaks: the registry is not bound to the token, and a token minted without an object binds to no object at all. The plan-to-object binding also collapses for any two objects with the same analysis state. §2 |
| Is the privacy boundary correct? | **The rule is right; the implementation leaks.** The clustering — the exact thing that must go out only as `group_N` — is published verbatim in `capabilities$contrast_levels` on the ordinary `FindClusters()` path, and the engine then accepts a cluster-vs-cluster contrast naming those labels. §3 |
| What is production-ready? | 4 subsystems of 15. Most of the layer is *works*; the provider-driven half is *works*, not *production-ready*, on 12–28 measured trials. §4 |
| Is the multi-round story implemented? | **Yes**, and correctly: one plan = one round, a round changes the fingerprint, the previous plan is refused, round 2 works. The gap is ergonomic, not semantic. §5 |

---

## 1. The real cost of adding one action

The registry is split across 11 `R/ai-registry-*.R` files (3,147 lines) and
assembled by `AIDefaultExecutionRegistry()` at `R/ai-execution.R:582-611`.
The live registry holds **34 actions** (measured, not grepped):

```r
Rscript -e 'devtools::load_all("."); length(AIDefaultExecutionRegistry(
  scuttle::mockSCE(ncells = 20, ngenes = 50), include = "all"))'   # 34
```

`NEWS.md:617` says `run_decontx` took the registry "33 → 34", so the most
recent addition is the best worked example. Here is every place that addition
had to touch, and what each one costs in judgement.

| # | What has to change | Where (file:line) | Mechanical? |
|---|---|---|---|
| 1 | The `AIAction()` entry: `name`, `description`, `handler`, `input_schema`, `prerequisites`, `output_schema`, `allowed_state_types`, `estimated_cost`, `idempotent` | `R/ai-registry-<group>.R` (11 files; e.g. `R/ai-registry-preprocess.R:133-236` for `run_decontx`) | **judgement** |
| 2 | Group added to the catalog allowlist — *only if the action opens a new group* | `R/ai-execution.R:587` (`allowed_groups`) | mechanical |
| 3 | Group factory called in the assembly | `R/ai-execution.R:599-609` | mechanical |
| 4 | Group added to the "which group holds this action" resolver — *new group only* | `R/ai-execution.R:626-643` (`sclet_ai_action_groups`) | mechanical |
| 5 | The roxygen `@param include` list must name the group | `R/ai-execution.R:571-575` | mechanical |
| 6 | The underlying sclet function must write the analysis **state record** the action claims | e.g. `sclet_set_analysis_state()` in the domain file | **judgement** |
| 7 | The action's `allowed_state_types` must name that state type, or execution fails | `R/ai-execution.R:1001-1010` (`sclet_ai_action_state_problems`) | mechanical |
| 8 | The **state-record id rule**, so a later step *in the same plan* can see the record before it exists | `R/ai-planning.R:259-298` (`sclet_ai_planned_state_records`, a `switch()` on the action name) | **judgement** |
| 9 | If the action produces a referenceable id, the deferred-promise marker list | `R/ai-planning.R:216-220` (`sclet_ai_is_deferred_promise`) | **judgement** |
| 10 | If the action produces a referenceable *output field*, the execution summary must observe it — today this is hard-coded per action name (`find_variable_features`, `run_de_test`, `run_cnv_inference`) | `R/ai-execution.R:1114-1143` | **judgement** |
| 11 | **Evidence recording**: a bounded evidence node + the recorder call. Values must be a flat named list of non-string scalars (the privacy layer rejects free text) | inside the handler; `sclet_ai_record_evidence_for_run()` / `RecordAIEvidence()` | **judgement** |
| 12 | **Report**: the action name must survive the "unmentioned suggestion" stem match | `R/ai-report.R:134-146` (automatic for `run_*` names) | mechanical |
| 13 | **Prompt surface A** — the goal→action keyword map | `R/ai-planning.R:1570-1596` (`sclet_ai_goal_action_map`) | **judgement** |
| 14 | **Prompt surface B** — the "unexplored" data-support suggestion that proposes it | `R/ai-data-support.R:440-500` | **judgement** |
| 15 | **Prompt surface C** — the data-support capability class the suggestion depends on (`supported("ambient_cleaning")` etc.) | `R/ai-data-support.R:306-321` | **judgement** |
| 16 | **Prompt surface D** — the round-progression remedy table (which warning the action answers) | `R/ai-data-support.R:604-626` | **judgement** |
| 17 | **Skills** — the action table and the ordering rules | `inst/skills/sclet-analysis/SKILL.md`, `inst/skills/sclet-driver/SKILL.md` | **judgement** |
| 18 | **Tests** — registration, schema, prerequisites, evidence, execution, state contract | `tests/testthat/test-ai-<name>.R` | **judgement** |
| 19 | **Benchmark coverage table** — a row | `../sclet-benchmark-ai/docs/action-coverage.md` | mechanical |
| 20 | **Benchmark script + a results record** (costs provider calls) | `../sclet-benchmark-ai/R/actions/<name>.R` + `results/<name>.jsonl` | **judgement** |
| 21 | **NEWS.md** release note | `NEWS.md:617` (the `run_decontx` entry) | mechanical |

### Answer: **21 places — 8 mechanical, 13 judgement**

The mechanical 8 are #2, #3, #4, #5, #7, #12, #19, #21 — and four of them
(#2–#5) are conditional on the action opening a *new* registry group. The 13
judgement calls are where the real thinking is: the description the model
reads, the input schema it must satisfy, the prerequisite that refuses a bad
plan, the evidence node the report will cite, the id rule that lets a later
step in the same plan see the record, the keyword that makes a goal map to
this action, and the skill text that tells the model how to order it.

### The important structural point

**There is no completeness pin.** The only pin in this area is one-directional:
`tests/testthat/test-ai-round-progression.R:262-281` asserts that every action
named in the *remedy table* exists in the registry. Nothing asserts the
reverse — that every registered action appears in the goal map, the remedy
table, the data-support suggestions, the skills, or the benchmark coverage
table. Measured on the live tree:

| Surface | Actions covered | Of 34, missing |
|---|---|---|
| `sclet_ai_goal_action_map()` (`R/ai-planning.R:1570`) | 13 | 21, incl. `run_decontx`, `compare_integration_routes` |
| `sclet_ai_round_progression_remedies()` (`R/ai-data-support.R:604`) | 7 | 27, incl. `run_milo`, `run_cci`, `run_cnv_inference` |
| `inst/skills/*/SKILL.md` (any) | 20 | 14, incl. `run_decontx`, `run_velocity`, all 6 spatial/multimodal |
| `../sclet-benchmark-ai/docs/action-coverage.md` | 32 | **`run_decontx`, `compare_integration_routes`** |

So an action can be registered, tested, and shipped while being invisible to
the goal matcher, to the round-progression judge, to the skills, and to the
coverage table. `run_decontx` is exactly that action today: it is in the
registry and in the remedy table, and it is in **neither** the goal map, the
skills, nor the coverage table.

### F1 — The coverage table under-reports the registry (a direct rule violation)

`../sclet-benchmark-ai/docs/action-coverage.md:1-10` states, in bold, that the
registry "returns **32** registered actions", re-counted 2026-10-19. The live
registry returns **34**. `run_decontx` and `compare_integration_routes` have no
row at all. `NEWS.md:617` records the 33→34 transition; the coverage table was
never updated. The project's own standing rule — *"Records and coverage update
in the same commit"* (`AGENTS.md`, AI-layer rules) — is broken, and the
consequence the rule was written to prevent is exactly what happened: the
table is now 2 rows short, and the two missing actions are invisible to
whoever plans the next round of work.

*One-line fix:* re-count the registry and add the two missing rows, and add a
test that reads `AIDefaultExecutionRegistry(include = "all")` and fails when
the count in the coverage table differs.

---

## 2. Is the validation-token flow airtight? — **No**

### What the flow is

`ValidateAIPlan()` (`R/ai-planning.R:659-1177`) returns an
`sclet_ai_plan_validation` object (`R/ai-planning.R:1144-1168`) carrying
`valid`, `errors`, `warnings`, `plan_id`, `plan_fingerprint`,
`context_fingerprint`, `requires_confirmation`, `confirmation_token`,
`claim_ceiling`, `plan`, `human_confirmations`, `checked_at`.

`ExecuteAIPlan()` (`R/ai-execution.R:1516-1552`) accepts a plan when
**all** of these hold:

1. `validation$valid` is TRUE (`R/ai-execution.R:1538`);
2. `validation$plan_id == plan$plan_id` (`R/ai-execution.R:1543`);
3. `validation$plan_fingerprint == sclet_ai_plan_fingerprint(plan)` —
   a SHA-256 over the canonicalised plan, minus `$validation`
   (`R/ai-planning.R:144-168`, `R/ai-execution.R:1544`);
4. `validation$context_fingerprint == GetAnalysisLedger(object)$fingerprint`
   (`R/ai-execution.R:1545-1546`) — **skipped entirely when the validation's
   `context_fingerprint` is `NULL`**;
5. the `confirmation` argument equals `validation$confirmation_token`
   (`R/ai-execution.R:1592-1598`).

`plan_fingerprint` binds the token to the **plan's content**. `context_fingerprint`
is the only thing that binds it to **an object**. The **registry is not bound
at all** — the validation object has no registry field (measured: the
validation's names are `valid, errors, warnings, plan_id, plan_fingerprint,
context_fingerprint, requires_confirmation, confirmation_token, claim_ceiling,
plan, human_confirmations, checked_at`).

### The attempts, as runnable R (all executed; full text in Appendix A)

| # | Attack | Result |
|---|---|---|
| A | Validate against object A, execute against object B (two independent `mockSCE()`) | **refused** — `sclet_ai_stale_plan`. The two objects have different fingerprints. |
| B | Validate with `object = NULL`, execute against any object | **BREAKS.** `status = "completed"`. |
| C | Validate against registry R1, execute against registry R2 with an impostor handler under the same action name | **BREAKS.** The impostor handler ran; `status = "completed"`. |
| D | Edit the plan's params after validation | **refused** — `sclet_ai_stale_plan`. |
| E | Swap the action name after validation | **refused** — `sclet_ai_stale_plan`. |
| F | Record a round between validation and execution | **refused** — `sclet_ai_stale_plan`. The fingerprint changed (`sclet-ai-046756e2` → `sclet-ai-2114cfe2`). |
| G | Validate against A, execute against B where B is A with its cells permuted | **BREAKS.** `status = "completed"`, and B's cell 1 carries different counts than A's cell 1. |
| H | Replay the same token twice with `record = FALSE` | Allowed. The token is **not single-use**; it is replayable until the object changes. |

### F2 — A token minted without an object binds to no object at all — **HIGH**

`ValidateAIPlan(plan, registry = reg)` with `object = NULL` is a legal call
(`object = NULL` is the documented default, `R/ai-planning.R:662`). It records
`context_fingerprint = current_fingerprint %||% plan$context_fingerprint %||% NULL`
(`R/ai-planning.R:1161`), i.e. `NULL`. `ExecuteAIPlan()` then skips the whole
object check, because the guard is
`(!is.null(validation$context_fingerprint) && !identical(...))`
(`R/ai-execution.R:1545-1546`).

The validation still warns — `"object was not supplied; object-specific
prerequisites were not checked"` (`R/ai-planning.R:1066`) — and still mints a
**confirmation token**. The token then authorises execution against *any*
object, with *no* prerequisite check having run.

This is not a theoretical path: `ValidateAIPlan()`'s own `@examples`
(`R/ai-planning.R:641-657`) shows the call with `object = sce`, but the
signature makes the object optional, and a caller who validates a plan for
review ("does this plan even parse?") and later executes it on the real object
will find the object check silently absent.

*One-line fix:* in `ExecuteAIPlan()`, refuse when `is.null(validation$context_fingerprint)`
and a non-dry run was requested — or have `ValidateAIPlan()` refuse to mint a
`confirmation_token` when no object was supplied.

### F3 — The validation token does not bind the registry — **HIGH**

`validation` carries no registry fingerprint, and `ExecuteAIPlan()` never
compares the registry it was handed with the registry the plan was validated
against. Demonstrated end to end: a validation minted against the real
registry was executed against a registry whose `normalize_data` action is a
one-line impostor. The impostor's handler ran and the execution reported
`status = "completed"`.

The practical consequence: the confirmation token does not certify *what will
run*, only *which plan text* will run. Any caller that rebuilds the registry
between validation and execution — `AIDefaultExecutionRegistry(object,
include = ...)` with a different `include`, a registry assembled per round, a
user-supplied extension registry — silently changes the code that the human
approved.

*One-line fix:* add a `registry_fingerprint` to the validation (hash of the
descriptors' `name`/`input_schema`/`requires_confirmation`/`mutates_object`)
and compare it in `ExecuteAIPlan()` alongside the other three checks.

### F4 — The object binding is a *state* digest, not an identity — **MEDIUM**

`GetAnalysisLedger(object)$fingerprint` is `sclet_ai_fingerprint(result)`
(`R/ai-context.R:148`), and `sclet_ai_fingerprint()` (`R/ai-context.R:1055-1063`)
is an **additive 31-bit checksum** of the serialised ledger:

```r
weights <- seq_along(raw) %% 100003
checksum <- sum((as.double(raw) + 1) * (weights + 1)) %% 2147483647
```

Two consequences, both measured:

1. **It is not a hash.** It is linear in the serialised bytes, and the weight
   sequence repeats every 100,003 bytes, so a two-byte edit at distance
   100,003 leaves the checksum unchanged. I did **not** construct such a
   collision end to end through the public API and am not claiming one; the
   property is a property of the formula. The space is 2³¹ ≈ 2.1×10⁹.
2. **It does not cover the data.** The ledger's `dataset` summary
   (`R/ai-context.R:696-720`) contains no expression values; the only
   data-derived quantities in the whole fingerprint are two aggregate scalars
   in `analysis_story$data_support$shapes$ambient`. So two objects with the
   same analysis state and the same two scalars share a fingerprint. Measured:
   `A` and `B = A[, sample(ncol(A))]` have the **identical** fingerprint
   (`sclet-ai-045fa69e`), their counts differ, and a plan validated against A
   executes against B with `status = "completed"`.

So the honest statement of what the token guarantees is: *the token binds the
plan to the analysis state, not to the data, and not to the registry.*

*One-line fix:* use `digest::digest()` (already used for the plan fingerprint
at `R/ai-planning.R:161-163`) over the ledger, and hash a registry descriptor
digest into the validation.

### What does work

Attempts D, E, F and H behave correctly. The `plan_fingerprint` is a real
SHA-256 over the canonicalised plan and catches every content edit, including
a swapped action name. The round-to-round refusal (F) is the designed
behaviour and is correct. The token is replayable only while the object is
unchanged (H) — a legitimate property for `record = FALSE` dry-run-style use,
but it should be documented.

---

## 3. The privacy boundary — the rule is right, the implementation leaks

### The rule, and where it is implemented

`AGENTS.md` states it precisely: *a design column's levels may be exposed
because the planner cannot invent them; the clustering (`ActiveIdent`) is the
result the analysis is meant to discover, so it goes out only as `group_N`
codes.*

The code implements the two halves separately:

* **Clustering → `group_N` codes**: `sclet_ai_build_identity_levels()`
  (`R/ai-context.R:939-967`) over `sclet_ai_ident_group_codes()`
  (`R/ai-context.R:969-991`). Measured on a 5-cluster object: `group_1..group_5`
  with `n_cells` each, `raw_values_included = FALSE`. **Correct.**
* **Design levels → verbatim**: `sclet_ai_build_contrast_levels()`
  (`R/ai-context.R:888-936`). **This is where the leak is.**
* The **privacy scanner** (`R/ai-privacy.R:61-265`) is a *name-driven*
  allow/deny/aggregate classifier (`sclet_ai_payload_field_kind()`,
  `R/ai-privacy.R:351-386`). It scans the **context only** — `object` is
  `invisible()` at `R/ai-privacy.R:67` — and the sanitised payload is what goes
  out (`R/ai-adapter.R:440-446`, serialised with `dput()` at
  `R/ai-adapter.R:499`).

### F5 — The clustering is published verbatim as a "design level" on the ordinary `FindClusters()` path — **CRITICAL**

`sclet_ai_build_contrast_levels()` excludes the identity column **by name**:

```r
ident_column <- tryCatch(ActiveIdent(object), error = function(e) NULL)   # R/ai-context.R:900
if (!is.null(ident_column) && !nzchar(ident_column)) ident_column <- NULL
...
if (!is.null(ident_column) && identical(column, ident_column)) next     # R/ai-context.R:904
```

But `FindClusters()` sets the identity to the **string** `"colLabels"`
(`R/clustering.R:130-131`), and `SingleCellExperiment` stores `colLabels` in
the **colData column literally named `label`**. So the comparison is
`"label" == "colLabels"` → **FALSE**, and the exclusion never fires on the
standard path.

Measured on the unmodified pipeline `mockSCE → NormalizeData → RunPCA →
FindNeighbors → FindClusters(resolution = 1.0)`:

```
colData$label (SCE stores colLabels here): 1,2,3,4,5
ActiveIdent(sce): colLabels   <-- the name used for the exclusion
  contrast column=Mutation_Status  levels=negative|positive
  contrast column=Cell_Cycle       levels=G0|G1|G2M|S
  contrast column=Treatment        levels=treat1|treat2
  contrast column=label            levels=1|2|3|4|5      <-- the clustering
```

The privacy scanner does not stop it. The field name is `levels`;
`sclet_ai_payload_field_kind("levels", ...)` returns `"unknown"`
(`R/ai-privacy.R:351-386`), and the "unknown → aggregate" reclassification at
`R/ai-privacy.R:122-128` promotes any character vector of length ≤ 12 to
`"aggregate"`, which `sclet_ai_payload_atomic()` then emits raw through
`sclet_ai_payload_bound_character()` (`R/ai-privacy.R:324-329`, `432-436`).
The sanitised payload at `privacy = "standard"` therefore contains the
clustering verbatim, with **zero redactions** and **zero warnings**.

Worse, the leak is not merely cosmetic — it is *actionable*. The engine
accepts a cluster-vs-cluster contrast that names the real cluster labels:

```
run_de_test(group='label', ident.1='1', ident.2='2')  ->  ACCEPTED
plan-level validation of that step: valid = TRUE
```

and with human-meaningful cluster names (the normal state of any annotated
object) the same path publishes them:

```
contrast column=label  levels=CLUSTER_T_CELL_SIGNATURE | CLUSTER_B_CELL_SIGNATURE | ...
sanitized payload: CLUSTER_T_CELL_SIGNATURE present: TRUE   (all five)
```

So the `group_N` anonymisation that `identity_levels` enforces is undone by
`contrast_levels`, and the model can ask for a differential expression between
two named clusters — precisely the analysis the boundary exists to prevent.

The same leak fires for **any** second column carrying the clustering under a
different name (`cell_type`, a `seurat_clusters` copy, an annotation column
derived from the clusters) — a state that is extremely common in real data.

*One-line fix:* exclude a column from `contrast_levels` when it **is** the
identity column **or** when its value set equals the identity's value set
(`identical(sort(unique(as.character(cd[[column]]))),
sort(unique(as.character(Idents(object)))))`), and treat the field name
`levels` as a `group` kind rather than letting it be promoted to `aggregate`.

### F6 — The privacy scanner cannot see the object — **MEDIUM**

`R/ai-privacy.R:67` is `invisible(object)`. The scanner is a pure function of
the *context*. That is a defensible design (the context is the only thing that
goes out), but it means the privacy layer has no way to notice that a context
field named `levels` is carrying the clustering: it cannot compare the values
against `ActiveIdent()`. F5 is therefore a *context-builder* defect that the
privacy layer structurally cannot catch. Any future context field can leak
the same way.

*One-line fix:* give the scanner the object (it already receives it) and add
a value-level check: any character vector that is a permutation of the
identity's levels is treated as a `group` vector, whatever it is called.

### F7 — The AI tool loop bypasses the privacy layer entirely — **MEDIUM**

`R/ai-adapter.R:440-446` sanitises the *initial* context. But the plan round
also hands the model a **tool registry** (`R/ai-tools.R:19-74`), and tool
results are returned to the model **verbatim** — `get_status`,
`get_ledger(detail = "full")`, `get_capabilities`, `get_analysis_record`,
`get_lineage`, `get_quality_checks`. The tool is handed to the provider with
`execute = spec$handler` (`R/ai-adapter.R:120`) — the raw handler, with no
sanitiser in the path — and the agent's own loop (`R/ai-adapter.R:564-585`)
returns the result to the model. None of the six tool outputs is passed
through `sclet_ai_sanitize_context()`.

Measured: `get_ledger` and `get_capabilities` return the raw
`capabilities$contrast_levels`, i.e. the same cluster labels that the
sanitised payload carries. The sanitiser is therefore not a boundary around
"what leaves the machine" — it is a filter on one path out of two.

*One-line fix:* wrap every tool result in `sclet_ai_sanitize_context()` before
it is returned to the provider, in `sclet_ai_call_tools()`.

### What the privacy layer does get right (measured)

* A bare design column's levels **do** go out: `CTRL_LEVEL` / `STIM_LEVEL`
  present in the sanitised payload at all three privacy levels
  (`strict`, `standard`, `local`) — which is what the design requires
  (a `group` contrast cannot be planned without them).
* The identity column is excluded from `contrast_levels` **when it is named
  exactly as the active ident** (e.g. `ActiveIdent(sce) <- "ClusterSecret"`),
  so the leak is specific to the `colLabels`/`label` mismatch and to duplicate
  columns.
* Matrix, S4, function, environment and factor values are all refused
  (`R/ai-privacy.R:137-148`, `334-337`); data frames are reduced to shape.
* `note`/`comment`, paths, barcodes, donor/patient/subject fields are on the
  deny list (`R/ai-privacy.R:356-362`).
* `privacy = "local"` + `user_consent = TRUE` does let unknown atomic values
  through raw (`R/ai-privacy.R:232-234`) — opt-in, consented, and `local` by
  name. Not a defect, but it should be documented in the `privacy` argument's
  `@param` (which currently advertises a `"off"` level that does not exist —
  `match.arg` at `R/ai-privacy.R:22` accepts only `strict`/`standard`/`local`).

---

## 4. What is actually production-ready

Verdicts are based on a runnable check, not on how the code reads. The
deterministic subsystems were checked by a 22-assertion smoke script
(Appendix A.5) — **22 pass, 0 fail** on the working tree. The provider-facing
subsystems are judged from the recorded benchmarks in
`../sclet-benchmark-ai/results/*.jsonl`, with the measurement date and the
`sclet` commit named.

**Legend.** `prototype` — runs, but the contract is not met. `works` — the
contract is met and measured, but the sample is small or the behaviour is not
yet load-bearing. `production-ready` — measured repeatedly on real data, with
a failure mode that is visible when it breaks. `unknown` — not verified by this
audit.

| Subsystem | Verdict | Evidence (runnable check or measured record) |
|---|---|---|
| Analysis context / ledger (`sclet_ai_context`, `GetAnalysisLedger`) | **works** | Smoke: fingerprint present, capabilities/story/actions present, no matrix in the payload. Bounded by construction, but it is the component that leaks the clustering (F5) — *works, with a known hole*. |
| Privacy / payload boundary (`R/ai-privacy.R`) | **prototype** | Boundary rule is right; the implementation leaks the clustering (F5) and the tool loop bypasses it (F7). Smoke script passes only because it tests the scanner, not the boundary. |
| Plan validation (`ValidateAIPlan`) | **production-ready** | Smoke: valid plan passes, empty plan refused, unknown action refused, unregistered-group action names its group. Measured `plan_guardrails.jsonl` 10/10 passed (2026-10-16). |
| Plan generation + repair (`AIPlanAnalysis`) | **works** | `ai_plan_repair.jsonl` (28 records, 2026-10-09→18): first-attempt valid **4/28 = 0.143**; after repair **16/28 = 0.571**; repair rescued **12/28 = 0.429**. i.e. *most plans are invalid on the first try and repair fixes fewer than half of them*. |
| Goal discrimination | **works** | `goal_discrimination.jsonl` (18 records, 2026-10-16→17): `plan_valid` 11/18 = 0.611, `hit_any` 14/18 = 0.778, but `hit_first` **2/18 = 0.111** — the right action is usually in the plan, rarely in the first attempt. |
| Execution engine (`ExecuteAIPlan`) | **production-ready** | `execute_plan.jsonl` (28 records, 2026-10-09→16): `ok` 28/28, `recorded` 23/23, `ledger_has_execution_record` 23/23, **`dry_run_mutated` 0/23** and **`dry_run_recorded` 0/23**, `n_missing_artifacts` max **0**. The dry-run guarantee and the output contract both hold. |
| Evidence authoring & re-run safety | **works** | `authorship_multinode.jsonl` (4 records, 2026-10-14): `ok` **3/4 = 0.750**. One in four runs did not produce a citable artifact. |
| Evidence independence / claim ceiling | **production-ready** | `derived_claim_levels.jsonl` (10 records, 2026-10-12→13): `derivation_available` 10/10, `ok` 10/10, `ceiling_status` `allowed` in 6/10 and `unsupported` in 4/10 — the ceiling discriminates. |
| Interpretive evidence references | **works** | `interpret_evidence_refs.jsonl` (4 records, 2026-10-13): findings-with-refs / findings = 0/14, 0/23, 3/7, 6/16. Two of four runs cited **nothing**. |
| Interpretation + report | **works** | `interpret_claims.jsonl` 2/2 ok (2026-10-12). The round report self-check (`sclet_ai_render_round_report`, `R/ai-functions.R:734`) is the thing that makes failure visible; the *content* of the findings is not reliably cited. |
| Clarification channel | **works** | `clarification_channel_bench.jsonl` (4 records, 2026-10-17→18): `used_decision` **2/4 = 0.500** — the model used a recorded user decision in half the runs. `decision_visible_in_context` 2/4, so the visibility and the use agree. |
| Round progression judge | **works** | `m13c_data_support_adoption.jsonl` (18 records, 2026-10-21): `injection_matches_arm` 18/18; the suggestion reached the model in 12/18 trials and was adopted in 5 → **adoption | suggestion = 0.42**. |
| Full real-data round (Kang IFN-β, 32,738 × 43,095) | **works** | `realdata_kang_ifnb_reliability.jsonl` (12 runs, 2026-10-21, commit `5947e4b`): `valid_plan_rate` **0.917 (11/12)**, `completed_rate` **0.917 (11/12)**, Wilson 95% CI **[0.646, 0.985]**. The one failure was a zero-action plan (`plan_shape`); a second run carried a `contrast_invalid` prerequisite message. `m12a-batch-5947e4b.jsonl` corroborates: `n_completed` = 10,10,**0**,10,10,10,10,10,10,11,10,10. |
| Scale probe (100k cells) | **unknown** | `scale_100k_probe.jsonl` is **one** record (2026-10-19) and has no boolean outcome field. It is a single observation, not a rate. |
| Domain actions (malignant / velocity / fate / spatial / multimodal) | **works**, unevenly | `run_cnv_inference` and `call_malignant_cells` are `run` on GSE144735 (`realdata_malignant_bench.jsonl`, 4 records, 2026-10-19) but the CNV step took **3,990 s and 22.6 GB RSS**. The six spatial/multimodal actions are *gated refusals* (`gate_readiness_bench.jsonl`, 6 records, 2026-10-18) — they have never executed, by design. `run_milo` is `none`: `milo` is not installed on this machine. |

### Which numbers are single observations, not rates

Stated plainly, because the coverage table itself warns about this:

* `scale_100k_probe.jsonl` — **1 record**, no outcome field. Not a rate.
* `run_annotation_bench.jsonl` — 6 records, but "accuracy is NOT measured"
  (the table says so); records 1–4 are defects.
* `run_trajectory_bench.jsonl` — **1 record** (2026-10-18).
* `run_velocity_chain.jsonl` — 2 records (2026-10-17).
* `fate_chain_bench.jsonl` / `fate_companion_bench.jsonl` — 3 + 3 records; the
  synthetic fixture produced a **single** CellRank macrostate, so the
  multi-fate path is unproven.
* `run_doublet_rare_bench.jsonl` — 5 records, and `run_rare_cell_detection`
  landed **0 evidence nodes** because the model ordered `run_de_test` last.
* `interpret_claims.jsonl` — 2 records; `interpret_null_signal.jsonl` — 2.
* `authorship_multinode.jsonl` — 3/4 is a rate over 4.
* The Kang 12-run batch is the **only** real-data reliability sample with a
  CI, and its CI is wide: [0.646, 0.985] for a point estimate of 0.917.
* The coverage table's own note — *"the 400-s 10/10 is still one draw rather
  than a rate"* — is still true of the 10/10 line, and the 12-run batch is the
  first thing that turns it into a rate.

---

## 5. The multi-round story — implemented, and correct

### The design, and whether the code implements it

> *One plan is one round; starting the next round is the caller's job. A
> multi-stage analysis is a sequence of plans, each validated against the
> fingerprint of the object it runs on; recording a round changes that
> fingerprint, so the previous plan is refused afterwards by design.*
> (`AGENTS.md`, AI-layer architecture)

Checked, one clause at a time, by running it:

| Clause | Verdict | Evidence |
|---|---|---|
| One plan is one round | **yes** | `ExecuteAIPlan()` walks `plan$actions` in order with per-step retry and nothing rewrites later steps (`R/ai-execution.R:1609-1720`); the repair loop in `AIPlanAnalysis()` runs against **validation errors and goal coverage only**, never against results (`R/ai-planning.R:1930-1975`). |
| A round changes the object's fingerprint | **yes** | Measured: round-1 `sclet-ai-42961c38` → round-2 `sclet-ai-23d61753` after the execution was recorded. The ledger is in `metadata`, so it is inside the hashed context. |
| The previous plan is refused afterwards | **yes** | Re-executing round 1's plan on round 1's own output: `ERROR: The supplied plan validation does not match the current object or plan. Revalidate before execution.` (`sclet_ai_stale_plan`) |
| The caller can start the next round | **yes, but by hand** | `RunAIAnalysis()` (`R/ai-functions.R:395-700`) is a one-round driver: plan → validate → dry-run → execute → record → report → optional interpretation. It returns `result$object`; the caller must call `AIPlanAnalysis()`/`RunAIAnalysis()` again on that object. There is **no** `RunAIAnalysisRounds()`-style helper and no "next round" driver anywhere in `R/` (the only `AIPlanAnalysis()` call site outside its own file is `R/ai-functions.R:428`). |
| Round N+1 is *planned from* round N | **partially** | `AIRecommendNextStep()` (`R/ai-tools.R:464-509`) proposes next actions and `sclet_ai_previous_round_progression()` compares round N with round N−1, but the **proposed follow-up is not itself a plan**: it is a free-text recommendation plus `proposed_actions`, with no validation, no success criteria, and no path into `ExecuteAIPlan()`. The roadmap's `proposed_followup_plan` field (`ai-product-roadmap.md:168`) is not implemented. |
| The data contradicting the plan mid-execution | **fails loudly, and the round is still recorded** | See below. |

### F8 — A contradiction mid-execution fails loudly, but the failed round is recorded and looks like a round — **MEDIUM**

`prerequisites` are checked at **validation** time only
(`R/ai-planning.R:1040-1051`). `ExecuteAIPlan()` does **not** re-run them; it
checks the parameter types, the return type, the *output* contract
(`sclet_ai_action_output_problems`, `R/ai-execution.R:946-1000`) and the
*state* contract (`R/ai-execution.R:1001-1010`). The deferral is explicit and
documented — `R/ai-registry-annotation.R:216-232` says the producing step's
"output contract and the execution-time check are the authority" — but the
execution-time check verifies what the step **produced**, not whether its
**input** is still valid.

Measured, with a plan whose step 2 depends on a cluster that the data will
not have:

```
valid at validation: TRUE
status: failed
    s1 find_clusters -> completed
    s2 run_de_test   -> failed | One or both of the identity groups have no cells
recorded: TRUE
ledger ai_execution records: 1
evidence on object: 0
```

Good: the failure is explicit, the plan stops, and no evidence is written.
Two things are not good:

* The failed round **is recorded** as an `ai_execution` ledger record and
  **is rendered as a round report** by `RunAIAnalysis()` (the report is
  rendered at `R/ai-functions.R:600-620` on every round, whatever the
  execution status). The `continue_blocked` flag
  (`R/ai-functions.R:620`) is set only when the **report self-check** fails —
  not when the **execution** failed. A caller reading
  `report$status == "failed"` gets the truth, but a caller reading
  `report$verification_passed == TRUE` gets a green light on a failed round.
* With `continue_on_error = TRUE` on the failing step the status becomes
  `completed_with_errors` (`R/ai-execution.R:1670-1674`) and the rest of the
  plan runs. That is the designed behaviour and is fine, but the report
  carries no distinct flag for it.

*One-line fix:* set `report$continue_blocked` when
`execution$status %in% c("failed", "completed_with_errors")`, not only when
the report self-check fails.

### The ergonomics gap

The design is honest ("starting the next round is the caller's job") but the
caller has to do three things by hand, in this order, and get the validation
object right:

```r
plan   <- AIPlanAnalysis(sce, goal = "...")            # plan$validation is the validation
exec   <- ExecuteAIPlan(sce, plan, registry,
                         validation = plan$validation,  # the OBJECT, not the token
                         dry_run = FALSE,
                         confirmation = plan$validation$confirmation_token)
sce2   <- exec$object                                   # round 2 runs on THIS
```

Getting that wrong is exactly the failure `AGENTS.md` already records
("passing only `confirmation =` makes the engine re-validate and mint a token
that can never match"). The package has no helper that performs a round and
hands back a *ready-to-use next-round* plan.

---

## 6. Ranked findings

| # | Severity | Finding | Where | One-line fix |
|---|---|---|---|---|
| F5 | **critical** | The clustering is published verbatim in `capabilities$contrast_levels` (the `colLabels`/`label` name mismatch defeats the exclusion), and the engine then accepts a cluster-vs-cluster contrast naming those labels. | `R/ai-context.R:900-904` (+ `919`, `R/ai-privacy.R:122-128`, `324-329`) | Exclude a column whose value set equals the identity's, and classify `levels` as `group`, not `aggregate`. |
| F2 | **high** | A validation minted with `object = NULL` mints a confirmation token that binds to no object; `ExecuteAIPlan()` then skips the object check entirely. | `R/ai-planning.R:1161`, `R/ai-execution.R:1545-1546` | Refuse a non-dry execution when `validation$context_fingerprint` is `NULL`. |
| F3 | **high** | The validation token does not bind the registry: a validation minted against R1 executes a different handler from R2 and reports `completed`. | `R/ai-execution.R:1144-1168`, `1516-1552` | Add a `registry_fingerprint` to the validation and compare it in `ExecuteAIPlan()`. |
| F1 | **high** | The benchmark coverage table says 32 actions; the registry has 34. `run_decontx` and `compare_integration_routes` have no row. | `../sclet-benchmark-ai/docs/action-coverage.md:1-10` | Re-count and add the two rows; add a test that fails on a count mismatch. |
| F7 | **medium** | The AI tool loop returns `get_ledger`/`get_capabilities` to the model verbatim, bypassing the privacy layer. | `R/ai-tools.R:19-74`, `R/ai-adapter.R:440-446` | Sanitize every tool result in `sclet_ai_call_tools()`. |
| F8 | **medium** | A round whose execution failed is still recorded and rendered as a report; `continue_blocked` is set only on a report self-check failure, not on an execution failure. | `R/ai-functions.R:600-625` | Set `continue_blocked` when the execution status is `failed`/`completed_with_errors`. |
| F4 | **medium** | The object binding is a 31-bit additive checksum over the ledger, and the ledger carries no expression values — two objects with the same state share a fingerprint, and a plan validated against one executes against the other. | `R/ai-context.R:1055-1063`, `148` | Use `digest::digest()` over the ledger, and include a data-level digest. |
| F9 | **medium** | No completeness pin: an action can ship while absent from the goal map, the remedy table, the skills and the coverage table. `run_decontx` is in that state today. | `R/ai-planning.R:1570`, `R/ai-data-support.R:604`, `inst/skills/` | Add a test asserting every registry action appears in the goal map or the remedy table. |
| F10 | **low** | `sclet_ai_call()`'s `@param privacy` advertises a `"off"` level that `match.arg` does not accept. | `R/ai-adapter.R:407-409` vs `R/ai-privacy.R:22` | Fix the `@param` text, or implement `"off"`. |
| F11 | **low** | The confirmation token is replayable while the object is unchanged (`record = FALSE` allows unlimited re-execution). | `R/ai-execution.R:1592-1598` | Document it, or bind the token to a one-shot nonce when `record = TRUE`. |
| F12 | **low** | `AIPlanAnalysis()` is the only entry point that produces a plan; there is no helper that runs a round and returns a validated next-round plan, so the multi-round contract is 4 hand-written calls. | `R/ai-functions.R:395-700` | Add `RunAIRound()` / `NextRoundPlan()` returning `{object, plan, validation}`. |

---

## 7. What I could not verify

* **No new provider call was made.** Every provider-dependent number is from a
  recorded benchmark; the most recent real-data batch is commit `5947e4b`
  (2026-10-21), and the working tree was **dirty** in that run
  (`sclet_dirty: true` in the records), so the tree has moved since.
* **The 100k scale probe is a single record with no outcome field.** It is
  not evidence that the layer works at scale.
* **`run_milo` cannot be exercised here** — `milo` is not installed on this
  machine; the action is `none` in the coverage table for that reason.
* **The 2-byte collision construction for the 31-bit checksum is a property of
  the formula, not a demonstrated end-to-end collision.** I did not construct
  one through the public API and do not claim one.
* **Spatial and multimodal actions have never executed** — all six are gated
  refusals by design, so their execution paths are unproven.
* I did not run `make check` (25 min) — the AI test directory was run instead
  (see §8); codoc and Rd checks were out of scope for this audit.

---

## 8. AI test suite

`devtools::test_dir("tests/testthat", filter = "^ai-", reporter = "summary")`
over the 69 `test-ai-*.R` files: see the summary printed at the end of this
file.

---

## Appendix A — the scripts

All scripts live in `/tmp/audit-ai/` (nothing was written into the repo
except this report).

| Script | What it does |
|---|---|
| `attack1.R` | Attempts A and B: validate against A, execute against B; validate with `object = NULL`, execute against B. |
| `attack2.R` | Attempts D (registry substitution) and the fingerprint sensitivity probes. |
| `attack3.R` | Attempts E, F, G: stale plan (params), stale plan (action swap), stale object (a round recorded in between), plus the cell-permutation fingerprint test. |
| `attack4.R` | End-to-end "validate against A, execute against the permuted B". |
| `replay.R` | Token replay: `record = FALSE` twice, then `record = TRUE`, then a replay on the recorded object. |
| `privacy1.R`–`privacy7.R` | The privacy sweep: secrets planted in colData / rowData / metadata / assay / layer names and in a duplicate clustering column; the sanitised payload and every tool output searched for them; the `FindClusters()` path; the accepted cluster-vs-cluster contrast. |
| `multiround.R`, `multiround2.R` | The mid-execution contradiction, `continue_on_error`, the stale-plan refusal, and the round-2 validation. |
| `maturity.R` | The 22-assertion deterministic smoke check behind the maturity table. |
| `skills.R` | The action-surface coverage computation (goal map / remedies / skills / coverage table). |
