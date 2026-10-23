# AIP Spec: the execution trace as a graph

> Status: **rev 2, post-review — awaiting sign-off** (2026-10-22)
> Rev 1 was reviewed by a factual audit (`/tmp/spec-review/facts.md`) and a
> design review (`/tmp/spec-review/design.md`). Rev 2 drops the parts the
> reviews showed were either already implemented or unverifiable; the revision
> note at the end lists exactly what changed and why.
> Baseline: `devel` @ `b3924ba`, clean tree, package version 1.1.0.
> Recorded baseline (not re-run for this spec): AI suite 72 files; `make check`
> 0/0/0 at M12 close (`NEWS.md`); frozen-tree reliability 12/12 at `706b6e6`.
> The factual audit replayed this spec's storage additions on a throwaway copy
> and ran 16 `test-ai-*` files: same results as baseline, 0 failures.

## 1. Problem

`GetAnalysisLedger()` exposes the execution records (`state$analyses`, read
through `R/ai-context.R`) and the state records
(`state$states$records[[type]][[id]]`). An AI reading them can list the steps of
the last round and list the evidence, but cannot answer the question a trace
exists for:

> **Which step produced which result, and which step consumed which earlier
> step's output?**

Three verified defects stand between the ledger and that answer.

- **D1 — the plan's topology is dropped at execution.** `depends_on` does not
  appear in `R/ai-execution.R` at all (`grep -c depends_on` = 0). The per-step
  rows written into `params$steps` carry only
  `id / action / status / duration_sec / params`
  (`R/ai-execution.R:1222-1230`). The plan was a DAG; the record is an array.
  The one place the edges are known — `${step_6.output.reduction}` resolution —
  keeps them in memory and then discards them.
- **D2 — no general step→record edge.** Record ids surface only where an
  action's own output schema happens to name one (`detest_id` at
  `R/ai-execution.R:1135`, `cnv_id` at `:1145`); there is no step→record link.
  The machinery for this already exists twice: the evidence diff
  (`before_evidence_ids` / `new_evidence_ids`, `R/ai-execution.R:1279`) and the
  state-type diff for contract checks (`sclet_ai_action_state_problems`,
  `R/ai-execution.R:1001-1010`). Neither is generalised to record ids.
- **D3 — both stores overwrite by id, with no history.**
  `sclet_set_state_record()` assigns `type_records[[id]] <- record`
  (`R/state.R:824`) and `sclet_set_analysis()` assigns
  `state$analyses[[name]] <- value` (`R/state.R:633`). A re-run with the same
  id replaces the node. This is out of scope to fix, but the graph must not
  pretend otherwise (see §3).

Two things the first draft claimed are **not** defects, and rev 2 corrects:

- The planned action set is **not** gone: `params$steps` carries every action
  with its resolved parameters, and `inputs` carries `plan_id`,
  `context_fingerprint` and `action_ids` (`R/ai-execution.R:1235-1239`). What is
  gone is the *topology*, the success criteria, the assumptions and the selected
  route.
- Round progression is **not** blocked and is **not** this spec's consumer:
  M13-b's `sclet_ai_previous_round_progression()` (`R/ai-data-support.R:938`)
  already rebuilds the previous plan from `params$steps`, pinned by
  `tests/testthat/test-ai-round-progression.R`.

So the value here is narrower and honest: **a round is explainable as a path,
not as a list.** That is also the precondition for every later question —
independence across rounds, "why is this claim allowed", an external agent
reading one round — without pretending the planning success rate changes.

## 2. Goals

1. After any non-dry `ExecuteAIPlan()`, the object alone reconstructs the
   executed round's topology: for every step, the upstream steps it declared and
   the records it wrote.
2. Record→record edges are recovered **where the record already declares the
   id** — no inference, no guessing.
3. One read-only projection, `AnalysisGraph(object)`, returns typed nodes and
   typed edges with stable ids, serialised with the object (no environments, no
   paths), plus three queries: `graph_neighbors()`, `graph_producers()`,
   `graph_consumers()`.
4. Legacy objects (records written before this change) still build a graph:
   nodes present, evidence edges present, per-round `complete = FALSE`, and
   every query returns `complete` alongside its answer.

## 3. Non-goals

- **No storage rewrite.** Records stay where they are; the graph is a read-time
  projection.
- **No graph database, no external service, no new dependency.**
- **No supersession / version history** (D3). The graph is a current-state DAG
  plus the execution trace; a re-run shows two producers for one node and the
  node holds the latest values. Stated, not fixed.
- **No plan snapshot and no `success_assessment` on the record.** Rev 1 proposed
  both; the review showed no live consumer (§ revision note).
- **No `state_after` fingerprint.** The fingerprint is option-dependent
  (`sclet_ai_fingerprint()` hashes the built context, `R/ai-context.R:1116`,
  operative call `:148`), so a value stored at step time is not comparable to a
  later read. Step order already carries what it was for.
- **No decision nodes.** Nothing links `ai_design_confirmation` to a step.
- **No `graph_path()`.** Cross-round shortest-path semantics are undefined and
  no consumer needs it.
- **No change to the planning prompt or the privacy gate.**
- **No report renderer change in this spec** (see §7, stretch).

## 4. Design

### 4.1 Node kinds

| kind | source | what it means |
|---|---|---|
| `round` | execution record (`state$analyses`, `type = "ai_execution"`) | one executed round; label is `plan_id` (the task is not stored) |
| `step` | record `params$steps[[i]]` | one action of that round, with `action`, `status`, `attempts`, `duration_sec` |
| `record` | `state$states$records[[type]][[id]]` | one analysis result written by a step |
| `evidence` | `state$states$records$ai_evidence[[id]]` | one evidence node |

All four kinds share one table, so the columns are the union, not per-kind:
`id`, `kind`, `label`, `type`, `status`, `created_at`, `complete`, `action`,
`duration_sec`, `attempts` (non-applicable cells are `NA`). Evidence is the one
kind whose label is its bare id; the rest use `<type>::<id>`.

Node ids: `round:<execution_id>`, `step:<execution_id>:<step_id>`,
`record:<type>::<id>`, `evidence:<id>`. An `ai_evidence` record is **always** the
node `evidence:<id>` — never also `record:ai_evidence::<id>`; the mapping lives
in one function so a step and an evidence-sourced edge cannot disagree.

Attributes are copied from an explicit whitelist, never from `artifacts` or
`params` wholesale. This is the privacy mechanism: `sclet_ai_evidence_value_ok()`
is a *value* whitelist (numerics and `group_N`-shaped strings,
`R/ai-evidence.R:9-21`) and would reject ids and timestamps, so it cannot be
applied to node attributes; a field whitelist can.

### 4.2 Edge kinds

| edge | from → to | source |
|---|---|---|
| `step_of` | step → round | `params$steps[[i]]` of that round |
| `depends_on` | step → step | `plan$actions[[i]]$depends_on`, stored by this change |
| `produced` | step → record \| evidence | ids present after the step and not before |
| `refreshed` | step → record \| evidence | ids present both sides whose `created_at` advanced |
| `consumed` | record → record | a value in `record$inputs` that equals an existing record id |
| `source_of` | evidence → record \| round | `evidence$source` (a round target means the source names an execution record) |
| `derived_from` | evidence → evidence \| record | `evidence$parents`, resolved by id |

`derived_from` may point at either kind: evidence `parents` are sometimes
**detest record ids** (`R/ai-registry-annotation.R:135`), so resolution decides
the target and the edge keeps a `parent_kind` attribute. Every edge whose target
does not exist — a dangling `depends_on`, a `record_ids` entry for a record that
is gone, an unresolvable evidence parent — is counted in
`meta$unresolved_edges`, never silently dropped.

`consumed` is derived only from declared inputs that resolve exactly
(`clustering.inputs$graph = "knn_graph"`, `enrichment.inputs$detest_id`); an
input that does not name an existing id produces no edge. A field ending in
`_id` resolves against that type first (`detest_id = "markers"` is a detest
record even when a graph record is also called `markers`), and provenance
bookkeeping fields (`task`, `plan_id`, `action_ids`, `*fingerprint*`,
`schema_version`) are never treated as references. Evidence nodes are excluded
from `consumed`: their `inputs` mirror `source` / `parents`, which `source_of`
and `derived_from` already express with a precise kind.

### 4.3 Storage changes (one writer: `ExecuteAIPlan()`)

In `params$steps[[i]]` of the execution record, add:

- `depends_on`: character vector, verbatim from the plan step.
- `attempts`: integer, copied from the step result (today it exists only in
  `artifacts$results`, not in `params$steps`).
- `record_ids`: `"<type>::<id>"` for records **created** by this step
  (`"analysis::<id>"` for the untyped store; evidence mapped to `evidence:<id>`
  by the reader).
- `refreshed_ids`: ids present before and after whose `created_at` advanced.

`schema_version` `0.1.2` → `0.1.3` (`R/state.R:3`). Nothing branches on it; it
is bumped so a future reader learns the field set changed. Every new field is
optional on read.

### 4.4 Capturing the ids (the only delicate part)

The handler is called by `do.call(...)` at `R/ai-execution.R:1683`, inside the
retry loop, and its value is committed at `:1739` only after the output and
state contract checks pass. Therefore:

- Capture `before_ids` from the object the step starts from (`before`), once per
  step.
- Compute `after_ids` at the commit point (`current <- value`), not around
  `do.call`: a value the contract gate discards must not produce edges.
- Guard: `if (!inherits(value, "SingleCellExperiment"))` → `produced` and
  `refreshed` are `character(0)`. Ten registered actions return `"value"`
  (`R/ai-registry-read.R:23,42,53`; six in `R/ai-registry-spatial-multimodal.R`;
  one in `R/ai-registry-trajectory-velocity-fate.R:458`) and by contract they do
  not change the object; calling `sclet_get_state()` on a returned list is an
  S4 error, which rev 1 would have hit on `RunAIAnalysis()`'s default read
  registry.
- `refreshed` must be `created_at`-based, not `intersect(before_ids,
  after_ids)`: every pre-existing record intersects, so the naive form reports
  the entire ledger as refreshed on every step (reproduced in the audit).
- `sclet_ai_record_ids(object)` walks `state$states$records` and
  `state$analyses`, returns namespaced keys, and guards against an empty
  `state$analyses` (`names(NULL)` would emit a phantom `analysis::` key).

Cost: O(records) per step (tens), the same order as the existing state-type
contract check at `R/ai-execution.R:1001`.

### 4.5 Projection: `AnalysisGraph()`

```r
AnalysisGraph(object)
#> list(
#>   meta  = list(schema = "sclet-trace-graph/1", n_rounds = 2,
#>                legacy_rounds = 1, unresolved_edges = 0,
#>                object_replaced = FALSE),
#>   nodes = data.frame(id, kind, label, type, status, created_at,
#>                      complete, action, duration_sec, attempts,
#>                      stringsAsFactors = FALSE),
#>   edges = data.frame(from, to, kind, parent_kind, declared_input,
#>                      stringsAsFactors = FALSE)
#> )
```

- No `detail` parameter: all rounds are always built. A round without stored
  topology has `complete = FALSE`; the evidence edges it does have are still
  present (rev 1 would have returned empty edge tables, deleting edges that work
  today). A round with no recorded steps is also `complete = FALSE` (an
  execution record can survive a plan that ran nothing), and an object with no
  rounds at all reports `complete = FALSE` — no trace is not a complete trace.
- Queries return the answer **and** the completeness of the round(s) involved:
  `graph_neighbors(graph, id, direction = c("in","out","both"))`,
  `graph_producers(graph, record_id)`, `graph_consumers(graph, step_id)`.
- `print()` is a five-line summary, not a drawing.
- New file `R/ai-trace-graph.R`; `AnalysisGraph` exported (the external-agent
  contract), the queries exported with examples.

## 5. Acceptance criteria

- **A1** Chain (`s1→s2→s3`), fan-out (`s2→s3,s4`) and fan-in (`s3,s4→s5`)
  rounds reconstruct their `depends_on` sets exactly from the record.
- **A2** A fixture handler that writes two records has both attributed to its
  own step, and to no other step of the round.
- **A3** Re-running the same action with the same record id yields `refreshed`
  (created_at advanced) and **not** `produced`; a re-run that writes a new id
  yields `produced`.
- **A4** A legacy execution record (no `depends_on`, no `record_ids`): round
  node present with `complete = FALSE`, evidence `source_of` / `derived_from`
  edges present, queries return `complete = FALSE`, no error.
- **A5** `AnalysisGraph()` is `identical()` after `saveRDS()`/`readRDS()`.
- **A6** Both tables carry only the documented columns. A record whose
  `artifacts`/`params`/`summary` contain a barcode-shaped string (e.g.
  `"AAACCTGAGA"`) and a `colData` value appear nowhere in either table. The
  guarantee is structural (nothing from `artifacts` is read); `status` and the
  round label are record metadata and are not sanitised.
- **A7** `dry_run = TRUE` and `record = FALSE` leave no `ai_execution` record
  and therefore no `step` nodes; assert on the record store, not on the return
  value.
- **A8** A `returns = "value"` step completes within a normal round without
  error and contributes a `step` node with no `produced` edges.
- **A9** Every test that passed before still passes; the audit's additive replay
  (16 `test-ai-*` files on a copy) is the baseline evidence, to be re-run after
  implementation.

## 6. Tests to pin

New file `tests/testthat/test-ai-trace-graph.R`:

| case | asserts |
|---|---|
| chain / fan-out / fan-in | A1 |
| two-record step | A2 |
| re-run, same id / new id | A3 |
| legacy round | A4 |
| `saveRDS` round-trip | A5 |
| whitelist + planted barcode | A6 |
| dry run, `record = FALSE` | A7 |
| `returns = "value"` step in a round | A8 |
| failed step mid-plan | step present, no `produced`, later steps unaffected |
| `continue_on_error` | failed step present, planner edges intact |
| retry (`attempts > 1`) | one node per step, `attempts` stored |
| `consumed` edges | clustering→graph id, enrichment→detest id |
| evidence parent that is a record id | `derived_from` targets the record, `parent_kind = "record"` |

Plus one regression assertion in `test-ai-execution.R`: `params$steps` contains
`depends_on` for every planned action — the failure mode is a later refactor
dropping it again, which is D1.

## 7. Stretch, not in this spec's acceptance

Two report columns (`depends_on`, `produced`) in
`sclet_ai_report_steps()` (`R/ai-report.R:342`) — it returns markdown rows, so
this is a row-format change, not a field addition. File it as a follow-up once
the data layer is in; if it ships early, A8's table gains two columns.

## 8. Risks

- **Fingerprint churn.** `sclet_ai_fingerprint()` hashes the built context, so
  the new fields change the fingerprint once, on the next write. Plans are not
  persisted, and `ValidateAIPlan()` already refuses a plan whose fingerprint no
  longer matches — verified: no pinned test compares a literal fingerprint or
  pins an execution record's exact shape. This is why A9's replay passed.
- **Record size.** Four small fields per step (ids are short strings). No
  artifacts, no params copied. Target: a round record stays in the low kB.
- **Overwrite (D3).** Documented in the graph's `meta` and in the Rd; a node
  whose values are newer than the step that produced it is possible. Detecting a
  rewrite needs the record's `created_at`, which made a missing stamp a silent
  attribution bug: a same-id re-run produced no `refreshed` edge at all. Every
  state-record writer must stamp (fixed for `R/celloracle.R` and the failed-run
  record in `R/malignancy-cnv.R`; `R/cellrank.R` and `R/regvelo.R` already did
  at runtime).
- **Object replacement.** A step that returns a different object (metacell
  aggregation) takes the earlier records with it. The step result records
  `replaces_object`, and the graph reports `meta$object_replaced = TRUE` so an
  absent round is explained by the action rather than read as a projection bug.
- **`consumed` false positives.** Only exact id matches; a string that equals an
  existing id is treated as an id. If this ever misfires, the edge carries
  `declared_input` so a reader can see which field produced it.
- **`returns = "value"` actions cannot contribute record edges.** Verified
  2026-10-23: none of the ten value-returning actions writes a state or analysis
  record, so the rule loses nothing. If one ever does, its write is invisible to
  the graph and that action needs the object-handle change this spec avoids.

## 9. Revision note (what rev 2 changed)

Dropped, because the review showed they were a second step or untestable:
plan snapshot and `success_assessment` (no live consumer — round progression
already rebuilds the plan from `params$steps`); `state_after` fingerprint
(option-dependent, unverifiable as an equality); `confirmed_by` decision edges
(no stored link); `derived_from evidence→evidence`-only typing (parents can be
record ids); `detail` parameter and empty-edge legacy tables (they deleted
working edges); the report consumer from acceptance (kept as a stretch).

Added, from the audit: `attempts` in `params$steps`; `refreshed` by
`created_at`; the commit-point placement of the diff; the `returns = "value"`
guard; evidence node id unified; unresolvable parents counted in `meta`;
whitelist-based node attributes instead of the evidence deny-list;
`consumed` edges from declared inputs (edges we already have).
