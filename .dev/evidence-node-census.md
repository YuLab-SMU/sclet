# Census: how many evidence nodes does a real pipeline produce?

Date: 2026-10-13. Offline only -- no provider call is needed to produce evidence
nodes, they are created by `AIAction` handlers during `ExecuteAIPlan()`.

Completed by the Lead from a partly-finished teammate run: the run script is
`.dev/census-run.R` (writes to `/tmp/sclet_census/`), the stage outputs it
produced are quoted below verbatim, and the two questions the script did not
reach are answered from the already-recorded benchmark logs instead. Nothing
here is inferred from reading the code; every count has a source line.

## 1. Measured result: the standard pipeline yields exactly ONE node

`sclet:::sclet_ai_evidence_get_all(object)` after each stage of the core
pbmc3k pipeline (`/tmp/sclet_census/dump-*.txt`, `ceiling-*.txt`):

| stage after | n_nodes | ceiling status |
|---|---|---|
| read fixture (`00_read`) | **0** | `unsupported` / "no evidence has been recorded on this object" |
| `normalize_data` | **0** | `unsupported` |
| `find_variable_features` | **0** | `unsupported` |
| `scale_data` | **0** | `unsupported` |
| `run_pca` | **0** | `unsupported` |
| `run_umap` | **0** | `unsupported` |

So the eight pipeline actions record no evidence at all. Nodes come only from
the analysis actions (`run_de_test`, `run_annotation`,
`run_rare_cell_detection`, `run_trajectory`, `run_velocity` -- the
`RecordAIEvidence()` call sites in `R/ai-execution.R`), and only when the
action names its output (`run_de_test(name = "markers")` -> `ev:findallmarkers`).

## 2. Cross-checked against every real-provider record

`n_evidence_nodes` from `sclet-benchmark-ai/results/*.jsonl`, all 17 records
that carry the field:

| log | records | nodes per record |
|---|---|---|
| `interpret_evidence_refs` | 4 | 1, 1, 1, 1 (`ev:findallmarkers`) |
| `derived_claim_levels` | 6 (`with_evidence`) + 4 (control) | 1 each; control 0 each |
| `execute_plan` | 7 executed runs | 1 (`evidence_ids = ev:findallmarkers`) |

The 7- and 8-step plans these benchmarks execute (normalize -> HVG -> scale ->
PCA -> [UMAP] -> neighbors -> clusters -> `run_de_test`) end with exactly one
node every time, on real pbmc3k, including the run that added `run_umap`.

## 3. What this does to the assumption in the authorship spec

`spec-ai-evidence-authorship.md` section 6 justified refusing on multi-node
objects by saying "a real analysis records several (`ev:qc`, `ev:clusters`,
`ev:findallmarkers`)". **That is false as measured.** No action in the pipeline
records `ev:qc` or `ev:clusters`; clustering and QC write state records, not
evidence nodes. The consequence is the opposite of what the spec assumed:

- the single-node case is the **common** case, so the attachment rule fires
  routinely rather than rarely;
- therefore the relevance defect (`analysis-authorship-attachment-risk.md`)
  bites in normal use, not at the margin;
- and the multi-node refusal -- the one behaviour the design treats as the
  safety property -- is the case that almost never occurs on a standard
  analysis, and has still never been produced by a real pipeline here.

A genuinely multi-node object requires two or more *analysis* actions
(e.g. `run_de_test` + `run_trajectory`, or a repeated named DE run), which the
benchmarks have never executed. Reaching >= 2 nodes offline therefore needs a
hand-built plan with two named analysis steps; the synthetic two-node objects
used in `analysis-authorship-attachment-risk.md` section 3 and in the
prose-name test were built exactly that way.

## 4. Node structure (from the recorders, verified against the one real node)

`ev:findallmarkers` on the benchmark objects: `kind = deterministic_summary`,
`claim_level = associated`, `status = completed`, `raw_values_included = false`,
`source` = the named DE run. `RecordAIEvidence()` defaults `claim_level` to
`measured` and derives `dependency_group` from a hash of (parents, scope), so
two nodes recorded without distinct `scope`/`parents`/`source` land in ONE
dependency group and count as one line of support -- measured: two such nodes
still give `n_independent_groups = 1`; adding a node with a distinct `scope`
gives 2 and lifts the ceiling to `consistent_with`.

## 5. Not measured

- No real-provider object with >= 2 evidence nodes exists in any log, so the
  multi-node ceiling and the multi-node refusal are untested outside the
  synthetic objects above.
- `run_annotation` / `run_trajectory` / `run_velocity` /
  `run_rare_cell_detection` node shapes were not exercised: the census script
  stopped before the analysis stages (annotation needs a celldex reference,
  velocity needs spliced+unspliced assays that the pbmc3k 10x gene-matrix
  fixture does not have). Their recorder code is the only source for their
  `kind`/`claim_level` values, quoted in section 4 as such, not as measurements.
