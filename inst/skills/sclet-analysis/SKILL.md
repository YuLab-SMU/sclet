---
name: sclet-analysis
description: >
  Single-cell RNA-seq analysis guidelines for planning and ordering sclet
  AIAction plans: the stage order that must not be reordered (normalize ->
  HVG -> scale -> PCA -> neighbors -> clusters -> markers), which steps
  produce evidence other steps need (put producers BEFORE consumers), which
  parameters are modelling choices the AI must never pick for the user
  (integration method, enrichment database, CCI engine, design roles), and
  what a conservative plan does when the object lacks annotation, batches or
  replicated samples. Use when proposing an AIPlanAnalysis plan, choosing
  between registered actions, ordering plan steps, filling action
  parameters, or deciding whether to ask the user instead of guessing.
---

# Single-cell analysis guidelines for sclet plans

These guidelines encode community best practice for single-cell RNA-seq
analysis (Heumos et al. 2023, *Best practices for single-cell analysis
across modalities*, Nat Rev Genet 24:550-572, DOI 10.1038/s41576-023-00586-w;
online book https://www.sc-best-practices.org) into rules an AI planner can
follow when proposing `AIPlanAnalysis` plans over sclet's registered
AIActions. They are guidelines, not laws of nature: every rule here says what
it protects, so you can tell when it does not apply. When a rule and the
engine's validation disagree, the engine is the authority -- a refused plan
with a reason is progress, not defeat.

## The standard stage order

For a counts-only object the pipeline is:

```
normalize_data -> find_variable_features -> scale_data -> run_pca
  -> find_neighbors -> find_clusters -> run_de_test
```

This order is not arbitrary. Each step consumes the previous step's output:
PCA needs scaled HVGs, the kNN graph needs PCs, clusters need the graph,
markers need an identity. A plan that runs a step before its input exists is
refused at validation -- but *ordering mistakes that validation cannot see*
are worse, because they execute fine and answer nothing:

- **Evidence producers must come BEFORE the steps that consume them.**
  `run_rare_cell_detection` grades each rare population on the independent
  signals available at the moment it runs -- a marker test, a doublet score,
  QC metrics, sample labels. Cluster size never counts as evidence. A plan
  that runs `run_de_test` after `run_rare_cell_detection` completes 8/8 steps
  and produces zero evidence, so no rare population can be claimed as real.
  If the goal mentions rare populations, order the marker test first.
- Same logic for `run_enrichment`: it interprets a marker table, so
  `run_de_test` must precede it. And `run_fate_analysis` needs recorded
  velocity, which needs the spliced/unspliced assays.
- `run_umap` is for visualisation. It is not a prerequisite of clustering,
  and including it in every plan adds cost without information. Include it
  when the goal is exploratory or the user wants to see the embedding; omit
  it when the goal is a concrete statistical question.

## Choosing the right action for the goal

Match the question to the analysis, not to the tool you reach for most
often:

| The user asks about... | The right action | Not the answer |
|---|---|---|
| which genes distinguish populations | `run_de_test` | |
| which processes/pathways the markers imply | `run_enrichment` (after `run_de_test`) | a second marker test |
| likely doublets | `run_doublet_detection` | |
| small populations the clustering averages over | `run_rare_cell_detection` (+ marker test first) | sub-clustering alone |
| cell types, against a reference | `run_annotation` | `run_de_test` |
| populations signalling to each other | `run_cci` (engine REQUIRED, see below) | `run_de_test` |
| neighbourhoods changing in abundance | `run_milo` (design REQUIRED, see below) | `run_de_test` |
| a trajectory and its terminal states | `run_trajectory` / `run_fate_analysis` | |
| batch effects to remove | `run_integration` (design REQUIRED, see below) | |

A differential-expression test is a legitimate *component* of many of these
analyses, so "the plan contains run_de_test" is not evidence of a wrong
choice. The failure to avoid is answering a question about signalling,
abundance or cell identity with a marker test *instead of* the analysis that
answers it.

## Modelling choices that belong to the user, never to the AI

Several actions deliberately have NO default for a parameter because picking
it is a biological modelling decision. When the goal does not supply the
value, the correct plan either asks via your rationale or omits the step --
it never guesses:

- `run_cci` requires `method` (CellChat / CellPhoneDB / NicheNet). The
  engines use different ligand-receptor databases and statistical tests, so
  "CellChat" is not the safe default, it is one choice among several. If
  `context$analysis_story$user_decisions` records an engine, the choice has
  already been made: pass that engine.
- `run_enrichment` requires `db` (GO or KEGG) and takes `orgDb`/`keyType`/
  `ont`. `db` takes the ontology *database* name, not the annotation
  package: `"GO"` or `"KEGG"`, never `"org.Hs.eg.db"` -- the OrgDb goes in
  `orgDb`.
- `run_integration` requires `batch` and needs `ConfirmAIDesignSemantics()`
  before it runs. It does not accept `dims` at all (the engine chooses
  internally) -- supplying one is refused.
- `run_milo` requires `sample_col`, `condition_col` and `contrasts`. Design
  roles are never guessed from column names that happen to look right: a
  column named `sample` is not necessarily the sampling unit.
- `run_trajectory` requires `start_cluster` to be named explicitly.

If the goal names none of these, say what is missing in the rationale. A
plan that asks is conservative; a plan that guesses is not.

## When the object cannot answer the question

Refusals are information, and they are the system working as designed:

- No annotation reference cached (`run_annotation`): the action needs a
  celldex dataset the user must choose. Do not plan around it with a marker
  test and present the clusters as cell types. When the user HAS named one,
  pass it literally: `ref = "celldex::MonacoImmuneData"`, `labels =
  "label.main"`. A bare dataset name (`MonacoImmuneData`) works too.
  `ref = "<celldex::MonacoImmuneData>"` is accepted and normalized, but a
  description of the object (`"<MonacoImmuneData SummarizedExperiment, Entrez
  rownames>"`) is not a value and is refused -- measured 2026-10-17, where two
  of four attempts died on exactly that distinction. If
  `context$analysis_story$user_decisions` records the reference (and labels
  column), the user has already chosen: use that value instead of asking
  again.
- No `spatialCoords` / no multimodal records: the spatial and multimodal
  read-only actions report readiness honestly. Do not plan a spatial
  analysis for a non-spatial object.
- No `spliced`/`unspliced` assays: velocity and fate actions are gated. Do
  not plan velocity for a counts-only object.
- No replicates (one sample per condition): `run_milo`'s sample replication
  signal is unavailable and condition-level claims are not supportable from
  one sample. Say so rather than planning the test anyway.

Report the gate to the user in your rationale, with what would unblock it.
Never propose a substitute analysis and let it read as the answer to the
original question.

## Parameter discipline

- Fill every parameter with a real value of the declared type. A type label,
  a placeholder (`"<cluster labels from step_6>"`), or a prose phrase where a
  value belongs is rejected -- and it costs a repair round.
- Reference earlier outputs with `${step_id.output.field}` -- all three
  segments, always. `${step_4.reduction}` and `${step_2.variable_features}`
  are the two most common malformations; both are corrected in the
  validation error, but the reference must be *usable* to be worth writing:
  check the producing action's `output_schema` exposes the field
  (`run_pca` -> `reduction`; `find_variable_features` -> `variable_features`;
  `run_de_test` -> `detest_id`).
- Every reference needs the referenced step listed in `depends_on`.
- `run_enrichment`'s `detest_id` must name the detest record an earlier step
  produced (or one already on the object); referencing
  `${step_N.output.detest_id}` of the `run_de_test` step is the intended
  form.
- Do not add `dims` to `run_integration` -- it does not accept one.
- `scale_data` without `features` scales all genes (slow on a whole-genome
  object); passing `${step_2.output.variable_features}` restricts it to the
  HVG set, which is the standard pipeline.

## What NOT to do

- Never claim a biological conclusion from QC alone. Clusters are
  data-driven groupings until markers and evidence support a type.
- Never remove cells the plan did not say it would remove: the doublet
  action's contract is to LABEL, not to subset.
- Never present a `hypothesis`-level statement as a finding; the claim
  ladder exists so recorded claims carry only what the evidence supports.
- Never invent a criterion unrelated to what the plan actually produces --
  `success_criteria` is a contract you will be graded against.
