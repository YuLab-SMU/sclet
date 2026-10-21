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
  replicated samples, and the four dimensions that separate a good plan from
  one that merely validates (goal coverage, data fitness, round progression,
  evidence increment) plus which data shape should decide the exploration
  order. Use when proposing an AIPlanAnalysis plan, choosing
  between registered actions, ordering plan steps, filling action
  parameters, deciding whether a plan is good enough to hand back, or
  deciding whether to ask the user instead of guessing.
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

Malignant-cell identification is a TWO-line call. Run
`run_cnv_inference(reference_col = "<normal-cell column>")` first -- the
reference must be exogenous (normal cells or normal samples the user named);
never use clusters the analysis itself produced as the reference. Then call
`call_malignant_cells(cnv_id = "<the run name>", signature_score = "<colData column>")`
or pass `signature`. The four states are malignant / likely_malignant /
undetermined / non_malignant, and `identity_genes` demotes a competing identity
programme (for a gastric tumour, the intestinal-metaplasia profile). A call
without a CNV record is refused, because signature-only is not the two-line
answer.

Reporting a four-state call means reporting what decided it, not just the
four counts. Three things belong in the same paragraph:

- **the signature cut and its source.** `signature_cut` defaults to the POOL
  MEDIAN (`signature_cut_is_median = TRUE` in the evidence), so "signature
  positive" then means "scores above half the cells" -- a relative statement
  that lands near 50% by construction (22,006 of 44,011 in a user report on
  GSE150290, 2026-10-20). A near-half malignant/signature fraction is the cut
  talking, not the biology.
- **the CNV threshold.** `cnv_reference_threshold` is the reference burden's
  quantile (default P99) widened by 3 x MAD. On 44,011 cells that can mean
  about 29 positives. `non_malignant` needs a CNV-negative line, so on a
  CNV-poor cohort "we cannot see a CNV" becomes "evidence of not malignant":
  in that report 90 of 201 DNA-hotspot-mutant cells came out non_malignant.
  State that limitation; do not let CNV-negative read as proof.
- **the pre-demote counts.** The identity demote is DATASET-LEVEL: one
  separation statistic above the threshold moves every supported/probable
  cell to undetermined, which can make `malignant` unreachable on that
  dataset. Read `<name>_call_pre_demote` (or `n_pre_demote_supported`) and
  report both states, so a reader sees what the evidence said before the axis
  check.

Also record whether the CNV run used the HMM stage (`parameter_hmm` in the
CNV evidence): `hmm = FALSE` gives a burden from the non-HMM profile, which is
a different measurement, not a faster version of the same one.

**The CNV run is the slowest action in the registry.** Measured 2026-10-19 on
1,400 cells with HMM on: 66 minutes, 22.6 GB of memory, 27 minutes of it HMM
Gibbs sampling. Do not add it because it sounds thorough; add it only when the
goal asks about copy number or malignant cells, and say in the rationale that
the user is starting a long step. Cluster-vs-cluster markers do not support a
malignant claim, so the alternative is to ask the user, never to skip the CNV
line silently.

Condition contrasts are NOT marker tests. To compare stimulated vs control
(any experimental group column), use
`run_de_test(group = "condition", ident.1 = "stim", ident.2 = "ctrl")`.
Add `covariates = "donor"` for the paired analysis "within donors" actually
means: the action then runs pseudobulk per donor x condition and limma with
design `~ condition + donor`. Without covariates it is a cell-level test and
the donor/batch structure is NOT adjusted -- do not present it as paired.
Never pass condition levels to ident.1/ident.2 without `group`: those are
identity levels, so the test would compare clusters.

**Where the two level sets come from is not the same place, and confusing them
is a measured failure** (2026-10-21, the only failure in 12 same-tree trials):
a plan wrote `ident.1/ident.2` as levels of `condition` that the object does
not have, and died with `contrast_invalid: ident.1/ident.2 must be levels of
'condition': ctrl, stim` -- the error listed the right levels, but the prompt
had told the model to take group values from `identity_levels$codes`, which is
the clustering. Read them from the right block:

- **with `group = <column>`**: the levels are that colData column's values.
  Read them from `context$capabilities$contrast_levels$columns`, take the entry
  whose `column` equals your group, and copy two of its `levels` verbatim. If
  the entry is missing, or `levels_available` is `FALSE`, the choice is the
  user's -- say so, never translate a label into a plausible neighbour.
- **without `group`** (cluster markers, `start_cluster`): the values are
  clustering codes, from `context$capabilities$identity_levels$codes`.

One symptom worth recognising: if the error's `(you wrote: X)` names a value
you did not see anywhere in the context, you invented it. Copy, do not infer.

| The user asks about... | The right action | Not the answer |
|---|---|---|
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
- `find_variable_features` takes `method` -- `scran` (the engine default),
  `scrapper` or `seurat` -- and the three select different genes from the
  same object. Write `nfeatures` always; write `method` when the choice
  matters to the question. Both are recorded in the object's `features$hvg`
  and surface in `GetAnalysisLedger()$features`, so the next round can see
  which selection is in force instead of assuming one.

## What counts as a good plan (four testable dimensions)

A plan that validates has only cleared the structural gate. These four
dimensions are what separates a plan that answers the question from one that
merely runs. Each is checkable against the ledger, so treat them as a
checklist before you hand the plan back, not as style advice.

1. **Goal coverage.** Every kind of work the goal names has an action, OR is
   explicitly declared out of scope with a reason in `success_criteria` /
   rationale. Silence about a named part of the goal is the failure mode: a
   goal that says "cluster and find markers" is not met by clustering alone,
   and not met by clustering plus an unstated decision to skip markers. What
   is refused by the data can be declared -- what is simply forgotten cannot.

2. **Data fitness.** The method and scale match the shape of THIS object, not
   the shape of the last object you saw. Check before planning:
   `context$dataset` (n_cells, assays, columns) and
   `context$analysis_story$data_support` (see the `sclet-provenance` skill).
   Two objects of different size get different answers; a 43k-cell object and
   a 3.7k-cell object are not the same planning problem simply because the
   action names are the same. When the data cannot support a step, the
   correct plan says so with the gate's own vocabulary (`counts_assay_missing`,
   `no_spatial_record`, ...) rather than substituting a nearby analysis.

3. **Round progression.** When a previous round's report exists, this plan
   must respond to it: the warnings and limitations it stated are either
   addressed or explicitly carried forward with a reason; recorded
   `evidence_gaps` get an action or a stated decision. Re-running an analysis
   that already ran, with unchanged parameters, is `no_progress` -- it costs
   a round and produces nothing. A recorded user decision (annotation
   reference, integration method, CCI engine) is settled: use it instead of
   asking again or choosing otherwise.

4. **Evidence increment.** Say what NEW evidence this round will add -- a new
   `ev:` node, or a new independent evidence line. A purely confirmatory
   round is legitimate, but then say that it is confirmatory. Do not let a
   round that adds no evidence describe itself as progress; that is the same
   defect as reporting a zero result as a positive one.

### Data shape decides the exploration order

The standard stage order above is for a counts-only object with no prior
state. Real objects differ, and the shape decides what comes first. Read the
shape first, then choose the route:

| Data shape (detectable from the object) | Plan first | Only then |
|---|---|---|
| Multiple samples x conditions, QC/doublet not yet recorded | QC and doublet detection, then confirm the design (`confirm_design`) | integration, then condition contrasts -- an unconfirmed design makes every integration step refusable and every condition claim unsafe |
| `counts` plus `spliced`/`unspliced` | velocity, then trajectory, then fate | fate needs recorded velocity; without all three assays the velocity chain is gated, and a counts-only object must not get one |
| Tumour epithelium with a reference population | CNV inference, then the four-state malignant call | subtype or identity claims about the epithelium wait until the CNV thresholds and their source are recorded |
| A cluster whose dominant label confidence is low (< ~0.5) | mark it mixed, do NOT give it a type | a type may only be assigned once independent evidence (markers, or a doublet explanation ruled out) exists |
| Many samples per condition | differential abundance (`run_milo`) is available | with one sample per condition the replication signal is absent; say so instead of planning the test |
| One sample, no condition column | the descriptive spine (QC -> clustering -> markers -> enrichment) | condition-level and abundance claims are not supportable; do not phrase them |
| No `spatialCoords`, no multimodal records | nothing spatial; the read-only gates will refuse | a spatial or multimodal question needs that data, not a proxy analysis |

Two rules that override the table: an evidence producer goes before its
consumer (see "The standard stage order"), and a modelling choice the user
owns never gets defaulted by you just because the data shape makes it
convenient.

## What NOT to do

- Never claim a biological conclusion from QC alone. Clusters are
  data-driven groupings until markers and evidence support a type.
- Never remove cells the plan did not say it would remove: the doublet
  action's contract is to LABEL, not to subset.
- Never present a `hypothesis`-level statement as a finding; the claim
  ladder exists so recorded claims carry only what the evidence supports.
- Never invent a criterion unrelated to what the plan actually produces --
  `success_criteria` is a contract you will be graded against.
