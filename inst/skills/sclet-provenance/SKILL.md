---
name: sclet-provenance
description: >
  Provenance and reproducibility guidelines for sclet AI plans and results:
  what a plan must state so an analysis can be reproduced (parameters that
  change the result, seeds, thresholds, the exact action names), what the
  ledger records automatically and what it cannot, why re-running a plan is
  the only honest reproduction check, and how to report an analysis chain so
  downstream readers can trace every number to the step and parameters that
  produced it, and how to read `analysis_story$data_support` (what the data
  can and cannot answer, and which analyses are unexplored) before planning.
  Use when writing a plan's rationale or success_criteria,
  interpreting the execution ledger, comparing two AIAnalyses, or reporting
  results back to the user.
---

# Provenance guidelines for sclet analyses

These guidelines align sclet's ledger with community provenance standards
for computational workflows (RO-Crate workflow-run provenance, Garcia-Bernardo
et al. 2024, PLOS ONE, DOI 10.1371/journal.pone.0309210; IEEE 2791
BioCompute Object for genomic pipelines; and the reproducibility checklist
for single-cell transcriptomics in "Twelve Tips for Reproducible Analysis of
Single-Cell Transcriptomics", ACM SQUR 2025, DOI 10.1145/3736731.3746138).
The engine already records most of this automatically; the guidelines are
about what the *plan* must state so the record is complete, and how to read
what the ledger gives back.

## What the engine records, and what only you can supply

Recorded automatically for every executed plan, keyed by execution id:

- every step's action name, resolved parameters, status and timing;
- the execution summary of each step (class, dimensions, assays,
  reductions, active ident/graph);
- one evidence node per evidence-producing action, with its own values,
  claim level and parent lineage;
- the plan's validation fingerprint, so a plan cannot be silently altered
  between validation and execution.

Read the per-step rows back with
`GetAnalysisLedger(object, detail = "full")$analyses[[execution_id]]$params$steps`:
one entry per step with `id`, `action`, `status`, `duration_sec` and
`params`. Scalars are stored verbatim; a longer vector keeps only its type
and length (a 2000-gene `features` vector is recorded as
`list(type = "character", length = 2000)`), so the state records are what
make that step reproducible.

Only the plan's `rationale` can record WHY: the biological motivation for
the chosen stages, the goal interpretation, and what the plan deliberately
does NOT do. A plan whose rationale is empty is reproducible but not
auditable -- a reader cannot tell whether a missing step was a mistake or a
decision.

## Parameters that change the result must be explicit

A plan is only reproducible if every result-affecting parameter is written
into it, not left to engine defaults that may change between versions:

- `resolution` and `seed` for `find_clusters` (the two most result-shaping
  choices; Louvain is stochastic, so the seed is what makes the clusters
  reproducible -- write both even when taking the defaults);
- `ncomponents` for `run_pca` and `dims` for `find_neighbors`/`run_umap`;
- `nfeatures` (or the selection method) for `find_variable_features`;
- `min.pct`, `logfc.threshold`, `pval_cutoff`, `logfc_cutoff` for the
  marker and enrichment layers;
- `sample_col`/`condition_col`/`contrasts` for `run_milo`, `batch` for
  `run_integration` -- design parameters ARE provenance.

State in the rationale which values are defaults and which are choices.
"A reader must be able to re-run this plan and get the same clusters" is the
test.

## The claim chain: every finding traces to evidence

The ledger enforces what the provenance standards ask for:

- findings cite `evidence_refs`; a finding citing nothing derives to
  `hypothesis`, which asserts nothing;
- an evidence node's `parents` tie derived measurements (enrichment) to the
  measurement they were derived from (the marker test), so the chain and
  the enrichment are one line of support, not two;
- the claim ceiling is computed from recorded evidence, not from prose.

When interpreting results, keep the chain intact: report which execution and
which evidence node a statement rests on. A statement like "cluster 3 is
responding to interferon" is only reportable at the level its cited evidence
supports -- if the evidence is one marker test on one sample, the statement
is a hypothesis, and writing it as a conclusion breaks the chain the engine
just built.

## Comparing analyses

`CompareAIAnalyses()` reports criterion-based differences and normalizes
`not_available` metrics. When comparing two plans:

- different parameters are different analyses, even with the same actions
  (resolution 0.4 vs 0.8 is not the same clustering);
- a plan that declared `success_criteria` is comparable on those criteria;
  one that declared none is comparable only on execution facts;
- never average or summarize across analyses whose fixtures or parameters
  differ -- that is the comparison-measures-nothing trap.

## Reading an anonymized group code back to its label

Every group that reaches an evidence node is anonymized to `cluster_N`, and the
codes are assigned by `paste0("cluster_", seq_along(unique(<column values>)))`
over the grouping column. They are positional, not stable: `cluster_7` means
"the seventh distinct value of that column on this object", and re-levelling
the factor, subsetting the object, or grouping by a different column all move
the mapping.

To read a code back, take the same vector the record used and index it:

```r
levels(factor(as.character(object[[the_group_column]])))[7]   # -> "cluster_7"'s label
```

`the_group_column` is the record's `inputs$group` (for a trajectory) or the
column named in the fate run's `cluster_key`. Measured on the pancreas object
(2026-10-20): `inputs$group = "cluster"`, whose seventh level is
`Ngn3 low EP`, which is the `start_group = cluster_7` the trajectory evidence
reported.

Two rules follow:

- **never carry a code across objects.** `cluster_7` on a filtered object is a
  different group than on the parent. Re-derive it from the object the record
  was written on.
- **a code is not a cell type.** Report the code, and add the label only when
  you resolved it this way from the same object -- and say it was resolved,
  because the evidence node itself carries only the code.

## Reading `analysis_story$data_support` before planning

`GetAnalysisLedger(object)$analysis_story` records what this object has *done*.
`analysis_story$data_support` records what this object can and cannot answer --
computed deterministically from the object, never from a model. Read it before
writing a plan; it is the difference between planning for this object and
planning for the last object you saw.

| Field | What it tells you | How to use it |
|---|---|---|
| `shapes` | Boolean facts about the data itself: counts present, spliced/unspliced present, several samples, replicates per condition, spatial coordinates, multimodal assays, candidate design columns | Do not infer these from the action list or from memory -- read them here. `has_counts = FALSE` means the velocity chain is gated, not "try it and see" |
| `supports` | Which question classes the data supports, each with a `reason` naming the gate and an `unblocked_by` | A `supported = FALSE` entry with its reason is the sentence to put in your rationale. Reuse the gate's own vocabulary so the record and the refusal agree |
| `blocked` | Actions currently gated on this object, with the gate name | These are not failures; they are the plan telling you what needs the user's input or different data |
| `unconfirmed_semantics` | Design choices still unrecorded (annotation reference, batch/condition roles, trajectory root) | These belong in a structured question to the user or in `confirm_design` -- never defaulted |
| `unexplored` | Analyses this object can support that the ledger has no result for: `what` + the `action` that would produce it | **This is the main input to "what should the next round do".** It is the list of data value not yet collected. An `evidence_gaps` entry says an analysis lacks evidence; an `unexplored` entry says the data was never analysed at all |

Two cautions. First, `unexplored` is a suggestion computed from the data, not a
mandate: an analysis can be unexplored and still not worth doing for the user's
question, and you say so rather than doing it to empty the list. Second, when a
field here disagrees with what you believe about the object, the field wins --
it is derived from the object, your belief is derived from an earlier message.

## What NOT to claim about provenance

- Never claim a step "was executed" unless it is in the ledger with a
  completed status. Dry-run and planned steps are not evidence.
- Never report a derived number (a fold change, an FDR quantile) without the
  thresholds it was computed under; they are parameters, and they change the
  result.
- Never present a re-analysis as a replication: same data, same plan, same
  engine is reproducibility, not replication -- an independent dataset is
  replication. The distinction belongs in the wording of the report.
