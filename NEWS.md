# sclet 1.1.0

## M12: making "reliable" measurable

M1–M11 proved the single-cell spine runs end to end on real data. M12 addresses the
next question: **is it reliable, and how would a reader know?** The answer was not a
number — it was the absence of one. Nine repeats produced "5 of 9" with no
denominator, no interval, and no way to tell *what* failed. Nothing could be
repaired on that evidence, and "reliable" was a feeling.

This milestone does not add analysis capability. It gives reliability a
denominator, gives "a good plan" four measurable dimensions, and closes the loop
between a run and the report the next round plans from. Plan and forensics:
`.dev/spec-m12-plan-quality-and-reliability.md` (§13 records what was actually
built and where it deviates from that plan).

### Reliability now has a denominator

- **Twelve same-tree trials** of the Kang goal on one commit (`5947e4b`),
  classified by a closed-set `failure_class` enum into eight shapes across five
  layers (`plan_empty`, `plan_shape`, `plan_knowledge`, `plan_refused`,
  `execute_failed`, `report_unverifiable`, `record_error`, `not_measured`).
  **11 of 12 completed (91.7%, Wilson 95% CI [64.6%, 98.5%]).** The historical
  42-trial aggregate across 19 commits stands at 66.7% [51.6%, 79.0%]. The two
  intervals overlap: the improvement is directional, **not a demonstrated
  difference**, and the cross-commit aggregate must not be quoted as the current
  tree's rate because it contains defects already fixed.
- The classifier (`sclet-benchmark-ai/R/actions/classify-reliability.R`) is
  mechanical — same record, same bucket, no hand grading — and self-checks that
  every code it recognises exists in the engine's own error messages
  (`classify-reliability-selfcheck.R`). The count that matters is the shape
  table, not the rate: fix one shape, re-run, and watch it leave the table.
- `validation_errors` holds a single error string per record, not a structured
  array, so the shapes are classified by rule rather than by parsing a code —
  the self-check is what makes that trustworthy, and it now covers the 12 codes
  the engine actually emits.

### What the data can support is no longer the planner's guess

- **`analysis_story$data_support`** (new, schema 1.1) records what this object
  can and cannot answer: `shapes` (counts present *and integral*,
  spliced/unspliced, replicates per condition, spatial, multimodal, candidate
  design columns), `supports`/`blocked` (per question class, with the **existing
  gate name** as the reason), `unconfirmed_semantics`, and `unexplored` — the
  analyses this object can support that the ledger has **never run**.
  `evidence_gaps` says an analysis lacks evidence; `unexplored` says the data was
  never analysed at all. That distinction is the input for choosing the next
  round.
- Every `reason` reuses a gate name the registry already uses, pinned by a test
  that calls the live prerequisites, so the block cannot drift into a private
  vocabulary that disagrees with the real gates.
- Verified on the six M1–M11 objects (`.dev/m12-question-bank.R`, 6/6): no object
  is told to redo finished work, and each empty suggestion list is *justified*
  rather than silently produced. The pancreas object, which has no `clustering`
  record because it went from the embedding straight to a trajectory, was
  wrongly reported `cluster_identity_missing` and had its whole suggestion list
  suppressed; grouping evidence is now read from any analysis that presupposes
  one, not from a record literally typed `clustering`.

### A round ends with a rendered report, not with execution returning

- **`RunAIAnalysis(render_report = TRUE, verify_before_continue = TRUE)`** (new
  arguments, both default) renders the round's report from the ledger and carries
  its self-check into `report$rendered` / `report$verification_status`. When the
  self-check fails, `report$continue_blocked` states it, so a caller cannot
  mistake "the analysis ran" for "the numbers can be carried forward".
  Rendering was previously an optional finishing move, which in practice meant
  it did not happen — and without it the self-check that makes a failure visible
  never ran, while the next round was left reading a ledger that grows with every
  round.
- Reliability is defined operationally in this milestone: **not "never wrong" but
  "always visible when wrong"** — the ledger catches it, the report says it, and
  an independent process recomputes it. Three failure classes now exist that
  used to be silent by construction.

### Validation errors say what to change

- `plan_empty`, `depends on unknown step(s)`, `references unknown step(s)`,
  `contrast_invalid` and `contrast_missing` now carry a `fix:` clause naming the
  repair. The 4 failures in 9 earlier trials were all correctly refused shapes
  that the messages did not tell the model how to clear.
- **`capabilities$contrast_levels`** (new): the levels of the colData columns a
  `group` contrast can name. The one failure in the 12 same-tree trials was
  `contrast_invalid: ident.1/ident.2 must be levels of 'condition': ctrl, stim` —
  not a validator defect. The prompt sent group-valued parameters to
  `identity_levels$codes`, which describes the **clustering**, while a condition
  contrast needs the levels of a **design column**, a set with no representation
  in the context at all. Both prompts now state which is which, and
  `contrast_invalid` names the invalid value you wrote alongside the valid ones.

### The design/result boundary, made explicit

- **Design semantics may be exposed; analysis results may not.** The privacy
  layer is not "no colData values leave": a `group` contrast cannot be planned
  without its column's levels, and withholding them is exactly what forced the
  model to invent one. `ActiveIdent` (the clustering) is the result the analysis
  exists to discover and stays anonymised as `group_N`; a design column
  (`condition`, `Treatment`, `donor`) is the user's own experiment and its levels
  must be present. `contrast_levels` therefore **excludes the active identity
  column by name** — a mock object whose clusters were named `treat1`/`treat2`
  had those labels re-emitted through `contrast_levels` until it did.
- Two traps found implementing it, both pinned: `levels` must be a character
  *vector* (an unnamed list of scalars is redacted as "unnamed nested values are
  not bounded"), and `note` is on the privacy deny list, so explanatory text
  belongs in `guidance`.

### Skills

- `sclet-analysis` gained the four plan-quality dimensions (goal coverage, data
  fitness, round progression, evidence increment) and a **data-shape → exploration
  order** table: which shape plans which route first, and what has to be recorded
  before the next stage is claimable. `sclet-provenance` documents reading
  `data_support`, with `unexplored` as the main input to "what should the next
  round do". `sclet-driver` gained the three measured failure shapes as
  right/wrong pairs, and the rule that a round ends with a rendered, self-checked
  report.
- Selection widened: a planning goal with no reporting vocabulary
  ("what else can this data support") now receives the provenance skill, which
  carries `data_support`. Every one of these sections is pinned by an assertion
  that it reaches the prompt, because `options(sclet.ai.use_skills)` is a switch
  and an edit that never reaches the model is invisible in every other test.

### Tests

- New files: `test-ai-data-support.R` (63 assertions; the determinism pin caught
  a `return()` inside `tryCatch` that collapsed an entire `shapes` block to a
  bare `TRUE`), `test-ai-round-report.R` (14), `test-ai-plan-error-guidance.R`
  (41). `test-ai-skills.R` grew to 68 and `test-ai-ident-groups.R` to 48. Full
  suite: **5951 passing, 0 failures**; `make check` 0 errors / 0 warnings /
  0 notes.

### Known gaps, stated rather than smoothed

- **Round progression (dimension 3 of 4) has no computed metric.** It needs two
  plans in view at once and this milestone plans from a single round's context;
  it exists in `sclet-analysis` as guidance a model can read, not as a number.
  Carried to M13.
- **`report_unverifiable` is currently empty** — M12-e came after every record in
  the bench set. Zero occurrences is not "self-checks never fail"; it is
  "the checks had not started running".
- Suggestion coverage in `unexplored` is deliberately conservative (9 classes).
  An omitted suggestion is silence, not an error, but silence is also not
  completeness.


Eleven closed loops (M1–M11) drove every single-cell action against real data.
Each loop produced a rendered report re-derived from the ledger alone, and the
defects below are what that surfaced — nearly all of them invisible to the
existing test suite because they were *first-execution* problems or semantics
that were wrong while every test passed. Per-defect forensics live in `notes/`
(backlog B1–B26, the per-milestone reports and verification records).

## Real-data loops

**M1–M11 ran the whole single-cell spine end to end**, and each action's first
real execution is what found the defects listed below.

- **PBMC 10k (33,538 × 11,769)**: normalize → HVG → scale → PCA → neighbors →
  clusters → markers → GO + KEGG, then annotation and doublet/rare-cell scoring.
- **Kang IFN-β (43,095 cells, 2 batches, 8 paired donors)**: design
  confirmation → fastMNN/Harmony routes → route comparison → paired pseudobulk
  DE → differential abundance.
- **GSE144735 CRC (1,400 cells)**: inferCNV copy-number inference → four-state
  malignant call.
- **scVelo pancreas (3,696 cells)**: velocity → slingshot trajectory → CellRank
  fate on real spliced/unspliced counts.
- All 33 registered actions now have an execution record; spatial and
  multimodal remain deliberate typed refusals.

## Report layer

- **`RenderAIReport()`** (new export) renders a recorded result from the ledger
  and **re-checks it**: every cited evidence id must resolve and every recorded
  claim level is re-derived from the evidence, so a report that cannot be
  verified says so in its own verification section. It calls no model and
  includes no raw values.
- The current milestone reports in `notes/` are renderer output, so `notes/` has
  one source of reports. Recorded analyses render method, design, cutoffs and
  size (clustering seed, detest thresholds, enrichment database), route
  comparisons render as a table, and each finding's stated limitations and
  uncertainty are printed beside its claim instead of staying in the record.
- `RecordAIResult()` now records the result's `warnings`; previously a
  ledger-rendered report could not show caveats the run itself raised.

## Evidence integrity

- **A re-run used to leave a stale evidence node behind, and every recorder did
  it.** An evidence id derives from the analysis name, so running the same
  action twice reproduces it; `RecordAIEvidence()` refuses the duplicate and
  every recorder wrapped that in `tryCatch(error = function(e) object)`, so the
  second run returned an object whose `colData` disagreed with the citable node
  — a node describing an older run, citable and rendering with no warning.
  `sclet_ai_record_evidence_for_run()` is now the single recorder: identical
  values are a no-op, different values without `replace = TRUE` raise
  `evidence_id_conflict`, and `replace = TRUE` refreshes. All ten recorders
  route through it and a test asserts each one's routing.
- **The clarification transcript reported intent, not outcome**:
  `recorded = TRUE` was written even when the recorder threw, so a later round
  could read a decision the ledger never stored. `recorded` now reports whether
  the write happened and carries `record_error`.
- **`ResolveAIClarifications()` collision fix**: a second answer to the same
  question silently overwrote nothing and was dropped; later answers now gain a
  distinct id instead of losing to the first.
- **Evidence nodes no longer carry fields that were never measured**
  (`top_up_ranks` was `1:10` by construction, `top_up_fractions` was ten `NA`s,
  `n_significant` restated the table's own pre-filter, two opaque integer codes
  had no consumer, `top_label_ranks` was positional).
- **Fate evidence distinguished two different populations**: the terminal-state
  share was among the cells CellRank *labels* terminal (≈1/k by construction)
  and read as "where cells go". The node now also carries
  `destination_fractions` (argmax over every cell's absorption probabilities,
  flagged `destination_fractions_are_argmax`), `n_terminal_state_cells` and
  `fate_dimension_mean_absorption`.
- **CNV and malignant nodes record the operating point** (`parameter_hmm`,
  `parameter_cutoff`, `signature_cut_value`, `signature_cut_is_median`,
  `cnv_reference_threshold`, `n_cnv_positive`, `n_pre_demote_supported`), so a
  four-state result can be re-derived or challenged from the ledger.
- **A write path that swallows its error yields a citable artifact describing an
  older run.** Measured: after a stricter re-run the record held 0 marker rows
  while the node still reported the previous `n_rows` and `min_pct`. Never wrap
  a write in `tryCatch(error = function(e) object)`;
  the audit of all 225 handlers in the package found three instances of this
  class and all are fixed.

## Claim discipline

- **The competing-hypothesis test was circular**: it defined the cells under
  discussion as those the signature had selected, so the competing programme was
  guaranteed to separate them, and its demotion swept `refuted` cells into
  `undetermined` — a refusal could be overturned by the check meant to
  challenge affirmative calls. Groups are now defined by the evidence
  (`in_scope` = at least one line positive), and `refuted` is never demoted.
- **The identity demote is dataset-level** and no longer overwrites the verdict:
  `<name>_call_pre_demote` keeps what the evidence said, so a report can show
  both. A single separation statistic could otherwise make `malignant`
  unreachable on a whole cohort.
- **`transcriptomic_cut` defaults to the pool median**, so "signature positive"
  is a relative statement that lands near half the cells by construction; the
  evidence now flags it (`signature_cut_is_median`) rather than letting a
  near-half fraction read as biology.
- **Independence rests on declaration, and is now enforced**: providers sharing
  a `source` count as one line, and re-registering an existing provider name is
  refused instead of silently replacing it.
- **The causal ceiling is enforced at record time**; `causal` is refused rather
  than silently clamped, and an unresolvable citation is an error for the writer
  and a flagged `unresolved_evidence_refs` for an audit-only caller — they used
  to disagree, so an audit-only path showed a quiet downgrade.
- `k-means` always returns two groups, so a degenerate split is reported as
  `split_degenerate` with a warning instead of returned as a classification, and
  a flipped k-means assignment is the same partition (agreement is
  `max(same, swapped)/n`).
- `AddModuleScore` control bins collapse every score to 0 when the scored matrix
  holds only signature genes, so scoring runs against the full feature set with
  both signatures excluded from each other's controls.

## Execution and validation

- **Documented spellings that could not execute.** `run_milo`'s own error
  message and the skill both said to write `contrasts = "stim - ctrl"`, but the
  string reached `miloR::testNhoods()`, which parses coefficient names; the
  contrast is now translated where the design is known and both spellings are
  recorded. The `run_milo` handler also assumed an SCE return where its callee
  returns a list, so every execution died after the expensive part.
- **Prerequisites now check what the handler actually needs**: `run_velocity`
  requires the `counts` assay its engine reads (`counts_assay_missing`; it used
  to validate and die inside Python), `run_enrichment` refuses KEGG with a
  SYMBOL table at validation instead of spending a run to return nothing
  (`key_type_incompatible`), and a non-count assay is refused rather than
  silently truncated.
- **`run_integration` could never execute** — the handler forwarded an argument
  the API does not have, and nothing had ever run the action through
  `ExecuteAIPlan()`.
- **Deferred projection across sequential steps**: a plan may reference the
  reduction, variable-feature set or identity its own earlier step produces,
  including list-valued outputs projected to multi-element vectors, which
  previously produced an unactionable "condition has length > 1".
- **The route comparison measured the wrong embedding** — it always read the
  uncorrected PCA, so three routes reported identical batch mixing while the
  per-route numbers differed 19-fold. Metrics now read each route's own
  embedding and record which one, and the comparison names a batch-mixing
  versus condition-signal trade-off instead of recommending the largest number.
- **Clustering was not reproducible**: `FindClusters()` neither seeded nor
  recorded the random stream, so two runs of one plan returned different cluster
  counts at the same resolution. It now takes `seed` (default `1L`), seeds
  locally, restores the caller's stream, and records the seed.
- **An internal `.dimred` reduction survived every run** — an unexplained
  duplicate of the PCA visible to every later step that enumerates reductions.
- **`run_enrichment` chose the annotation database for the user** (`%||% "GO"`)
  although every sibling action refuses to fill in such a choice.
- **An empty enrichment run reported `completed`** with `n_enriched = 0` and no
  result; it now fails with `enrichment_no_result` and registers no evidence.
- **Marker evidence recorded thresholds the run was not given** (constants
  instead of the plan's own `min.pct`/`logfc.threshold`), so the artifact
  described a different analysis than the one that executed.
- Pinning `FindAllMarkers` cutoffs into the detest record, so `get_detest()`
  alone says what produced its table.
- `ValidateAIPlan()` gained stop conditions, a data-scale check, human
  confirmations, success criteria and a repair loop that keeps the plan with the
  fewest validation errors and stops when a round stops improving.

## Ledger and privacy

- The AI-facing ledger carries `sclet_version` and the install timestamp, so
  "which build produced this" is answerable from the ledger; a stale install had
  previously been indistinguishable from a current one at the same version.
- The ledger records each executed step's parameters and duration, so a run is
  reproducible from its record rather than only from the driver script.
- `analysis_story$conflicts` no longer fires on every multi-round analysis: the
  append-only AI log types are excluded while genuine state conflicts still
  report.
- Payload allowlist extended for the new structural fields; the classifier's
  scope guard now also exempts `output_schema` paths, so schema metadata is not
  mistaken for group-label data.

## New and documented API

Features that landed between 1.0.2 and this release, each with the commit that
introduced it (this section is what issue #29 asked the changelog to carry).

- **Multi-route integration** — `RunIntegrationRoutes()` is the documented entry
  point: it registers the raw baseline, runs the requested strategies through
  `RunIntegration()`, and records a uniform metric contract per route (batch
  mixing, biological preservation, cluster stability, runtime), storing a
  missing metric as `status = "not_available"` with a reason instead of
  fabricating it. `CompareAIAnalyses()` orders routes by the requested
  criterion, detects batch-versus-biological trade-offs, and never
  auto-selects a best route. (`92e9351`)
- **Read-only diagnostics**, one per analysis stage, each reporting whether its
  prerequisites are present and naming what is missing: `check_enrichment_readiness()`
  (`0638704`), `check_communication_readiness()` (`7eb59ee`),
  `check_fate_readiness()` (`ff2073e`), `check_fate_drivers_readiness()`
  (`d972242`), plus `check_spatial_readiness()` and
  `check_multimodal_readiness()` (`b5398c6`).
- **Read-only comparison layers**: `check_claim_ceiling()` answers "how high may
  a claim citing these evidence ids go" without recording anything
  (`118b699`); `compare_trajectory_roots()` compares completed trajectory runs
  (`dc576af`); `summarize_spatial_evidence()`, `compare_spatial_analyses()`,
  `summarize_multimodal_evidence()` and `compare_multimodal_evidence()` cover
  the spatial and multimodal catalogs (`b5398c6`).
- **Interpretation auditing** reaches the read-only facades: `AIInvestigate()`,
  `AIExplainAnalysis()` and `AIRecommendNextStep()` audit claims the same way
  the write gate does, so advice cannot assert more than the ledger supports.
- **Advanced chains are plannable end to end**: `ValidateAIPlan()` checks each
  step against the outputs earlier steps declare, not only against the object as
  it exists before execution.


## Correctness outside the AI layer

- `RunSCENIC()` called `prune2df()` with the raw database paths and the
  GRNBoost2 adjacency in place of the ranking databases and modules; it now
  builds modules from the adjacency, wraps the ranking databases, and passes
  them correctly. It also passes `num_workers` into GRNBoost2, gained
  checkpoint/resume support so a late failure no longer reruns the multi-hour
  step, pins a complete tested dependency set, and restores `np.object` on
  numpy ≥ 1.24.
- `source = "auto"` now works across program/regulon accessors instead of
  failing with a misleading message after the real cause was swallowed.
- inferCNV integration: a gene order with no id overlap is a named error rather
  than a zero-row crash, chromosome coverage is checked (a window slides along
  a chromosome, so genes piled onto one chromosome are not ready), duplicated
  gene rows merge into one interval, and named provider states are matched by
  name rather than positionally.
- `n_cnv_regions` was written to `colData` and the state record but never
  computed, so a promised measurement did not exist; the column is gone.
- Sparse-assay gene curves and heatmaps no longer omit structural zeros, and
  `ScaleData` computes true per-row variance instead of one scalar for the whole
  matrix.

# sclet 1.0.2

+ **R-native AI integration — Phase 1**: Added a bounded, versioned `GetAnalysisLedger()` / `sclet_ai_context()` view for `SingleCellExperiment` analysis state, together with the provider-neutral `sclet_ai_call()` adapter and structured `sclet_ai_result` contract. AI responses now retain findings, evidence references, warnings, recommendations, and conservative claim levels; causal claims require explicit evidence references. The existing `sclet_copilot()` and `AuditAnalysisChain()` entry points remain available through the unified adapter.
+
+  + **aisdk native integration**: Uses YuLab-SMU `aisdk` native models, tool descriptors, `z_*` schemas, and structured-output support without requiring users to leave R. Full expression matrices are not included in the AI-facing context by default.
+
++ **R-native AI integration — Phase 2**: Added a read-only AI tool registry for deterministic status, ledger, capability, analysis-record, lineage, and quality-check queries. Added `AIStatus()`, `AIReviewQC()`, `AIExplainAnalysis()`, `AIRecommendNextStep()`, `AIAuditAnalysisChain()`, and explicit `RecordAIResult()` ledger recording without changing the active computational state.
+
+  + **Dual-stage structured summaries**: High-level AI functions can run an `aisdk::create_agent()` tool loop followed by `aisdk::generate_object(mode = "tool")` and `z_*` schema validation. If the structured-summary stage is unavailable, the tool-loop response is retained with an explicit fallback warning; strict callers can disable fallback.
+  + **Safety and maintenance**: Execution tools remain disabled in this phase; the registry is read-only, AI results are not written back implicitly, and Makefile-based package checks are documented in `AGENTS.md`.
+  + Added mock coverage for native tool conversion, structured summaries, causal-claim validation, fallback behavior, and ledger recording.
+
 # sclet 1.0.1

+ **`source = "auto"` now works across program/regulon accessors** (reported in [#28](https://github.com/YuLab-SMU/sclet/issues/28)): `sclet_resolve_program_activity()` never handled `source = "auto"`, so the documented `plot_program_dotplot(..., source = "auto")` call failed with a misleading `No program activity data could be resolved.` message (the underlying `Unsupported program activity source 'auto'` was swallowed by a `tryCatch`). `source = "auto"` now tries `geneset_scoring` then `scenic` and returns the first source that supplies the requested program, so `get_program()`, `plot_program_dotplot()`, `plot_program_heatmap()`, and `has_program()` all honor it. `plot_program_dotplot()` also collects per-program resolution failures and reports them instead of hiding the real cause when nothing resolves.

+ **`RunSCENIC()` correctness, dependency pinning, and resume support** (reported in [#27](https://github.com/YuLab-SMU/sclet/issues/27)):
  + **Fixed the incorrect `prune2df()` call.** pySCENIC expects `prune2df(rnkdbs, modules, motif_annotations_fname, ...)`, but sclet was calling `prune2df(dbs, motif_f, adj, ...)` (passing the raw database paths and the GRNBoost2 adjacency in place of the ranking databases and the modules). `RunSCENIC()` now builds co-expression modules from the adjacency via `pyscenic.utils.modules_from_adjacencies()`, wraps each `.feather` cisTarget ranking database as a `ctxcore.rnkdb.FeatherRankingDatabase`, and calls `prune2df(db_objects, modules, motif_annotations_path, ...)` before `df2regulons()`.
  + **Passed `num_workers` into GRNBoost2.** GRNBoost2 is driven by a thread-based `dask.distributed.LocalCluster` when `num_workers > 1` (and falls back to arboreto's default scheduler otherwise), so GRNBoost2 and motif pruning/AUCell no longer run with mismatched parallelism. Motif pruning uses `client_or_address = "custom_multiprocessing"` to bypass the Dask scheduler entirely, avoiding the `Nanny failed to start` / `diagnostics_port` regressions.
  + **Checked-point / resume support.** `RunSCENIC()` gained `output_dir`, `save_intermediate`, and `resume` arguments. With `save_intermediate = TRUE` the GRNBoost2 adjacency and the pruned motif data frame are persisted as parquet checkpoints; with `resume = TRUE` completed steps are reloaded from `output_dir`, so a failure in motif pruning or AUCell no longer forces GRNBoost2 (the multi-hour step) to rerun. AUCell is always recomputed against the SCE.
  + **Pinned a full, tested SCENIC dependency set** instead of only the top-level packages. `sclet_scenic_env` now also pins `setuptools==80.9.0` (keeps `pkg_resources` importable; removed in setuptools >= 81), `arboreto==0.1.6`, `dask==2024.2.1`, `distributed==2024.2.1`, `dask-expr==0.5.3`, and `pyarrow==20.0.0` (used for the parquet checkpoints). This resolves the `No module named 'pkg_resources'` and the `Nanny failed to start` / `Server.__init__()` / `diagnostics_port` failures.
  + **Restored `np.object` for pySCENIC 0.12.1.** numpy >= 1.24 removed `np.object`, which pySCENIC 0.12.1 still references (`AttributeError: module 'numpy' has no attribute 'object'`). Every basilisk run now re-adds `np.object = object` in the launched process before any callback imports pySCENIC (a no-op on numpy versions that still ship it).
  + **Actionable environment-setup errors.** When `sclet_scenic_env` (or any basilisk environment) fails to create/activate, sclet now raises a typed `sclet_env_setup_failed` error explaining that basilisk deletes the partially-created environment and pointing at the documented external-Python repair (e.g. `Sys.setenv(BASILISK_CUSTOM_PYTHON_3_10 = ..., BASILISK_EXTERNAL_DIR = ..., BASILISK_NO_PYENV = '1')` plus a manually created virtualenv) instead of a generic R-side message.

+ **BiocNeighbors Compatibility Fix**: `RunKNNPredict()` (and thus `RunReferenceMapping(method = "KNN")`) no longer fails with `unable to find an inherited method for function 'queryKNN' for signature 'X = "dgCMatrix"'` on BiocNeighbors 2.0.x-2.4.x (Bioconductor 3.20-3.22, R 4.4/4.5), as reported in [#25](https://github.com/YuLab-SMU/sclet/issues/25). Matrices reaching `BiocNeighbors::queryKNN()`/`findKNN()` now pass through a dispatch-aware helper that keeps sparse inputs sparse on BiocNeighbors >= 2.6 and transparently falls back to ordinary matrices where sparse dispatch is unavailable, so no user-side version management is needed. Also fixed the `RunSymphonyMapping()` install hint: symphony has been removed from CRAN, so the message now points to `remotes::install_github('immunogenomics/symphony')`.
+ **fastMNN / `BatchRemover` semantic alignment**: `BatchRemover()` (and thus `RunIntegration(method = "fastMNN")`) now graft the correction artifacts back onto the *original* object instead of replacing it with the stripped batchelor output. This fixes three issues reported in [#24](https://github.com/YuLab-SMU/sclet/issues/24):
  + The `corrected` embedding (fastMNN's intended clustering/UMAP artifact) is kept in `reducedDim("corrected")`, registered as a `reduction` state, and made the active reduction, so `RunUMAP` / `FindNeighbors` consume it directly.
  + The original `counts`/`logcounts` and **all genes** are preserved, so marker detection and differential expression work again after integration.
  + The "corrected" layer is no longer backed by the low-rank `reconstructed` matrix (which batchelor marks *not for quantitative analysis* and which spans only the HVGs when `correct.all = FALSE`). It now points to the full gene-resolution `reconstructed` assay when `correct.all = TRUE`, otherwise to the original source expression; the integration record reports `corrected_expression` and `embedding_only` flags so downstream tools can distinguish an expression-corrected integration from an embedding-only one.
+ **In-Silico Perturbation State Contract**: `RunInSilicoPerturbation()` now registers results through the unified `analysis-state` contract: each run creates a typed `perturbation` state record (id `celloracle_<gene>`) with inputs/artifacts/params plus a command log entry, while the shift reduction stays in `reducedDim`. `get_perturbation()` reads typed records by id, keeps the most recent record active across multiple target genes, and falls back to legacy metadata for old objects.
+ **RegVelo GPU Backend**: Added `RunRegVelo()` as an RNA velocity backend that runs RegVelo through either the packaged Python stack or a user-managed `reticulate` Python, stores the inferred velocity as a `SingleCellExperiment` assay, and registers it in the standard velocity state layer.
+ **RegVelo to CellRank Workflow**: Extended `RunCellRank()` so CellRank can consume velocity assays produced by `RunRegVelo()`, including the GPU/SLURM-friendly `reticulate` path. Added CellRank diagnostics via `CellRankSummary()`, `VelocityFateCorrelation()`, and `plot_velocity_fate_correlation()`.
+ **RegVelo Smoke-Test Kit**: Added a self-contained `inst/regvelo-smoke-test` workflow for pyenv-based server setup, small h5ad preparation, SLURM RegVelo execution, and CellRank-after-RegVelo validation without conda/mamba.
+ **Basilisk Python Environment Reliability**:
  + Fixed `RunCellFate()`/`RunCellRank()` failing with an unsatisfiable conda dependency: `sclet_cellrank_env` pinned `python=3.12` together with `scvelo=0.3.2`, whose conda build requires `scikit-learn <1.2.0`, but no such scikit-learn build exists for Python 3.12. Bumped `scvelo=0.3.4`, which relaxed the `scikit-learn` upper bound.
  + **`GLIBCXX_... / CXXABI_... not found` — clarified host requirement.** A `GLIBCXX_3.4.31 not found` (or the `CXXABI_...` sibling, e.g. `CXXABI_1.3.15`) import failure is a **host/toolchain** problem, not a sclet bug: a compiled Python extension was built against a newer C++ runtime than the executing process provides. sclet cannot upgrade a system `libstdc++` (that needs root), so it makes every basilisk run (a) **fully activate** the environment (`full.activation = TRUE`) and (b) prepend the environment's own `lib`/`lib64` to `LD_LIBRARY_PATH` for the launched Python process. When the environment actually *bundles* a newer `libstdc++` (conda-backed basilisk environments ship `libstdc++-ng`), this makes the extension resolve the environment's library instead of an outdated system one — including when the import runs inside a basilisk worker (fork/PSOCK cluster). It is a harmless no-op for the default pip/venv backend, which bundles no `libstdc++`. Running these backends therefore still requires the host to satisfy **either** of: (a) a system C++ runtime providing the required symbol — e.g. Ubuntu ≥ 22.04 with GCC ≥ 12 (Ubuntu **20.04 and 22.04 do not** provide `GLIBCXX_3.4.31`/`CXXABI_1.3.15` by default), **or** (b) a conda-backed basilisk environment, which bundles its own newer `libstdc++`. A failed Python import now raises an actionable `sclet_glibcxx_not_found` error (matching both symbol families) pointing to these two options instead of a raw `ImportError`.
  + Fixed `RunSCENIC()` failing with `PackagesNotFoundError` for `pyscenic=0.12.1` and `loompy=3.0.7`, as reported in [#23](https://github.com/YuLab-SMU/sclet/issues/23). Neither is resolvable from conda-forge (pyscenic is a bioconda/PyPI package; `loompy=3.0.7` only exists on PyPI), so `sclet_scenic_env` now installs both via `pip`. This lets them resolve on both the conda and the pip/venv basilisk backends.
  + Added an internal `.env()` accessor so the documented troubleshooting command `basilisk::obtainEnvironmentPath(sclet:::.env("sclet_cellrank_env"))` actually works (the earlier guidance referenced a helper that did not exist).
  + Fixed `RunCellRank()`/`RunCellFate()` failing with `'Categorical' object has no attribute 'get_values'` on pandas >= 2.2, which removed the legacy `get_values()` API from `Series`/`Categorical`/`Index` (reported in [#23](https://github.com/YuLab-SMU/sclet/issues/23)). Every basilisk run now restores `get_values()` as an alias for `to_numpy()` in the launched Python process — a no-op on older pandas — so any dependency in the sclet Python stack that still calls it keeps working on current pandas.
+ **Sparse-assay gene curve / heatmap fix**: `plot_genecurve()` (`genecurve_plot`) and `pseudo_heatmap()` no longer fail with `t.default(...): argument is not a matrix` (or `cannot coerce class 'dgCMatrix' to a data.frame`) when the expression assay is a sparse `dgCMatrix` (the default `logcounts` sclet produces on large datasets), as reported in [#26](https://github.com/YuLab-SMU/sclet/issues/26). The requested feature block is now materialized with `as.matrix()` before the transposition, which keeps zero-expression cells in the data (so the smoothed/GAM gene curves are not biased upward) while preserving the fast positional `cbind` + `melt` reshape. The same latent bug in `pseudo_heatmap()` is fixed alongside it.

# sclet 1.0.0

+ **State Priority & Perturbation Sensitivity Mainline**: Introduced `RunPerturbationPriority()` wrapping the Augur framework for ranking cell types by perturbation responsiveness, `RunRareCellDetection()` as a reserved interface for rare-cell identification, and `RunStatePriorityWorkflow()` as the semantic shell over the state-priority analysis mainline. Added `get_perturbation_priority()` / `has_perturbation_priority()` and `get_rare_cells()` / `has_rare_cells()` accessors.
+ **Spatial Context & Niche Mainline**: Modernized `RunSpatialDeconvolution()` with `name` parameter, typed state record, and command logging. Introduced `RunSpatialColocalization()` and `RunSpatialNiche()` as reserved interfaces for colocalization and niche detection, and `RunSpatialWorkflow()` as the semantic shell over the spatial-context mainline. Updated `get_spatial()` to support id-based retrieval across deconvolution, colocalization, and niche records. Added `get_colocalization()` / `has_colocalization()` and `get_spatial_niche()` / `has_spatial_niche()` accessors.
+ **Multimodal Expansion Mainline**: Introduced `has_multimodal_adt()` / `has_multimodal_atac()` as modality detection helpers, `register_multimodal_adt()` / `register_multimodal_atac()` to record altExp ADT/ATAC metadata in the typed state layer, and `RunMultimodalWorkflow()` as the semantic shell that auto-detects and registers available modalities. This mainline is intentionally reserved — the object contract is established now so that modality-specific analysis backends can be plugged in later without breaking the navigation model.
+ **Visualization Layer Additions**: Added `plot_perturbation_ranking()` to visualize Augur perturbation priority AUCs as a ranked bar chart. Added `plot_program_dotplot()` for program activity summarized across groups (dot size = percent expressed, color = mean activity). Added `get_program()` / `has_program()` as a unified program-level accessor that transparently resolves signatures, pathways, and regulons from gene set scoring and SCENIC sources.
+ **Unified GRN Entry Point**: Added `RunGRN()` as a semantic frontend for gene regulatory network analysis. Currently dispatches to `RunSCENIC()`, keeping "GRN analysis" as the user-facing concept and leaving room for future backends (decoupleR, FigR, etc.) through the same interface.
+ **Spatial Visualization Layer**: Added `plot_spatial_deconvolution()` (ranked bar chart of mean cell type proportions from cell2location) and `plot_spatial_composition()` (heatmap of cell type proportions across spatial spots). Both consume the typed spatial state records produced by `RunSpatialDeconvolution()`.
+ **Velocity Latent Time Visualization**: Added `plot_velocity_latent_time()` to color cells on an embedding by velocity pseudotime (latent time), showing the predicted temporal ordering inferred by scVelo's dynamical mode.
+ **Symphony Reference Mapping**: Added `RunSymphonyMapping()` as a full Symphony backend for `RunReferenceMapping(method = "Symphony")`. The workflow builds a harmonized reference atlas from the reference SCE with `symphony::buildReference()`, maps query cells with `symphony::mapQuery()`, and predicts labels with `symphony::knnPredict()`. Supports batch variable integration (`vars` parameter) and confidence scoring. Predicted labels are stored in `colData$symphony_predicted` and confidence scores in `colData$symphony_confidence`, fully integrated with the analysis-state contract.
+ **Rare Cell Detection with Density**: Replaced the `RunRareCellDetection()` RareQ/Seurat backend with a fully SCE-native `"density"` method. The algorithm builds a kNN graph via `BiocNeighbors`, computes mutual connectivity Q-values, identifies high-density seeds, propagates labels, and flags small clusters as rare. No external package dependencies beyond `BiocNeighbors`. The RareQ backend has been removed.
+ **Spatial Colocalization with SVP**: Replaced the `RunSpatialColocalization()` stub with a full SVP backend. `SVP::runGLOBALBV()` computes global bivariate spatial autocorrelation between feature columns (e.g., deconvolution-derived cell type proportions), producing a pairwise correlation matrix. Features are auto-detected from prior deconvolution results if available, or specified explicitly. The correlation matrix is stored in the colocalization state record.
+ **Spatial Niche Detection with SVP LISA**: Replaced the `RunSpatialNiche()` stub with a full SVP backend. `SVP::runLISA()` computes Local Indicators of Spatial Association for each feature, categorizing spots as High-High (hotspot/niche), Low-Low, High-Low, Low-High, or not significant. LISA categories are stored in `colData` with a `lisa_` prefix. Features are auto-detected from prior deconvolution results.


# sclet 0.99.7

+ **Data Cleaning**: Introduced `RunDoubletFinder()` to detect doublets using `scDblFinder`, and `RunDecontX()` to remove ambient RNA contamination using `celda::decontX`. Both functions natively support the `SingleCellExperiment` object and integrate seamlessly with the `analysis-state` architecture, storing clean matrices safely as new assays (e.g., `decontXcounts`) without overwriting raw data.
+ **Imputation**: Added `RunImputation()` to perform zero-preserving imputation using the `ALRA` method. The function handles sparse matrices and integrates with the analysis-state machine by storing the imputed results in a separate assay (`alra_imputed`).
+ **Diffusion Map**: Added `RunDiffusionMap()` as a robust dimensional reduction alternative for complex trajectory inference. It uses the `destiny` package to calculate diffusion components, storing them in `reducedDim(sce, "DM")`.
+ **Integration**: Refined `RunIntegration(method = "scVI")` to avoid dense matrix materialization, converting sparse matrix inputs directly to `scipy.sparse.csr_matrix` via `reticulate` for improved performance and reduced memory usage on large datasets.
+ **Memory Safety**: Audited other large-matrix workflows and now keep sparse representations as long as possible in `RunSpatialDeconvolution()`, `RunSCENIC()`, `RunCellRank()`, `RunInSilicoPerturbation()`, `RunKNNPredict()`, while also avoiding full-assay transpose in `pseudo_heatmap()` and `genecurve_plot()` unless only the requested features are needed.
+ **Phenotype Association**: Added `RunPhenotypeAssociation()` to directly link single-cell transcriptomics with bulk clinical data (e.g., survival analysis) using the `scPAS` package. Supports Cox regression, binomial, and gaussian models, tracking the association results and parameters directly in the `analysis-state`.
+ **Fate Mapping Workflow**: Refined the `velocity -> fate` chain around `SingleCellExperiment` state records. `RunVelocity()` now registers velocity outputs through the standard `analysis-state` contract, `RunCellRank()` stores terminal states, named fate probabilities, and lineage driver statistics in a structured way, and `RunCellFate()` provides a thin user-facing entry point that currently dispatches to `CellRank`.
+ **Fate Visualization**: Added native fate plotting helpers following the `plot_<module>_<type>` convention: `plot_fate_probability()`, `plot_fate_terminal_states()`, and `plot_fate_driver_trends()`. These functions consume the stored CellRank artifacts directly and expose fate probabilities, terminal state assignments, and lineage driver expression trends on top of the tracked reductions/layers.
+ **Mainline Workflow Shells**: Introduced three high-level `Run<Mainline>Workflow()` entry points as semantic shells over existing atomic interfaces: `RunTrajectoryWorkflow()` (diffusion map + Slingshot + velocity + CellRank/CellFate), `RunProgramWorkflow()` (gene-set scoring + SCENIC regulon inference), and `RunReferenceWorkflow()` (SingleR/KNN reference mapping). These shells are designed to make the recommended analysis path easier to follow and document without hiding the underlying methods.
+ **Reference Mapping Layer**: Added `RunReferenceMapping()` as a unified entry point for reference-based annotation, dispatching to `SingleR` or lightweight `KNN` label transfer with a consistent state contract and accessor layer. Added `plot_reference_label_transfer_heatmap()` and `plot_reference_label_confidence()` for downstream inspection of label transfer quality across query groups.
+ **Program & Regulon Activity Visualization**: Added `plot_program_activity()` and `plot_program_heatmap()` to unify the visualization layer across `RunGeneSetScoring()` and `RunSCENIC()` outputs. Both consumers resolve program-level activity from the typed state records, making the source (geneset scoring vs. regulon) transparent to the plotting call. Added `plot_scenic_activity()` for regulon-level embedding visualization, and `plot_trajectory_overview()` as the first native overview plot for the trajectory mainline.
+ **State Contract Refinements**: `RunGeneSetScoring()` now writes a typed `geneset_scoring` state record; `RunSCENIC()`, `RunVelocity()`, `RunSingleR()`, and `RunCellRank()` all accept a `name` parameter and register results through the unified analysis-state layer instead of ad hoc metadata slots. `get_geneset_scoring()` now reads both new typed records and legacy formats transparently.

# sclet 0.99.6

+ **Pathway Scoring & Marker Detection**: Enhanced `RunGeneSetScoring()` with `as_altExp = TRUE` to store pathway scores as an alternative experiment. Added `FindMarkerPathways()` to detect cell-type specific pathways robustly using `FindAllMarkers`.
+ **KEGG Integration**: Added `RunKEGG()` as a one-stop wrapper to download KEGG pathways via `clusterProfiler`, score cells, and perform differential marker pathway analysis seamlessly.
+ **Cell-Cell Communication**: Introduced `RunCCI()` as a unified interface for cell-cell communication analysis, seamlessly supporting `CellChat`, `CellPhoneDB`, and `NicheNet` backends. This eliminates the need for manual data formatting and environment setup (using `basilisk` for Python backends), storing results in a standardized format within the analysis state.
+ **Cell-Cell Communication Visualization**: Added two native CCI visualization methods: `plot_cci_sigmoid` (S-curve connection plot) and `plot_cci_arrow` (bidirectional arrow plot with dynamic expression mapping). These native functions eliminate the need for heavy third-party visualization dependencies.
+ **Unified Plotting Names**: Established a `plot_<module>_<type>` naming convention for visualizations (e.g., `plot_cci_bubble`, `plot_velocity`, `plot_lineage`), while retaining classic aliases like `DimPlot` for Seurat users.

# sclet 0.99.5

+ **AI Copilot & Evidence Governance**: Introduced `sclet_copilot()` as an intelligent diagnostic agent. Powered by the `aisdk` framework, it reads the complete `analysis-state` provenance of the `SingleCellExperiment` object to provide context-aware diagnostics.
+ **Cross-Chain Auditing**: Added `AuditAnalysisChain()` to perform rigorous cross-chain error control. It calculates a "State-Dependency Confidence Score" to help distinguish true biological signals from computational artifacts (e.g., over-alignment during integration), actively combating "Transcriptomic Overload".
+ **Spatial Deconvolution**: Added `RunSpatialDeconvolution()` to leverage the `cell2location` deep learning model via `basilisk`. It maps single-cell references onto spatial transcriptomics spots for high-resolution abundance inference.
+ **In-Silico Perturbation**: Added `RunInSilicoPerturbation()` to wrap the `celloracle` framework. This allows for predictive simulation of cell fate trajectory shifts following virtual gene knockouts.
+ **Advanced Fate Mapping**: Added `RunCellRank()` to combine directed RNA velocity with transcriptomic similarity, enabling robust inference of terminal states and fate probabilities for complex developmental processes.

# sclet 0.99.4

+ **Python Interoperability**: Introduced `basilisk` infrastructure (`sclet_env`) to seamlessly manage Python environments and dependencies, enabling robust integration with the Python single-cell ecosystem without polluting user's local environments.
+ **RNA Velocity**: Added `RunVelocity()` and `VelocityPlot()` to support RNA velocity analysis via `velociraptor` (`scVelo`), with results seamlessly integrated into the `velocity` analysis state.
+ **Gene Regulatory Networks**: Added `RunSCENIC()` to run the complete `pySCENIC` pipeline (GRN inference, motif enrichment, and AUCell scoring) via `basilisk`. The resulting AUC matrix is stored as an `altExp` and tracked in the `scenic` analysis state.
+ **Gene Set Scoring**: Added `RunGeneSetScoring()` to calculate gene set or pathway activity scores using `UCell`, `AUCell`, or `GSVA`. Scores are added to `colData` and registered in the `geneset_scoring` analysis state.
+ **State Management**: Updated `Status()` and accessors (`get_velocity()`, `get_scenic()`, `get_geneset_scoring()`) to fully support the newly added high-level analysis modules.

# sclet 0.99.3

+ `analysis-state contract`: Introduced a unified, lightweight analysis-state registry to manage downstream analysis states (annotation, integration, mapping, trajectory, communication, detest, enrichment) and their dependencies.
+ **Integration Workflow**: Added `RunIntegration()` as a unified entry point for multi-backend data integration, currently supporting `fastMNN` and `Harmony`. Downstream functions like `FindNeighbors` and `RunUMAP` can now automatically perceive and consume the integration states.
+ **Reference Mapping**: Added `RunKNNPredict()` and `RunReferenceMapping()` for lightweight label transfer based on `BiocNeighbors::queryKNN`. Added `ProjectionPlot()` for query-reference visualization.
+ **DE & Enrichment**: Added `RunDEtest()` (wrapping `FindMarkers`/`FindAllMarkers`) and `RunEnrichment()` (wrapping `clusterProfiler`). Included corresponding visualization tools: `DEtestPlot()` (volcano plot) and `EnrichmentPlot()` (dotplot).
+ **Visualization Ecosystem**: Enhanced visualization capabilities by leveraging the `ggsc` package. Added `CellDimPlot()`, `FeatureDimPlot()`, `GroupHeatmap()`, and `CellStatPlot()`.
+ **Technical Debt & Bug Fixes**:
  - Migrated to `scrapper::aggregateAcrossCells.se` for pseudobulk aggregation and `scrapper` for HVG selection to eliminate deprecated warnings from `scuttle` and `scran`.
  - Fixed missing dependencies (`enrichplot`, `scran`) and non-ASCII character warnings, achieving a clean `R CMD check` (0 errors, 0 warnings, 0 notes).

# sclet 0.99.2

+ Performance improvement: `FindAllMarkers()` uses `presto` for faster calculation (2026-02-06, Fri)
+ Bug fix: `BatchRemover()` robustness improvements (2026-02-06, Fri)
+ Internal: Renamed source files from `.r` to `.R` (2026-02-06, Fri)
+ `RunMilo` and `refit_milo` for Milo differential abundance analysis (2026-01-30, Fri)

# sclet 0.99.1

+ `RunSuperCell` for metacell construction (2026-01-16, Fri)
+ `RunSlingshot` for trajectory inference
+ `RunEnrichment` for functional enrichment analysis
+ `RunExplorer` for interactive data exploration using iSEE

# sclet 0.99.0

+ `ReadH5AD` and `WriteH5AD` for AnnData interoperability (2025-06-04, Wed)
+ `as.Seurat` and `as.SCE` for object conversion (2025-06-04, Wed)
+ Performance improvements: support DelayedMatrix in normalization (2025-06-04, Wed)

# sclet 0.1.0

+ `AggregateExpression` for pseudobulk analysis (2025-06-03, Tue)
+ `FindMarkers` supports pseudobulk DE via `DESeq2` (2025-06-03, Tue)
+ `BatchRemover` now records parameters in metadata (2025-06-03, Tue)

# sclet 0.0.9

+ `DimPlot` wrapper for `ggsc::sc_dim` (2025-06-02, Mon)
+ `RunSingleR` wrapper for automatic cell type annotation (2025-06-02, Mon)

# sclet 0.0.8

+ Fix `Idents` and `FindClusters` consistency (2025-05-30, Fri)
+ Fix `subset_cell` bug (2025-05-30, Fri)
+ Fix `gene_summary_table` bugs (2025-05-30, Fri)
+ Add Suggests package checks (2025-05-30, Fri)
+ Improve `NormalizeData` and `FindVariableFeatures` (2025-05-30, Fri)

# sclet 0.0.7

+ `RunSuperCell()` to merge similar cells to metacell (2025-05-28, Wed)

# sclet 0.0.6

+ `sc_merge()` to merge a named list of SingleCellExperiment objects to a single one (2025-04-17, Thu)
+ `BatchRemover()` for batch correction (2025-03-20, Thu)
+ `FindVariableFeatures()` and `VariableFeatures()` supports `method = 'scran'` (2025-03-17, Mon)

# sclet 0.0.5

+ Slingshot supports (2025-03-12, Wed)
    - genecurve_plot
    - RunSlingshot
    - lineage_plot
    - pseudo_heatmap
    - pseudo_plot

# sclet 0.0.4

+ re-export `scuttle::logNormCounts()` and `scuttle::perCellQCMetrics(()` (2024-10-08, Tue)
+ Import 'Matrix' for calculate stats of SparseMatrix (2024-10-07, Mon)
+ re-export `aplot::plot_list()` (2024-10-06, Sun)

# sclet 0.0.3

+ Adding cell-cell communication module. (2024-10-06, Sun)
    - `RunCellChat()` wrapper function for running CellChat.
+ `gene_summary_table()` to add gene summary information to marker gene table (2024-10-02, Wed)

# sclet 0.0.2

+ `FindMarkers()` and `FindAllMarkers()` (2024-10-02, Wed)
+ remove vignettes and host online book (2024-10-01, Tue)
    - <https://yulab-smu.top/sclet>
+ re-export `scater::plotColData()` 
+ `QCMetrics()` function to add seurat-like QC metrics ('nFeature_RNA' and 'nCount_RNA' to 'SingleCellExperiment' object)

# sclet 0.0.1

+ Seurat-like functions for SingleCellExperiment (2024-09-20, Fri)
    - `Read10X()`
    - `PercentageFeatureSet()`
    - `VlnPlot()`
    - `FeatureScatter()`
    - `subset()`
    - `NormalizeData()`
    - `FindVariableFeatures()`
    - `VariableFeatures()`
    - `VariableFeaturePlot()`
    - `ScaleData()`
    - `ElbowPlot()`
    - `FindNeighbors()`
    - `FindClusters()`
    - `Idents()`
    - `RenameIdents()`
    - `RunUMAP()`
    - `RunPCA()`
+ re-export:
    - `SummarizedExperiment::colData()`
    - `SummarizedExperiment::rowData()`
    - `SummarizedExperiment::assay()`

## M13: what a dirty dataset asks of the planner

M12 gave "reliable" a denominator and gave "a good plan" four dimensions. All of
that was measured on data this package had already analysed. M13 asks the next
question — **given a dataset nobody has seen, does the planner make a judgement it
could not make before?** — and pays a debt M12 left: it never proved the model
*uses* `analysis_story$data_support`, only that the field was correct.

Spec and forensics: `.dev/spec-m13-dirty-data-judgement.md` (§10–12 record what was
actually built and where it deviates). Reports: `notes/m13-ambient-report.md`,
`notes/m13-data-support-adoption.md`.

### Round progression now has a number (closing M12's B18)

- **`sclet_ai_round_progression()`** judges two rounds against each other and
  returns `advanced` / `confirmatory` / `no_progress` / **`regressed`**, plus
  `first_round` when there is no previous round. Four judges, all relations between
  the previous plan, the previous report's own fields, and the current plan:
  warnings answered (by an action that remedies them, or explicitly carried with a
  reason), unchanged steps repeated, evidence gaps closed, recorded user decisions
  actually used.
- **`regressed` is the class that had to exist.** A dimension that can only report
  "good" or "average" cannot report "worse", and re-deriving an answer the previous
  round already discarded is the most expensive move a plan can make. On the real
  M7 object the judge returns `regressed`, and correctly: that report said cluster
  12's label came from a doublet-dominated population, and a plan that re-annotates
  without decontaminating is re-deriving a discarded answer. Two independent paths
  — the engine's calculation and M7's own human-written warning — agree on the
  shape.
- `report$round_progression` carries it into `RunAIAnalysis()`'s result.
- A previous round is **two records, not one**: the plan comes from the newest
  `ai_execution` record's `params$steps`, the report from the newest `ai_report`
  record whose warnings live in `artifacts` (not `summary`). Looking for either
  half in the other place returns `not_comparable` on every real object, which
  reads exactly like "D3 was never computed here".

### Ambient RNA: the correct answer was not a plannable action

- **`run_decontx` is a new action (33 → 34).** `RunDecontX()` has been exported
  while having no AIAction, which meant the one correct response to ambient
  contamination — clean before interpreting — was not a step a planner could
  choose. M13 measured that gap and closes it. This is the first new action since
  M11, taken deliberately and recorded as an exception.
- It wraps `celda::decontX`, writes `decontXcounts` and `decontX_contamination`,
  and **removes no cell and switches no active assay** — downstream steps must name
  the cleaned matrix explicitly. A non-integer or negative `counts` assay is
  refused (`counts_assay_not_count`) rather than silently estimated from.
- The evidence node records a distribution (q25/median/q75/q95/max), never
  per-cell values, and states `contamination_call_is_model_estimate = TRUE`.
- **decontX is deterministic** on a fixed input (two runs, identical vectors, max
  abs diff 0), so an identical re-run is an idempotent no-op — the recorder's
  no-op path is the correct outcome here, not a missed conflict. A *changed* input
  under the same node name is a conflict, and that branch is pinned separately.
  Worth recording because the temptation runs the other way: a stochastic-looking
  estimator is assumed to vary, and a test asserting "a re-run conflicts" would be
  pinned to a false premise.

### The data now says how dirty it is

- **`sclet_ai_ambient_pressure()`**: two screens reported **separately and not
  multiplied** — a low-UMI tail share and the detection uniformity of the top 200
  genes — each with its own threshold, plus a coarse `label`. Multiplying two soft
  signals produces a figure with two decimal places and no meaning, and hides
  which half moved. The function states its own limits in its return value: a high
  detection rate is also what housekeeping genes look like, so `high` means
  "cleaning is worth considering", not "contamination was measured".
- On the GSE158055 sample the screen reads **high** (low-UMI tail 0.0335,
  uniformity 0.9139, both thresholds crossed), and `data_support` reports
  `ambient_cleaning: TRUE` with **`run_decontx` as its single `unexplored`
  suggestion**, carrying the screen's own numbers as the reason.
- The fixture is one sample of GSE158055 (27,943 × 8,129), not the 1,462,702-cell
  cohort: the full matrix makes decontX unrunnable here, and pooling samples would
  confound ambient signal with between-sample differences so the screen would
  measure the wrong thing. Extraction scans 12% of the file by exploiting the fact
  that barcodes are grouped by sample and the mtx is column-sorted, and it
  **fails loudly** if a sample's columns are ever not contiguous rather than
  returning a partial matrix.

### Does the planner use `data_support`? Directionally yes, not demonstrably

Two arms, same fixture, same goal (deliberately **not** mentioning contamination,
or the test would measure compliance with the goal text instead of use of the
context): 6 planning trials each.

| arm | planned `run_decontx` | rate | Wilson 95% |
|---|---|---|---|
| control (no `data_support`) | 0 / 6 | 0% | [0%, 39%] |
| treatment (with) | 2 / 6 | 33.3% | [9.7%, 70%] |

- **The zero baseline is the more persuasive half**: `run_decontx` is a registered
  action, so the planner could have reached for it unprompted, and never did in six
  tries. It is not the model's default habit.
- **Fisher exact p = 0.455.** The intervals overlap heavily. The honest sentence
  is "directionally consistent, not distinguishable from chance at 12 trials" —
  **not** "`data_support` improved planning by 33 points".
- Injection was verified 12/12 against the arm, because otherwise `adoption = 0`
  could just be a switch that never fired.
- **The finding that matters more than the rate:** of the four trials that did not
  adopt, **none mentioned ambient anywhere in its rationale**. They are not
  "considered and deliberately declined" (which the spec treats as legitimate —
  `unexplored` is a claim about the data, the goal is the user's intent). The
  suggestion was in the context and did not reach the reasoning. That is a
  visibility problem of the same kind as M12's skill-topic gap, and it is recorded
  as B20 rather than papered over. One adopting trial also changed `run_decontx`
  without saying why (B21).

### Packaging defects found while closing out

- **The committed `DESCRIPTION` could not build**: `Collate` listed 84 of 95 R
  files, missing all eleven `ai-registry-*.R`. `make check` had been passing only
  because `devtools::check()` regenerates documentation first and that run added
  them. `make check` is now **`devtools::check(document = FALSE)`** — documentation
  regeneration is `make rd`, which `make all` already runs.
- Two exported objects (`ValidateAIEvidenceRefs`, `sclet_ai_evidence_independence`)
  had roxygen blocks but no committed Rd; their documentation is now tracked.
- A stray 21 KB file named `--args` had been left at the package root and is gone.

### Tests

New: `test-ai-decontx.R` (42), `test-ai-round-progression.R` (58).
AI suite **6074 passing, 0 failures**. `make check`: **0 errors, 0 warnings, 1
note** — the note is the intentional `unlockBinding()` in `R/zzz.R` that patches a
basilisk internal, unchanged by this milestone.

### Known gaps

- **`unexplored` still does not reliably reach the reasoning** (B20). Until that is
  fixed, `data_support`'s product value is smaller than its correctness suggests.
- **12 trials cannot separate direction from chance.** More trials or a stronger
  design would be needed to claim an effect.
- Ambient pressure was screened on **one sample**, so mixed-sample ambient is
  untested.
