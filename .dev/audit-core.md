# Core-layer audit — sclet (branch `devel`, 2026-10-22)

Independent audit of the **CORE (non-AI) analysis layer** of `sclet`.
Nothing in `R/`, `man/`, `tests/`, `NAMESPACE` or `DESCRIPTION` was modified.
New files written by this audit (all read-only helpers, plus this report):

| file | what it does |
|---|---|
| `.dev/tools/audit-duplicate-definitions.R` | every top-level name in `R/`, the files that assign it, and which one wins under `Collate` |
| `.dev/tools/audit-signature-docs.R` | `formals()` vs the Rd `\arguments` / `\usage`, per Rd topic; also dumps every exported signature with `file:line` |
| `.dev/tools/audit-state-contract.R` | per module: how many exported analysis functions register a state record, call a readiness gate, or `stop()` |
| `.dev/tools/audit-core-runtime-evidence.R` | the runtime claims below (E1–E5), re-runnable in one command |

Reproduce everything:

```bash
export R_LIBS_USER=~/Rlib
cd /home/wang/data/source/omics/sclet
Rscript .dev/tools/audit-duplicate-definitions.R .
Rscript .dev/tools/audit-signature-docs.R .
Rscript .dev/tools/audit-state-contract.R .
Rscript .dev/tools/audit-core-runtime-evidence.R .
```

---

## 0. The five findings that matter most

| # | finding | evidence | one-line fix |
|---|---|---|---|
| **1** | **11 exported accessors return the *wrong* record** when asked for an id that does not exist — silently, with a fabricated `method` | `R/analysis-accessors.R:514` (`sclet_normalize_analysis_record`), used at `:28 :83 :138 :193 :248 :303 :570 :823 :879 :937 :1001 :1170`; runtime proof in E1 | in `sclet_normalize_analysis_record()`, only use the legacy `analyses` slot when **no `id` was requested**; return `NULL` when a requested id is absent |
| **2** | **13 of the 14 exported `check_*_readiness()` gates are never called by the analysis they gate** — the analysis layer never enforces its own prerequisites | `R/malignancy-cnv.R:303` is the *only* core call site; E3 lists the rest | call the matching gate at the top of each `Run*` and abort on `status != "ready"`/`not_ready` |
| **3** | **`check_de_readiness()` checks the wrong package**: it gates on `limma`, but the default marker path needs `presto` | `R/ai-interpretation-gaps.R:275` (gate) vs `R/markers.R:387` (`FindMarkers_Presto` stops on missing `presto`) | gate on `presto` (and `limma` only for the covariate/pseudobulk path) |
| **4** | **Six different analyses share one state bucket** (`type = "annotation"`), so one "active" slot serves all of them | `R/annotation.R:193`, `R/phenotype-evidence.R:557`, `R/malignancy.R:735`, `R/malignancy-cnv.R:427`, `:673`, `R/malignancy-combine.R:173`; E6 | add `cnv`, `malignancy`, `phenotype` state types to `sclet_state_types()` and give each module its own |
| **5** | **`sclet_resolve_group_column()` turns a missing active identity into base R's `"argument is of length zero"`** — an undiagnosable error in `RunCellChat` / `RunCellPhoneDB` / `RunNicheNet` | `R/utils-cci.R:24` (`!group %in% colnames(meta)` with `group = NULL`); E2 | guard `is.null(group)` and `stop("No group column ...")` before the `%in%` test |

Ranked 6–15 are in §2–§5 below.

---

## 1. DUPLICATE DEFINITIONS OF EXPORTED NAMES

### 1.1 Verdict: **no duplicate top-level definitions remain** — verified two independent ways

* **Parser-based** (`.dev/tools/audit-duplicate-definitions.R`): all 95 files in `R/` are parsed with `keep.source = TRUE`; every top-level `<-` / `=` assignment is collected with its file, line and `Collate` rank.
  * total top-level names: **757**
  * names assigned **more than once anywhere in the package: 0**
  * exported names defined in more than one file: **0**
  * exported names assigned twice in one file: **0**
* **Text-based cross-check**: `grep -nE '^[A-Za-z._][A-Za-z0-9._]*\s*(<-|=)[^=]' R/*.R` → 764 column-0 assignments; `cut -f1 | sort | uniq -d` → **empty**.
* `sclet_copilot()` now exists exactly once: `R/copilot.R:43`. The dead second definition is gone and the explanatory comment at `R/ai-functions.R:87-98` now states the truth.

### 1.2 `Collate` is roxygen2-generated and in sync — checked, not assumed

`DESCRIPTION`'s `Collate` is **not** in C-collation order (`state.R` sits at position 52, not 86). That looked like a hand-edit, so it was tested: `roxygen2:::generate_collate()` reproduces the committed order exactly, because `R/data_cleaning.R:1` and `R/imputation.R:1` both carry `#' @include state.R`, and roxygen2's topological sort pulls `state.R` ahead of them. Running `roxygen2::roxygenise(".", roclets = "collate")` on a copy of the tree (`/tmp/sclet-collate-probe`) produced **no change** to `DESCRIPTION`.

* `Collate`: 95 entries / 95 files in `R/`, `setdiff()` in both directions empty, no duplicates.
* Consequence: **no exported name's runtime definition depends on collation order any more.** The `sclet_copilot` bug class is gone.

### 1.3 Other duplicate-definition channels: clean

* No `setGeneric()` / `setClass()` anywhere in `R/`; the only S4 method is `setMethod("subset", "SingleCellExperiment", ...)` at `R/subset.R:4` — one definition, no conflict with another package's method for that generic.
* No duplicated `S3method()` registrations in `NAMESPACE`.
* The 8 names re-exported from other packages in `R/reexport.R` (`plot_list`, `colData`, `rowData`, `colData<-`, `rowData<-`, `assay`, `plotColData`, `logNormCounts`) are **not** defined locally — no shadowing.

### 1.4 Residual risk found: **8 roxygen tags that roxygen2 cannot see**

Plain `#` comments where `#'` was intended. roxygen2 silently ignores them, so whatever they were meant to do simply does not happen:

| file:line | tag | consequence today |
|---|---|---|
| `R/clustering.R:92` | `# @importFrom SingleCellExperiment 'colLabels<-'` | inert; the import survives only because another file also imports it (`NAMESPACE:298`). Remove that other import and `colLabels<-` breaks. |
| `R/sling_plot.R:305` | `# @export` (+ `# @noRd`) | inert. The function (`add_plot_layer`) is internal and unexported — correct today. **Fixing the apostrophe would export an internal helper and suppress its Rd**, i.e. an "undocumented exported object" `R CMD check` failure. |
| `R/sling_plot.R:340` | `# @export` (+ `# @noRd`) | same (`base_plot`) |
| `R/sling_plot.R:380` | `# @export` (+ `# @noRd`) | same (`gam_analysis_SCE`) |
| `R/sling_plot.R:429` | `# @export` (+ `# @noRd`) | same (`compute_gradient`) |
| `R/ai-planning.R:1205,1206,1213` | `# @param`, `# @return`, `# @noRd` | AI layer, out of this audit's scope, but same defect class |

**Fix**: either restore the apostrophes or delete the dead tags — but for the four `sling_plot.R` helpers, *delete* them, do not repair them.

### 1.5 One stale comment (cosmetic)

`R/ai-functions.R:89-90` says "DESCRIPTION's Collate lists this file at 164 and `copilot.R` at 199". The Collate has 95 entries; `ai-functions.R` is #16 and `copilot.R` is #51. The conclusion is still right, the positions are not.

---

## 2. ARGUMENT-NAME AND SIGNATURE DRIFT

### 2.1 `@param` vs `formals()` — **0 mismatches** (after accounting for two legitimate shapes)

Checked all 254 Rd topics that carry an exported function (`.dev/tools/audit-signature-docs.R`):

* **A. a formal with no `\item` in `\arguments`: 0**
* **B. an `\item` for an argument that does not exist: 0**

Two shapes had to be handled or they produce false positives:
* a getter and its setter share one topic (`ActiveIdent` / `ActiveIdent<-` in `man/ActiveIdent.Rd`) — compared against the **union** of their formals;
* one `\item` may name several arguments (`\item{n_ctrl, n_bins}` in `man/RunCNVInference.Rd`) — item names are split on commas.

This matches the known-good state: the roxygen/Rd drift was already repaired.

### 2.2 `Run*` entry points: the *same object* has three different argument names

48 exported `Run*` functions. First argument:

| name | count | examples |
|---|---|---|
| `object` | **19** | `RunPCA` (`R/dimred.R:261`), `RunUMAP` (`R/dimred.R:344`), `RunDEtest` (`R/markers.R:227`), `RunIntegration` (`R/integration.R:31`), `RunSingleR` (`R/annotation.R:105`) |
| `sce` | **28** | `RunCellChat` (`R/cci.R:70`), `RunMilo` (`R/runMilo.R:72`), `RunVelocity` (`R/velocity.R:25`), `RunSCENIC` (`R/scenic.R:41`), `RunSlingshot` (`R/runSlingshot.R:25`) |
| `sce_spatial` | **1** | `RunSpatialDeconvolution` (`R/spatial.R:26`) |

The split is not module-driven — it cuts across modules: `R/cci.R` alone has `RunCCI(object=)` at `:21`, `RunCellChat(sce=)` at `:70`, `RunCellPhoneDB(sce=)` at `:262`, `RunNicheNet(sce=)` at `:549`. `R/dimred.R` has `RunPCA(object=)`, `RunUMAP(object=)` but `RunDiffusionMap(sce=)` (`R/dimred.R:452`).

**Fix**: pick one (`sce` is the majority in the domain modules; `object` in the accessor-facing modules) and rename; `sce_spatial` should become `sce` + `ref`.

### 2.3 The analysis-record id is called `name` when written and `id` when read

| written as `name` | read as `id` |
|---|---|
| `RunCNVInference(..., name = "cnv")` `R/malignancy-cnv.R:236` | `get_detest(object, id)` `R/markers.R:306` |
| `RunCellChat(..., name = "cellchat")` `R/cci.R:70` | `get_enrichment(object, id)` `R/Enrichment.R:145` |
| `RunDEtest(..., name = "detest")` `R/markers.R:227` | `get_velocity(object, element, id)` `R/analysis-accessors.R:812` |
| `RunEnrichment(..., name = ...)` `R/Enrichment.R:28` | `get_milo(object, element, id)` `R/analysis-accessors.R:67` |
| `RunGRN(..., name = ...)` `R/scenic.R:419` | `get_cellchat(object, element, id)` `R/analysis-accessors.R:177` |
| `RunMilo(..., name = ...)` `R/runMilo.R:72` | `get_scenic(object, element, id)` `R/analysis-accessors.R:868` |
| `RunVelocity(..., name = "velocity")` `R/velocity.R:25` | `get_cellrank(object, element, id)` `R/analysis-accessors.R:984` |
| `RunSpatialNiche(..., name = ...)` `R/spatial.R:484` | `get_spatial(object, element, id)` `R/analysis-accessors.R:1045` |
| `RunSCENIC(..., name = ...)` `R/scenic.R:41` | `get_supercell(object, element, id)` `R/analysis-accessors.R:122` |
| `RunEnrichment(object, id = NULL, ...)` `R/Enrichment.R:28` — **the odd one out: a `Run*` that already uses `id`** | `get_graph(object, name, element)` `R/analysis-accessors.R:711` — **the odd one out: an accessor that uses `name`** |
| | `has_graph(object, name)` `R/analysis-accessors.R:741` |

Two exceptions inside the same family make it worse, not better: `RunEnrichment()` already takes `id` (`R/Enrichment.R:28`) while every other `Run*` takes `name`; `get_graph()`/`has_graph()` take `name` while every other accessor takes `id`.

**Fix**: `name` → `id` in the `Run*` family and `name` → `id` in `get_graph`/`has_graph` (or the reverse, but pick one and apply it to all 48 + 55).

### 2.4 Assay / layer selection is spelled five different ways

| spelling | functions |
|---|---|
| `assay_name` | `RunCellChat` `R/cci.R:70`, `RunCellPhoneDB` `R/cci.R:262`, `RunNicheNet` `R/cci.R:549` |
| `assay_use` | `RunGRN` `R/scenic.R:419`, `RunGeneSetScoring` `R/geneset_scoring.R:30`, `RunProgramWorkflow` `R/workflows.R:181`, `RunSCENIC` `R/scenic.R:41` |
| `assay` | `RunDecontX` `R/data_cleaning.R:88`, `RunImputation` `R/imputation.R:22`, `RunPerturbationPriority` `R/priority.R:54`, `RunCNVInference` `R/malignancy-cnv.R:236`, `RunSuperCell` `R/SuperCell.R:26`, `NormalizeData` `R/preprocessing.R:74` |
| `layer` | `RunPCA` `R/dimred.R:261`, `RunIntegration` `R/integration.R:31`, `RunKNNPredict` `R/annotation.R:348`, `RunPerturbationPriority` `R/priority.R:54`, `RunSuperCell` `R/SuperCell.R:26` |
| `assay.type` | `RunSingleR` `R/annotation.R:105`, `RunKNNPredict` `R/annotation.R:348`, `RunReferenceMapping` `R/annotation.R:275` |

Note `assay.type` also uses a **dot** where the rest of the package uses `_`. `RunPCA`/`RunIntegration` use `layer` for the assay and have no `assay` argument at all, while `RunSuperCell` takes *both* `assay` and `layer`.

### 2.5 The grouping column is spelled six different ways

| spelling | functions |
|---|---|
| `group.by` | `RunKEGG` `R/pathways.R:84` |
| `group_col` | `RunPerturbationPriority` `R/priority.R:54` |
| `sample_col` / `condition_col` | `RunMilo` `R/runMilo.R:72`, `RunPerturbationPriority` `R/priority.R:54` |
| `cluster_key` | `RunCellFate` `R/cellrank.R:500`, `RunCellRank` `R/cellrank.R:30` |
| `group` | `RunCellChat` `R/cci.R:70`, `RunCellPhoneDB` `R/cci.R:262`, `RunNicheNet` `R/cci.R:549`, `RunDEtest` `R/markers.R:227` |
| `ident.1` / `ident.2` | `RunDEtest` `R/markers.R:227`, `FindMarkers` `R/markers.R:172` |

The `group` vs `group_col` split is the dangerous one: `RunPerturbationPriority(condition_col = , group_col = )` and `RunCellChat(group = )` both name a colData column, but a caller who learned one has no reason to guess the other.

### 2.6 `RunVelocity()` accepts **two spellings of the same argument**

`R/velocity.R:25-26`:

```r
RunVelocity <- function(sce, mode = ..., use.dimred = NULL,
                         subset.row = NULL, subset_row = NULL, name = "velocity", ...)
```

`R/velocity.R:45-50` then has to arbitrate:

```r
if (!is.null(subset.row) && !is.null(subset_row)) {
    stop("Please supply only one of 'subset.row' or 'subset_row'.")
}
if (is.null(subset.row) && !is.null(subset_row)) {
    subset.row <- subset_row
}
```

`subset_row` is documented as "Legacy alias of `subset.row`" (`R/velocity.R:11`). This is the only function in the package with a duplicated argument (checked mechanically: normalising `.`/`_` in every exported signature gives exactly one hit). Every other function uses one spelling — and note `RunPCA` uses `subset_row` (`R/dimred.R:261`), so the two velocity/PCA entry points spell the same idea differently.

**Fix**: drop `subset_row` from `RunVelocity()` (keep it as a deprecated wrapper only if external callers exist).

### 2.7 `use.dimred` vs `reduction`

`RunVelocity(sce, use.dimred = NULL, ...)` (`R/velocity.R:25`) vs `RunSlingshot(sce, reduction, ...)` (`R/runSlingshot.R:25`), `RunCellRank(sce, reduction, ...)` (`R/cellrank.R:30`), `RunPCA(object, layer, ...)`. Same concept, two names, one of them with a dot.

### 2.8 Two parallel state-writing APIs with confusingly similar names

`R/state.R` exposes **three** closely-named writers:

| function | line | what it is |
|---|---|---|
| `sclet_set_analysis(object, name, value)` | `R/state.R:629` | writes the **legacy** `state$analyses[[name]]` slot, arbitrary shape |
| `sclet_set_state_record(object, type, value, id, active)` | `R/state.R:784` | the real record writer |
| `sclet_set_analysis_state(object, type, id, method, inputs, artifacts, params, summary, active)` | `R/state.R:825` | a thin wrapper over the previous one that fills in the canonical five fields |

27 core functions use the wrapper, 10 use the record writer directly, and `sclet_set_analysis` is used by 30+ call sites for the legacy slot. `sclet_set_analysis` vs `sclet_set_analysis_state` is a one-token difference for two completely different things.

---

## 3. SILENT ERROR PATHS

Ranked by how badly a user can be misled.

### 3.1 (rank 1) Eleven exported accessors return the **wrong record** for an id that does not exist

`R/analysis-accessors.R:514-552`:

```r
sclet_normalize_analysis_record <- function(type, id, record = NULL, analysis = NULL, default_method = NULL) {
    if (is.null(record) && is.null(analysis)) {
        return(NULL)
    }
    result <- record
    if (is.null(result)) { result <- list() }
    if (!is.null(analysis)) { for (nm in names(analysis)) { ... } }   # <-- legacy slot wins
    if (is.null(result$method) && !is.null(default_method)) {
        result$method <- default_method                               # <-- method invented
    }
    ...
    if (is.null(result$inputs))   { result$inputs   <- list() }        # <-- empty is indistinguishable
    if (is.null(result$artifacts)){ result$artifacts<- list() }
    if (is.null(result$summary))  { result$summary  <- list() }
```

Every accessor passes the legacy slot as the fallback, so when the requested record id is absent the accessor silently substitutes the legacy analysis and stamps it with a **hard-coded `default_method`** — the name of a backend that may never have run.

Runtime proof (E1, `.dev/tools/audit-core-runtime-evidence.R`): on an object that has a trajectory record `run1` and a legacy `analyses$trajectory` slot, `get_trajectory(sce, id = "typo")` returns the legacy record with `method = "slingshot"`; `has_trajectory(sce, id = "typo")` returns `TRUE`.

Eleven exported accessors are affected — the fallback call sites and the `default_method` each one stamps in:

| accessor | fallback site | `default_method` it reports |
|---|---|---|
| `get_trajectory` | `R/analysis-accessors.R:28` | `slingshot` |
| `get_milo` | `:83` | `miloR` |
| `get_supercell` | `:138` | `SuperCell::SCimplify` |
| `get_cellchat` | `:193` | `CellChat` |
| `get_cellphonedb` | `:248` | `CellPhoneDB` |
| `get_nichenet` | `:303` | `NicheNet` |
| `get_batch` | `:570` | (no default, still returns the wrong record) |
| `get_velocity` | `:823` | `velociraptor::scvelo` |
| `get_scenic` | `:879` | `pySCENIC` |
| `get_geneset_scoring` | `:937` | (none) |
| `get_cellrank` | `:1001` | `CellRank` |
| `get_perturbation` | `:1170` | (none) |

These are not theoretical: the legacy slot is written by the corresponding `Run*` functions (`R/velocity.R:93`, `R/scenic.R:272`, `R/cellrank.R:149`, `R/runMilo.R:231`, `R/runSlingshot.R:86`, `R/cci.R:158/467/693`, `R/SuperCell.R:101`, `R/geneset_scoring.R:136`), so after a real run the slot is populated and the fallback fires on every typo.

A user comparing two runs — `get_velocity(sce1, id = "run2")` vs `get_velocity(sce2, id = "run2")` — gets *the same run-1 record from both objects*, with `method = "velociraptor::scvelo"`, and no indication that the id was wrong.

**Fix**: in `sclet_normalize_analysis_record()`, take the legacy `analysis` only when no `id` was requested; when an `id` was requested and no record matched, return `NULL`.

### 3.2 (rank 2) A missing active identity produces base R's `"argument is of length zero"`

`R/utils-cci.R:12-28`:

```r
sclet_resolve_group_column <- function(object, meta, group = NULL) {
    if (is.null(group)) { group <- ActiveIdent(object) }      # NULL when no active identity
    if (identical(group, "colLabels")) { ... }
    else if (!group %in% colnames(meta)) {                     # <-- R/utils-cci.R:24
        stop(sprintf("Group column '%s' not found in colData.", group))
    }
```

With `group = NULL`, `!NULL %in% colnames(meta)` is `logical(0)`, and `if (logical(0))` is base R's *argument is of length zero*. The carefully written message on line 25 is never reached.

Runtime proof (E2): `RunCellChat(sce)` and `RunCellPhoneDB(sce)` on a plain `mockSCE()` both fail with `Error: argument is of length zero` — no mention of groups, identities or `colData`. The same resolver is used by all three CCI entry points (`R/cci.R:107`, `:288`, `:584`).

Note `check_communication_readiness()` would have said exactly the right thing ("no group column is resolvable … set an active identity or add a cluster column"), but nothing calls it (§5).

**Fix**: `if (is.null(group)) stop("No group column: set an active identity (ActiveIdent) or pass `group = <colData column>`.")`.

### 3.3 (rank 3) `RunCCI()` forwards an argument one of its three targets cannot accept

`R/cci.R:21-35`:

```r
RunCCI <- function(object, method = c("CellChat", "CellPhoneDB", "NicheNet"), keep_object = FALSE, ...) {
    method <- match.arg(method)
    if (method == "CellChat") {
        args <- list(...)
        if (is.null(args$return)) { args$return <- "sce" }
        return(do.call(RunCellChat, c(list(sce = object, keep_object = keep_object), args)))
    } else if (method == "CellPhoneDB") {
        return(RunCellPhoneDB(object, ...))      # keep_object NOT forwarded
    } else if (method == "NicheNet") {
        return(RunNicheNet(object, ...))          # keep_object NOT forwarded
    }
}
```

`keep_object` is forwarded **only** to the CellChat branch, but it is a parameter of `RunCCI` itself, so `RunCCI(sce, method = "NicheNet", keep_object = TRUE)` passes it into `RunNicheNet()`. `formals(RunNicheNet)` has **no `keep_object` and no `...`** (`R/cci.R:549`), so R raises `unused argument (keep_object = TRUE)` — after the analysis has already started. (`RunCellPhoneDB` does have `...` (`R/cci.R:262`), so there the stray flag is swallowed by `...` and forwarded to the Python backend instead.)

The two branches also dispatch on different names: `RunCCI` says `method`, `RunCellChat` says `db_item`/`type`; `RunCCI` takes `object`, `RunCellChat` takes `sce`.

**Fix**: forward `keep_object` to all three branches and add `keep_object` to `RunNicheNet()` (or drop it from `RunCCI()`).

### 3.4 (rank 4) An analysis that produced nothing is indistinguishable from one that never ran

`sclet_normalize_analysis_record()` fills `inputs`, `artifacts` and `summary` with **empty lists** when they are missing (`R/analysis-accessors.R:541-549`). A record for an analysis that ran and found nothing, and a record for an analysis that ran but recorded no summary, are byte-identical to the caller:

* `get_trajectory(sce)$summary$n_clusters` → `NULL` in both cases (runtime proof in E2-style check: `summary` has length 0).
* `get_enrichment(sce)$artifacts$result` → `NULL` whether the enrichment produced an empty table or produced no record at all.
* `has_*()` is `!is.null(get_*)` (`R/analysis-accessors.R:441` etc.), so a record with an empty summary reports `has_* == TRUE`.

The package has the vocabulary to fix this (`summary$status = "unavailable"` is used by the CNV path, `R/malignancy-cnv.R:681`) but the accessors do not surface it.

**Fix**: keep `NULL` for `summary`/`artifacts` instead of substituting `list()`, or standardise on an explicit `status` field that the accessors expose.

### 3.5 (rank 5) `RunCNVInference()` degrades to "unavailable" — this one is **correct** (a positive example)

`R/malignancy-cnv.R:303-312` is the only core call site of a readiness gate:

```r
readiness <- check_cnv_readiness(object, reference_cells = ..., reference_col = ...,
                                 gene_order = gene_order, species = species)
if (!isTRUE(readiness$ready)) {
    return(sclet_cnv_record_unavailable(object, name, readiness))
}
```

`sclet_cnv_record_unavailable()` (`R/malignancy-cnv.R:648`) writes NA columns, a `status = "unavailable"` record with the reasons, **and emits a `warning()`** (`R/malignancy-cnv.R:687-693`). This is the model the other 13 gates should follow.

### 3.6 Lower-severity silent paths (listed, not ranked)

| file:line | path | note |
|---|---|---|
| `R/pseudobulk.R:44-47` | `scrapper::aggregateAcrossCells.se()` error → `NULL` → silent fallback to `scuttle::aggregateAcrossCells()` | legitimate fallback; the two backends can aggregate differently, and nothing tells the user which ran |
| `R/pseudobulk.R:90-100` | `DESeq2::DESeqDataSet()` error → `NULL` → rebuild via `DESeqDataSetFromMatrix` | reasonable fallback, but the original error is discarded |
| `R/regvelo.R:435-441` | Python call error → stored in `last_error` and `NULL` returned | the error is captured, then presumably re-raised — acceptable |
| `R/cellrank.R:211-214` | `sclet_velocity_assay_name()` error → `NULL` → the "missing assay" check at `:216-225` is **skipped** | the guard that would have caught a missing assay is bypassed when the lookup itself fails |
| `R/plot-cci.R:228-231` | `strwidth()` error → `0.2` | benign (a layout default) |
| `R/annotation.R:31-32` | `getExportedValue("celldex", name)` error → `NULL` → `return(NULL)` | a reference dataset that cannot be loaded becomes "no reference"; callers cannot tell it from "reference not requested" |
| `R/pipeline_summary.R:16` | `return(invisible(NULL))` when nothing is recorded | fine — the function prints a message first |
| `R/basilisk.R:362` | `return(invisible(NULL))` when there is no bundled lib dir | fine — documented no-op |

---

## 4. STATE-CONTRACT ADOPTION

### 4.1 Method

For each of the **65 exported core analysis functions** (`Run*`, `Find*`, `Call*`, `Classify*`, `Aggregate*`, `Compare*`, `Normalize*`, `Scale*`, `subset_*`, `Read*`/`Write*`, `Batch*`), the body was parsed and searched for:

* `sclet_set_analysis_state(` — the canonical writer (`R/state.R:825`);
* `sclet_set_state_record(` — the underlying writer (`R/state.R:784`), used directly by 10 functions;
* `check_*_readiness(` — any readiness gate;
* `stop(` / `cli::cli_abort(` — does the function fail loudly.

Both record writers produce the same five-field record (`method` / `inputs` / `artifacts` / `params` / `summary` — `R/state.R:825-849` wraps `R/state.R:784`), so a function is counted as "adopting the contract" if it calls either.

### 4.2 Headline numbers

| | count |
|---|---|
| exported core analysis functions | **65** |
| register a state record (`sclet_set_analysis_state`) | 27 |
| register a state record (`sclet_set_state_record` directly) | 10 |
| **register a state record (either way)** | **37 (57 %)** |
| **write a result onto the object with NO state record** | **28 (43 %)** |
| call a readiness gate | **1** |
| `stop()` somewhere in the body | 52 |

### 4.3 Per-module table (inconsistent modules marked)

| module | n | state | no state | gate |
|---|---|---|---|---|
| annotation | 4 | 3 | 1 | 0 |
| basic-workflow | 1 | 0 | **1** | 0 |
| batch | 1 | 1 | 0 | 0 |
| cci | 4 | 3 | 1 | 0 |
| celloracle | 1 | 0 | 1 | 0 |
| cellrank | 2 | 0 | **2** | 0 |
| clustering | 2 | 2 | 0 | 0 |
| data_cleaning | 2 | 2 | 0 | 0 |
| dimred | 3 | 2 | 1 | 0 |
| Enrichment | 1 | 1 | 0 | 0 |
| geneset_scoring | 1 | 0 | 1 | 0 |
| imputation | 1 | 1 | 0 | 0 |
| integration | 1 | 1 | 0 | 0 |
| interactive | 1 | 0 | **1** | 0 |
| interop | 2 | 0 | **2** | 0 |
| io | 1 | 0 | **1** | 0 |
| malignancy | 2 | 0 | **2** | 0 |
| malignancy-cnv | 1 | 0 | 1 | **1** |
| malignancy-combine | 1 | 0 | 1 | 0 |
| markers | 3 | 1 | **2** | 0 |
| multimodal | 1 | 0 | **1** | 0 |
| pathways | 2 | 0 | **2** | 0 |
| phenotype | 1 | 1 | 0 | 0 |
| phenotype-evidence | 1 | 0 | 1 | 0 |
| preprocessing | 3 | 1 | **2** | 0 |
| priority | 3 | 2 | 1 | 0 |
| pseudobulk | 2 | 0 | **2** | 0 |
| regvelo | 1 | 0 | 1 | 0 |
| runMilo | 1 | 1 | 0 | 0 |
| runSlingshot | 1 | 1 | 0 | 0 |
| scenic | 2 | 0 | **2** | 0 |
| spatial | 4 | 3 | 1 | 0 |
| subset | 2 | 0 | **2** | 0 |
| SuperCell | 1 | 1 | 0 | 0 |
| velocity | 1 | 0 | 1 | 0 |
| workflow | 1 | 0 | **1** | 0 |
| workflows | 3 | 0 | **3** | 0 |

**Modules that never adopt the contract** (0 of every function registers a record): `basic-workflow`, `cellrank`, `interactive`, `interop`, `io`, `malignancy`, `multimodal`, `pathways`, `pseudobulk`, `scenic`, `subset`, `velocity`, `workflow`, `workflows`.

**Modules that always adopt it**: `batch`, `clustering`, `data_cleaning`, `Enrichment`, `imputation`, `integration`, `phenotype`, `runMilo`, `runSlingshot`, `SuperCell`.

**Modules that are internally inconsistent** (some functions record, some do not): `annotation` (3/4), `cci` (3/4), `cellrank` (0/2), `dimred` (2/3), `markers` (1/3), `malignancy` (0/2), `malignancy-combine` (0/1), `preprocessing` (1/3), `priority` (2/3), `scenic` (0/2), `spatial` (3/4), `velocity` (0/1).

### 4.4 The 28 functions that write a result onto the object with no state record

The consequential ones (the rest are thin dispatchers or I/O helpers that delegate to functions which *do* record):

| function | file:line | writes onto the object? |
|---|---|---|
| `RunDiffusionMap` | `R/dimred.R:452` | yes — a reduction |
| `FindVariableFeatures` | `R/preprocessing.R:165` | yes — the HVG selection |
| `ScaleData` | `R/preprocessing.R:729` | yes — a scaled assay |
| `RunTrajectoryWorkflow` | `R/workflows.R:32` | yes — trajectory result |
| `RunProgramWorkflow` | `R/workflows.R:181` | yes — program scoring |
| `RunReferenceWorkflow` | `R/workflows.R:304` | yes — mapping result |
| `RunSpatialWorkflow` | `R/spatial.R:628` | yes |
| `RunMultimodalWorkflow` | `R/multimodal.R:175` | yes |
| `RunStatePriorityWorkflow` | `R/priority.R:540` | yes |
| `RunGRN` | `R/scenic.R:419` | yes — regulatory network |
| `RunCellFate` | `R/cellrank.R:500` | yes |
| `RunStandardPipeline` | `R/workflow.R:19` | yes (delegates to functions that record) |
| `FindAllMarkers` | `R/markers.R:30` | returns a table — not an object mutation, but no record either |
| `FindMarkers` | `R/markers.R:172` | same |
| `FindMarkerPathways` | `R/pathways.R:25` | same |
| `RunKEGG` | `R/pathways.R:84` | same |
| `AggregateExpression` | `R/pseudobulk.R:26` | returns an SE — no record |
| `FindMarkers_pseudobulk` | `R/pseudobulk.R:85` | returns a table |
| `FindMalignantSeedGenes` | `R/malignancy.R:81` | returns a table |
| `ClassifyMalignantEpithelium` | `R/malignancy.R:374` | yes (writes `type = "annotation"`, see §4.5) |
| `CallMalignantCells` | `R/malignancy-combine.R:84` | yes (writes `type = "annotation"`) |
| `CallPhenotypeCells` | `R/phenotype-evidence.R:260` | yes (writes `type = "annotation"`) |
| `subset_cell` / `subset_feature` | `R/subset.R:61` / `:26` | return a subset — arguably out of scope for a state record |
| `Read10X` / `ReadH5AD` / `WriteH5AD` | `R/io.R:18`, `R/interop.R:16`, `:40` | I/O — out of scope |

**The pattern**: `FindAllMarkers()`/`FindMarkers()`/`FindMarkerPathways()`/`RunKEGG()` return a *table* and record nothing, so a marker table computed once is indistinguishable from one that was never computed — and `RunEnrichment()` reads it back through `get_detest()`, which only exists for the `RunDEtest()` path. `RunDEtest()` **does** record (`R/markers.R:266`), so the record exists; but `FindAllMarkers()` called directly does not.

### 4.5 Six analyses share one state bucket: `type = "annotation"`

`sclet_state_types()` (`R/state.R:686`) has 23 registered types. There is **no** `cnv`, `malignancy` or `phenotype` type, so three unrelated modules write into the `annotation` bucket:

| writer | file:line | what it is |
|---|---|---|
| `RunSingleR` / `RunKNNPredict` / `RunSymphonyMapping` | `R/annotation.R:193` | label transfer |
| `CallPhenotypeCells` | `R/phenotype-evidence.R:557` | phenotype assignment |
| `ClassifyMalignantEpithelium` | `R/malignancy.R:735` | malignant-cell classification |
| `RunCNVInference` | `R/malignancy-cnv.R:427` and `:673` | CNV inference |
| `CallMalignantCells` | `R/malignancy-combine.R:173` | malignant-cell calling |

Because there is exactly one "active" slot per type (`sclet_get_active_state(object, "annotation")`), **whichever of the six ran last becomes the active annotation**, and `get_annotation()` (`R/analysis-accessors.R:402`) returns it.

Runtime proof (E6 in the evidence script):

```
get_annotation(sce) -> method: inferCNV per-cell copy-number burden  id: cnv
get_annotation(sce, id='cnv')$summary$n_issues: 2
has_annotation(sce): TRUE
sclet_get_state_records(sce,'annotation') ids: cnv, annotation
```

So a user who ran a CNV analysis and then asked `get_annotation(sce)` gets a CNV record, and `has_annotation(sce)` answers `TRUE`.

**Fix**: register `cnv`, `malignancy` and `phenotype` in `sclet_state_types()` (`R/state.R:686`) and give each module its own bucket.

### 4.6 Two writers, one contract — the canonical five fields are not always present

27 functions use `sclet_set_analysis_state()` (which guarantees `method`/`inputs`/`artifacts`/`params`/`summary`). 10 use `sclet_set_state_record()` directly, and those must supply the shape themselves. `RunVelocity()` (`R/velocity.R:93-127`) does supply all five — but it also adds `results`, `reductions` and `coldata_fields` as siblings, so the record has eight top-level fields, three of which are outside the contract. `RunCellRank()`, `RunSCENIC()`, `RunSuperCell()`, `RunGeneSetScoring()`, `RunRegVelo()`, `RunCellChat()`, `RunCellPhoneDB()`, `RunNicheNet()`, `RunTrajectoryWorkflow()` and `RunProgramWorkflow()` are the ten.

---

## 5. THE READINESS GATES

### 5.1 Inventory — 14 exported gates, all reachable from the exported API

All 14 are in `NAMESPACE` (`NAMESPACE:142-157`):

```
check_annotation_readiness      check_fate_readiness
check_cnv_readiness             check_integration_readiness
check_communication_readiness   check_milo_readiness
check_de_readiness              check_multimodal_readiness
check_enrichment_readiness      check_rare_cell_readiness
check_fate_drivers_readiness    check_spatial_readiness
                                check_trajectory_readiness
                                check_velocity_readiness
```

(plus `check_claim_ceiling` and `check_phenotype_provider`, which are not readiness gates.)

### 5.2 Call graph — **only one of them is enforced by the analysis layer**

From E3 (all of `R/`, roxygen `@examples` and the definitions excluded):

| gate | call sites |
|---|---|
| `check_cnv_readiness` | **`R/malignancy-cnv.R:303`** — inside `RunCNVInference()` ✅ |
| `check_communication_readiness` | `R/ai-communication.R:162` (AI layer only) |
| `check_integration_readiness` | `R/ai-investigate.R:26`, `R/ai-investigate.R:55` (AI layer only) |
| `check_multimodal_readiness` | `R/ai-registry-spatial-multimodal.R:166` (AI layer only) |
| `check_rare_cell_readiness` | `R/ai-context.R:556` (AI layer only) |
| `check_spatial_readiness` | `R/ai-registry-spatial-multimodal.R:28` (AI layer only) |
| `check_annotation_readiness` | **no call site anywhere** |
| `check_de_readiness` | **no call site anywhere** |
| `check_enrichment_readiness` | **no call site anywhere** |
| `check_fate_drivers_readiness` | **no call site anywhere** |
| `check_fate_readiness` | **no call site anywhere** |
| `check_milo_readiness` | **no call site anywhere** |
| `check_trajectory_readiness` | **no call site anywhere** |
| `check_velocity_readiness` | **no call site anywhere** |

**13 of 14 gates are advisory-only.** Every analysis function bypasses its own gate:

| gate | the analysis it should gate |
|---|---|
| `check_velocity_readiness` | `RunVelocity()` `R/velocity.R:25` |
| `check_trajectory_readiness` | `RunSlingshot()` `R/runSlingshot.R:25`, `RunTrajectoryWorkflow()` `R/workflows.R:32` |
| `check_fate_readiness` / `check_fate_drivers_readiness` | `RunCellFate()` `R/cellrank.R:500`, `RunCellRank()` `R/cellrank.R:30` |
| `check_de_readiness` | `RunDEtest()` `R/markers.R:227` |
| `check_enrichment_readiness` | `RunEnrichment()` `R/Enrichment.R:28` |
| `check_integration_readiness` | `RunIntegration()` `R/integration.R:31` |
| `check_milo_readiness` | `RunMilo()` `R/runMilo.R:72` |
| `check_spatial_readiness` | `RunSpatialDeconvolution()` `R/spatial.R:26`, `RunSpatialColocalization()` `:193`, `RunSpatialNiche()` `:484` |
| `check_communication_readiness` | `RunCellChat()` `R/cci.R:70`, `RunCellPhoneDB()` `:262`, `RunNicheNet()` `:549` |
| `check_annotation_readiness` | `RunSingleR()` `R/annotation.R:105`, `RunKNNPredict()` `:348`, `RunSymphonyMapping()` `:509` |
| `check_rare_cell_readiness` | `RunRareCellDetection()` `R/priority.R:288` |
| `check_multimodal_readiness` | `RunMultimodalWorkflow()` `R/multimodal.R:175` |

The gates are well-written diagnostics — they are simply never consulted by the code they describe. `check_communication_readiness()` would have produced "no group column is resolvable" for the exact case that currently dies with "argument is of length zero" (§3.2).

### 5.3 Do the gates check the thing they claim?

**Package-dependency coverage** (E4):

| gate | packages it checks | what the gated analysis actually requires |
|---|---|---|
| `check_de_readiness` | **limma** (`R/ai-interpretation-gaps.R:275`) | **`presto`** (`R/markers.R:387`) for the default marker path; limma only for the covariate/pseudobulk path (`R/markers.R:476`) |
| `check_enrichment_readiness` | clusterProfiler + `org.Hs.eg.db` | clusterProfiler ✅ |
| `check_milo_readiness` | miloR | miloR ✅ |
| `check_cnv_readiness` | infercnv + rjags | infercnv ✅ |
| `check_communication_readiness` | CellChat / basilisk / nichenetr | CellChat ✅ |
| `check_velocity_readiness` | — | velociraptor (`R/velocity.R:27`) |
| `check_trajectory_readiness` | — | slingshot (`R/runSlingshot.R:36`) |
| `check_integration_readiness` | — | harmony / basilisk / zellkonverter (`R/integration.R`) |
| `check_spatial_readiness` | — | basilisk (`R/spatial.R:26`), SVP (`R/spatial.R:193`, `:484`) |
| `check_annotation_readiness` | — | SingleR + celldex (`R/annotation.R`), BiocNeighbors (`:348`), symphony (`:509`) |
| `check_rare_cell_readiness` | — | BiocNeighbors (`R/priority.R:288`) |
| `check_fate_readiness` | — | basilisk + reticulate (`R/cellrank.R:30`) |

**9 of 14 gates check no package at all**, although the analysis they gate hard-`stop()`s with an install message when the Suggests package is missing. The gate would report "ready" and the analysis would then refuse to run — the gate's whole value is in the dependency check it omits.

### 5.4 The `check_de_readiness` package mismatch, in detail

`check_de_readiness()` (`R/ai-interpretation-gaps.R:234`) computes:

```r
de_package_available <- base::requireNamespace("limma", quietly = TRUE)   # :275
...
if (!de_package_available) {
    questions <- c(questions, "The 'limma' package is not installed, which the marker test relies on.")
}
ready <- group_resolvable && enough_groups && enough_cells && de_package_available   # :309
```

But the marker test it describes is `RunDEtest()` → `FindAllMarkers()`/`FindMarkers()` → `FindMarkers_Presto()`, which is:

```r
if (!requireNamespace("presto", quietly = TRUE)) {          # R/markers.R:387
    stop("Package 'presto' is required for marker testing. ...")
}
```

`limma` is only reached by the *covariate* path (`sclet_run_condition_de()`, `R/markers.R:476`). So on a machine with limma but no presto the gate says `ready_for_diagnostic` and the analysis dies; on a machine with presto but no limma the gate says `not_ready` and blocks a marker test that would have worked.

Two smaller points in the same gate:
* `counts_available` is computed and reported as a question (`R/ai-interpretation-gaps.R:239`) but is **not** part of `ready` (`:309`) — a `counts` assay is asked for but never required, while `RunDEtest()` resolves the assay through `resolve_marker_assay()` which accepts `logcounts`.
* `check_de_readiness(object, group = NULL)` takes a `group` argument, but `RunDEtest()`'s grouping argument is also called `group` and means something different (a `colData` column vs. an identity level). The gate's `group` is a column name; `RunDEtest(group = ...)` is a column name too — but `ident.1`/`ident.2` are identity *levels*. The gate never considers `ident.1`/`ident.2`.

### 5.5 Gates that are genuinely good

For contrast, `check_enrichment_readiness()` (`R/ai-enrichment.R:50`) is the model: it checks clusterProfiler **and** the org.Db, checks that the marker table actually carries a `cluster` column (because `RunEnrichment()` groups markers by cluster and a pooled table would compare one set against itself), checks that a cluster column is resolvable, and states the reason in plain language. It is simply never called.

---

## 6. Looked for and did NOT find

Recording these because a negative result is information.

* **No duplicate definitions of exported names anywhere.** Two independent methods (R parser over all 95 files, 757 top-level names; and a column-0 assignment grep, 764 hits) both find zero. No exported name is defined twice in one file, none in two files. `sclet_copilot()` exists once, at `R/copilot.R:43`.
* **`DESCRIPTION`'s `Collate` is in sync with roxygen2**, even though it is not in C-collation order — `state.R` is pulled forward by `#' @include state.R` in `R/data_cleaning.R:1` and `R/imputation.R:1`. Running roxygen2's collate roclet on a copy produced no change.
* **No duplicated `S3method()` registrations**, no `setGeneric()`/`setClass()` in `R/`, and no name in `R/reexport.R` that is also defined locally.
* **No `@param` / `formals()` mismatch** in any of the 254 Rd topics that carry an exported function — 0 undocumented formals, 0 phantom documented arguments (after handling shared getter/setter topics and multi-name `\item{}`s).
* **No `Rd` file is missing for an exported name**, and no exported name lacks an Rd alias (category E of `audit-signature-docs.R` = 0).
* **No analysis function returns `invisible(NULL)` on a failure path** other than the two benign cases (`R/basilisk.R:362`, a documented no-op; `R/pipeline_summary.R:16`, a print-only function that messages first).
* **No `tryCatch` handler in the core layer swallows an error on a result-producing path without recording it.** All 29 `tryCatch` calls outside `R/ai-*.R` were read; the only ones that discard an error are the `scrapper`/`DESeq2` fallbacks in `R/pseudobulk.R` (both re-raise via a different route) and `R/cellrank.R:211` (whose guard is then skipped — noted in §3.6).
* **The one enforced gate degrades loudly.** `RunCNVInference()` → `sclet_cnv_record_unavailable()` writes NA columns, a `status = "unavailable"` record with the reasons, and a `warning()` (`R/malignancy-cnv.R:683`). This is the behaviour the other 13 should have.
* **No S4 method conflicts** with other packages: `setMethod("subset", "SingleCellExperiment", ...)` (`R/subset.R:4`) is the package's only S4 method.
* **Not verified** (stated as unverified rather than assumed): the behaviour of `RunNicheNet()` and `RunCellPhoneDB()` under a real nichenetr / CellPhoneDB run — `nichenetr` is not installed on this machine, so the `unused argument` finding in §3.3 rests on `formals(RunNicheNet)` having neither `keep_object` nor `...` (verified) rather than on an executed run.

---

## 7. Suggested order of work

1. `sclet_normalize_analysis_record()` — one function, fixes 11 accessors (§3.1). Highest value per line changed.
2. `sclet_resolve_group_column()` — one `if` clause, fixes 3 entry points (§3.2).
3. `R/utils-cci.R` + the gates: call each `check_*_readiness()` at the top of its `Run*` and abort on `not_ready` (§5.2); fix `check_de_readiness`'s package check (§5.4).
4. Add `cnv` / `malignancy` / `phenotype` state types and move the three modules out of the `annotation` bucket (§4.5).
5. Rename the drifted argument families: `object`/`sce` (§2.2), `name`/`id` (§2.3), `assay_name`/`assay_use`/`assay`/`layer`/`assay.type` (§2.4), the six group-column spellings (§2.5). This is a breaking API change — do it in one release with a deprecation shim.
6. Register the 28 missing state records (§4.4).
7. Delete the 8 inert roxygen tags (§1.4) and fix the stale Collate positions in the comment at `R/ai-functions.R:89`.

---

## 附录 B — 整改状态（2026-10-22，Lead 复核后）

本报告的 5 条结论已由 Lead 独立复现；其中 4 条已修并提交，1 条记录待决。

| 发现 | 状态 | commit |
|---|---|---|
| #1 11 个 accessor 对不存在的 id 返回别的记录 | **已修** | `42bfed0` |
| #5 缺 active ident 时报 "argument is of length zero" | **已修** | `42bfed0` |
| #3 `check_de_readiness` 卡 limma 而非 presto | **已修** | `ea1468d` |
| 状态契约采用率 37/65 | 记录待决 | — |
| #4 六个分析共用 `type="annotation"` 状态桶 | **记录待决**（见下） | — |
| #2 14 个闸门里 13 个从不被它所守的分析调用 | **记录待决**（见下） | — |

### #4 为什么没当场修

`get_annotation()` 会返回一条 CNV 记录（运行时确认：`id="cnv_run_1",
method="infercnv::inferCNV"`）。修法是把 `cnv` / `malignancy` / `phenotype`
加进 `sclet_state_types()` 并改 6 个写入点。

没有当场做的原因：这是**状态语义变更**，不是局部修复。改完之后，
已经用 `type="annotation"` 记下 CNV 记录的对象将不再被 `get_annotation()`
看到——这是对的，但同时也意味着 `required_states = list(annotation = ...)`
的 AI action 在那些对象上会判定为缺状态。需要先决定：旧的 annotation 类型
记录要不要迁移，以及 AI 层的 required_states 要不要跟着扩。

### #2 为什么没当场修

14 个闸门里 13 个只被 `R/ai-*.R` 调用，分析函数自己不执行自己的闸门。
修法是"每个 `Run*` 开头调自己的闸门，not_ready 就 abort"。

没有当场做的原因：这会**改变已导出函数的失败行为**——现在能跑通的调用会
开始报错。哪些闸门的约束是硬性的（cnv 的 infercnv 依赖）、哪些只是建议性的
（de 的细胞数下限），必须先逐个定性，否则会把本来能出结果的调用挡掉。
