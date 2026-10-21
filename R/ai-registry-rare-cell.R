# sclet_ai_registry_rare_cell(): built-in AI actions for doublet and rare-cell detection.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# Rare-cell results reach the evidence layer as counts, not identities.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_rare_cell <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("rare_cell" %in% include) {
    actions <- c(actions, list(
        run_doublet_detection = AIAction(
            name = "run_doublet_detection",
            description = paste(
                "Score doublets with scDblFinder, store per-cell doublet class/score in colData, and",
                "register one bounded evidence node (aggregate counts and score quantiles) so the call",
                "is citable from the ledger. The class label is one classifier's estimate, not a",
                "measurement of doubleness, and no cell is removed or filtered.",
                "No per-cell score or class vector enters the evidence node."
            ),
            handler = function(object, params) {
                scored <- RunDoubletFinder(object)
                # The preprocess record id is fixed by RunDoubletFinder, so
                # the evidence id is stable and a re-run degrades to a no-op
                # inside the recorder rather than duplicating or erroring.
                sclet_ai_record_doublet_evidence(scored, name = "doublet_finder")
            },
            input_schema = list(),
            prerequisites = function(object, params, planned = NULL) {
                if (!"counts" %in% SummarizedExperiment::assayNames(object)) {
                    return("counts_assay_missing: run_doublet_detection requires a 'counts' assay in the object")
                }
                if (!base::requireNamespace("scDblFinder", quietly = TRUE)) {
                    return("optional_package_missing: scDblFinder is required for run_doublet_detection; install via BiocManager::install('scDblFinder')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(preprocess = "any", ai_evidence = "any"),
                note = paste(
                    "Per-cell doublet calls are added to colData. No cell is removed or filtered by",
                    "this action. The registered evidence carries only the number of cells, the number",
                    "called doublet, that fraction, and score quantiles; a single classifier's call is",
                    "recorded at claim level 'associated' because nothing independent corroborates it."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("preprocess", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE
        ),
        run_rare_cell_detection = AIAction(
            name = "run_rare_cell_detection",
            description = "Label density-based rare populations and register one bounded evidence node per small population, graded by the number of independent signals available.",
            handler = function(object, params) {
                reduction <- params$reduction %||% "PCA"
                name <- params$name %||% "rareq"
                dims <- if (is.null(params$dims)) 1:20 else as.integer(unlist(params$dims, use.names = FALSE))
                rare_threshold <- params$rare_threshold %||% 10
                result <- RunRareCellDetection(
                    object,
                    method = "density",
                    reduction = reduction,
                    dims = dims,
                    k = params$k %||% 20,
                    q_threshold = params$q_threshold %||% 0.25,
                    rare_threshold = rare_threshold,
                    min_cluster_size = params$min_cluster_size %||% 3L,
                    name = name
                )
                sclet_ai_record_rare_cell_evidence(result, name, rare_threshold)
            },
            input_schema = list(
                reduction = "character",
                dims = "integer_vector",
                k = "integer",
                q_threshold = "number",
                rare_threshold = "integer",
                # Density artefacts below this size are dissolved into
                # their densest neighbour instead of being reported as
                # populations (the 2419-for-2700-cells defect).
                min_cluster_size = "integer",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                reduction <- params$reduction %||% "PCA"
                if (!is.character(reduction) || length(reduction) != 1L || !nzchar(reduction)) {
                    return("reduction_invalid: 'reduction' must be a single non-empty reduction name")
                }
                # planned$reductions first (deferred projection): a valid
                # plan runs run_pca in an earlier step and the reduction
                # only exists in the planned view at validation time.
                # Reading reducedDimNames() alone made any plan putting
                # run_rare_cell_detection after run_pca fail with
                # "reduction 'PCA' is not available" -- found offline by
                # the doublet/rare benchmark dry run 2026-10-14 before a
                # billed call, with nothing wrong in the plan.
                available <- unique(c(
                    as.character(planned$reductions %||% character()),
                    tryCatch(SingleCellExperiment::reducedDimNames(object),
                        error = function(e) character())))
                if (!reduction %in% available) {
                    return(paste0(
                        "reduction_missing: reduction '", reduction,
                        "' is not available; run the dimensionality reduction first. Available reductions: ",
                        paste(available, collapse = ", ")
                    ))
                }
                if (!base::requireNamespace("BiocNeighbors", quietly = TRUE)) {
                    return("optional_package_missing: BiocNeighbors is required for density-based rare-cell detection; install via BiocManager::install('BiocNeighbors')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(rare_cells = "any"),
                note = paste(
                    "Rare populations are only labeled (colData rare_cluster). No cell is removed, filtered or merged.",
                    "Evidence claim_level is graded by the number of independent signals: zero signals records no",
                    "evidence at all, one signal records 'associated' with low_confidence = TRUE, and two or more",
                    "signals record 'consistent_with'. A signal counts only when it is informative for that",
                    "population: markers must be attributable to the population itself, QC columns must vary and be",
                    "observed on both sides, doublet calls must cover the population and the rest, and sample",
                    "replication needs at least two sample labels. Cluster size alone never produces evidence."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("rare_cells", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE
        )
    ))
    }

    actions
}
