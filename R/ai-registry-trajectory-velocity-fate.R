# sclet_ai_registry_trajectory_velocity_fate(): built-in AI actions for trajectory, velocity and fate analysis.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# Three groups, one file, in the order the registry has always emitted them.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_trajectory_velocity_fate <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("trajectory" %in% include) {
    actions <- c(actions, list(
        run_trajectory = AIAction(
            name = "run_trajectory",
            description = "Infer a slingshot trajectory from an explicitly supplied start cluster and record a bounded aggregate evidence node.",
            handler = function(object, params) {
                name <- params$name %||% "slingshot"
                start_cluster <- sclet_ai_resolve_ident_value(object, params$start_cluster)
                result <- RunSlingshot(
                    object,
                    group = params$group,
                    reduction = params$reduction %||% "UMAP",
                    start_cluster = start_cluster,
                    end_cluster = params$end_cluster,
                    reverse = isTRUE(params$reverse),
                    align_start = isTRUE(params$align_start),
                    seed = params$seed %||% 2025,
                    name = name
                )
                sclet_ai_record_trajectory_evidence(
                    result,
                    name = name,
                    group = params$group,
                    start_cluster = start_cluster,
                    reduction = params$reduction %||% "UMAP"
                )
            },
            input_schema = list(
                group = list(type = "character", required = TRUE),
                start_cluster = list(type = "character", required = TRUE),
                reduction = "character",
                end_cluster = "character",
                reverse = "logical",
                align_start = "logical",
                seed = "integer",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                start <- sclet_ai_resolve_ident_value(object, params$start_cluster)
                if (is.null(start) || !is.character(start) || length(start) != 1L || !nzchar(trimws(start))) {
                    return(paste0(
                        "start_cluster_missing: trajectory root must be an explicit cluster label supplied by the user. ",
                        "Do not rely on slingshot's automatic root selection: inspect the current cluster identities and ",
                        "state which cluster represents the origin."
                    ))
                }
                group <- params$group
                if (is.null(group) || !is.character(group) || length(group) != 1L || !nzchar(trimws(group))) {
                    return("group_missing: run_trajectory requires an explicit 'group' naming the cluster column in colData")
                }
                columns <- colnames(SummarizedExperiment::colData(object))
                if (!group %in% columns) {
                    return(paste0(
                        "group_column_missing: cluster column '", group, "' is not present in colData. ",
                        "Available columns: ", paste(columns, collapse = ", ")
                    ))
                }
                idents <- tryCatch(Idents(object), error = function(e) NULL)
                has_idents <- !is.null(idents) && length(idents) == ncol(object)
                has_cluster_col <- has_idents ||
                    any(c("cluster", "ident", "seurat_clusters") %in% columns)
                if (!has_cluster_col) {
                    return("cluster_identity_missing: run FindClusters() or provide a cluster column before inferring a trajectory")
                }
                labels <- if (group %in% columns) {
                    as.character(SummarizedExperiment::colData(object)[[group]])
                } else {
                    as.character(idents)
                }
                if (!start %in% unique(stats::na.omit(labels))) {
                    return(paste0(
                        "start_cluster_unknown: start_cluster '", start,
                        "' is not one of the current cluster identities in column '", group, "'"
                    ))
                }
                reduction <- params$reduction %||% "UMAP"
                # Same deferred-projection treatment as
                # run_rare_cell_detection's prerequisite (both found the
                # same way): at validation time the reduction exists only
                # in planned$reductions. This trajectory branch is fixed by
                # the same pattern but NOT end-to-end measured in the
                # 2026-10-14 wave (its engine package monocle3 is not
                # installed here); run_doublet_rare_bench.R covers the rare
                # half.
                available <- unique(c(
                    as.character(planned$reductions %||% character()),
                    tryCatch(SingleCellExperiment::reducedDimNames(object),
                        error = function(e) character())))
                if (!reduction %in% available) {
                    return(paste0(
                        "reduction_missing: reduction '", reduction,
                        "' is not available; run the corresponding dimensionality reduction first. ",
                        "Available reductions: ", paste(available, collapse = ", ")
                    ))
                }
                if (!base::requireNamespace("slingshot", quietly = TRUE)) {
                    return("optional_package_missing: slingshot is required for run_trajectory; install via BiocManager::install('slingshot')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(trajectory = "any"),
                note = paste(
                    "Pseudotime is a relative ordering computed from the supplied start_cluster and reduction;",
                    "it is not absolute time. Cluster order is never interpreted as temporal order and no time-point",
                    "labels are generated. No cluster identity is overwritten."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("trajectory", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = TRUE
        )
    ))
    }

    if ("velocity" %in% include) {
    actions <- c(actions, list(
        run_velocity = AIAction(
            name = "run_velocity",
            description = "Estimate RNA velocity with velociraptor::scvelo and record a bounded aggregate evidence node.",
            handler = function(object, params) {
                name <- params$name %||% "velocity"
                result <- RunVelocity(
                    object,
                    mode = params$mode %||% "deterministic",
                    use.dimred = params$use.dimred %||% NULL,
                    subset.row = if (is.null(params$subset.row)) NULL else sclet_ai_vector_param(params$subset.row, "character"),
                    name = name
                )
                sclet_ai_record_velocity_evidence(result, name = name)
            },
            input_schema = list(
                mode = "character",
                use.dimred = "character",
                subset.row = "character_vector",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                if (!all(c("spliced", "unspliced") %in% available)) {
                    return(paste0(
                        "spliced_unspliced_missing: run_velocity requires 'spliced' and 'unspliced' assays; available: ",
                        paste(available, collapse = ", ")
                    ))
                }
                # velociraptor's scvelo method also reads a `counts` assay
                # (for its normalisation) and fails inside the Python call
                # with "'counts' not in names(assays(...))" when it is
                # absent -- a gate that passed validation and died at
                # execution, which is what the prerequisite is for. The
                # spliced/unspliced velocity exports in this repository
                # (hermann sperm, regvelo smoke) carry only logcounts, so
                # this is the ordinary case, not an exotic one (measured
                # 2026-10-20).
                if (!"counts" %in% available) {
                    return(paste0(
                        "counts_assay_missing: run_velocity also needs a 'counts' assay (velociraptor normalises from it); available: ",
                        paste(available, collapse = ", "),
                        ". Supply the raw counts, or use a velocity object that keeps them."
                    ))
                }
                reduction <- params$use.dimred %||% planned$active_reduction %||%
                    tryCatch(DefaultReduction(object), error = function(e) NULL)
                available_reductions <- planned$reductions %||% tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
                if (is.null(reduction) || !nzchar(reduction) || !reduction %in% available_reductions) {
                    return(paste0(
                        "reduction_missing: a dimensionality reduction is required for run_velocity; available: ",
                        paste(available_reductions, collapse = ", ")
                    ))
                }
                if (!base::requireNamespace("velociraptor", quietly = TRUE)) {
                    return("optional_package_missing: velociraptor is required for run_velocity; install via BiocManager::install('velociraptor')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(velocity = "any"),
                note = paste(
                    "Velocity direction and velocity_pseudotime are model-fitted estimates from scVelo,",
                    "not measured biological facts; they depend on the chosen mode and the spliced/unspliced",
                    "ratio and must not be read as ground-truth cell-state transitions. No cell is removed",
                    "or filtered by this action."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("velocity", "ai_evidence"),
            estimated_cost = "high",
            idempotent = FALSE
        )
    ))
    }

    if ("fate" %in% include) {
    actions <- c(actions, list(
        run_fate_analysis = AIAction(
            name = "run_fate_analysis",
            description = "Estimate terminal states and fate probabilities with CellRank from recorded RNA velocity and record a bounded aggregate evidence node.",
            handler = function(object, params) {
                name <- params$name %||% "cellfate"
                cluster_key <- params$cluster_key %||%
                    tryCatch(ActiveIdent(object), error = function(e) NULL)
                result <- RunCellFate(
                    object,
                    method = params$method %||% "CellRank",
                    reduction = params$reduction %||% "PCA",
                    cluster_key = cluster_key,
                    velocity_id = params$velocity_id,
                    backend = params$backend %||% "basilisk",
                    python = params$python,
                    name = name
                )
                # The output schema only requires that *some* trajectory
                # record exists, which an earlier step could have supplied.
                # Verify this action's own result through the same resolver
                # the recorder uses, so a run that produced nothing cannot
                # be reported as completed.
                resolved <- sclet_ai_resolve_cellrank_record(result, cellrank_id = name)
                if (is.null(resolved$id)) {
                    stop(
                        "fate_no_result: RunCellFate returned without a completed CellRank trajectory record.",
                        call. = FALSE
                    )
                }
                result <- sclet_ai_record_fate_evidence(
                    result,
                    name = name,
                    cluster_key = cluster_key
                )
                evidence_id <- sclet_ai_fate_evidence_id(name)
                evidence <- tryCatch(
                    sclet_ai_evidence_get_all(result)[[evidence_id]],
                    error = function(e) NULL
                )
                if (is.null(evidence)) {
                    stop(
                        "fate_evidence_missing: CellRank completed but its bounded evidence node was not recorded.",
                        call. = FALSE
                    )
                }
                result
            },
            input_schema = list(
                method = "character",
                reduction = "character",
                cluster_key = "character",
                velocity_id = "character",
                backend = "character",
                python = "character",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                records <- sclet_get_state_records(object, "velocity")
                completed <- vapply(records, function(record) {
                    identical(as.character(record$status)[[1L]], "completed")
                }, logical(1L))
                completed_ids <- names(records)[completed]
                # A plan chaining run_velocity -> run_fate_analysis has no
                # velocity record at validation time, because nothing has
                # executed yet. The projected id of the earlier step is the
                # promise; RunCellRank() still stops if it turns out broken,
                # and that surfaces as this step failing rather than as a
                # silently empty result.
                planned_velocity <- as.character(
                    planned$state_records$velocity %||% character())
                usable <- unique(c(completed_ids, planned_velocity))
                if (!length(usable)) {
                    return(paste0(
                        "velocity_record_missing: run_fate_analysis requires a completed velocity record; ",
                        "run run_velocity first. Available velocity records: ",
                        paste(names(records), collapse = ", ")
                    ))
                }
                velocity_id <- params$velocity_id
                if (!is.null(velocity_id) && !velocity_id %in% usable) {
                    return(paste0(
                        "velocity_id_unknown: velocity record '", velocity_id,
                        "' is neither a completed velocity record nor produced by an earlier step of ",
                        "this plan. Available completed velocity records: ",
                        paste(completed_ids, collapse = ", "),
                        if (length(planned_velocity)) paste0(
                            "; planned: ", paste(planned_velocity, collapse = ", ")) else ""
                    ))
                }
                cluster_key <- params$cluster_key %||%
                    tryCatch(ActiveIdent(object), error = function(e) NULL)
                columns <- colnames(SummarizedExperiment::colData(object))
                if (is.null(cluster_key) || !is.character(cluster_key) ||
                    length(cluster_key) != 1L || !nzchar(trimws(cluster_key)) ||
                    !cluster_key %in% columns) {
                    return(paste0(
                        "cluster_key_missing: run_fate_analysis requires a cluster column in colData; ",
                        "supply 'cluster_key' or run FindClusters() first. Available columns: ",
                        paste(columns, collapse = ", ")
                    ))
                }
                reduction <- params$reduction %||% "PCA"
                available <- unique(c(
                    as.character(planned$reductions %||% character()),
                    tryCatch(SingleCellExperiment::reducedDimNames(object),
                        error = function(e) character())))
                if (!reduction %in% available) {
                    return(paste0(
                        "reduction_missing: reduction '", reduction,
                        "' is not available; run the corresponding dimensionality reduction first. ",
                        "Available reductions: ", paste(available, collapse = ", ")
                    ))
                }
                if (!base::requireNamespace("reticulate", quietly = TRUE)) {
                    return("fate_backend_unavailable: reticulate is required for run_fate_analysis; install via install.packages('reticulate')")
                }
                backend <- params$backend %||% "basilisk"
                if (identical(backend, "basilisk") &&
                    !base::requireNamespace("basilisk", quietly = TRUE)) {
                    return("optional_package_missing: basilisk is required for backend 'basilisk'; install via BiocManager::install('basilisk')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(trajectory = "any", ai_evidence = "any"),
                note = paste(
                    "Terminal states and absorption probabilities are CellRank model estimates computed from the",
                    "velocity kernel and the supplied cluster labels; they are not measured cell fate and are not",
                    "temporal ground truth. They depend on the chosen method, reduction, cluster annotation and",
                    "velocity estimates. No cell is removed or filtered by this action."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("trajectory", "ai_evidence"),
            estimated_cost = "high",
            idempotent = FALSE
        ),
        confirm_fate_terminal_states = AIAction(
            name = "confirm_fate_terminal_states",
            description = paste(
                "Record the user's own explicit terminal-state interpretation for a CellRank result.",
                "The terminal states named here are a user-chosen interpretation layer over the model output,",
                "supplied by the user and merely checked against the codes actually present in the stored",
                "terminal-state column. Confirming them records a decision; it does not make them measured",
                "biology. This action never picks, suggests, ranks or guesses a terminal state itself."
            ),
            handler = function(object, params) {
                resolved <- sclet_ai_resolve_cellrank_record(object, params$cellrank_id)
                sclet_ai_record_fate_terminal_state_confirmation(
                    object,
                    cellrank_id = resolved$id,
                    terminal_state_codes = sclet_ai_validated_terminal_state_codes(
                        object, resolved$record, params$terminal_state_codes
                    )
                )
            },
            input_schema = list(
                cellrank_id = "character",
                terminal_state_codes = list(type = "character_vector", required = TRUE),
                rationale_code = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                resolved <- sclet_ai_resolve_cellrank_record(object, params$cellrank_id)
                if (!is.null(resolved$message)) return(resolved$message)
                codes <- params$terminal_state_codes
                if (is.null(codes) || !length(codes)) {
                    return(paste0(
                        "terminal_state_codes_missing: confirm_fate_terminal_states requires the user to name the",
                        " terminal states explicitly, e.g. terminal_state_codes = c(\"cluster_1\"). This action never",
                        " picks a terminal state for the user."
                    ))
                }
                codes <- sclet_ai_vector_param(codes, "character")
                if (!length(codes) || any(is.na(codes)) || any(!nzchar(trimws(codes)))) {
                    return(paste0(
                        "terminal_state_codes_missing: terminal_state_codes must be a non-empty character vector",
                        " with no blank entry."
                    ))
                }
                valid <- sclet_ai_terminal_state_codes(object, resolved$record)
                unknown <- setdiff(codes, valid)
                if (length(unknown)) {
                    return(paste0(
                        "terminal_state_unknown: terminal state code(s) ", paste(unknown, collapse = ", "),
                        " are not present in the stored terminal-state column. Valid codes: ",
                        if (length(valid)) paste(valid, collapse = ", ") else "(none)"
                    ))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(trajectory = "any", ai_evidence = "any"),
                note = paste(
                    "The confirmed terminal states are the USER's interpretation layer over CellRank's model",
                    "output, not a measurement. Confirming them records a decision and does not make them",
                    "measured biology or temporal ground truth. They are validated against the codes present in",
                    "the stored terminal-state column, and this action never picks, suggests, ranks or guesses",
                    "which states are terminal. No cell is removed or filtered."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("trajectory", "ai_evidence"),
            estimated_cost = "low",
            idempotent = TRUE
        ),
        summarize_fate_drivers = AIAction(
            name = "summarize_fate_drivers",
            description = paste(
                "Read-only bounded summary of the CellRank lineage drivers for one user-specified lineage.",
                "It reports aggregate statistics only and withholds gene identity; it never ranks lineages,",
                "never nominates a \"best\" lineage, and never writes an evidence node or mutates the object."
            ),
            handler = function(object, params) {
                summarize_fate_lineage_drivers(
                    object,
                    id = params$cellrank_id,
                    lineage = params$lineage,
                    top_n = params$top_n %||% 10L
                )
            },
            input_schema = list(
                cellrank_id = "character",
                lineage = "character",
                top_n = "integer"
            ),
            prerequisites = function(object, params, planned = NULL) {
                resolved <- sclet_ai_resolve_cellrank_record(object, params$cellrank_id)
                if (!is.null(resolved$message)) return(resolved$message)
                drivers <- resolved$record$artifacts$lineage_drivers
                if (is.null(drivers) || !is.data.frame(drivers) || !nrow(drivers)) {
                    return(paste0(
                        "lineage_driver_table_missing: the CellRank record '", resolved$id,
                        "' has no stored lineage driver table; rerun the fate analysis with driver computation",
                        " available before summarizing drivers."
                    ))
                }
                TRUE
            },
            returns = "value",
            requires_confirmation = FALSE,
            output_schema = list(
                required_states = list(trajectory = "any", ai_evidence = "any"),
                note = paste(
                    "Read-only: reports bounded aggregate statistics for one lineage the user named, and",
                    "deliberately withholds gene identity. It never ranks lineages, never nominates a",
                    "\"best\" lineage, and writes no evidence node; the object is not mutated."
                )
            ),
            mutates_object = FALSE,
            allowed_state_types = c("trajectory", "ai_evidence"),
            estimated_cost = "low",
            idempotent = TRUE
        )
    ))
    }

    actions
}
