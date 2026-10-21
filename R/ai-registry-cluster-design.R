# sclet_ai_registry_cluster_design(): built-in AI actions for clustering and design confirmation.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# Cluster actions come first, then confirm_design, which sits behind a
#   composite predicate: it appears when any of integration, annotation or
#   malignant is requested, because every design-aware group needs a confirmed
#   batch/condition column. Keeping both here preserves confirm_design's
#   position between find_clusters and run_integration.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_cluster_design <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("cluster" %in% include) {
    actions <- c(actions, list(
        find_clusters = AIAction(
            name = "find_clusters",
            description = paste(
                "Find Louvain clusters from the registered KNN graph.",
                "Louvain consumes random numbers, so `seed` (default 1) decides",
                "the labels; it is recorded with the clustering state so the",
                "same resolution reproduces the same clusters."
            ),
            handler = function(object, params) FindClusters(
                object,
                resolution = params$resolution %||% 0.5,
                seed = params$seed %||% 1L
            ),
            input_schema = list(resolution = "number", seed = "integer"),
            prerequisites = function(object, params, planned = NULL) {
                graph <- planned$active_graph %||% DefaultGraph(object) %||% "knn_graph"
                available <- planned$graphs %||% names(sclet_get_state(object)$graphs %||% list())
                if (!graph %in% available) {
                    return("knn_graph is required; run find_neighbors first")
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(
                state = "clustering:louvain_clusters",
                required_states = list(clustering = "louvain_clusters"),
                active_ident = "colLabels"
            ),
            allowed_state_types = "clustering",
            estimated_cost = "medium",
            idempotent = TRUE
        )
    ))
    }

    # Design semantics, as an action.
    #
    # `run_integration` refuses to run until the batch column has been confirmed
    # with ConfirmAIDesignSemantics(), and no registered action could perform
    # that confirmation. On full Kang (2026-10-18) the goal asked to "integrate
    # the two sequencing batches" and every plan -- including the goal-coverage
    # repair round -- was rejected with `design_semantics_not_confirmed`, because
    # the only way to clear it was a function call the planner cannot make.
    #
    # This action says "this colData column plays this design role", which is a
    # statement about the data, so it mutates the object and is recorded like
    # any other analysis step. It deliberately does not guess: every role must be
    # named in the plan, and an unknown column is refused.
    # Design confirmation mutates the object, so it must not live in the
    # read-only default group (that group is asserted to be non-mutating).
    # Integration needs it; annotation/malignant claims use design columns too.
    if (any(c("integration", "annotation", "malignant") %in% include)) {
    actions <- c(actions, list(            confirm_design = AIAction(
            name = "confirm_design",
            description = paste0(
                "Confirm which colData columns play the design roles before ",
                "integration or design-aware analysis. Exposed roles: batch, ",
                "condition, individual (maps to the subject role) and ",
                "specimen (maps to the sample role), e.g. ",
                "batch = \"batch\", condition = \"condition\", individual = \"donor\"."
            ),
            handler = function(object, params) {
                roles <- list()
                exposed <- c(batch = "batch", condition = "condition",
                    individual = "subject", specimen = "sample")
                for (key in names(exposed)) {
                    if (!is.null(params[[key]])) {
                        roles[[exposed[[key]]]] <-
                            sclet_ai_vector_param(params[[key]], "character")[[1L]]
                    }
                }
                if (!length(roles)) {
                    stop("confirm_design needs at least one role column")
                }
                ConfirmAIDesignSemantics(object, design = roles)
            },
            input_schema = list(
                # Names avoid "subject"/"sample": the privacy sanitizer
                # redacts fields matching identifier patterns, so a schema
                # key called `subject` was stripped out of available_actions
                # (with a warning) before the planner ever saw it. These map
                # onto the design roles internally.
                batch = "character",
                condition = "character",
                individual = "character",
                specimen = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                cd_names <- planned$colData_columns %||%
                    colnames(SummarizedExperiment::colData(object))
                roles <- c("batch", "condition", "individual", "specimen")
                given <- Filter(function(role) !is.null(params[[role]]), roles)
                if (!length(given)) {
                    return(paste0("design_role_missing: name at least one of ",
                        paste(roles, collapse = ", "), " as a colData column"))
                }
                for (role in given) {
                    value <- params[[role]]
                    if (!is.character(value) || length(value) != 1L || !nzchar(value)) {
                        return(paste0("design_role_invalid: '", role,
                            "' must be one colData column name"))
                    }
                    if (!value %in% cd_names) {
                        return(paste0("design_role_invalid: '", role,
                            "' points to '", value,
                            "', which is not in colData; available: ",
                            paste(cd_names, collapse = ", ")))
                    }
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(required_states = list(ai_design_confirmation = "any")),
            allowed_state_types = c("ai_design_confirmation"),
            mutates_object = TRUE,
            estimated_cost = "low",
            idempotent = TRUE
        )
    ))
    }

    actions
}
