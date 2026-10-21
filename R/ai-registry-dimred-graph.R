# sclet_ai_registry_dimred_graph(): built-in AI actions for dimensional reduction and the neighbour graph.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# Two groups, one file: they are the PCA/UMAP-then-neighbours pair, and the
#   registry emits dimred before graph in this order.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_dimred_graph <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("dimred" %in% include) {
    actions <- c(actions, list(
        run_pca = AIAction(
            name = "run_pca",
            description = "Compute and register a PCA reduction.",
            handler = function(object, params) RunPCA(
                object,
                subset_row = sclet_ai_gene_selector(
                    object,
                    sclet_ai_resolve_deferred_selector(object, params$subset_row, "subset_row"),
                    arg = "subset_row"
                ),
                exprs_values = params$exprs_values %||% NULL,
                layer = params$layer %||% NULL,
                ncomponents = params$ncomponents %||% 50
            ),
            input_schema = list(
                subset_row = "character_vector",
                exprs_values = "character",
                layer = "character",
                ncomponents = "integer"
            ),
            prerequisites = function(object, params, planned = NULL) {
                available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                source <- params$exprs_values %||% params$layer %||%
                    planned$active_assay %||% tryCatch(DefaultLayer(object), error = function(e) NULL)
                if (!is.null(source) && !source %in% available) {
                    return(paste0("PCA source assay/layer is not available: ", source))
                }
                if (!length(available)) return("at least one assay is required")
                if (sclet_ai_is_deferred_gene_selector(params$subset_row)) {
                    return(TRUE)
                }
                sclet_ai_selector_problem(object, params$subset_row, "subset_row")
            },
            returns = "sce",
            output_schema = list(
                reduction = "PCA",
                required_reductions = "PCA",
                required_states = list(reduction = "pca"),
                active_reduction = "PCA"
            ),
            allowed_state_types = "reduction",
            estimated_cost = "medium",
            idempotent = TRUE
        ),
        run_umap = AIAction(
            name = "run_umap",
            description = "Compute and register a UMAP reduction from an existing reduction.",
            handler = function(object, params) RunUMAP(
                object,
                dims = if (is.null(params$dims)) NULL else sclet_ai_vector_param(params$dims, "integer"),
                reduction = params$reduction %||% NULL,
                layer = params$layer %||% NULL
            ),
            input_schema = list(
                dims = "integer_vector",
                reduction = "character",
                layer = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                reduction <- params$reduction %||% planned$active_reduction %||% DefaultReduction(object)
                # Scalar guard first: `${step_N.output.<field>}` naming a
                # multi-valued schema key resolves (after projection) to a
                # vector, and `!reduction %in% available` on length > 1
                # raises a bare condition error instead of returning a
                # validation message. Found 2026-10-14 via the integration
                # benchmark (a plan referencing run_integration's
                # three-value allowed_reductions).
                if (length(reduction) > 1L || (!is.null(reduction) && !is.character(reduction))) {
                    return("reduction must be a single reduction name")
                }
                available <- planned$reductions %||% SingleCellExperiment::reducedDimNames(object)
                if (is.null(reduction) || !reduction %in% available) {
                    return("a valid source reduction is required; run run_pca first")
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(
                reduction = "UMAP",
                required_reductions = "UMAP",
                required_states = list(reduction = "umap"),
                active_reduction = "UMAP"
            ),
            allowed_state_types = "reduction",
            estimated_cost = "medium",
            idempotent = TRUE
        )
    ))
    }

    if ("graph" %in% include) {
    actions <- c(actions, list(
        find_neighbors = AIAction(
            name = "find_neighbors",
            description = "Build and register a KNN/SNN graph from a reduction.",
            handler = function(object, params) FindNeighbors(
                object,
                dims = sclet_ai_vector_param(params$dims, "integer"),
                reduction = params$reduction %||% NULL,
                k = params$k %||% 10
            ),
            input_schema = list(
                dims = list(type = "integer_vector", required = TRUE),
                reduction = "character",
                k = "integer"
            ),
            prerequisites = function(object, params, planned = NULL) {
                reduction <- params$reduction %||% planned$active_reduction %||% DefaultReduction(object) %||% "PCA"
                # Same scalar guard as run_umap: a projected multi-valued
                # reference arrives as a vector and made the %in% test
                # below raise "condition has length > 1" (measured
                # 2026-10-14 via the integration benchmark).
                if (length(reduction) > 1L || !is.character(reduction)) {
                    return("reduction must be a single reduction name")
                }
                available <- planned$reductions %||% SingleCellExperiment::reducedDimNames(object)
                if (!reduction %in% available) {
                    return(paste0("reduction is not available: ", reduction))
                }
                dims <- if (is.null(params$dims)) NULL else sclet_ai_vector_param(params$dims, "integer")
                # When the reduction is only PLANNED (registered by an
                # upstream step, e.g. HARMONY from run_integration), the
                # object does not carry it yet and neither does
                # `planned$reduction_dims` -- and `[[` on a missing name
                # THROWS rather than returning NULL, so even the `%||%`
                # fallback could not save it (the `invalid subscript 'type'`
                # error was from inside reducedDim after the first throw was
                # bypassed; both layers must be guarded). The third layer of
                # the same deferred-projection gap found 2026-10-14. Skip
                # the upper bound when the dimension count is not knowable
                # at planning time; execution stays the authority.
                rd <- planned$reduction_dims
                max_dim <- if (!is.null(rd) && reduction %in% names(rd)) {
                    rd[[reduction]]
                } else {
                    tryCatch(ncol(SingleCellExperiment::reducedDim(object, reduction)),
                        error = function(e) NULL)
                }
                if (is.null(dims) || any(dims < 1L) ||
                    (!is.null(max_dim) && any(dims > max_dim))) {
                    return("dims must be within the selected reduction")
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(
                graph = "knn_graph",
                active_graph = "knn_graph",
                required_graphs = "knn_graph",
                required_states = list(graph = "knn_graph")
            ),
            allowed_state_types = "graph",
            estimated_cost = "medium",
            idempotent = TRUE
        )
    ))
    }

    actions
}
