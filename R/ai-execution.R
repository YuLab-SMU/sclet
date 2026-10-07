#' Define one explicitly registered AI execution action
#'
#' @param name Stable action name used in plans.
#' @param handler Function receiving `(object, params)` and returning an updated
#'   `SingleCellExperiment` by default.
#' @param description Human-readable action description.
#' @param prerequisites Function receiving `(object, params)` and returning
#'   `TRUE`, `FALSE`, or a character explanation.
#' @param returns Either `"sce"` or `"value"`.
#' @param requires_confirmation Logical. Defaults to `TRUE`.
#' @param input_schema Optional named parameter schema. Character vectors of
#'   length greater than one are treated as enumerations; scalar values name
#'   primitive types.
#' @param output_schema Optional bounded output schema metadata.
#' @param mutates_object Logical. Whether the handler may modify the SCE.
#' @param allowed_state_types Optional state types the handler may write.
#' @param estimated_cost Optional cost label such as `"low"`, `"medium"`, or `"high"`.
#' @param idempotent Logical. Whether repeating the action is expected to be safe.
#' @return An action descriptor.
#' @export
AIAction <- function(
    name,
    handler,
    description = NULL,
    prerequisites = function(object, params) TRUE,
    returns = c("sce", "value"),
    requires_confirmation = TRUE,
    input_schema = list(),
    output_schema = list(),
    mutates_object = NULL,
    allowed_state_types = character(),
    estimated_cost = "low",
    idempotent = FALSE
) {
    returns <- match.arg(returns)
    estimated_cost <- match.arg(estimated_cost, c("low", "medium", "high"))
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
        stop("`name` must be a single non-empty character string.")
    }
    if (!is.function(handler)) {
        stop("`handler` must be a function.")
    }
    if (!is.function(prerequisites)) {
        stop("`prerequisites` must be a function.")
    }
    if (!is.list(input_schema) || is.null(names(input_schema)) && length(input_schema)) {
        stop("`input_schema` must be a named list.")
    }
    if (!is.list(output_schema)) {
        stop("`output_schema` must be a list.")
    }
    if (is.null(mutates_object)) {
        mutates_object <- identical(returns, "sce")
    }
    if (!is.logical(mutates_object) || length(mutates_object) != 1L || is.na(mutates_object)) {
        stop("`mutates_object` must be a single non-missing logical value.")
    }
    if (!is.logical(requires_confirmation) || length(requires_confirmation) != 1L || is.na(requires_confirmation)) {
        stop("`requires_confirmation` must be a single non-missing logical value.")
    }
    if (!is.logical(idempotent) || length(idempotent) != 1L || is.na(idempotent)) {
        stop("`idempotent` must be a single non-missing logical value.")
    }
    allowed_state_types <- as.character(allowed_state_types)
    descriptor <- list(
        name = name,
        description = description %||% paste("Execute registered action", name),
        handler = handler,
        prerequisites = prerequisites,
        returns = returns,
        requires_confirmation = isTRUE(requires_confirmation),
        input_schema = input_schema,
        output_schema = output_schema,
        mutates_object = isTRUE(mutates_object),
        allowed_state_types = allowed_state_types,
        estimated_cost = estimated_cost,
        idempotent = isTRUE(idempotent)
    )
    class(descriptor) <- c("sclet_ai_action", "list")
    descriptor
}

sclet_ai_vector_param <- function(value, type = c("character", "integer")) {
    type <- match.arg(type)
    if (is.list(value)) value <- unlist(value, use.names = FALSE)
    if (type == "character") return(as.character(value))
    as.integer(value)
}
sclet_ai_action_param_problems <- function(params, schema) {
    if (is.null(params)) {
        params <- list()
    }
    if (!is.list(params)) {
        return("action params must be a list")
    }
    if (!length(schema)) {
        return(character())
    }
    schema_names <- names(schema)
    if (is.null(schema_names) || any(!nzchar(schema_names))) {
        return("action input_schema must be named")
    }
    problems <- character()
    required_names <- names(schema)[vapply(schema, function(spec) {
        is.list(spec) && isTRUE(spec$required)
    }, logical(1))]
    missing_required <- setdiff(required_names, names(params))
    if (length(missing_required)) {
        problems <- c(problems, paste0("missing required parameter(s): ", paste(missing_required, collapse = ", ")))
    }
    unknown <- setdiff(names(params), schema_names)
    if (length(unknown)) {
        problems <- c(problems, paste0("unknown parameter(s): ", paste(unknown, collapse = ", ")))
    }
    for (name in intersect(names(params), schema_names)) {
        value <- params[[name]]
        spec <- schema[[name]]
        required <- FALSE
        if (is.list(spec) && !is.null(spec$type)) {
            required <- isTRUE(spec$required)
            spec <- spec$type
        }
        if (is.null(value)) {
            if (required) {
                problems <- c(problems, paste0("missing required parameter: ", name))
            }
            next
        }
        if (is.character(spec) && length(spec) > 1L) {
            if (!as.character(value)[[1L]] %in% spec) {
                problems <- c(problems, paste0("parameter ", name, " must be one of: ", paste(spec, collapse = ", ")))
            }
            next
        }
        type <- as.character(spec)[[1L]]
        valid <- switch(
            type,
            character_vector = (is.character(value) && length(value) >= 1L) ||
                (is.list(value) && length(value) >= 1L && all(vapply(value, function(x) is.character(x) && length(x) == 1L, logical(1)))),
            integer_vector = ((is.numeric(value) && length(value) >= 1L && all(value == as.integer(value))) ||
                (is.list(value) && length(value) >= 1L && all(vapply(value, function(x) is.numeric(x) && length(x) == 1L && x == as.integer(x), logical(1))))),
            numeric_vector = is.numeric(value) && length(value) >= 1L,
            character = is.character(value) && length(value) == 1L,
            string = is.character(value) && length(value) == 1L,
            logical = is.logical(value) && length(value) == 1L,
            boolean = is.logical(value) && length(value) == 1L,
            numeric = is.numeric(value) && length(value) == 1L,
            number = is.numeric(value) && length(value) == 1L,
            integer = is.numeric(value) && length(value) == 1L && value == as.integer(value),
            list = is.list(value),
            object = (is.atomic(value) && length(value) > 0L) || is.list(value) || inherits(value, "SummarizedExperiment"),
            TRUE
        )
        if (!isTRUE(valid)) {
            problems <- c(problems, paste0("parameter ", name, " has invalid type"))
        }
    }
    unique(problems)
}

#' Build the safe built-in action catalog
#'
#' The `read` group exposes bounded deterministic ledger views. The optional
#' analysis groups wrap existing sclet functions and therefore require the
#' normal plan validation and confirmation safeguards.
#'
#' @param object A `SingleCellExperiment` object.
#' @param include Character vector of action groups. Supported groups are
#'   `"read"`, `"preprocess"`, `"dimred"`, `"graph"`, `"cluster"`,
#'   `"integration"`, `"annotation"`, `"rare_cell"`, `"trajectory"`, and
#'   `"all"`. Defaults to `"read"`.
#' @return A `sclet_ai_execution_registry` containing the requested actions.
#' @export
AIDefaultExecutionRegistry <- function(object, include = "read") {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    include <- unique(as.character(include))
    allowed_groups <- c("read", "preprocess", "dimred", "graph", "cluster", "integration", "annotation", "rare_cell", "trajectory", "all")
    if (!length(include) || any(!include %in% allowed_groups)) {
        stop("`include` must contain only: ", paste(allowed_groups, collapse = ", "))
    }
    if ("all" %in% include) {
        include <- setdiff(allowed_groups, "all")
    }
    actions <- list()
    if ("read" %in% include) {
        actions <- c(actions, list(
            inspect_status = AIAction(
                name = "inspect_status",
                description = "Inspect deterministic sclet status without modifying the object.",
                handler = function(object, params) Status(object),
                returns = "value",
                requires_confirmation = FALSE,
                mutates_object = FALSE,
                output_schema = list(type = "object"),
                estimated_cost = "low",
                idempotent = TRUE
            ),
            inspect_ledger = AIAction(
                name = "inspect_ledger",
                description = "Inspect the bounded AI-facing analysis ledger.",
                handler = function(object, params) {
                    detail <- params$detail %||% "summary"
                    target <- params$target %||% NULL
                    GetAnalysisLedger(object, detail = detail, target = target)
                },
                input_schema = list(
                    detail = c("summary", "full"),
                    target = "character"
                ),
                returns = "value",
                requires_confirmation = FALSE,
                mutates_object = FALSE,
                output_schema = list(type = "object"),
                estimated_cost = "low",
                idempotent = TRUE
            ),
            check_qc = AIAction(
                name = "check_qc",
                description = "Inspect deterministic quality checks without modifying the object.",
                handler = function(object, params) GetAnalysisLedger(object)$quality_checks,
                returns = "value",
                requires_confirmation = FALSE,
                mutates_object = FALSE,
                output_schema = list(type = "object"),
                estimated_cost = "low",
                idempotent = TRUE
            )
        ))
    }
    if ("preprocess" %in% include) {
        actions <- c(actions, list(
            normalize_data = AIAction(
                name = "normalize_data",
                description = "Normalize an assay into the sclet logcounts layer.",
                handler = function(object, params) NormalizeData(
                    object,
                    scale.factor = params$scale.factor %||% 10000,
                    assay = params$assay %||% "counts"
                ),
                input_schema = list(
                    scale.factor = "number",
                    assay = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    assay <- params$assay %||% "counts"
                    available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                    if (!assay %in% available) {
                        return(paste0("assay is not available: ", assay))
                    }
                    TRUE
                },
                returns = "sce",
                output_schema = list(
                    assay = "logcounts",
                    required_assays = "logcounts",
                    required_states = list(preprocess = "normalize_logcounts"),
                    active_assay = "logcounts"
                ),
                allowed_state_types = "preprocess",
                estimated_cost = "medium",
                idempotent = TRUE
            ),
            find_variable_features = AIAction(
                name = "find_variable_features",
                description = "Identify highly variable genes using sclet's variance model.",
                handler = function(object, params) FindVariableFeatures(
                    object,
                    nfeatures = params$nfeatures %||% 2000,
                    method = params$method %||% "scran"
                ),
                input_schema = list(
                    nfeatures = "integer",
                    method = c("scran", "scrapper", "seurat")
                ),
                prerequisites = function(object, params, planned = NULL) {
                    available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                    if (!"logcounts" %in% available) {
                        return("logcounts assay is required; run normalize_data first")
                    }
                    TRUE
                },
                returns = "sce",
                output_schema = list(
                    required_commands = "FindVariableFeatures",
                    required_hvg = TRUE
                ),
                allowed_state_types = character(),
                estimated_cost = "medium",
                idempotent = TRUE
            ),
            scale_data = AIAction(
                name = "scale_data",
                description = "Scale an assay and register the scaled layer.",
                handler = function(object, params) ScaleData(
                    object,
                    features = if (is.null(params$features)) NULL else sclet_ai_vector_param(params$features, "character"),
                    assay = params$assay %||% "logcounts"
                ),
                input_schema = list(
                    features = "character_vector",
                    assay = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    assay <- params$assay %||% "logcounts"
                    available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                    if (!assay %in% available) {
                        return(paste0("assay is not available: ", assay))
                    }
                    TRUE
                },
                returns = "sce",
                output_schema = list(
                    required_assays = "scaled",
                    active_assay = "scaled"
                ),
                allowed_state_types = character(),
                estimated_cost = "medium",
                idempotent = TRUE
            )
        ))
    }
    if ("dimred" %in% include) {
        actions <- c(actions, list(
            run_pca = AIAction(
                name = "run_pca",
                description = "Compute and register a PCA reduction.",
                handler = function(object, params) RunPCA(
                    object,
                    subset_row = if (is.null(params$subset_row)) NULL else sclet_ai_vector_param(params$subset_row, "character"),
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
                    TRUE
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
                    available <- planned$reductions %||% SingleCellExperiment::reducedDimNames(object)
                    if (!reduction %in% available) {
                        return(paste0("reduction is not available: ", reduction))
                    }
                    dims <- if (is.null(params$dims)) NULL else sclet_ai_vector_param(params$dims, "integer")
                    max_dim <- planned$reduction_dims[[reduction]] %||%
                        ncol(SingleCellExperiment::reducedDim(object, reduction))
                    if (is.null(dims) || any(dims < 1L) || any(dims > max_dim)) {
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
    if ("cluster" %in% include) {
        actions <- c(actions, list(
            find_clusters = AIAction(
                name = "find_clusters",
                description = "Find Louvain clusters from the registered KNN graph.",
                handler = function(object, params) FindClusters(
                    object,
                    resolution = params$resolution %||% 0.5
                ),
                input_schema = list(resolution = "number"),
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
    if ("integration" %in% include) {
        actions <- c(actions, list(
            run_integration = AIAction(
                name = "run_integration",
                description = "Run batch-correction integration via fastMNN, Harmony, or scVI and register the result.",
                handler = function(object, params) {
                    method <- match.arg(params$method %||% "fastMNN", c("fastMNN", "Harmony", "scVI"))
                    batch <- params$batch
                    features <- if (!is.null(params$features)) sclet_ai_vector_param(params$features, "character") else NULL
                    dims <- if (!is.null(params$dims)) sclet_ai_vector_param(params$dims, "integer") else NULL
                    RunIntegration(
                        object,
                        method = method,
                        batch = batch,
                        features = features,
                        layer = params$layer %||% NULL,
                        reduction = params$reduction %||% NULL,
                        name = params$name %||% tolower(method),
                        dims = dims
                    )
                },
                input_schema = list(
                    batch = list(type = "character", required = TRUE),
                    method = c("fastMNN", "Harmony", "scVI"),
                    name = "character",
                    dims = "integer_vector",
                    features = "character_vector",
                    layer = "character",
                    reduction = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    method <- match.arg(params$method %||% "fastMNN", c("fastMNN", "Harmony", "scVI"))
                    batch <- params$batch
                    if (is.null(batch) || !is.character(batch) || length(batch) != 1L || !nzchar(batch)) {
                        return("missing required parameter: batch (a single colData column name)")
                    }
                    cd_names <- planned$colData_columns %||% colnames(SummarizedExperiment::colData(object))
                    if (!batch %in% cd_names) {
                        return(paste0("batch column is not available in colData: ", batch))
                    }
                    confirmation <- sclet_ai_find_design_confirmation(object, batch = batch)
                    if (is.null(confirmation)) {
                        return(paste0(
                            "design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, design = list(batch = '",
                            batch, "', ...)) before running integration"
                        ))
                    }
                    if (identical(method, "Harmony")) {
                        if (!requireNamespace("harmony", quietly = TRUE)) {
                            return("optional_package_missing: harmony is required for method Harmony; install.packages('harmony')")
                        }
                    }
                    if (identical(method, "scVI")) {
                        if (!requireNamespace("basilisk", quietly = TRUE)) {
                            return("optional_package_missing: basilisk is required for method scVI")
                        }
                        if (!requireNamespace("zellkonverter", quietly = TRUE)) {
                            return("optional_package_missing: zellkonverter is required for method scVI")
                        }
                    }
                    TRUE
                },
                returns = "sce",
                requires_confirmation = TRUE,
                output_schema = list(
                    required_states = list(integration = "any_qualified"),
                    allowed_reductions = c("HARMONY", "scVI")
                ),
                mutates_object = TRUE,
                allowed_state_types = c("integration", "reduction"),
                estimated_cost = "medium",
                idempotent = FALSE
            )
        ))
    }
    if ("annotation" %in% include) {
        actions <- c(actions, list(
            run_de_test = AIAction(
                name = "run_de_test",
                description = "Run a differential-expression or all-cluster marker analysis and register evidence.",
                handler = function(object, params) {
                    ident.1 <- if (!is.null(params$ident.1)) sclet_ai_vector_param(params$ident.1, "character") else NULL
                    ident.2 <- if (!is.null(params$ident.2)) sclet_ai_vector_param(params$ident.2, "character") else NULL
                    result <- RunDEtest(
                        object,
                        ident.1 = ident.1,
                        ident.2 = ident.2,
                        name = params$name %||% if (is.null(ident.1)) "findallmarkers" else "detest",
                        min.pct = params$min.pct %||% 0.01,
                        logfc.threshold = params$logfc.threshold %||% 0.1
                    )
                    rec_name <- params$name %||% if (is.null(ident.1)) "findallmarkers" else "detest"
                    de_records <- sclet_get_state_records(result, "detest")
                    rec <- de_records[[rec_name]]
                    artifact <- if (!is.null(rec)) {
                        tryCatch(rec$artifacts$result, error = function(e) NULL)
                    } else NULL
                    if (!is.null(artifact) && is.data.frame(artifact)) {
                        logfc_col <- grep("logFC|avg_log2FC|avg_logFC|log2FC", colnames(artifact), value = TRUE, ignore.case = TRUE)[[1L]]
                        if (!is.na(logfc_col)) {
                            p_col <- grep("^p_val_adj$|^adj.P.Val$|padj|FDR|p_adj", colnames(artifact), value = TRUE, ignore.case = TRUE)[[1L]]
                            p_col <- if (!is.na(p_col)) p_col else grep("^p_val$|P.Value|^pvalue$", colnames(artifact), value = TRUE, ignore.case = TRUE)[[1L]]
                            p_thresh <- 0.05
                            fc_thresh <- 0.1
                            sig <- if (is.na(p_col)) {
                                rep(TRUE, nrow(artifact))
                            } else {
                                artifact[[p_col]] < p_thresh
                            }
                            sig[is.na(sig)] <- FALSE
                            up <- sig & artifact[[logfc_col]] >= fc_thresh
                            top_n <- 10L
                            ord <- order(ifelse(up, -artifact[[logfc_col]], NA_real_), na.last = NA)
                            top_genes <- if (length(ord)) as.character(rownames(artifact)[utils::head(ord[seq_len(min(top_n, length(ord)))], top_n)]) else character()
                            idents <- Idents(result)
                            groups <- unique(as.character(idents))
                            aggregate_values <- list(
                                n_rows = as.integer(nrow(artifact)),
                                n_significant = as.integer(sum(sig)),
                                n_upregulated = as.integer(sum(up)),
                                n_groups = as.integer(length(groups)),
                                smallest_group_n = as.integer(min(as.integer(table(as.character(idents), useNA = "no")))),
                                pvalue_threshold = p_thresh,
                                logfc_threshold = fc_thresh,
                                raw_values_included = FALSE
                            )
                            top <- rep(NA_integer_, 10L)
                            if (length(top_genes) > 0L) top[seq_len(min(10L, length(top_genes)))] <- seq_len(min(10L, length(top_genes)))
                            aggregate_values$n_top_up <- as.integer(sum(!is.na(top)))
                            aggregate_values$top_up_ranks <- as.integer(top)
                            aggregate_values$top_up_fractions <- as.numeric(rep_len(NA_real_, 10L))
                            aggregate_values$annotation_method_code <- 1L
                            aggregate_values$annotation_scope_code <- 2L
                            evidence <- list(
                                id = if (is.null(ident.1)) "ev:findallmarkers" else paste0("ev:detest_", paste(ident.1, collapse = "_"), if (length(ident.2)) paste0("_vs_", paste(ident.2, collapse = "_"))),
                                kind = "deterministic_summary",
                                values = aggregate_values,
                                claim_level = "associated"
                            )
                            tryCatch(RecordAIEvidence(result, evidence,
                                source = rec_name,
                                parents = if (!is.null(rec$id)) as.character(rec$id) else character(),
                                scope = NULL),
                                error = function(e) result)
                        } else result
                    } else result
                },
                input_schema = list(
                    ident.1 = "character_vector",
                    ident.2 = "character_vector",
                    name = "character",
                    min.pct = "number",
                    logfc.threshold = "number"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    ident.1 <- if (is.null(params$ident.1)) NULL else sclet_ai_vector_param(params$ident.1, "character")
                    ident.2 <- if (is.null(params$ident.2)) NULL else sclet_ai_vector_param(params$ident.2, "character")
                    active_ident <- ActiveIdent(object)
                    if (is.null(active_ident)) {
                        return("cluster identity is not set: run FindClusters() or set ActiveIdent(object) <- 'colname' before running marker/DE tests")
                    }
                    idents <- Idents(object)
                    if (is.null(idents)) {
                        return("cluster identity is available at slot but has no values")
                    }
                    min_cells_per_group <- 2L
                    ident_vec <- as.character(idents)
                    grp_n <- table(ident_vec, useNA = "no")
                    if (length(grp_n) < 1L || min(as.integer(grp_n)) < min_cells_per_group) {
                        return(paste0("at least one group/identity is too small; minimum required per group: ", min_cells_per_group))
                    }
                    if (is.null(ident.1) && length(grp_n) < 2L) {
                        return("FindAllMarkers requires at least two distinct identity groups; object has 1")
                    }
                    if (!is.null(ident.1) && !all(ident.1 %in% names(grp_n))) {
                        return(paste0("ident.1 contains unknown group labels: ", paste(setdiff(ident.1, names(grp_n)), collapse = ", ")))
                    }
                    if (!is.null(ident.2) && !all(ident.2 %in% names(grp_n))) {
                        return(paste0("ident.2 contains unknown group labels: ", paste(setdiff(ident.2, names(grp_n)), collapse = ", ")))
                    }
                    if (!is.null(ident.1) && !is.null(ident.2) && length(intersect(ident.1, ident.2))) {
                        return("ident.1 and ident.2 overlap; DE comparison groups must be disjoint")
                    }
                    TRUE
                },
                returns = "sce",
                requires_confirmation = TRUE,
                output_schema = list(
                    required_states = list(detest = "any")
                ),
                mutates_object = TRUE,
                allowed_state_types = c("detest", "ai_evidence"),
                estimated_cost = "medium",
                idempotent = FALSE
            ),
            run_annotation = AIAction(
                name = "run_annotation",
                description = "Run SingleR reference-based annotation with an explicitly supplied reference and record aggregate evidence.",
                handler = function(object, params) {
                    ref_arg <- params$ref
                    if (is.character(ref_arg) && length(ref_arg) == 1L && nzchar(ref_arg)) {
                        if (grepl("Data$", ref_arg) && requireNamespace("celldex", quietly = TRUE)) {
                            fun <- tryCatch(getExportedValue("celldex", ref_arg), error = function(e) NULL)
                            if (!is.function(fun)) stop(paste0("reference '", ref_arg, "' is not a function in celldex"), call. = FALSE)
                            ref <- do.call(fun, list())
                        } else {
                            stop(paste0("reference '", ref_arg, "' is not a supported celldex dataset name (must end with 'Data'). Provide a SummarizedExperiment ref directly or a celldex function name."), call. = FALSE)
                        }
                    } else {
                        ref <- ref_arg
                    }
                    labels_arg <- params$labels
                    result <- RunSingleR(
                        object,
                        ref = ref,
                        labels = labels_arg,
                        layer = params$layer %||% NULL,
                        name = params$name %||% "singler"
                    )
                    name <- params$name %||% "singler"
                    ann_records <- sclet_get_state_records(result, "annotation")
                    ann_rec <- ann_records[[name]]
                    if (!is.null(ann_rec)) {
                        artifacts <- if (!is.null(ann_rec)) tryCatch(ann_rec$artifacts, error = function(e) NULL) else NULL
                        labels_col <- if (!is.null(artifacts)) as.character(artifacts$labels_col) else NULL
                        score_col <- if (!is.null(artifacts)) as.character(artifacts$score_col) else NULL
                        cd <- SummarizedExperiment::colData(result)
                        lab_vec <- if (!is.null(labels_col) && labels_col %in% colnames(cd)) cd[[labels_col]] else NULL
                        score_vec <- if (!is.null(score_col) && score_col %in% colnames(cd)) cd[[score_col]] else NULL
                        idents <- if (!is.null(ActiveIdent(result))) as.character(Idents(result)) else NULL
                        n_cells <- if (!is.null(idents)) length(idents) else ncol(result)
                        idents_vec <- if (is.null(idents)) rep("group_0", n_cells) else {
                            u <- unique(as.character(idents))
                            map <- paste0("group_", seq_along(u))
                            names(map) <- u
                            as.character(map[as.character(idents)])
                        }
                        aggregate_values <- list(n_cells = as.integer(n_cells))
                        if (!is.null(lab_vec)) {
                            lab_vec_ch <- as.character(lab_vec)
                            lab_vec_ch[!nzchar(lab_vec_ch) | is.na(lab_vec_ch)] <- "NA"
                            tb <- table(lab_vec_ch, useNA = "no")
                            top_n <- 15L
                            ord <- order(as.integer(tb), decreasing = TRUE)
                            top_labels <- names(tb)[utils::head(ord, top_n)]
                            top_fractions <- round(as.integer(tb[top_labels]) / n_cells, 4L)
                            aggregate_values$n_labels <- as.integer(length(tb))
                            tl <- rep(NA_integer_, 15L)
                            if (length(top_labels) > 0L) tl[seq_len(min(15L, length(top_labels)))] <- seq_len(min(15L, length(top_labels)))
                            tf <- rep_len(NA_real_, 15L)
                            if (length(top_fractions) > 0L) tf[seq_len(min(15L, length(top_fractions)))] <- top_fractions
                            aggregate_values$n_top_labels <- as.integer(sum(!is.na(tl)))
                            aggregate_values$top_label_ranks <- as.integer(tl)
                            aggregate_values$top_label_fractions <- as.numeric(tf)
                            if (!is.null(idents) && length(idents) == n_cells) {
                                tb2 <- table(idents_vec, lab_vec_ch, useNA = "no")
                                cluster_totals <- rowSums(tb2); cluster_totals[cluster_totals == 0L] <- NA_integer_
                                prop <- as.data.frame(tb2 / cluster_totals[row(tb2)])
                                colnames(prop) <- c("cluster", "label", "fraction")
                                prop <- prop[prop$fraction > 0.05 & !is.na(prop$fraction), , drop = FALSE]
                                prop <- prop[order(-prop$fraction), , drop = FALSE]
                                keep_n <- 30L; prop <- utils::head(prop, keep_n)
                                n_rows <- as.integer(nrow(prop))
                                aggregate_values$cluster_summary_rows <- n_rows
                                groups_idx <- integer(keep_n)
                                fractions_vec <- numeric(keep_n)
                                group_levels <- unique(idents_vec)
                                idx_map <- seq_along(group_levels)
                                names(idx_map) <- group_levels
                                if (n_rows > 0L) {
                                  ok_clusters <- as.character(prop$cluster[seq_len(n_rows)]) %in% names(idx_map)
                                  groups_idx[seq_len(n_rows)] <- ifelse(ok_clusters, idx_map[as.character(prop$cluster[seq_len(n_rows)])], NA_integer_)
                                  fractions_vec[seq_len(n_rows)] <- round(as.numeric(prop$fraction[seq_len(n_rows)]), 4L)
                                }
                                if (keep_n > n_rows) {
                                  groups_idx[seq(n_rows + 1L, keep_n)] <- NA_integer_
                                  fractions_vec[seq(n_rows + 1L, keep_n)] <- NA_real_
                                }
                                aggregate_values$cluster_summary_group_index <- as.integer(groups_idx)
                                aggregate_values$cluster_summary_fractions <- as.numeric(fractions_vec)
                            }
                        }
                        if (!is.null(score_vec)) {
                            score_num <- as.numeric(score_vec)
                            aggregate_values$mean_score <- round(mean(score_num, na.rm = TRUE), 4L)
                            qs <- round(stats::quantile(score_num, c(0.25, 0.5, 0.75), na.rm = TRUE), 4L)
                            aggregate_values$score_quantile_25 <- as.numeric(qs[[1L]])
                            aggregate_values$score_quantile_50 <- as.numeric(qs[[2L]])
                            aggregate_values$score_quantile_75 <- as.numeric(qs[[3L]])
                        }
                        aggregate_values$annotation_method_code <- 2L
                        aggregate_values$annotation_scope_code <- 3L
                        aggregate_values$raw_values_included <- FALSE
                        evidence <- list(
                            id = paste0("ev:annotation_", name),
                            kind = "deterministic_summary",
                            values = aggregate_values,
                            claim_level = "consistent_with"
                        )
                        result <- tryCatch(RecordAIEvidence(
                            result,
                            evidence,
                            source = name,
                            parents = if (!is.null(ann_rec$id)) as.character(ann_rec$id) else character(),
                            scope = NULL
                        ), error = function(e) result)
                    }
                    result
                },
                input_schema = list(
                    ref = list(type = "object", required = TRUE),
                    labels = list(type = "object", required = TRUE),
                    name = "character",
                    layer = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    ref_arg <- params$ref
                    labels_arg <- params$labels
                    if (missing(ref_arg) || is.null(ref_arg)) {
                        return("reference_missing: annotation requires an explicit 'ref' parameter: either a celldex dataset name (e.g. HumanPrimaryCellAtlasData) or a SummarizedExperiment reference. Do not rely on the package default of HumanPrimaryCellAtlasData.")
                    }
                    ref_is_name <- is.character(ref_arg) && length(ref_arg) == 1L && nzchar(ref_arg)
                    ref_is_object <- inherits(ref_arg, "SummarizedExperiment") ||
                        (is.matrix(ref_arg) && !is.null(rownames(ref_arg)) && !is.null(colnames(ref_arg)))
                    if (!ref_is_name && !ref_is_object) {
                        return("reference_invalid: 'ref' must be either a single non-empty string naming a celldex dataset or a SummarizedExperiment/matrix with dimnames.")
                    }
                    if (missing(labels_arg) || is.null(labels_arg)) {
                        return("labels_missing: annotation requires an explicit 'labels' parameter selecting which labels column of the reference to use (e.g. label.main, label.fine) or a vector of per-reference-sample labels.")
                    }
                    labels_is_name <- is.character(labels_arg) && length(labels_arg) == 1L && nzchar(labels_arg)
                    labels_is_vec <- is.atomic(labels_arg) && length(labels_arg) > 0L &&
                        (is.character(labels_arg) || is.factor(labels_arg))
                    if (!labels_is_name && !labels_is_vec) {
                        return("labels_invalid: 'labels' must be either a single non-empty string identifying a reference colData column or one label per reference sample (character or factor vector).")
                    }
                    if (labels_is_vec && ref_is_object) {
                        ref_n <- if (is.matrix(ref_arg)) ncol(ref_arg) else ncol(ref_arg)
                        if (length(labels_arg) != ref_n) {
                            return(paste0("labels_mismatch: length(labels)=", length(labels_arg), " does not match ncol(ref)=", ref_n))
                        }
                    }
                    real_req_ns <- base::requireNamespace
                    if (!real_req_ns("SingleR", quietly = TRUE)) {
                        return("optional_package_missing: SingleR is required for run_annotation; install via BiocManager::install('SingleR')")
                    }
                    if (ref_is_name && grepl("Data$", ref_arg) && !real_req_ns("celldex", quietly = TRUE)) {
                        return(paste0("optional_package_missing: celldex is required to resolve ref='", ref_arg, "'; install via BiocManager::install('celldex')"))
                    }
                    TRUE
                },
                returns = "sce",
                requires_confirmation = TRUE,
                output_schema = list(
                    required_states = list(annotation = "any"),
                    note = "Predicted labels are written into colData columns <name>_labels / <name>_pruned.labels and are NEVER used to overwrite Idents(object) or any existing cluster identity."
                ),
                mutates_object = TRUE,
                allowed_state_types = c("annotation", "mapping", "ai_evidence"),
                estimated_cost = "medium",
                idempotent = FALSE
            )
        ))
    }
    if ("rare_cell" %in% include) {
        actions <- c(actions, list(
            run_doublet_detection = AIAction(
                name = "run_doublet_detection",
                description = "Score doublets with scDblFinder and store per-cell doublet class/score in colData.",
                handler = function(object, params) {
                    RunDoubletFinder(object)
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
                    required_states = list(preprocess = "any"),
                    note = "Per-cell doublet calls are added to colData. No cell is removed or filtered by this action."
                ),
                mutates_object = TRUE,
                allowed_state_types = c("preprocess"),
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
                    name = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    reduction <- params$reduction %||% "PCA"
                    if (!is.character(reduction) || length(reduction) != 1L || !nzchar(reduction)) {
                        return("reduction_invalid: 'reduction' must be a single non-empty reduction name")
                    }
                    available <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
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
    if ("trajectory" %in% include) {
        actions <- c(actions, list(
            run_trajectory = AIAction(
                name = "run_trajectory",
                description = "Infer a slingshot trajectory from an explicitly supplied start cluster and record a bounded aggregate evidence node.",
                handler = function(object, params) {
                    name <- params$name %||% "slingshot"
                    result <- RunSlingshot(
                        object,
                        group = params$group,
                        reduction = params$reduction %||% "UMAP",
                        start_cluster = params$start_cluster,
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
                        start_cluster = params$start_cluster,
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
                    start <- params$start_cluster
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
                    available <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
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
    AIExecutionRegistry(actions)
}
sclet_ai_record_trajectory_evidence <- function(object, name, group, start_cluster, reduction) {
    record <- sclet_get_state_record(object, "trajectory", name)
    if (is.null(record)) return(object)
    groups <- as.character(SummarizedExperiment::colData(object)[[group]])
    group_levels <- unique(stats::na.omit(groups))
    n_groups <- length(group_levels)
    start_index <- match(as.character(start_cluster), group_levels)
    # Only the numeric per-lineage columns (slingPseudotime_1, slingPseudotime_2, ...)
    # are summarized here. The plain "slingPseudotime" column is a per-cell
    # data.frame and is deliberately never read into evidence.
    pseudotime_columns <- grep("^slingPseudotime_[0-9]+$",
        colnames(SummarizedExperiment::colData(object)), value = TRUE)
    lineage_summaries <- list()
    n_lineages <- as.integer(record$summary$n_lineages %||% 0L)
    for (column in pseudotime_columns) {
        values <- as.numeric(SummarizedExperiment::colData(object)[[column]])
        values <- values[is.finite(values)]
        if (!length(values)) next
        quantiles <- round(unname(stats::quantile(values, c(0.25, 0.5, 0.75))), 4L)
        lineage_summaries[[column]] <- list(
            n_observed = as.integer(length(values)),
            mean = round(mean(values), 4L),
            quantile_25 = as.numeric(quantiles[[1L]]),
            median = as.numeric(quantiles[[2L]]),
            quantile_75 = as.numeric(quantiles[[3L]]),
            max = round(max(values), 4L)
        )
    }
    values <- list(
        n_lineages = as.integer(n_lineages),
        n_groups = as.integer(n_groups),
        start_group_code = as.integer(start_index),
        group_codes = paste0("cluster_", seq_len(n_groups)),
        start_group = paste0("cluster_", start_index),
        relative_ordering_only = TRUE,
        pseudotime_is_absolute_time = FALSE,
        raw_values_included = FALSE
    )
    if (length(lineage_summaries)) {
        first <- lineage_summaries[[1L]]
        values$n_summarized_lineages <- as.integer(length(lineage_summaries))
        values$first_lineage_median <- first$median
        values$first_lineage_iqr <- round(first$quantile_75 - first$quantile_25, 4L)
        values$first_lineage_max <- first$max
    }
    evidence <- list(
        id = paste0("ev:trajectory_", name),
        kind = "deterministic_summary",
        values = values,
        claim_level = "consistent_with"
    )
    tryCatch(
        RecordAIEvidence(object, evidence, source = name, parents = character(), scope = NULL),
        error = function(e) object
    )
}

sclet_ai_record_rare_cell_evidence <- function(object, name, rare_threshold) {
    summary <- summarize_small_cluster_evidence(object, cluster = "rare_cluster", size_threshold = rare_threshold)
    if (!identical(summary$status, "available") || !length(summary$clusters)) {
        return(object)
    }
    notes <- character()
    for (entry in summary$clusters) {
        label <- entry$cluster
        n_signals <- as.integer(entry$n_independent_signals_available)
        if (n_signals < 1L) {
            notes <- c(notes, paste0(
                "no_evidence_recorded: population ", label, " (size ", entry$size,
                ") has no independent signal available; cluster size alone cannot establish that this population is ",
                "real. Run doublet detection, marker tests, or add sample metadata before making any claim."
            ))
            next
        }
        signal <- entry$independent_signals
        values <- list(
            population_label = label,
            population_size = as.integer(entry$size),
            population_fraction = as.numeric(entry$fraction_of_total),
            n_independent_signals = n_signals,
            qc_available = isTRUE(signal$qc$available),
            doublet_available = isTRUE(signal$doublet$available),
            marker_available = isTRUE(signal$marker$available),
            sample_replication_available = isTRUE(signal$sample_replication$available),
            low_confidence = n_signals < 2L,
            raw_values_included = FALSE
        )
        if (isTRUE(signal$doublet$available)) {
            values$doublet_fraction <- as.numeric(signal$doublet$doublet_fraction)
            values$doublet_fraction_other_cells <- as.numeric(signal$doublet$doublet_fraction_other_cells)
        }
        if (isTRUE(signal$qc$available)) {
            values$qc_metric_count <- as.integer(signal$qc$n_qc_metrics)
            values$qc_max_eta_squared <- as.numeric(signal$qc$max_eta_squared)
        }
        if (isTRUE(signal$sample_replication$available)) {
            values$replicated_across_groups <- isTRUE(signal$sample_replication$present_in_multiple_samples)
            values$group_presence_count <- as.integer(signal$sample_replication$n_samples_present)
        }
        if (isTRUE(signal$marker$available)) {
            values$marker_evidence_present <- TRUE
            values$marker_up_gene_count <- as.integer(signal$marker$n_significant_up_genes %||% 0L)
        }
        evidence <- list(
            id = paste0("ev:rare_", name, "_", label),
            kind = "deterministic_summary",
            values = values,
            claim_level = if (n_signals >= 2L) "consistent_with" else "associated"
        )
        recorded <- tryCatch(
            RecordAIEvidence(object, evidence, source = name, parents = character(), scope = NULL),
            error = function(e) e
        )
        if (inherits(recorded, "error")) {
            notes <- c(notes, paste0(
                "evidence_record_failed: population ", label, " was not registered as evidence (",
                base::conditionMessage(recorded), "); no claim was recorded for it."
            ))
        } else {
            object <- recorded
        }
    }
    if (length(notes)) {
        attr(object, "sclet_ai_note") <- notes
    }
    object
}

#' Build an allowlisted AI execution registry
#'
#' An empty registry is the default. Actions are never discovered from the
#' global environment or invoked by name unless explicitly registered here.
#'
#' @param actions Named list of `AIAction()` descriptors.
#' @return An object of class `sclet_ai_execution_registry`.
#' @export
AIExecutionRegistry <- function(actions = list()) {
    if (is.null(actions)) {
        actions <- list()
    }
    if (!is.list(actions)) {
        stop("`actions` must be a list of AIAction descriptors.")
    }
    normalized <- list()
    for (i in seq_along(actions)) {
        action <- actions[[i]]
        if (!inherits(action, "sclet_ai_action")) {
            if (is.list(action) && is.function(action$handler)) {
                action <- do.call(
                    AIAction,
                    c(
                        list(
                            name = action$name %||% names(actions)[[i]],
                            handler = action$handler
                        ),
                        action[setdiff(names(action), c("name", "handler"))]
                    )
                )
            } else {
                stop("Every execution registry entry must be an AIAction descriptor.")
            }
        }
        name <- action$name
        if (is.null(name) || !nzchar(name)) {
            stop("Every execution action must have a name.")
        }
        if (!is.null(normalized[[name]])) {
            stop("Duplicate execution action: ", name)
        }
        normalized[[name]] <- action
    }
    class(normalized) <- c("sclet_ai_execution_registry", "list")
    normalized
}

sclet_ai_action_output_problems <- function(object, descriptor) {
    schema <- descriptor$output_schema %||% list()
    problems <- character()
    if (length(schema$required_assays)) {
        missing <- setdiff(as.character(schema$required_assays), SummarizedExperiment::assayNames(object))
        if (length(missing)) problems <- c(problems, paste0("missing output assay(s): ", paste(missing, collapse = ", ")))
    }
    if (length(schema$required_reductions)) {
        missing <- setdiff(as.character(schema$required_reductions), SingleCellExperiment::reducedDimNames(object))
        if (length(missing)) problems <- c(problems, paste0("missing output reduction(s): ", paste(missing, collapse = ", ")))
    }
    if (length(schema$required_graphs)) {
        missing <- vapply(as.character(schema$required_graphs), function(name) {
            is.null(sclet_get_graph(object, name))
        }, logical(1))
        if (any(missing)) problems <- c(problems, paste0("missing output graph(s): ", paste(as.character(schema$required_graphs)[missing], collapse = ", ")))
    }
    if (isTRUE(schema$required_hvg) && is.null(sclet_get_hvg_nfeatures(object))) {
        problems <- c(problems, "highly variable feature state was not registered")
    }
    if (length(schema$active_assay) && !identical(sclet_get_active_assay(object), schema$active_assay)) {
        problems <- c(problems, paste0("active assay is not ", schema$active_assay))
    }
    if (length(schema$active_reduction) && !identical(DefaultReduction(object), schema$active_reduction)) {
        problems <- c(problems, paste0("active reduction is not ", schema$active_reduction))
    }
    if (length(schema$active_ident) && !identical(ActiveIdent(object), schema$active_ident)) {
        problems <- c(problems, paste0("active identity is not ", schema$active_ident))
    }
    if (length(schema$required_commands)) {
        commands <- sclet_get_commands(object)
        command_names <- vapply(commands, function(x) x$command %||% "", character(1))
        missing <- setdiff(as.character(schema$required_commands), command_names)
        if (length(missing)) problems <- c(problems, paste0("missing output command(s): ", paste(missing, collapse = ", ")))
    }
    if (length(schema$required_states)) {
        for (type in names(schema$required_states)) {
            expected <- as.character(schema$required_states[[type]])
            records <- tryCatch(sclet_get_state_records(object, type), error = function(e) list())
            if (identical(expected, "any")) {
                if (!length(records)) {
                    problems <- c(problems, paste0("missing output state(s) for ", type, ": any"))
                }
            } else if (identical(expected, "any_qualified")) {
                if (!length(records)) {
                    problems <- c(problems, paste0("missing output state(s) for ", type, ": any_qualified"))
                }
            } else {
                missing <- setdiff(expected, names(records))
                if (length(missing)) problems <- c(problems, paste0("missing output state(s) for ", type, ": ", paste(missing, collapse = ", ")))
            }
        }
    }
    unique(problems)
}
sclet_ai_action_state_problems <- function(before, after, descriptor) {
    before_types <- names(sclet_get_state(before)$states$records %||% list())
    after_types <- names(sclet_get_state(after)$states$records %||% list())
    new_types <- setdiff(after_types, before_types)
    allowed <- descriptor$allowed_state_types %||% character()
    if (length(new_types) && !all(new_types %in% allowed)) {
        return(paste0("action registered unexpected state type(s): ", paste(setdiff(new_types, allowed), collapse = ", ")))
    }
    character()
}

sclet_ai_execution_summary <- function(value, descriptor = NULL) {
    if (inherits(value, "SingleCellExperiment")) {
        summary <- list(
            class = class(value),
            n_features = nrow(value),
            n_cells = ncol(value),
            assays = SummarizedExperiment::assayNames(value),
            reductions = SingleCellExperiment::reducedDimNames(value),
            active_assay = sclet_get_active_assay(value),
            active_reduction = tryCatch(DefaultReduction(value), error = function(e) NULL),
            fingerprint = tryCatch(GetAnalysisLedger(value)$fingerprint, error = function(e) NULL)
        )
        if (!is.null(descriptor$output_schema)) {
            contract <- descriptor$output_schema
            for (name in c("assay", "reduction", "graph", "state")) {
                if (!is.null(contract[[name]])) summary[[name]] <- contract[[name]]
            }
        }
        return(summary)
    }
    if (is.null(value) || is.atomic(value) && length(value) <= 50L) {
        return(value)
    }
    list(class = class(value), length = length(value))
}

sclet_ai_execution_id <- function(plan_id) {
    safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(plan_id))
    paste0("ai_execution_", safe_id, "_", format(Sys.time(), "%Y%m%d%H%M%S"))
}

sclet_ai_record_execution <- function(object, plan, results, status, dry_run) {
    execution_id <- sclet_ai_execution_id(plan$plan_id)
    record <- list(
        id = execution_id,
        type = "ai_execution",
        method = "ExecuteAIPlan",
        inputs = list(
            plan_id = plan$plan_id,
            context_fingerprint = plan$context_fingerprint,
            action_ids = vapply(plan$actions, function(x) x$id, character(1))
        ),
        params = list(
            dry_run = isTRUE(dry_run),
            requires_confirmation = isTRUE(plan$requires_confirmation)
        ),
        artifacts = list(results = results),
        summary = list(
            status = status,
            n_actions = length(plan$actions),
            n_completed = sum(vapply(results, function(x) identical(x$status, "completed"), logical(1))),
            n_failed = sum(vapply(results, function(x) identical(x$status, "failed"), logical(1)))
        ),
        created_at = Sys.time()
    )
    list(object = sclet_set_analysis(object, execution_id, record), id = execution_id)
}

sclet_ai_assess_success_criteria <- function(plan, outputs, object_after, before_evidence_ids) {
    criteria <- sclet_ai_normalize_success_criteria(plan$success_criteria %||% list())
    if (!length(criteria)) return(list())
    after_evidence <- tryCatch(
        sclet_ai_evidence_get_all(object_after),
        error = function(e) list()
    )
    new_evidence_ids <- setdiff(names(after_evidence), before_evidence_ids)
    step_ids <- vapply(plan$actions %||% list(), function(s) as.character(s$id %||% ""), character(1))
    lapply(criteria, function(crit) {
        if (identical(crit$source, "action_output")) {
            step_output <- outputs[[crit$step_id]]
            if (is.null(step_output)) {
                return(list(
                    id = crit$id,
                    description = crit$description,
                    status = "not_available",
                    reason = paste("step", crit$step_id, "did not complete"),
                    observed_value = NULL
                ))
            }
            if (is.list(step_output) && identical(step_output$status, "failed")) {
                return(list(
                    id = crit$id,
                    description = crit$description,
                    status = "not_available",
                    reason = paste("step", crit$step_id, "did not complete"),
                    observed_value = NULL
                ))
            }
            if (identical(crit$check, "exists")) {
                if (!is.null(crit$field) && nzchar(crit$field)) {
                    extracted <- sclet_ai_extract_dotted_field(step_output, crit$field)
                    if (is.null(extracted)) {
                        return(list(
                            id = crit$id,
                            description = crit$description,
                            status = "not_met",
                            reason = paste("field", crit$field, "not found in output of step", crit$step_id),
                            observed_value = NULL
                        ))
                    }
                    return(list(
                        id = crit$id,
                        description = crit$description,
                        status = "met",
                        reason = paste("field", crit$field, "exists in output of step", crit$step_id),
                        observed_value = extracted
                    ))
                }
                return(list(
                    id = crit$id,
                    description = crit$description,
                    status = "met",
                    reason = paste("output of step", crit$step_id, "exists"),
                    observed_value = step_output
                ))
            }
            extracted <- sclet_ai_extract_dotted_field(step_output, crit$field)
            if (is.null(extracted) || length(extracted) == 0L) {
                return(list(
                    id = crit$id,
                    description = crit$description,
                    status = "not_available",
                    reason = paste("field", crit$field, "not found in output of step", crit$step_id),
                    observed_value = NULL
                ))
            }
            comparison <- tryCatch(
                sclet_ai_compare_criterion(extracted, crit$check, crit$value),
                error = function(e) list(status = "not_available", reason = conditionMessage(e))
            )
            if (is.na(comparison$status) || is.null(comparison$status)) {
                comparison$status <- "not_available"
                comparison$reason <- "comparison produced NA"
            }
            list(
                id = crit$id,
                description = crit$description,
                status = comparison$status,
                reason = comparison$reason,
                observed_value = extracted
            )
        } else if (identical(crit$source, "evidence")) {
            matched_id <- intersect(crit$evidence_id, new_evidence_ids)
            if (!length(matched_id)) {
                return(list(
                    id = crit$id,
                    description = crit$description,
                    status = "not_available",
                    reason = paste("no evidence node found for", crit$evidence_id),
                    observed_value = NULL
                ))
            }
            evidence_node <- after_evidence[[matched_id[[1L]]]]
            if (identical(crit$check, "exists")) {
                if (!is.null(crit$field) && nzchar(crit$field)) {
                    extracted <- sclet_ai_extract_dotted_field(evidence_node, crit$field)
                    if (is.null(extracted)) {
                        return(list(
                            id = crit$id,
                            description = crit$description,
                            status = "not_met",
                            reason = paste("field", crit$field, "not found in evidence", crit$evidence_id),
                            observed_value = NULL
                        ))
                    }
                    return(list(
                        id = crit$id,
                        description = crit$description,
                        status = "met",
                        reason = paste("field", crit$field, "exists in evidence", crit$evidence_id),
                        observed_value = extracted
                    ))
                }
                return(list(
                    id = crit$id,
                    description = crit$description,
                    status = "met",
                    reason = paste("evidence node", crit$evidence_id, "exists"),
                    observed_value = evidence_node
                ))
            }
            extracted <- sclet_ai_extract_dotted_field(evidence_node$values %||% evidence_node, crit$field)
            if (is.null(extracted) || length(extracted) == 0L) {
                return(list(
                    id = crit$id,
                    description = crit$description,
                    status = "not_available",
                    reason = paste("field", crit$field, "not found in evidence", crit$evidence_id),
                    observed_value = NULL
                ))
            }
            comparison <- tryCatch(
                sclet_ai_compare_criterion(extracted, crit$check, crit$value),
                error = function(e) list(status = "not_available", reason = conditionMessage(e))
            )
            if (is.na(comparison$status) || is.null(comparison$status)) {
                comparison$status <- "not_available"
                comparison$reason <- "comparison produced NA"
            }
            list(
                id = crit$id,
                description = crit$description,
                status = comparison$status,
                reason = comparison$reason,
                observed_value = extracted
            )
        } else {
            list(
                id = crit$id,
                description = crit$description,
                status = "not_available",
                reason = paste("unknown source type:", crit$source),
                observed_value = NULL
            )
        }
    })
}

sclet_ai_compare_criterion <- function(observed, check, expected) {
    if (identical(check, "equals")) {
        if (identical(observed, expected)) {
            return(list(status = "met", reason = paste("observed value equals expected")))
        }
        return(list(status = "not_met", reason = paste("observed value does not equal expected")))
    }
    if (identical(check, "gte")) {
        if (is.numeric(observed) && is.numeric(expected) && !is.na(observed) && !is.na(expected)) {
            if (observed >= expected) {
                return(list(status = "met", reason = paste("observed", observed, ">= expected", expected)))
            }
            return(list(status = "not_met", reason = paste("observed", observed, "< expected", expected)))
        }
        return(list(status = "not_available", reason = "values are not numeric or contain NA"))
    }
    if (identical(check, "lte")) {
        if (is.numeric(observed) && is.numeric(expected) && !is.na(observed) && !is.na(expected)) {
            if (observed <= expected) {
                return(list(status = "met", reason = paste("observed", observed, "<= expected", expected)))
            }
            return(list(status = "not_met", reason = paste("observed", observed, "> expected", expected)))
        }
        return(list(status = "not_available", reason = "values are not numeric or contain NA"))
    }
    if (identical(check, "in")) {
        if (is.null(expected) || !length(expected)) {
            return(list(status = "not_available", reason = "expected set is empty"))
        }
        if (any(is.na(observed))) {
            return(list(status = "not_available", reason = "observed value contains NA"))
        }
        if (all(observed %in% expected)) {
            return(list(status = "met", reason = "observed value is in expected set"))
        }
        return(list(status = "not_met", reason = "observed value is not in expected set"))
    }
    list(status = "not_available", reason = paste("unknown check type:", check))
}

#' Execute a validated AI analysis plan under explicit safety controls
#'
#' The default is a side-effect-free dry run. Non-dry execution requires the
#' issued confirmation token returned by `ValidateAIPlan()` when the plan
#' contains a mutating or confirmation-required action, and invokes only
#' actions present in the supplied `AIExecutionRegistry()`.
#'
#' @param object A `SingleCellExperiment` object.
#' @param plan An `sclet_ai_plan` or compatible plan list.
#' @param registry An `AIExecutionRegistry`.
#' @param validation Optional result from `ValidateAIPlan()`.
#' @param dry_run Logical. If `TRUE`, do not invoke handlers or modify the object.
#' @param confirmation Confirmation token returned by `ValidateAIPlan()`.
#' @param record Logical. Record execution results in the analysis ledger.
#' @return An object of class `sclet_ai_execution` containing the resulting
#'   object, action results, execution status, and `success_assessment` (a list
#'   of per-criterion evaluation results; empty for dry runs). The
#'   `success_assessment` field is purely observational and does not alter
#'   `status`.
#' @export
ExecuteAIPlan <- function(
    object,
    plan,
    registry,
    validation = NULL,
    dry_run = TRUE,
    confirmation = NULL,
    record = TRUE
) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    if (missing(registry)) {
        registry <- AIExecutionRegistry()
    }
    if (!inherits(registry, "sclet_ai_execution_registry")) {
        registry <- AIExecutionRegistry(registry)
    }
    plan$actions <- sclet_ai_normalize_plan_actions(plan$actions %||% list())
    if (is.null(validation)) {
        validation <- ValidateAIPlan(plan, object = object, registry = registry, strict = TRUE)
    } else {
        if (!isTRUE(validation$valid)) {
            sclet_ai_error("sclet_ai_invalid_plan", "The supplied plan validation is not valid.")
        }
        current_fingerprint <- GetAnalysisLedger(object)$fingerprint
        if (!identical(as.character(validation$plan_id), as.character(plan$plan_id)) ||
            (!is.null(validation$context_fingerprint) &&
                !identical(as.character(validation$context_fingerprint), as.character(current_fingerprint)))) {
            sclet_ai_error(
                "sclet_ai_stale_plan",
                "The supplied plan validation does not match the current object or plan. Revalidate before execution."
            )
        }
    }
    if (isTRUE(dry_run)) {
        results <- lapply(plan$actions, function(step) {
            list(
                id = step$id,
                action = step$action,
                status = "planned",
                params = step$params,
                requires_confirmation = isTRUE(step$requires_confirmation)
            )
        })
        return(structure(
            list(
                object = object,
                plan = plan,
                validation = validation,
                results = results,
                status = "dry_run",
                dry_run = TRUE,
                recorded = FALSE,
                execution_id = NULL,
                success_assessment = list()
            ),
            class = c("sclet_ai_execution", "list")
        ))
    }
    if (isTRUE(validation$requires_confirmation) &&
        !identical(as.character(confirmation), as.character(validation$confirmation_token))) {
        sclet_ai_error(
            "sclet_ai_confirmation_required",
            "Non-dry execution requires the confirmation token returned by ValidateAIPlan()."
        )
    }

    current <- object
    outputs <- list()
    results <- list()
    status <- "completed"
    continued_failure <- FALSE
    before_evidence_ids <- tryCatch(
        names(sclet_ai_evidence_get_all(object)),
        error = function(e) character()
    )
    for (step in plan$actions) {
        descriptor <- registry[[step$action]]
        before <- current
        started <- Sys.time()
        resolved_params <- tryCatch(
            sclet_ai_resolve_plan_value(step$params, outputs),
            error = function(e) e
        )
        if (inherits(resolved_params, "error")) {
            status <- "failed"
            results[[length(results) + 1L]] <- list(
                id = step$id,
                action = step$action,
                status = "failed",
                error = conditionMessage(resolved_params),
                duration_sec = 0
            )
            break
        }
        bound_problems <- sclet_ai_action_param_problems(
            resolved_params,
            descriptor$input_schema %||% list()
        )
        if (length(bound_problems)) {
            status <- "failed"
            results[[length(results) + 1L]] <- list(
                id = step$id,
                action = step$action,
                status = "failed",
                error = paste("resolved parameters are invalid:", paste(bound_problems, collapse = "; ")),
                attempts = 0L,
                duration_sec = 0
            )
            break
        }
        attempts <- 0L
        value <- NULL
        while (attempts <= step$max_retries) {
            attempts <- attempts + 1L
            value <- tryCatch(
                do.call(descriptor$handler, list(current, resolved_params)),
                error = function(e) e
            )
            if (!inherits(value, "error")) break
        }
        elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        if (inherits(value, "error")) {
            status <- "failed"
            failure <- list(
                id = step$id,
                action = step$action,
                status = "failed",
                error = conditionMessage(value),
                attempts = attempts,
                duration_sec = elapsed
            )
            results[[length(results) + 1L]] <- failure
            outputs[[step$id]] <- list(status = "failed", error = conditionMessage(value))
            if (isTRUE(step$continue_on_error)) {
                continued_failure <- TRUE
                status <- "completed_with_errors"
                next
            }
            break
        }
        if (identical(descriptor$returns, "sce")) {
            if (!inherits(value, "SingleCellExperiment")) {
                status <- "failed"
                results[[length(results) + 1L]] <- list(
                    id = step$id,
                    action = step$action,
                    status = "failed",
                    error = "registered SCE action did not return a SingleCellExperiment",
                    attempts = attempts,
                    duration_sec = elapsed
                )
                break
            }
            output_problems <- sclet_ai_action_output_problems(value, descriptor)
            state_problems <- sclet_ai_action_state_problems(before, value, descriptor)
            contract_problems <- c(output_problems, state_problems)
            if (length(contract_problems)) {
                status <- "failed"
                results[[length(results) + 1L]] <- list(
                    id = step$id,
                    action = step$action,
                    status = "failed",
                    error = paste(contract_problems, collapse = "; "),
                    attempts = attempts,
                    duration_sec = elapsed
                )
                break
            }
            current <- value
        }
        output_summary <- sclet_ai_execution_summary(value, descriptor)
        outputs[[step$id]] <- output_summary
        results[[length(results) + 1L]] <- list(
            id = step$id,
            action = step$action,
            status = "completed",
            output = output_summary,
            attempts = attempts,
            duration_sec = elapsed
        )
    }
    success_assessment <- sclet_ai_assess_success_criteria(
        plan, outputs, current, before_evidence_ids
    )
    execution_id <- NULL
    recorded <- FALSE
    if (isTRUE(record)) {
        recorded_result <- sclet_ai_record_execution(current, plan, results, status, dry_run = FALSE)
        current <- recorded_result$object
        execution_id <- recorded_result$id
        current <- sclet_log_command(
            current,
            command = "ExecuteAIPlan",
            params = list(plan_id = plan$plan_id, execution_id = execution_id),
            outputs = list(status = status, execution_id = execution_id)
        )
        recorded <- TRUE
    }
    structure(
        list(
            object = current,
            plan = plan,
            validation = validation,
            results = results,
            status = status,
            dry_run = FALSE,
            recorded = recorded,
            execution_id = execution_id,
            success_assessment = success_assessment
        ),
        class = c("sclet_ai_execution", "list")
    )
}

#' Validate and run an AI analysis plan through one controlled workflow
#'
#' @param object A `SingleCellExperiment` object.
#' @param plan An `sclet_ai_plan` or compatible plan list.
#' @param registry An `AIExecutionRegistry`.
#' @param dry_run Logical. Defaults to `TRUE`.
#' @param confirm Logical or character. If `TRUE`, use the token generated by
#'   validation; a character value is used as an explicit token.
#' @param record Logical. Record non-dry execution results in the ledger.
#' @return An `sclet_ai_execution` object.
#' @export
RunAIPlan <- function(
    object,
    plan,
    registry,
    dry_run = TRUE,
    confirm = FALSE,
    record = TRUE
) {
    validation <- ValidateAIPlan(
        plan,
        object = object,
        registry = registry,
        strict = TRUE
    )
    token <- if (isTRUE(confirm)) {
        validation$confirmation_token
    } else if (is.character(confirm) && length(confirm) == 1L) {
        confirm
    } else {
        NULL
    }
    ExecuteAIPlan(
        object = object,
        plan = plan,
        registry = registry,
        validation = validation,
        dry_run = dry_run,
        confirmation = token,
        record = record
    )
}
#' @export
print.sclet_ai_execution <- function(x, ...) {
    cat("sclet AI execution [", x$status, "]\n", sep = "")
    cat("Actions:", length(x$results), " Dry run:", isTRUE(x$dry_run), "\n", sep = "")
    if (!is.null(x$execution_id)) {
        cat("Execution record:", x$execution_id, "\n", sep = "")
    }
    invisible(x)
}
