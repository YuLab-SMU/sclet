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
            object = is.list(value),
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
#'   `"read"`, `"preprocess"`, `"dimred"`, `"graph"`, `"cluster"`, and
#'   `"all"`. Defaults to `"read"`.
#' @return A `sclet_ai_execution_registry` containing the requested actions.
#' @export
AIDefaultExecutionRegistry <- function(object, include = "read") {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    include <- unique(as.character(include))
    allowed_groups <- c("read", "preprocess", "dimred", "graph", "cluster", "all")
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
    AIExecutionRegistry(actions)
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
            missing <- setdiff(expected, names(records))
            if (length(missing)) problems <- c(problems, paste0("missing output state(s) for ", type, ": ", paste(missing, collapse = ", ")))
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
#'   object, action results, and execution status.
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
                execution_id = NULL
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
            execution_id = execution_id
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
