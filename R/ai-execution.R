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

#' Build the safe built-in read-only action catalog
#'
#' The returned actions never mutate the supplied SCE and do not require a
#' confirmation token. They expose only bounded, deterministic ledger views.
#'
#' @param object A `SingleCellExperiment` object.
#' @return A `sclet_ai_execution_registry` containing read-only actions.
#' @export
AIDefaultExecutionRegistry <- function(object) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    AIExecutionRegistry(list(
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

sclet_ai_execution_summary <- function(value) {
    if (inherits(value, "SingleCellExperiment")) {
        return(list(
            class = class(value),
            n_features = nrow(value),
            n_cells = ncol(value),
            assays = SummarizedExperiment::assayNames(value),
            fingerprint = tryCatch(GetAnalysisLedger(value)$fingerprint, error = function(e) NULL)
        ))
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
    results <- list()
    status <- "completed"
    for (step in plan$actions) {
        descriptor <- registry[[step$action]]
        started <- Sys.time()
        value <- tryCatch(
            do.call(descriptor$handler, list(current, step$params)),
            error = function(e) e
        )
        elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        if (inherits(value, "error")) {
            status <- "failed"
            results[[length(results) + 1L]] <- list(
                id = step$id,
                action = step$action,
                status = "failed",
                error = conditionMessage(value),
                duration_sec = elapsed
            )
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
                    duration_sec = elapsed
                )
                break
            }
            current <- value
        }
        results[[length(results) + 1L]] <- list(
            id = step$id,
            action = step$action,
            status = "completed",
            output = sclet_ai_execution_summary(value),
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

#' @export
print.sclet_ai_execution <- function(x, ...) {
    cat("sclet AI execution [", x$status, "]\n", sep = "")
    cat("Actions:", length(x$results), " Dry run:", isTRUE(x$dry_run), "\n", sep = "")
    if (!is.null(x$execution_id)) {
        cat("Execution record:", x$execution_id, "\n", sep = "")
    }
    invisible(x)
}
