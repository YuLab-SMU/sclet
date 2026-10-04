sclet_ai_collect_plan_references <- function(value) {
    if (is.character(value) && length(value) == 1L && grepl("^\\$\\{[^}]+\\}$", value)) {
        token <- sub("^\\$\\{", "", sub("\\}$", "", value))
        return(strsplit(token, "\\.", fixed = FALSE)[[1L]])
    }
    if (is.list(value)) {
        return(unlist(lapply(value, sclet_ai_collect_plan_references), use.names = FALSE))
    }
    character()
}

sclet_ai_resolve_plan_value <- function(value, outputs) {
    if (is.character(value) && length(value) == 1L && grepl("^\\$\\{[^}]+\\}$", value)) {
        token <- sub("^\\$\\{", "", sub("\\}$", "", value))
        path <- strsplit(token, "\\.", fixed = FALSE)[[1L]]
        if (length(path) < 3L || !identical(path[[2L]], "output")) {
            stop("Invalid plan output reference: ", value)
        }
        current <- outputs[[path[[1L]]]]
        if (is.null(current)) stop("Plan output is not available: ", value)
        for (field in path[-c(1L, 2L)]) {
            if (!is.list(current) || is.null(current[[field]])) {
                stop("Plan output field is not available: ", value)
            }
            current <- current[[field]]
        }
        return(current)
    }
    if (is.list(value)) {
        return(lapply(value, sclet_ai_resolve_plan_value, outputs = outputs))
    }
    value
}
sclet_ai_plan_id <- function(task = "analysis_plan") {
    stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    suffix <- paste(sample(c(letters, 0:9), 8L, replace = TRUE), collapse = "")
    paste(c("plan", task, stamp, suffix), collapse = "_")
}

sclet_ai_project_plan_value <- function(value, planned_outputs) {
    if (is.character(value) && length(value) == 1L && grepl("^\\$\\{[^}]+\\}$", value)) {
        token <- sub("^\\$\\{", "", sub("\\}$", "", value))
        path <- strsplit(token, "\\.", fixed = FALSE)[[1L]]
        if (length(path) < 3L || !identical(path[[2L]], "output")) return(value)
        current <- planned_outputs[[path[[1L]]]]
        if (is.null(current)) return(value)
        for (field in path[-c(1L, 2L)]) {
            if (!is.list(current) || is.null(current[[field]])) return(value)
            current <- current[[field]]
        }
        return(current)
    }
    if (is.list(value)) return(lapply(value, sclet_ai_project_plan_value, planned_outputs = planned_outputs))
    value
}
sclet_ai_plan_capabilities <- function(object) {
    state <- sclet_get_state(object)
    list(
        assays = SummarizedExperiment::assayNames(object),
        reductions = SingleCellExperiment::reducedDimNames(object),
        reduction_dims = setNames(
            lapply(SingleCellExperiment::reducedDimNames(object), function(name) {
                ncol(SingleCellExperiment::reducedDim(object, name))
            }),
            SingleCellExperiment::reducedDimNames(object)
        ),
        graphs = names(state$graphs %||% list()),
        hvg = !is.null(sclet_get_hvg_nfeatures(object)),
        active_assay = sclet_get_active_assay(object),
        active_reduction = tryCatch(DefaultReduction(object), error = function(e) NULL),
        active_graph = tryCatch(DefaultGraph(object), error = function(e) NULL),
        active_ident = tryCatch(ActiveIdent(object), error = function(e) NULL)
    )
}

sclet_ai_apply_planned_output <- function(planned, descriptor, params) {
    schema <- descriptor$output_schema %||% list()
    planned$assays <- unique(c(planned$assays, as.character(schema$required_assays %||% character())))
    planned$reductions <- unique(c(planned$reductions, as.character(schema$required_reductions %||% character())))
    planned$graphs <- unique(c(planned$graphs, as.character(schema$required_graphs %||% character())))
    if (isTRUE(schema$required_hvg)) planned$hvg <- TRUE
    if (length(schema$active_assay)) planned$active_assay <- schema$active_assay
    if (length(schema$active_reduction)) planned$active_reduction <- schema$active_reduction
    if (length(schema$active_ident)) planned$active_ident <- schema$active_ident
    if (length(schema$active_graph)) planned$active_graph <- schema$active_graph
    if (identical(descriptor$name, "run_pca")) {
        planned$reduction_dims$PCA <- as.integer(params$ncomponents %||% 50)
        planned$active_reduction <- "PCA"
    }
    planned
}

sclet_ai_call_prerequisites <- function(descriptor, object, params, planned) {
    fn <- descriptor$prerequisites
    fn_names <- names(formals(fn))
    if ("planned" %in% fn_names || "..." %in% fn_names) {
        return(fn(object, params, planned = planned))
    }
    fn(object, params)
}
sclet_ai_normalize_plan_actions <- function(actions) {
    if (is.null(actions)) {
        return(list())
    }
    if (!is.list(actions)) {
        actions <- as.list(actions)
    }
    lapply(seq_along(actions), function(i) {
        action <- actions[[i]]
        if (is.character(action) && length(action) == 1L) {
            action <- list(action = action)
        }
        if (!is.list(action)) {
            return(list(
                id = paste0("step_", i),
                action = NULL,
                params = list(),
                depends_on = character(),
                description = as.character(action)[[1L]]
            ))
        }
        name <- action$action %||% action$name %||% action$command
        params <- action$params %||% list()
        if (!is.list(params)) {
            params <- list(value = params)
        }
        depends_on <- action$depends_on %||% action$depends %||% character()
        list(
            id = as.character(action$id %||% paste0("step_", i))[[1L]],
            action = if (is.null(name)) NULL else as.character(name)[[1L]],
            params = params,
            depends_on = as.character(depends_on),
            description = action$description %||% action$rationale %||% NULL,
            expected_outputs = action$expected_outputs %||% list(),
            requires_confirmation = if (is.null(action$requires_confirmation)) TRUE else isTRUE(action$requires_confirmation),
            max_retries = as.integer(action$max_retries %||% 0L),
            continue_on_error = if (is.null(action$continue_on_error)) FALSE else isTRUE(action$continue_on_error)
        )
    })
}

#' Construct an analysis plan returned by the R-native AI planner
#'
#' @param task Short planning task identifier.
#' @param actions A list of normalized action steps.
#' @param context_fingerprint Fingerprint of the ledger used to make the plan.
#' @param rationale Optional natural-language planning rationale.
#' @param requires_confirmation Logical. Whether execution requires explicit confirmation.
#' @param ai_result Optional `sclet_ai_result` that produced the plan.
#' @param metadata Additional plan metadata.
#' @return An object of class `sclet_ai_plan`.
new_sclet_ai_plan <- function(
    task = "analysis_plan",
    actions = list(),
    context_fingerprint = NULL,
    rationale = NULL,
    requires_confirmation = TRUE,
    ai_result = NULL,
    metadata = list()
) {
    if (!is.character(task) || length(task) != 1L || is.na(task) || !nzchar(task)) {
        stop("`task` must be a single non-empty character string.")
    }
    plan <- list(
        plan_id = sclet_ai_plan_id(task),
        task = task,
        rationale = rationale,
        context_fingerprint = context_fingerprint,
        actions = sclet_ai_normalize_plan_actions(actions),
        requires_confirmation = isTRUE(requires_confirmation),
        ai_result = ai_result,
        metadata = metadata
    )
    class(plan) <- c("sclet_ai_plan", "list")
    plan
}

#' Validate an AI-generated analysis plan
#'
#' @param plan An `sclet_ai_plan` or compatible list.
#' @param object Optional `SingleCellExperiment` used for fingerprint and prerequisite checks.
#' @param registry An `AIExecutionRegistry`.
#' @param strict Logical. If `TRUE`, stop on invalid plans; otherwise return validation details.
#' @return An object of class `sclet_ai_plan_validation` with `valid`, `errors`,
#'   `warnings`, and a one-time confirmation token when valid.
#' @export
ValidateAIPlan <- function(
    plan,
    object = NULL,
    registry = AIExecutionRegistry(),
    strict = FALSE
) {
    errors <- character()
    warnings <- character()
    if (!is.list(plan)) {
        errors <- c(errors, "plan must be a list")
        plan <- list(plan_id = NULL, actions = list())
    }
    actions <- plan$actions %||% list()
    if (!is.list(actions)) {
        errors <- c(errors, "plan$actions must be a list")
        actions <- list()
    }
    if (is.null(plan$plan_id) || !is.character(plan$plan_id) ||
        length(plan$plan_id) != 1L || !nzchar(plan$plan_id)) {
        errors <- c(errors, "plan_id must be a non-empty character value")
    }
    if (!inherits(registry, "sclet_ai_execution_registry")) {
        registry <- tryCatch(AIExecutionRegistry(registry), error = function(e) e)
        if (inherits(registry, "error")) {
            errors <- c(errors, paste("invalid execution registry:", conditionMessage(registry)))
            registry <- AIExecutionRegistry()
        }
    }
    normalized <- sclet_ai_normalize_plan_actions(actions)
    ids <- vapply(normalized, function(x) x$id %||% "", character(1))
    if (any(!nzchar(ids))) {
        errors <- c(errors, "every plan action must have a non-empty id")
    }
    if (anyDuplicated(ids)) {
        errors <- c(errors, "plan action ids must be unique")
    }
    if (length(normalized)) {
        known_ids <- character()
        for (i in seq_along(normalized)) {
            step <- normalized[[i]]
            if (is.null(step$action) || !nzchar(step$action)) {
                errors <- c(errors, paste0("action ", step$id, " has no action name"))
            } else if (is.null(registry[[step$action]])) {
                errors <- c(errors, paste0("action is not registered: ", step$action))
            } else {
                descriptor <- registry[[step$action]]
                parameter_problems <- sclet_ai_action_param_problems(
                    step$params,
                    descriptor$input_schema %||% list()
                )
                if (length(parameter_problems)) {
                    errors <- c(errors, paste0(
                        "invalid params for action ", step$id, ": ",
                        paste(parameter_problems, collapse = "; ")
                    ))
                }
                if (is.na(step$max_retries) || step$max_retries < 0L || step$max_retries > 3L) {
                    errors <- c(errors, paste0("action ", step$id, " max_retries must be an integer from 0 to 3"))
                }
                if (step$max_retries > 0L && !isTRUE(descriptor$idempotent)) {
                    errors <- c(errors, paste0("action ", step$id, " is not idempotent and cannot be retried"))
                }
                if (isTRUE(descriptor$mutates_object) &&
                    !isTRUE(descriptor$requires_confirmation)) {
                    errors <- c(errors, paste0(
                        "mutating action ", step$id,
                        " must require confirmation"
                    ))
                }
            }
            missing_dependencies <- setdiff(step$depends_on, ids)
            if (length(missing_dependencies)) {
                errors <- c(errors, paste0(
                    "action ", step$id, " depends on unknown step(s): ",
                    paste(missing_dependencies, collapse = ", ")
                ))
            }
            forward_dependencies <- intersect(step$depends_on, ids[seq.int(i, length(ids))])
            if (length(forward_dependencies)) {
                errors <- c(errors, paste0(
                    "action ", step$id, " depends on a later step: ",
                    paste(forward_dependencies, collapse = ", ")
                ))
            }
            reference_path <- sclet_ai_collect_plan_references(step$params)
            if (length(reference_path)) {
                reference_steps <- unique(reference_path[seq(1L, length(reference_path), by = 3L)])
                unknown_reference_steps <- setdiff(reference_steps, ids)
                if (length(unknown_reference_steps)) {
                    errors <- c(errors, paste0(
                        "action ", step$id, " references unknown step(s): ",
                        paste(unknown_reference_steps, collapse = ", ")
                    ))
                }
                undeclared <- setdiff(reference_steps, step$depends_on)
                if (length(undeclared)) {
                    errors <- c(errors, paste0(
                        "action ", step$id, " references undeclared dependency step(s): ",
                        paste(undeclared, collapse = ", ")
                    ))
                }
            }
            known_ids <- c(known_ids, step$id)
        }
    }

    current_fingerprint <- NULL
    if (!is.null(object)) {
        if (!inherits(object, "SingleCellExperiment")) {
            errors <- c(errors, "object must be a SingleCellExperiment")
        } else {
            current_fingerprint <- tryCatch(
                GetAnalysisLedger(object)$fingerprint,
                error = function(e) NULL
            )
            full_fingerprint <- tryCatch(
                GetAnalysisLedger(object, detail = "full")$fingerprint,
                error = function(e) NULL
            )
            planned_fingerprint <- plan$context_fingerprint
            if (!is.null(planned_fingerprint) && !any(
                identical(as.character(planned_fingerprint), as.character(current_fingerprint)),
                identical(as.character(planned_fingerprint), as.character(full_fingerprint))
            )) {
                errors <- c(errors, "plan context fingerprint does not match the current object")
            }
            planned <- sclet_ai_plan_capabilities(object)
            planned_outputs <- list()
            for (step in normalized) {
                descriptor <- if (!is.null(step$action) && nzchar(step$action)) registry[[step$action]] else NULL
                planned_params <- sclet_ai_project_plan_value(step$params, planned_outputs)
                if (is.null(descriptor) || !is.function(descriptor$prerequisites)) {
                    next
                }
                check <- tryCatch(
                    sclet_ai_call_prerequisites(descriptor, object, planned_params, planned),
                    error = function(e) e
                )
                if (inherits(check, "error")) {
                    errors <- c(errors, paste0("prerequisite check failed for ", step$id, ": ", conditionMessage(check)))
                } else if (isFALSE(check)) {
                    errors <- c(errors, paste0("prerequisites not met for action ", step$id))
                } else if (is.character(check) && length(check)) {
                    errors <- c(errors, paste0("prerequisites not met for action ", step$id, ": ", paste(check, collapse = "; ")))
                }
                planned <- sclet_ai_apply_planned_output(planned, descriptor, planned_params)
                planned_outputs[[step$id]] <- descriptor$output_schema %||% list()
            }
        }
    } else if (length(normalized)) {
        warnings <- c(warnings, "object was not supplied; object-specific prerequisites were not checked")
    }

    errors <- unique(errors[nzchar(errors)])
    valid <- !length(errors)
    requires_confirmation <- any(vapply(normalized, function(step) {
        descriptor <- if (!is.null(step$action) && nzchar(step$action)) {
            registry[[step$action]]
        } else {
            NULL
        }
        !is.null(descriptor) && (isTRUE(descriptor$requires_confirmation) ||
            isTRUE(descriptor$mutates_object))
    }, logical(1)))
    confirmation_token <- if (valid && requires_confirmation) {
        paste0("sclet-confirm-", paste(sample(c(letters, 0:9), 32L, replace = TRUE), collapse = ""))
    } else {
        NULL
    }
    validation <- list(
        valid = valid,
        errors = errors,
        warnings = unique(warnings[nzchar(warnings)]),
        plan_id = plan$plan_id %||% NULL,
        context_fingerprint = current_fingerprint %||% plan$context_fingerprint %||% NULL,
        requires_confirmation = requires_confirmation,
        confirmation_token = confirmation_token,
        plan = plan,
        checked_at = Sys.time()
    )
    class(validation) <- c("sclet_ai_plan_validation", "list")
    if (isTRUE(strict) && !valid) {
        sclet_ai_error(
            "sclet_ai_invalid_plan",
            paste("Invalid AI plan:", paste(errors, collapse = "; "))
        )
    }
    validation
}

#' Ask the AI to propose a conservative, executable analysis plan
#'
#' The function only creates a plan. It never executes proposed actions.
#'
#' @param object A `SingleCellExperiment` object.
#' @param model Optional aisdk model.
#' @param goal Optional natural-language user goal to include in the planning context.
#' @param structured_output Logical. Use the tool-loop plus structured-output summary stage.
#' @param fallback_on_structure_error Logical. Retain the tool-loop response if structured output fails.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return An object of class `sclet_ai_plan`.
#' @export
AIPlanAnalysis <- function(
    object,
    goal = NULL,
    model = NULL,
    structured_output = TRUE,
    fallback_on_structure_error = TRUE,
    ...
) {
    context <- sclet_ai_context(object, target = "planning", detail = "full")
    result <- sclet_ai_call(
        task = "analysis_plan",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        fallback_on_structure_error = fallback_on_structure_error,
        system_prompt = paste(
            "Propose a conservative analysis plan from the supplied ledger.",
            if (!is.null(goal) && length(goal) == 1L && nzchar(trimws(goal))) {
                paste("The user's analysis goal is:", goal)
            } else {
                "No additional user goal was supplied."
            },
            "Return only actions that can be represented by a registered action name",
            "and parameter list. List prerequisites, dependencies, expected outputs,",
            "and uncertainty. Never execute an action or claim that it was executed.",
            sep = "\n"
        ),
        ...
    )
    plan <- new_sclet_ai_plan(
        task = result$task,
        actions = result$proposed_actions,
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        rationale = result$answer,
        requires_confirmation = TRUE,
        ai_result = result,
        metadata = list(
            provider = result$metadata$provider %||% NULL,
            structured_output = result$metadata$structured_output %||% NULL,
            goal = goal %||% NULL
        )
    )
    plan
}

#' @export
print.sclet_ai_plan <- function(x, ...) {
    cat("sclet AI plan [", x$plan_id, "]\n", sep = "")
    cat("Actions:", length(x$actions), " Confirmation required:", isTRUE(x$requires_confirmation), "\n", sep = "")
    if (!is.null(x$rationale) && nzchar(x$rationale)) {
        cat(x$rationale, "\n", sep = "")
    }
    invisible(x)
}
