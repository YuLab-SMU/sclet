sclet_ai_plan_id <- function(task = "analysis_plan") {
    stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    suffix <- paste(sample(c(letters, 0:9), 8L, replace = TRUE), collapse = "")
    paste(c("plan", task, stamp, suffix), collapse = "_")
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
            requires_confirmation = if (is.null(action$requires_confirmation)) TRUE else isTRUE(action$requires_confirmation)
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
            planned_fingerprint <- plan$context_fingerprint
            if (!is.null(planned_fingerprint) &&
                !identical(as.character(planned_fingerprint), as.character(current_fingerprint))) {
                errors <- c(errors, "plan context fingerprint does not match the current object")
            }
            for (step in normalized) {
                descriptor <- registry[[step$action]]
                if (is.null(descriptor) || !is.function(descriptor$prerequisites)) {
                    next
                }
                check <- tryCatch(
                    descriptor$prerequisites(object, step$params),
                    error = function(e) e
                )
                if (inherits(check, "error")) {
                    errors <- c(errors, paste0("prerequisite check failed for ", step$id, ": ", conditionMessage(check)))
                } else if (isFALSE(check)) {
                    errors <- c(errors, paste0("prerequisites not met for action ", step$id))
                } else if (is.character(check) && length(check)) {
                    errors <- c(errors, paste0("prerequisites not met for action ", step$id, ": ", paste(check, collapse = "; ")))
                }
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
#' @param structured_output Logical. Use the tool-loop plus structured-output summary stage.
#' @param fallback_on_structure_error Logical. Retain the tool-loop response if structured output fails.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return An object of class `sclet_ai_plan`.
#' @export
AIPlanAnalysis <- function(
    object,
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
        context_fingerprint = result$context$fingerprint %||% context$fingerprint,
        rationale = result$answer,
        requires_confirmation = TRUE,
        ai_result = result,
        metadata = list(
            provider = result$metadata$provider %||% NULL,
            structured_output = result$metadata$structured_output %||% NULL
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
