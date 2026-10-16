#' Split one `${...}` token into its path segments
#'
#' @param value A parameter value, possibly nested.
#' @return A list with one character vector of segments per reference found.
#'   Returning the segments *per reference* rather than one flat vector
#'   matters: a flat vector only aligns to references if every reference
#'   happens to have the same number of segments, and callers that walked it
#'   with a fixed stride then read a field name as a step id (reporting
#'   "references unknown step(s): y" for `${step_4.output.y}`) or read an
#'   empty first segment and reported "references unknown step(s): " with no
#'   name at all.
sclet_ai_collect_plan_references <- function(value) {
    if (is.character(value) && length(value) == 1L && grepl("^\\$\\{[^}]+\\}$", value)) {
        token <- sub("^\\$\\{", "", sub("\\}$", "", value))
        return(list(strsplit(token, "\\.", fixed = FALSE)[[1L]]))
    }
    if (is.list(value)) {
        return(unlist(lapply(value, sclet_ai_collect_plan_references), recursive = FALSE))
    }
    list()
}

#' Is a reference path well formed, i.e. `${step_id.output.field}`?
#'
#' The one definition of the grammar, shared by the malformed-reference
#' detector and by the validator's dependency analysis so the two cannot
#' disagree about which references count.
#'
#' @param path Character vector of path segments.
#' @return A single logical value.
sclet_ai_reference_is_well_formed <- function(path) {
    length(path) >= 3L &&
        nzchar(path[[1L]]) &&
        identical(path[[2L]], "output")
}

#' Find value references that are syntactically malformed
#'
#' The plan value-reference syntax is `${step_id.output.field}`. Only
#' [sclet_ai_collect_plan_references()] is permissive: it splits any
#' `${a.b.c}` token, so a malformed reference still parses into segments.
#' [sclet_ai_project_plan_value()] cannot resolve one and returns it
#' unchanged, and [sclet_ai_resolve_plan_value()] refuses it outright at
#' execution time. Detecting the malformed form during validation turns a
#' silent pass-through -- which surfaces later as an unrelated "reduction is
#' not available: ${step_4.reduction}" from some action's prerequisite --
#' into an error that names the actual mistake.
#'
#' @param value A parameter value, possibly nested.
#' @return Character vector of malformed reference strings (empty if none).
sclet_ai_malformed_plan_references <- function(value) {
    if (is.character(value) && length(value) == 1L && grepl("^\\$\\{[^}]+\\}$", value)) {
        token <- sub("^\\$\\{", "", sub("\\}$", "", value))
        path <- strsplit(token, "\\.", fixed = FALSE)[[1L]]
        if (!sclet_ai_reference_is_well_formed(path)) {
            return(value)
        }
        return(character())
    }
    if (is.list(value)) {
        return(unlist(lapply(value, sclet_ai_malformed_plan_references), use.names = FALSE))
    }
    character()
}

#' Suggest the corrected spelling of a malformed value reference
#'
#' Pairing a complaint with the corrected string is what makes the error
#' usable as feedback: the model can copy the suggestion instead of
#' re-deriving the grammar. `${step_4.reduction}` becomes
#' `${step_4.output.reduction}`.
#'
#' @param value A malformed reference string.
#' @return A suggested corrected reference, or `NA_character_` if the token
#'   has no step segment to build on.
sclet_ai_suggest_reference_fix <- function(value) {
    token <- sub("^\\$\\{", "", sub("\\}$", "", value))
    path <- strsplit(token, "\\.", fixed = FALSE)[[1L]]
    if (!length(path) || !nzchar(path[[1L]])) {
        return(NA_character_)
    }
    step <- path[[1L]]
    rest <- path[-1L]
    if (length(rest) && identical(rest[[1L]], "output")) {
        rest <- rest[-1L]
    }
    if (!length(rest)) {
        rest <- "<field>"
    }
    paste0("${", step, ".output.", paste(rest, collapse = "."), "}")
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

# Stable digest of every plan field that can affect validation or execution.
# `digest` is suggested rather than imported, so retain a deterministic base-R
# fallback for installations without it. Sorting named lists makes harmless
# element-order differences irrelevant while preserving action/list order.
sclet_ai_plan_fingerprint <- function(plan) {
    canonical <- function(value) {
        if (is.list(value)) {
            nms <- names(value)
            if (!is.null(nms)) {
                ord <- order(nms)
                value <- value[ord]
                names(value) <- nms[ord]
            }
            return(lapply(value, canonical))
        }
        if (is.atomic(value)) return(unname(value))
        value
    }
    serialized <- serialize(canonical(plan), NULL, version = 2L)
    if (requireNamespace("digest", quietly = TRUE)) {
        return(paste0("sclet-ai-plan-", digest::digest(serialized, algo = "sha256", serialize = FALSE)))
    }
    path <- tempfile("sclet-ai-plan-")
    on.exit(unlink(path), add = TRUE)
    writeBin(as.raw(serialized), path)
    paste0("sclet-ai-plan-", unname(tools::md5sum(path)))
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
        # A deferred promise is not projected. The static value of such a field
        # is a marker ("hvg") that stands for a set which only exists after the
        # producing step runs, so substituting it here would hand validation a
        # meaningless literal -- `scale_data(features = "hvg")` -- that neither
        # the gene selector nor the handler could act on. Passing the reference
        # through instead lets `prerequisites` recognise it as a promise, and
        # the handler resolves it against the real object at execution time.
        if (sclet_ai_is_deferred_promise(current)) return(value)
        return(current)
    }
    if (is.list(value)) return(lapply(value, sclet_ai_project_plan_value, planned_outputs = planned_outputs))
    value
}

# Is this a static output_schema value standing in for something that only
# exists after execution? The marker is deliberately a scalar string so it
# cannot be mistaken for a real gene vector or a real record id. `hvg` stands
# for the variable-feature set; `detest_record_id` stands for the id
# run_de_test registers (written into the execution summary as an observation,
# so projecting the marker here would hand the enrichment gate a literal that
# no planned-record list can contain -- measured as the reason a chained
# DE -> enrichment plan still failed after detest_id became referenceable).
sclet_ai_is_deferred_promise <- function(value) {
    is.character(value) && length(value) == 1L &&
        (identical(value, "hvg") || identical(value, "detest_record_id"))
}
sclet_ai_plan_capabilities <- function(object) {
    state <- sclet_get_state(object)
    records <- state$states$records %||% list()
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
        active_ident = tryCatch(ActiveIdent(object), error = function(e) NULL),
        # Record ids already on the object, keyed by state type. A step that
        # consumes an upstream record reads this, and steps that produce one
        # add to it as the walk proceeds, so a dependency inside a single plan
        # is visible to validation without executing anything. Without it a
        # chain such as run_de_test -> run_enrichment can never validate,
        # because at validation time no detest record exists yet -- measured by
        # test-ai-single-cell-chain.R.
        state_records = lapply(records, function(by_id) names(by_id %||% list()))
    )
}

# Which analysis-state records a step will register, keyed by state type.
#
# A step that consumes an upstream record needs to see records that earlier
# steps in the SAME plan will produce, because validation runs before anything
# executes. These id rules mirror the corresponding handlers -- run_de_test
# registers "findallmarkers" for an all-cluster test and "detest" for a
# targeted one unless a name is given -- and test-ai-single-cell-chain.R asserts
# the projection equals what the action really registers, so a handler whose id
# rule drifts fails there rather than silently disabling this projection.
sclet_ai_planned_state_records <- function(descriptor, params) {
    if (is.null(descriptor)) return(list())
    name <- descriptor$name
    provided <- params$name
    single <- function(type, default) {
        id <- if (is.character(provided) && length(provided) == 1L && nzchar(provided)) {
            provided
        } else {
            default
        }
        stats::setNames(list(id), type)
    }
    switch(name,
        run_de_test = {
            # Mirrors the handler: an explicit name wins, otherwise the id
            # depends on whether this is an all-cluster marker test.
            default <- if (is.null(params$ident.1)) "findallmarkers" else "detest"
            single("detest", default)
        },
        run_enrichment = single("enrichment", "enrichment"),
        run_velocity = single("velocity", "velocity"),
        run_fate_analysis = single("trajectory", "cellfate"),
        run_doublet_detection = list(preprocess = "doublet_finder"),
        list()
    )
}

sclet_ai_apply_planned_output <- function(planned, descriptor, params) {
    schema <- descriptor$output_schema %||% list()
    planned$assays <- unique(c(planned$assays, as.character(schema$required_assays %||% character())))
    planned$reductions <- unique(c(planned$reductions, as.character(schema$required_reductions %||% character())))
    # `allowed_reductions` lists names the step MAY register (method-dependent
    # -- run_integration registers corrected, HARMONY or the scVI reduction).
    # Planning must know those exist: measured 2026-10-14,
    # `find_neighbors(reduction = "HARMONY")` after a `run_integration` step
    # failed validation with "reduction is not available: HARMONY" although
    # execution really does register it -- nothing read `allowed_reductions`.
    # This is deferred projection, the same treatment as `planned$active_ident`:
    # validation accepts the promise because the step declares it, and the
    # execution-time contract (the prerequisite's real reducedDimNames check)
    # remains the authority if the plan picked a method whose reduction is
    # absent.
    planned$reductions <- unique(c(planned$reductions, as.character(schema$allowed_reductions %||% character())))
    planned$graphs <- unique(c(planned$graphs, as.character(schema$required_graphs %||% character())))
    if (isTRUE(schema$required_hvg)) planned$hvg <- TRUE
    if (length(schema$active_assay)) planned$active_assay <- schema$active_assay
    if (length(schema$active_reduction)) planned$active_reduction <- schema$active_reduction
    if (length(schema$active_ident)) planned$active_ident <- schema$active_ident
    if (length(schema$active_graph)) planned$active_graph <- schema$active_graph
    if (identical(descriptor$name, "run_pca")) {
        planned$reduction_dims$PCA <- sclet_ai_coerce_integer(params$ncomponents, default = 50L)
        planned$active_reduction <- "PCA"
    }
    # Add the records this step will register, so a later step in the same plan
    # can depend on them. Deferred projection again: validation accepts the
    # promise, and if the step really fails to register, its own execution-time
    # contract (or the next step's execution) is still the authority.
    produced <- sclet_ai_planned_state_records(descriptor, params)
    for (type in names(produced)) {
        existing <- planned$state_records[[type]] %||% character()
        planned$state_records[[type]] <- unique(c(existing, as.character(produced[[type]])))
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
#' Coerce a possibly-malformed model-supplied value to a single integer
#'
#' Model output is untrusted and may put a non-numeric string, `NA`, or a
#' length-0 value where an integer is expected (e.g. `max_retries`,
#' `ncomponents`). Plain `as.integer()` on such input still returns a
#' value (`NA`) but raises an "NAs introduced by coercion" warning, which
#' pollutes logs during normalization/validation of a plan that is
#' expected to be checked for exactly this kind of malformation. This
#' coerces quietly and substitutes `default` when the result is `NA`.
#'
#' @param value The raw value to coerce.
#' @param default Integer to use when `value` is `NULL`, length-0, or not
#'   coercible to a finite integer.
#' @return A single integer, never `NA`.
sclet_ai_coerce_integer <- function(value, default = 0L) {
    if (is.null(value) || !length(value)) {
        return(as.integer(default))
    }
    coerced <- suppressWarnings(as.integer(value[[1L]]))
    if (is.na(coerced)) {
        return(as.integer(default))
    }
    coerced
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
            max_retries = sclet_ai_coerce_integer(action$max_retries, default = 0L),
            continue_on_error = if (is.null(action$continue_on_error)) FALSE else isTRUE(action$continue_on_error)
        )
    })
}

sclet_ai_normalize_success_criteria <- function(success_criteria) {
    if (is.null(success_criteria)) {
        return(list())
    }
    if (!is.list(success_criteria)) {
        success_criteria <- as.list(success_criteria)
    }
    lapply(seq_along(success_criteria), function(i) {
        criterion <- success_criteria[[i]]
        if (is.character(criterion) && length(criterion) == 1L) {
            criterion <- list(description = criterion)
        }
        if (!is.list(criterion)) {
            return(list(
                id = paste0("criterion_", i),
                description = as.character(criterion)[[1L]],
                source = "action_output",
                step_id = NULL,
                evidence_id = NULL,
                check = "exists",
                field = NULL,
                value = NULL
            ))
        }
        list(
            id = as.character(criterion$id %||% paste0("criterion_", i))[[1L]],
            description = as.character(criterion$description %||% "")[[1L]],
            source = as.character(criterion$source %||% "action_output")[[1L]],
            step_id = if (is.null(criterion$step_id)) NULL else as.character(criterion$step_id)[[1L]],
            evidence_id = if (is.null(criterion$evidence_id)) NULL else as.character(criterion$evidence_id)[[1L]],
            check = as.character(criterion$check %||% "exists")[[1L]],
            field = if (is.null(criterion$field)) NULL else as.character(criterion$field)[[1L]],
            value = criterion$value
        )
    })
}

sclet_ai_normalize_stop_conditions <- function(stop_conditions) {
    if (is.null(stop_conditions)) {
        return(list())
    }
    if (!is.list(stop_conditions)) {
        stop_conditions <- as.list(stop_conditions)
    }
    lapply(seq_along(stop_conditions), function(i) {
        condition <- stop_conditions[[i]]
        if (is.character(condition) && length(condition) == 1L) {
            condition <- list(description = condition)
        }
        if (!is.list(condition)) {
            return(list(
                id = paste0("stop_condition_", i),
                description = as.character(condition)[[1L]],
                source = "action_output",
                step_id = NULL,
                evidence_id = NULL,
                check = "exists",
                field = NULL,
                value = NULL
            ))
        }
        list(
            id = as.character(condition$id %||% paste0("stop_condition_", i))[[1L]],
            description = as.character(condition$description %||% "")[[1L]],
            source = as.character(condition$source %||% "action_output")[[1L]],
            step_id = if (is.null(condition$step_id)) NULL else as.character(condition$step_id)[[1L]],
            evidence_id = if (is.null(condition$evidence_id)) NULL else as.character(condition$evidence_id)[[1L]],
            check = as.character(condition$check %||% "exists")[[1L]],
            field = if (is.null(condition$field)) NULL else as.character(condition$field)[[1L]],
            value = condition$value
        )
    })
}

sclet_ai_extract_dotted_field <- function(obj, field) {
    if (is.null(field) || !nzchar(field)) return(NULL)
    parts <- strsplit(field, ".", fixed = TRUE)[[1L]]
    current <- obj
    for (part in parts) {
        if (!is.list(current) || is.null(current[[part]])) return(NULL)
        current <- current[[part]]
    }
    current
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
#' @param assumptions Optional character vector or list of planning assumptions.
#' @param candidate_routes Optional list of named route descriptors (each with at
#'   minimum `id` and `description`).
#' @param selected_route Optional character naming one element of
#'   `candidate_routes`, or free text when `candidate_routes` is empty.
#' @param success_criteria Optional list of named criterion descriptors.
#' @param risks Optional character vector or list of known risks.
#' @param stop_conditions Optional list of named stop-condition descriptors.
#' @return An object of class `sclet_ai_plan`.
new_sclet_ai_plan <- function(
    task = "analysis_plan",
    actions = list(),
    context_fingerprint = NULL,
    rationale = NULL,
    requires_confirmation = TRUE,
    ai_result = NULL,
    metadata = list(),
    assumptions = list(),
    candidate_routes = list(),
    selected_route = NULL,
    success_criteria = list(),
    risks = list(),
    stop_conditions = list()
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
        metadata = metadata,
        assumptions = as.list(assumptions),
        candidate_routes = if (is.null(candidate_routes)) list() else candidate_routes,
        selected_route = selected_route,
        success_criteria = sclet_ai_normalize_success_criteria(success_criteria),
        risks = as.list(risks),
        stop_conditions = sclet_ai_normalize_stop_conditions(stop_conditions)
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
#'   `warnings`, and a one-time confirmation token when valid. Validation errors
#'   include `success_criterion_unresolvable` (criterion references unknown step
#'   or has empty evidence_id), `success_criterion_incomplete` (non-exists check
#'   missing field or value), `stop_condition_unresolvable` (stop condition
#'   references unknown step or has empty evidence_id), `stop_condition_incomplete`
#'   (non-exists check missing field or value), `data_scale_incompatible`
#'   (`run_pca` step requests more components than the object supports),
#'   `selected_route_unknown` (selected_route not in candidate_routes). The
#'   `human_confirmations` field is a read-only derived summary listing design
#'   confirmations already on record (status `"confirmed"`) and design
#'   confirmations the plan still requires (status `"required_not_confirmed"`).
#'   The `claim_ceiling` field is a read-only, informational pre-check: the
#'   strongest claim the evidence already recorded on the object could carry if
#'   the analysis stopped right now, plus the independent lines of support
#'   behind that cap. It is a ceiling and never a recommendation, it adds no
#'   suggestion or ranking, and it never affects `valid`, `errors`, `warnings`,
#'   `requires_confirmation` or `confirmation_token`. It degrades to
#'   `list(status = "unavailable", reason = ...)` when no object was supplied,
#'   when the object is not a `SingleCellExperiment`, or when the ceiling
#'   cannot be evaluated.
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
    # A plan with no actions is not a plan. ValidateAIPlan() used to accept it
    # (nothing inside can fail) and ExecuteAIPlan() then returned
    # status="completed" with zero results, so a repair attempt that gave up
    # and returned nothing was indistinguishable, in every consumer of status,
    # from an analysis that ran. Measured on a real provider 2026-10-16
    # (goal_discrimination, the signalling case): the repair round returned a
    # zero-action plan after the first attempt put an engine placeholder in
    # run_cci's method, and it validated with validation_error_counts "1,0".
    if (!length(actions)) {
        errors <- c(errors, paste0(
            "plan_empty: the plan has no actions. A plan with zero steps ",
            "performs no analysis and cannot be executed; re-propose the ",
            "plan, or if the goal genuinely needs no computation, answer ",
            "directly instead of planning."
        ))
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
    unusable_param_steps <- character()
    # Step id -> action name, so a reference can be resolved to the descriptor
    # of the step it names (references name step ids, not actions).
    action_by_step_id <- list()
    for (s in normalized) {
        action_by_step_id[[as.character(s$id)]] <- s$action
    }
    if (length(normalized)) {
        known_ids <- character()
        for (i in seq_along(normalized)) {
            step <- normalized[[i]]
            if (is.null(step$action) || !nzchar(step$action)) {
                errors <- c(errors, paste0("action ", step$id, " has no action name"))
            } else if (is.null(registry[[step$action]])) {
                errors <- c(errors, paste0(
                    "action is not registered: ", step$action,
                    "; the `action` field must be an exact registry key, and a",
                    " placeholder such as \"none\" is not a way to drop a step",
                    " -- remove the step from the plan instead"
                ))
            } else {
                descriptor <- registry[[step$action]]
                input_schema <- descriptor$input_schema %||% list()
                parameter_problems <- sclet_ai_action_param_problems(
                    step$params,
                    input_schema
                )
                if (length(parameter_problems)) {
                    # Only a *present* parameter of the wrong type makes the
                    # step's parameters unusable for a prerequisite verdict. A
                    # merely absent parameter must not: the prerequisite is
                    # often what produces the best message about it (e.g.
                    # start_cluster_missing).
                    if (length(sclet_ai_action_param_type_problems(step$params, input_schema))) {
                        unusable_param_steps <- c(unusable_param_steps, step$id)
                    }
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
            # Dependency analysis uses only well-formed references. A malformed
            # one is reported on its own (below); deriving dependencies from a
            # token we cannot parse would only add a second, less precise error
            # for the model to chase.
            references <- sclet_ai_collect_plan_references(step$params)
            well_formed <- Filter(sclet_ai_reference_is_well_formed, references)
            reference_steps <- unique(vapply(
                well_formed,
                function(path) path[[1L]],
                character(1)
            ))
            if (length(reference_steps)) {
                unknown_reference_steps <- setdiff(reference_steps, ids)
                if (length(unknown_reference_steps)) {
                    errors <- c(errors, paste0(
                        "action ", step$id, " references unknown step(s): ",
                        paste(unknown_reference_steps, collapse = ", "),
                        "; a reference may only name a step of this plan"
                    ))
                }
                undeclared <- setdiff(reference_steps, step$depends_on)
                if (length(undeclared)) {
                    errors <- c(errors, paste0(
                        "action ", step$id, " references undeclared dependency step(s): ",
                        paste(undeclared, collapse = ", "),
                        "; add them to this step's depends_on"
                    ))
                }
                # A reference is only useful if execution can resolve it, and
                # execution resolves against the runtime output summary, not
                # against the declared output_schema. Checking the field name
                # here is what keeps validation from accepting a plan that then
                # dies at execution with "Plan output field is not available".
                unresolvable <- character()
                for (path in well_formed) {
                    referenced_id <- path[[1L]]
                    referenced_action <- action_by_step_id[[referenced_id]]
                    descriptor_of <- if (!is.null(referenced_action) && !is.na(referenced_action)) {
                        registry[[referenced_action]]
                    } else {
                        NULL
                    }
                    if (is.null(descriptor_of)) {
                        next
                    }
                    allowed <- sclet_ai_referenceable_output_fields(descriptor_of)
                    if (!path[[3L]] %in% allowed) {
                        unresolvable <- c(unresolvable, paste0(
                            "${", paste(path, collapse = "."), "}",
                            " (step ", referenced_id, "=", referenced_action,
                            " exposes: ", paste(allowed, collapse = ", "), ")"
                        ))
                    }
                }
                if (length(unresolvable)) {
                    unusable_param_steps <- c(unusable_param_steps, step$id)
                    errors <- c(errors, paste0(
                        "value_reference_unresolvable: action ", step$id,
                        " references ", paste(unresolvable, collapse = "; ")
                    ))
                }
            }
            malformed_references <- sclet_ai_malformed_plan_references(step$params)
            if (length(malformed_references)) {
                # Report the syntax mistake itself. Left undetected, the
                # template is passed through as a literal string and the step's
                # prerequisite reports something unrelated to the real cause.
                unusable_param_steps <- c(unusable_param_steps, step$id)
                suggestions <- vapply(
                    malformed_references,
                    sclet_ai_suggest_reference_fix,
                    character(1),
                    USE.NAMES = FALSE
                )
                fixed <- suggestions[!is.na(suggestions)]
                errors <- c(errors, paste0(
                    "value_reference_malformed: action ", step$id, " uses ",
                    paste(malformed_references, collapse = ", "),
                    if (length(fixed)) {
                        paste0("; write ", paste(fixed, collapse = ", "), " instead")
                    } else {
                        ""
                    },
                    " (the syntax is ${step_id.output.field})"
                ))
            }
            known_ids <- c(known_ids, step$id)
        }
    }

    normalized_criteria <- sclet_ai_normalize_success_criteria(plan$success_criteria %||% list())
    if (length(normalized_criteria)) {
        criterion_ids <- vapply(normalized_criteria, function(x) x$id %||% "", character(1))
        if (anyDuplicated(criterion_ids)) {
            errors <- c(errors, "success_criterion ids must be unique")
        }
        for (i in seq_along(normalized_criteria)) {
            crit <- normalized_criteria[[i]]
            cid <- crit$id
            if (identical(crit$source, "action_output")) {
                if (is.null(crit$step_id) || !nzchar(crit$step_id)) {
                    errors <- c(errors, paste0("success_criterion_unresolvable: ", cid, " has empty step_id"))
                } else if (!crit$step_id %in% ids) {
                    errors <- c(errors, paste0("success_criterion_unresolvable: ", cid, " references unknown step ", crit$step_id))
                }
            } else if (identical(crit$source, "evidence")) {
                if (is.null(crit$evidence_id) || !nzchar(crit$evidence_id)) {
                    errors <- c(errors, paste0("success_criterion_unresolvable: ", cid, " has empty evidence_id"))
                }
            }
            if (!identical(crit$check, "exists")) {
                if (is.null(crit$field) || is.null(crit$value)) {
                    errors <- c(errors, paste0("success_criterion_incomplete: ", cid, " requires field and value for check '", crit$check, "'"))
                }
            }
        }
    }
    normalized_stop_conditions <- sclet_ai_normalize_stop_conditions(plan$stop_conditions %||% list())
    if (length(normalized_stop_conditions)) {
        stop_ids <- vapply(normalized_stop_conditions, function(x) x$id %||% "", character(1))
        if (anyDuplicated(stop_ids)) {
            errors <- c(errors, "stop_condition ids must be unique")
        }
        for (i in seq_along(normalized_stop_conditions)) {
            cond <- normalized_stop_conditions[[i]]
            cid <- cond$id
            if (identical(cond$source, "action_output")) {
                if (is.null(cond$step_id) || !nzchar(cond$step_id)) {
                    errors <- c(errors, paste0("stop_condition_unresolvable: ", cid, " has empty step_id"))
                } else if (!cond$step_id %in% ids) {
                    errors <- c(errors, paste0("stop_condition_unresolvable: ", cid, " references unknown step ", cond$step_id))
                }
            } else if (identical(cond$source, "evidence")) {
                if (is.null(cond$evidence_id) || !nzchar(cond$evidence_id)) {
                    errors <- c(errors, paste0("stop_condition_unresolvable: ", cid, " has empty evidence_id"))
                }
            }
            if (!identical(cond$check, "exists")) {
                if (is.null(cond$field) || is.null(cond$value)) {
                    errors <- c(errors, paste0("stop_condition_incomplete: ", cid, " requires field and value for check '", cond$check, "'"))
                }
            }
        }
    }
    candidate_routes <- plan$candidate_routes %||% list()
    selected_route <- plan$selected_route
    if (length(candidate_routes) && !is.null(selected_route)) {
        route_ids <- vapply(candidate_routes, function(r) as.character(r$id %||% ""), character(1))
        if (!as.character(selected_route) %in% route_ids) {
            errors <- c(errors, paste0("selected_route_unknown: '", selected_route, "' is not in candidate_routes"))
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
                if (!is.null(step$action) && identical(step$action, "run_pca")) {
                    requested_ncomp <- step$params$ncomponents
                    if (!is.null(requested_ncomp)) {
                        requested_ncomp <- suppressWarnings(as.integer(requested_ncomp))
                        max_possible <- min(ncol(object), nrow(object))
                        if (!is.na(requested_ncomp) && requested_ncomp > max_possible) {
                            errors <- c(errors, paste0(
                                "data_scale_incompatible: step ", step$id,
                                " requests ncomponents=", requested_ncomp,
                                " but the object only has min(ncol, nrow)=", max_possible,
                                " (ncol=", ncol(object), ", nrow=", nrow(object), ")"
                            ))
                        }
                    }
                }
                if (is.null(descriptor) || !is.function(descriptor$prerequisites)) {
                    next
                }
                if (step$id %in% unusable_param_steps) {
                    # This step has a parameter that is present but unusable --
                    # either of a type the schema rejects, or a malformed value
                    # reference that cannot be resolved. Probing the prerequisite
                    # anyway derives a verdict from unusable input: it may raise
                    # an opaque error ("missing value where TRUE/FALSE needed"),
                    # or -- worse -- return a confident but meaningless answer,
                    # e.g. a numeric comparison silently becoming a lexicographic
                    # one on a character value. Either way the model gets a
                    # second, misleading error to chase instead of the one
                    # actionable error it already has. Skip the probe, keep the
                    # real error, and say so.
                    #
                    # Note this must not extend to an *absent* parameter: there
                    # the prerequisite is usually the best diagnostic available.
                    warnings <- c(warnings, paste0(
                        "prerequisite check skipped for action ", step$id,
                        ": parameters are unusable, so a prerequisite ",
                        "verdict would be derived from invalid input"
                    ))
                } else {
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
    human_confirmations <- list()
    if (!is.null(object) && inherits(object, "SingleCellExperiment")) {
        design_confirms <- tryCatch(
            GetAnalysisLedger(object)$analysis_story$design_confirmations,
            error = function(e) list()
        )
        if (is.list(design_confirms) && length(design_confirms)) {
            for (dc in design_confirms) {
                design_map <- dc$design
                if (!is.list(design_map)) next
                confirmed_at <- dc$created_at
                for (role in names(design_map)) {
                    human_confirmations[[length(human_confirmations) + 1L]] <- list(
                        role = as.character(role),
                        column = as.character(design_map[[role]]),
                        status = "confirmed",
                        confirmed_at = confirmed_at
                    )
                }
            }
        }
    }
    unconfirmed_errors <- errors[grepl("design_semantics_not_confirmed", errors, fixed = TRUE)]
    if (length(unconfirmed_errors)) {
        for (ue in unconfirmed_errors) {
            human_confirmations[[length(human_confirmations) + 1L]] <- list(
                role = NA_character_,
                column = NA_character_,
                status = "required_not_confirmed",
                confirmed_at = NULL
            )
        }
    }
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

    # Read-only pre-check: the strongest claim the evidence already recorded on
    # this object could carry if the analysis stopped right now. Computed only
    # after `valid`, `errors`, `requires_confirmation` and `confirmation_token`
    # are final, so it cannot influence any of them (nor the RNG the token
    # draws from). `ValidateAIPlan()` runs on the hot path of every plan, so
    # every failure mode -- no object, a non-`SingleCellExperiment` object, an
    # unavailable `check_claim_ceiling()`, an error inside the ceiling
    # evaluation -- degrades to a short `unavailable` status rather than
    # propagating. A planning failure caused by an informational field would be
    # a regression.
    claim_ceiling <- tryCatch(
        sclet_ai_plan_claim_ceiling(object),
        error = function(e) {
            reason <- conditionMessage(e)
            if (is.null(reason) || length(reason) != 1L) {
                reason <- NA_character_
            } else {
                reason <- as.character(reason)
            }
            if (is.na(reason) || !nzchar(reason)) {
                reason <- "claim ceiling pre-check could not be evaluated"
            }
            list(status = "unavailable", reason = substr(reason, 1L, 200L))
        }
    )
    validation <- list(
        valid = valid,
        errors = errors,
        warnings = unique(warnings[nzchar(warnings)]),
        plan_id = plan$plan_id %||% NULL,
        # Bind validation (and its confirmation token) to the complete
        # normalized plan, not just the caller-controlled plan_id.
        # ExecuteAIPlan() normalizes plan$actions before it re-fingerprints,
        # so the digest here must be computed over the SAME normalized shape
        # (default-filled max_retries / requires_confirmation / params ...),
        # otherwise every not-yet-normalized caller plan is misrejected as
        # sclet_ai_stale_plan.
        plan_fingerprint = {
            fingerprinted <- plan
            fingerprinted$actions <- normalized
            sclet_ai_plan_fingerprint(fingerprinted)
        },
        context_fingerprint = current_fingerprint %||% plan$context_fingerprint %||% NULL,
        requires_confirmation = requires_confirmation,
        confirmation_token = confirmation_token,
        claim_ceiling = claim_ceiling,
        plan = plan,
        human_confirmations = human_confirmations,
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

# Read-only claim-ceiling pre-check for a plan (internal)
#
# Answers one question and nothing else: if the analysis stopped right now and
# a report were written from what is already recorded on the object, what is
# the strongest claim that evidence could carry? It delegates to
# `check_claim_ceiling()` with no `claims` and no `evidence_ids`, so the answer
# is derived from the object's currently recorded evidence and from nothing
# else -- a genuine pre-check, because a planner learns the ceiling before
# executing a single step.
#
# The result is a cap, never a hint: it reports what the evidence can support,
# never what should be claimed, and it carries no suggestion, recommendation or
# ranking. Only structural facts are returned -- a claim level, counts of lines
# of support, and the reasons behind them -- so `raw_values_included` is always
# `FALSE` and no data value ever reaches this field.
#
# Informational only: this function never mutates the object, never records
# evidence, and never feeds `valid`, `errors` or `warnings`. `ValidateAIPlan()`
# runs on the hot path of every plan, so every failure mode degrades to
# `list(status = "unavailable", reason = ...)` instead of propagating -- a
# planning failure caused by an informational field would be a regression. The
# failure modes are a missing object, an object that is not a
# `SingleCellExperiment`, an unavailable `check_claim_ceiling()`, an error
# raised by the ceiling evaluation, and a result that is not a list or carries
# no status.
#
# @param object The `SingleCellExperiment` the plan was built for, or `NULL`.
# @return A nested list with `status` (the ceiling evaluation's verdict,
#   `"unsupported"` or `"allowed"`, or `"unavailable"` when this pre-check
#   could not run), `ceiling_claim_level` (`NULL` when nothing is supported),
#   `n_independent_groups`, `n_supporting_nodes`, `reasons` and
#   `raw_values_included = FALSE`; or
#   `list(status = "unavailable", reason = ...)` when the pre-check cannot be
#   evaluated.
# @noRd
sclet_ai_plan_claim_ceiling <- function(object) {
    unavailable <- function(reason) {
        if (is.null(reason) || length(reason) != 1L) {
            reason <- NA_character_
        } else {
            reason <- as.character(reason)
        }
        if (is.na(reason) || !nzchar(reason)) {
            reason <- "claim ceiling pre-check could not be evaluated"
        }
        # A short reason only: this field is informational and stays small.
        list(status = "unavailable", reason = substr(reason, 1L, 200L))
    }
    if (is.null(object)) {
        return(unavailable("object was not supplied"))
    }
    if (!inherits(object, "SingleCellExperiment")) {
        return(unavailable("object is not a SingleCellExperiment"))
    }
    if (!exists("check_claim_ceiling", mode = "function")) {
        return(unavailable("check_claim_ceiling is not available"))
    }
    ceiling <- check_claim_ceiling(object = object, claims = NULL, evidence_ids = NULL)
    if (!is.list(ceiling)) {
        return(unavailable("check_claim_ceiling returned no result"))
    }
    # `status` mirrors the ceiling evaluation's own verdict vocabulary --
    # `"unsupported"` when the recorded evidence cannot support any claim at
    # all, `"allowed"` when a ceiling was computed -- so the plan layer speaks
    # the same language as the write gate in `RecordAIResult()`. `"unavailable"`
    # is reserved for this pre-check failing, and is the only value a caller
    # ever needs to distinguish: the others describe the ceiling itself.
    status <- if (is.character(ceiling$ceiling_status) && length(ceiling$ceiling_status) == 1L &&
        !is.na(ceiling$ceiling_status) && nzchar(ceiling$ceiling_status)) {
        ceiling$ceiling_status
    } else if (is.character(ceiling$status) && length(ceiling$status) == 1L &&
        !is.na(ceiling$status) && nzchar(ceiling$status)) {
        ceiling$status
    } else {
        character()
    }
    if (!length(status)) {
        return(unavailable("check_claim_ceiling reported no status"))
    }
    # Counts are reported as they are read. A missing or non-numeric count is
    # `NA_integer_` -- an honest "unknown" -- never a zero, which would read as
    # a real verdict about how much support exists.
    count <- function(value, default = NA_integer_) {
        if (is.null(value) || !is.numeric(value) || length(value) != 1L || is.na(value)) {
            return(default)
        }
        as.integer(value)
    }
    reasons <- ceiling$ceiling_reasons
    if (is.null(reasons)) reasons <- ceiling$reasons
    reasons <- if (is.character(reasons)) {
        reasons
    } else if (is.list(reasons)) {
        unlist(reasons, use.names = FALSE)
    } else {
        character()
    }
    reasons <- as.character(reasons)
    reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
    list(
        status = status,
        ceiling_claim_level = if (is.character(ceiling$ceiling_claim_level) &&
            length(ceiling$ceiling_claim_level) == 1L &&
            !is.na(ceiling$ceiling_claim_level) && nzchar(ceiling$ceiling_claim_level)) {
            ceiling$ceiling_claim_level
        } else {
            NULL
        },
        n_independent_groups = count(ceiling$n_independent_groups),
        n_supporting_nodes = count(ceiling$n_supporting_nodes),
        reasons = reasons,
        raw_values_included = FALSE
    )
}

#' Build the repair prompt block from validation errors
#'
#' Pure helper: character vector of validation errors -> character vector of
#' prompt lines to append to the system prompt. Returns `character(0)` for
#' empty/NULL input so callers can append unconditionally.
#'
#' @param errors Character vector of validation error strings (verbatim).
#' @return Character vector of prompt lines, or `character(0)`.
#' @noRd
sclet_ai_repair_prompt_block <- function(errors) {
    if (is.null(errors) || !length(errors) || !any(nzchar(errors))) {
        return(character(0))
    }
    errors <- as.character(errors)
    errors <- errors[nzchar(errors)]
    c(
        "The previous plan you proposed was rejected by deterministic",
        "validation against the registry and the object; it was NOT executed",
        "and must not be treated as having run. Fix ONLY the problems listed",
        "below; preserve the rest of the plan unchanged.",
        "Validation errors (verbatim):",
        paste0("- ", errors),
        "If a step cannot be fixed using only the registered actions and",
        "the parameter types shown in `context$available_actions`, DELETE",
        "that step from the plan. Do NOT invent a new action name, do NOT",
        "rephrase an action name, and do NOT substitute a similar-sounding",
        "name: the `action` field must be an exact key from",
        "`context$available_actions`.",
        "Important: the values shown in each action's `input_schema` are",
        "TYPE DESCRIPTIONS (e.g. \"integer\", \"character\", \"numeric\"),",
        "NOT literal values to pass. Pass real values of the declared type",
        "(e.g. an actual assay name for an `assay` parameter, an actual",
        "integer for `ncomponents`), not the type label itself.",
        "Pass numbers as bare JSON numbers (`\"ncomponents\": 50`), never as",
        "quoted strings (`\"50\"`), and never as a placeholder or an example.",
        "The errors above state what each rejected parameter actually",
        "received; use that to correct the shape."
    )
}

#' One-shot plan proposal (internal)
#'
#' Extracts the body of [AIPlanAnalysis()] so the repair loop and the
#' single-call path cannot drift apart. Returns a `sclet_ai_plan` with
#' `$validation` attached; never throws on an invalid plan.
#'
#' @inheritParams AIPlanAnalysis
#' @param repair_errors Optional character vector of validation errors from
#'   a previous attempt. When non-empty, a repair block is appended to the
#'   system prompt.
#' @return An `sclet_ai_plan` with `$validation` attached.
#' @noRd
sclet_ai_plan_once <- function(
    object,
    goal,
    model,
    structured_output,
    fallback_on_structure_error,
    registry,
    repair_errors = NULL,
    ...
) {
    context <- sclet_ai_context(object, target = "planning", detail = "full", registry = registry)
    base_prompt <- paste(
        "Propose a conservative analysis plan from the supplied ledger.",
        if (!is.null(goal) && length(goal) == 1L && nzchar(trimws(goal))) {
            paste("The user's analysis goal is:", goal)
        } else {
            "No additional user goal was supplied."
        },
        "Return only actions whose `action` field is an exact, verbatim key",
        "from `context$available_actions` (e.g. \"run_de_test\"), with a",
        "`params` list using only the parameter names listed in that action's",
        "input_schema. Do not invent action names. If no listed action fits",
        "a needed step, say so in your rationale instead of inventing one.",
        "Pass each parameter a real value of the declared type: a bare JSON",
        "number for integer/numeric (ncomponents = 50, NOT \"50\", NOT",
        "\"<integer>\", NOT \"e.g. 50\"), a single JSON string for character, and",
        "a JSON boolean for logical. Never pass a type label, a placeholder, or",
        "an example where the value belongs.",
        "For group-valued parameters (ident.1, ident.2, start_cluster), use",
        "only the `code` values from `context$capabilities$identity_levels$codes`.",
        "When `all_levels_available` is FALSE or no code fits, omit the",
        "parameter if the action supports omission (e.g. run_de_test without",
        "ident.1 runs the all-cluster marker path) rather than guessing.",
        "When a parameter must carry a value that an earlier step produces,",
        "take it from the producing action's `output_schema` (for example",
        "run_pca's output_schema declares reduction = \"PCA\"), or reference it",
        "with the plan's ${step_id.output.field} syntax. Never put a prose",
        "description -- such as \"PCA name from step_4\" -- where a concrete",
        "value or reference is required.",
        "Whenever a step uses another step's output, also list that step's id",
        "in this step's `depends_on` (e.g. depends_on = [\"step_2\"]). A",
        "reference alone is not enough: an undeclared dependency is rejected.",
        "Every step's `depends_on` must name only earlier steps of this plan.",
        "List prerequisites, dependencies, expected outputs, and uncertainty.",
        "Never execute an action or claim that it was executed.",
        "Two OPTIONAL top-level fields are also accepted: `success_criteria`",
        "and `stop_conditions`. Each is a list of criterion objects, each with",
        "fields `id`, `description`, `source`, `step_id` (when source =",
        "\"action_output\"), `evidence_id` (when source = \"evidence\"),",
        "`check`, and `field`/`value` as needed. Both are optional -- omit",
        "them entirely if no criterion corresponds to something the plan will",
        "actually produce. Do not invent criteria unrelated to the proposed",
        "actions.",
        "Source rules: source = \"action_output\" is the safe default and",
        "names a step of THIS plan via its `step_id` (use the step's id from",
        "the proposed_actions list). source = \"evidence\" names an",
        "`evidence_id` the plan is expected to write.",
        paste0(
            "Check must be one of: ",
            paste(sprintf('"%s"', sclet_ai_criterion_checks()), collapse = ", "),
            "."
        ),
        "Any check other than \"exists\" also requires BOTH `field` and",
        "`value`. Prefer a `field` the referenced step's execution output",
        "actually carries, so the criterion can be evaluated: n_features,",
        "n_cells, assays, reductions, active_assay, active_reduction,",
        "fingerprint, class, plus assay, reduction, graph, state where the",
        "referenced step's output_schema declares them. A criterion whose",
        "field is not found is reported as `not_available` at execution",
        "rather than failing the plan, so a wrong field name costs you a",
        "useless criterion, not an error.",
        sep = "\n"
    )
    repair_block <- sclet_ai_repair_prompt_block(repair_errors)
    system_prompt <- if (length(repair_block)) {
        paste(c(base_prompt, "", repair_block), collapse = "\n")
    } else {
        base_prompt
    }
    result <- sclet_ai_call(
        task = "analysis_plan",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        fallback_on_structure_error = fallback_on_structure_error,
        system_prompt = system_prompt,
        ...
    )
    plan <- new_sclet_ai_plan(
        task = result$task,
        actions = result$proposed_actions,
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        rationale = result$answer,
        requires_confirmation = TRUE,
        success_criteria = result$success_criteria %||% list(),
        stop_conditions = result$stop_conditions %||% list(),
        ai_result = result,
        metadata = list(
            provider = result$metadata$provider %||% NULL,
            structured_output = result$metadata$structured_output %||% NULL,
            goal = goal %||% NULL
        )
    )
    plan$validation <- ValidateAIPlan(
        plan,
        object = object,
        registry = registry,
        strict = FALSE
    )
    plan
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
#' @param registry Optional `AIExecutionRegistry`. Its action names,
#'   descriptions, and input schemas are exposed to the model via
#'   `context$available_actions` so proposed actions can reference real,
#'   registered action names. Also used to validate the returned plan
#'   (see `$validation` below). Defaults to
#'   `AIDefaultExecutionRegistry(object, include = "all")`.
#' @param max_repair_rounds Integer in `0:3`. When the first proposal fails
#'   validation, the validation errors are fed back to the model and it is
#'   asked to propose a corrected plan, up to this many times. `0L` preserves
#'   the pre-repair behaviour (exactly one provider call). Each repair round
#'   roughly doubles the provider calls of a plain invocation (one agent
#'   round-trip plus a second structured-output call), so raising this above
#'   the default `1L` has a real cost; the actual number of rounds used is
#'   recorded in `plan$metadata$repair_rounds_used`. Non-integer, `NA`, or
#'   out-of-range values are an error.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return An object of class `sclet_ai_plan`. The model's proposed actions
#'   are not guaranteed to be well-formed: it may omit an action name or
#'   name an action that is not registered. The plan's `$validation` field
#'   (a non-strict `ValidateAIPlan()` result) reports this; check
#'   `plan$validation$valid` before treating `plan$actions` as executable,
#'   or re-run `ValidateAIPlan()` explicitly later if the object has
#'   changed. When repair is enabled, the returned plan is the one with the
#'   fewest validation errors across all attempts (ties go to the earliest
#'   attempt); `plan$metadata$repair_rounds_used` records how many repair
#'   calls were made, `plan$metadata$validation_error_counts` records
#'   the error count per attempt (first attempt first), and
#'   `plan$metadata$first_attempt_errors` keeps the first attempt's error
#'   text -- the only record of what the model got wrong before any repair
#'   feedback, since the returned plan's own errors are usually empty once
#'   repair succeeds. This function itself
#'   never executes proposed actions.
#' @export
AIPlanAnalysis <- function(
    object,
    goal = NULL,
    model = NULL,
    structured_output = TRUE,
    fallback_on_structure_error = TRUE,
    registry = NULL,
    max_repair_rounds = 1L,
    ...
) {
    if (!is.numeric(max_repair_rounds) || length(max_repair_rounds) != 1L ||
        is.na(max_repair_rounds) || !is.finite(max_repair_rounds) ||
        max_repair_rounds != as.integer(max_repair_rounds) ||
        as.integer(max_repair_rounds) < 0L || as.integer(max_repair_rounds) > 3L) {
        sclet_ai_error(
            "sclet_ai_invalid_arg",
            "`max_repair_rounds` must be an integer in 0:3."
        )
    }
    max_repair_rounds <- as.integer(max_repair_rounds)
    if (is.null(registry)) {
        registry <- AIDefaultExecutionRegistry(object, include = "all")
    }

    plan <- sclet_ai_plan_once(
        object = object, goal = goal, model = model,
        structured_output = structured_output,
        fallback_on_structure_error = fallback_on_structure_error,
        registry = registry, repair_errors = NULL, ...
    )
    error_counts <- length(plan$validation$errors %||% character())
    first_attempt_errors <- plan$validation$errors %||% character()
    best <- plan
    best_count <- error_counts
    rounds_used <- 0L

    while (!isTRUE(best$validation$valid) && rounds_used < max_repair_rounds) {
        rounds_used <- rounds_used + 1L
        candidate <- sclet_ai_plan_once(
            object = object, goal = goal, model = model,
            structured_output = structured_output,
            fallback_on_structure_error = fallback_on_structure_error,
            registry = registry,
            repair_errors = best$validation$errors, ...
        )
        candidate_count <- length(candidate$validation$errors %||% character())
        error_counts <- c(error_counts, candidate_count)
        if (candidate_count < best_count) {
            best <- candidate
            best_count <- candidate_count
        } else {
            break
        }
    }

    best$metadata$repair_rounds_used <- rounds_used
    best$metadata$validation_error_counts <- as.integer(error_counts)
    # Keep the first attempt's error text, not just its count. It is the only
    # record of what the model actually got wrong before any repair feedback,
    # and it is what distinguishes "the model never proposed a valid plan"
    # from "the model proposed one and the loop discarded it". Counts alone
    # cannot tell those apart, and the returned plan's errors are usually
    # empty by the time anyone looks.
    best$metadata$first_attempt_errors <- as.character(first_attempt_errors)
    best
}

#' @export
print.sclet_ai_plan <- function(x, ...) {
    cat("sclet AI plan [", x$plan_id, "]\n", sep = "")
    cat("Actions:", length(x$actions), " Confirmation required:", isTRUE(x$requires_confirmation), "\n", sep = "")
    if (!is.null(x$validation)) {
        cat("Validation: ", if (isTRUE(x$validation$valid)) "valid" else "INVALID", "\n", sep = "")
        if (!isTRUE(x$validation$valid) && length(x$validation$errors)) {
            cat("  - ", paste(x$validation$errors, collapse = "\n  - "), "\n", sep = "")
        }
    }
    if (!is.null(x$rationale) && nzchar(x$rationale)) {
        cat(x$rationale, "\n", sep = "")
    }
    invisible(x)
}
