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
#' not available: `$\{step_4.reduction\}`" from some action's prerequisite --
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
#' re-deriving the grammar. `$\{step_4.reduction\}` becomes
#' `$\{step_4.output.reduction\}`.
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


#' Resolve a plan output reference against the outputs produced so far
#'
#' Turns a `${step_id.output.field}` reference into the value the referenced
#' step produced. A reference that cannot be resolved is an error: a plan that
#' points at an output nobody produced would otherwise run against `NULL`.
#'
#' @param value The reference, or any other value (returned unchanged).
#' @param outputs Named list of the outputs produced by the steps run so far.
#' @return The resolved value, or `value` itself when it is not a reference.
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
#
# `$validation` is EXCLUDED. It is the output of ValidateAIPlan(), and that
# function stores the fingerprint on the validation itself -- so hashing it made
# every fingerprint depend on whether the caller had already attached the
# validation to the plan object. Measured 2026-10-21: the same plan hashed two
# different values before and after `plan$validation <- v`, and ExecuteAIPlan()
# then refused with sclet_ai_stale_plan on a plan that had just validated. A
# self-reference in a content hash is a defect in the hash, not a caller error;
# the plan's identity is its task, actions, params and metadata.
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
    plan_for_hash <- plan
    plan_for_hash$validation <- NULL
    serialized <- serialize(canonical(plan_for_hash), NULL, version = 2L)
    if (requireNamespace("digest", quietly = TRUE)) {
        return(paste0("sclet-ai-plan-", digest::digest(serialized, algo = "sha256", serialize = FALSE)))
    }
    path <- tempfile("sclet-ai-plan-")
    on.exit(unlink(path), add = TRUE)
    writeBin(as.raw(serialized), path)
    paste0("sclet-ai-plan-", unname(tools::md5sum(path)))
}

# A digest of what a registry WILL RUN.
#
# The confirmation token bound the plan TEXT, not the code: a plan validated
# against the real registry could be executed against one whose `normalize_data`
# handler was a one-line impostor, and the execution reported
# `status = "completed"`. Measured 2026-10-22. The practical shape of that
# mistake is ordinary -- a caller rebuilds `AIDefaultExecutionRegistry(object,
# include = ...)` between the two calls, or assembles a registry per round -- so
# the token certified a plan whose approved code had been replaced.
#
# The descriptor METADATA is not enough to catch it: the impostor copied the
# description, the input schema and the prerequisites, and was still a different
# function. So the handler's own formals and body are part of the digest.
#
# Known limit, stated rather than implied: `deparse()` captures the handler's
# source, not the environment it closes over. Two closures with identical bodies
# that capture different values would digest the same. That is a much narrower
# door than "any registry at all", and closing it would mean digesting the
# environment, which is both large and not stable across sessions.
sclet_ai_registry_fingerprint <- function(registry) {
    if (is.null(registry) || !length(registry)) return(NA_character_)
    nms <- sort(names(registry))
    parts <- vapply(nms, function(nm) {
        descriptor <- registry[[nm]]
        handler <- descriptor$handler
        handler_formals <- if (is.function(handler)) {
            paste(names(formals(handler)), collapse = ",")
        } else {
            ""
        }
        handler_body <- if (is.function(handler)) {
            paste(deparse(body(handler)), collapse = " ")
        } else {
            paste(deparse(handler), collapse = " ")
        }
        paste(
            nm,
            as.character(descriptor$description %||% "")[[1L]],
            paste(deparse(descriptor$input_schema), collapse = " "),
            paste(deparse(descriptor$requires_confirmation), collapse = " "),
            paste(deparse(descriptor$mutates_object), collapse = " "),
            handler_formals,
            handler_body,
            sep = "\u001f"
        )
    }, character(1))
    canonical <- paste(parts, collapse = "\u001e")
    if (requireNamespace("digest", quietly = TRUE)) {
        return(paste0("sclet-ai-registry-",
                      digest::digest(canonical, algo = "sha256", serialize = FALSE)))
    }
    path <- tempfile("sclet-ai-registry-")
    on.exit(unlink(path), add = TRUE)
    writeLines(canonical, path)
    paste0("sclet-ai-registry-", unname(tools::md5sum(path)))
}

#' Project a plan output reference onto the outputs a plan PROMISES
#'
#' Resolves a `${step_id.output.field}` reference against the outputs the plan
#' declares it will produce, rather than against the outputs that exist yet.
#' Used when validating a plan, where the referenced step has not run.
#'
#' @param value The reference, or any other value (returned unchanged).
#' @param planned_outputs Named list of the outputs the plan declares.
#' @return The projected value, or `value` itself when it is not a reference.
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
# run_de_test registers; `cnv_record_id` stands for the id run_cnv_inference
# registers (measured 2026-10-19: a model wrote
# `cnv_id = "$\{step_1.output.cnv_id\}"`, the projection substituted this
# marker, and the call step was refused with `cnv_record_missing` naming
# "cnv_record_id" -- a literal no planned CNV run could ever have). Each is
# written into the execution summary as an observation, so projecting the
# marker here would hand a prerequisite a literal that no planned-record list
# can contain.
sclet_ai_is_deferred_promise <- function(value) {
    is.character(value) && length(value) == 1L &&
        (identical(value, "hvg") || identical(value, "detest_record_id") ||
            identical(value, "cnv_record_id"))
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
        confirm_design = {
            # Mirror ConfirmAIDesignSemantics()'s record id closely enough for a
            # later step to see that the confirmation will exist. Integration
            # checks that (see run_integration prerequisites), so a plan with
            # confirm_design before run_integration validates -- measured
            # 2026-10-18, where the goal to "integrate the two batches" was
            # unplannable because validation only looked at the object.
            roles <- c("batch", "condition", "subject", "sample")
            given <- Filter(function(role) !is.null(params[[role]]), roles)
            values <- vapply(given, function(role)
                as.character(params[[role]])[[1L]], character(1))
            id <- paste0("design;", paste0(given, "=", values, collapse = ","))
            list(ai_design_confirmation = id)
        },
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
    # `run_integration` registers its reduction under a method-dependent name,
    # so its schema cannot state one statically. Without this, the natural
    # chain `find_neighbors(reduction = "$\{integration.output.reduction\}")`
    # cannot project (the reference stays a literal) and its prerequisite
    # refuses the plan with "reduction is not available:
    # $\{step_6.output.reduction\}" -- measured 2026-10-19 on the Kang goal, where
    # it made the standard integrate -> cluster -> annotate -> DE pipeline
    # unplannable. Same deferred-projection treatment as run_pca above.
    if (identical(descriptor$name, "run_integration")) {
        planned$active_reduction <- sclet_ai_integration_reduction(params$method %||% "fastMNN")
    }
    # Add the records this step will register, so a later step in the same plan
    # can depend on them. Deferred projection again: validation accepts the
    # promise, and if the step really fails to register, its own execution-time
    # contract (or the next step's execution) is still the authority.
    # Records that are not analysis-state records but are still chainable:
    # `call_malignant_cells` consumes a CNV run by id, and the run does not exist
    # when the plan is validated. Project the id so a two-step plan validates;
    # execution re-reads the real record.
    # A design confirmation THIS plan will perform, as a boolean. Using
    # `planned$state_records$ai_design_confirmation` was wrong: it also contains
    # records already on the object, so a STALE confirmation made integration
    # validate (measured immediately: two stale-design tests failed). Only a
    # confirm_design step in the plan sets this.
    if (identical(descriptor$name, "confirm_design")) {
        planned$design_confirmed <- TRUE
    }
    if (identical(descriptor$name, "run_cnv_inference")) {
        planned$cnv_ids <- unique(c(planned$cnv_ids %||% character(),
            as.character(params$name %||% "cnv")))
    }
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
#'
#' @details
#' Exported so that an EXTERNAL agent can drive sclet without the `aisdk`
#' package: the agent reads the shipped skills (`AISkills()` /
#' `ReadAISkill()`), writes the plan itself with this constructor, and hands
#' it to `ValidateAIPlan()` / `ExecuteAIPlan()` -- the same validation,
#' fingerprint binding, ledger recording and evidence registration a
#' provider-written plan gets. The AI provider is only required when the
#' planner itself is a model (`AIPlanAnalysis()`); the engine never is.
#' The plan's `context_fingerprint` must come from
#' `GetAnalysisLedger(object)$fingerprint` for the object the plan will run
#' against, or validation refuses it as stale.
#'
#' @examples
#' \dontrun{
#' sce <- your_SingleCellExperiment()
#' plan <- new_sclet_ai_plan(
#'     task = "pbmc_markers",
#'     context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
#'     actions = list(
#'         list(id = "s1", action = "normalize_data"),
#'         list(id = "s2", action = "find_variable_features",
#'              params = list(nfeatures = 2000), depends_on = "s1"),
#'         list(id = "s3", action = "run_de_test",
#'              params = list(name = "markers"), depends_on = "s2")
#'     )
#' )
#' v <- ValidateAIPlan(plan, object = sce,
#'     registry = AIDefaultExecutionRegistry(sce, include = "all"))
#' executed <- RunAIPlan(sce, plan,
#'     AIDefaultExecutionRegistry(sce, include = "all"),
#'     dry_run = FALSE, confirm = TRUE)
#' }
#' @export
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' registry <- AIDefaultExecutionRegistry(sce, include = "all")
#' plan <- list(
#'     plan_id = "plan_demo",
#'     task = "demo",
#'     actions = list(list(id = "s1", action = "normalize_data"))
#' )
#' validation <- ValidateAIPlan(plan, object = sce, registry = registry)
#' validation$valid
#' # The `action` field must be an exact registry key.
#' bad <- list(
#'     plan_id = "plan_bad",
#'     task = "demo",
#'     actions = list(list(id = "s1", action = "no_such_action"))
#' )
#' ValidateAIPlan(bad, object = sce, registry = registry)$valid
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
        # Measured failure shape (2026-10-19 Kang reliability repeats): the model
        # returned a zero-action plan for a 43k-cell goal and the error said
        # only "re-propose the plan". It has no actions to repair and no
        # statement of what a non-empty plan must contain, so the same shape
        # recurred. The fix states the floor a plan has to clear.
        errors <- c(errors, paste0(
            "plan_empty: the plan has no actions. A plan with zero steps ",
            "performs no analysis and cannot be executed. ",
            "fix: a goal that names N kinds of work needs at least one action ",
            "per step of that work; list those actions, or state in ",
            "`success_criteria` (or the rationale) which part of the goal is ",
            "deliberately not done and why. Returning zero steps is not a way ",
            "to record 'nothing to do'"
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
    # An action that returns a DIFFERENT object (a metacell object, say) takes
    # the ledger with it: every earlier step's record is gone from the object the
    # round is recorded on. Allowing one inside a multi-step plan would report
    # "completed" while silently discarding the earlier steps, so it is only
    # valid as a plan of its own.
    replacement_steps <- Filter(function(step) {
        descriptor <- registry[[as.character(step$action %||% "")]]
        !is.null(descriptor) && isTRUE(descriptor$replaces_object)
    }, normalized)
    if (length(replacement_steps) && length(normalized) > 1L) {
        errors <- c(errors, paste0(
            "object_replacement_not_standalone: ",
            paste(vapply(replacement_steps, function(step) {
                as.character(step$id %||% "")
            }, character(1L)), collapse = ", "),
            " replaces the SingleCellExperiment with a different object and does not",
            " carry the ledger over, so it cannot share a plan with other steps.",
            " fix: run it as the only action of its own plan (or on a copy), and",
            " keep the analyses that produce this object in the earlier plan"
        ))
    }
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
                # Distinguish "no such action" from "exists, but you did not ask
                # for its group". The second case is by far the more common and
                # the old text sent the reader looking for a typo: `run_de_test`
                # is a real action in the `annotation` group, and a plan that
                # named it after a default `include` was told only that it was
                # "not registered" (measured 2026-10-21, running the book).
                elsewhere <- tryCatch(
                    sclet_ai_action_groups(object, step$action),
                    error = function(e) character()
                )
                errors <- c(errors, if (length(elsewhere)) {
                    paste0(
                        "action is not registered: ", step$action,
                        "; it exists, but in the ", paste(shQuote(elsewhere), collapse = "/"),
                        " registry group, which this registry did not request",
                        " -- add ", paste(shQuote(elsewhere), collapse = " or "),
                        " to `include` (or pass a registry that contains it)"
                    )
                } else {
                    paste0(
                        "action is not registered: ", step$action,
                        "; the `action` field must be an exact registry key, and a",
                        " placeholder such as \"none\" is not a way to drop a step",
                        " -- remove the step from the plan instead"
                    )
                })
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
                # Measured failure shape (2026-10-19, 2 of 9 Kang trials): the
                # redundant-DE dedupe removed a step that a later `depends_on`
                # still named, and the error said only which id was unknown.
                # The engine could not know whether the model meant "remove that
                # dependency too" or "put the step back", so it named neither.
                # State both repairs; the choice stays the model's, but it must
                # be a choice rather than a guess.
                errors <- c(errors, paste0(
                    "action ", step$id, " depends on unknown step(s): ",
                    paste(missing_dependencies, collapse = ", "),
                    ". fix: either remove those ids from this step's ",
                    "`depends_on` (and from any ${<id>.output.<field>} ",
                    "reference or success criterion that names them), or add ",
                    "the step back to `actions`. Do not leave a dependency on ",
                    "a step that is not in the plan"
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
                        "; a reference may only name a step of this plan. ",
                        "fix: point the reference at a step that exists in ",
                        "`actions`, or drop the reference and pass a literal ",
                        "value"
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
                # `run_integration` cannot declare its reduction statically (the
                # name depends on the method), but a later step may reference it.
                # Fill the projected output with the name this method registers,
                # so `$\{integration.output.reduction\}` resolves instead of
                # passing the template through to a prerequisite.
                if (identical(step$action, "run_integration")) {
                    planned_outputs[[step$id]]$reduction <- sclet_ai_integration_reduction(
                        planned_params$method %||% "fastMNN")
                }
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
    # The token authorises EXECUTION, and execution's object check is guarded on
    # a non-NULL `context_fingerprint`:
    #
    #     if (!is.null(validation$context_fingerprint) && !identical(...))
    #
    # With `object = NULL` there is nothing to bind the token to, so a token
    # minted here would authorise the plan against ANY object with no
    # prerequisite check having run -- and `object = NULL` is a documented,
    # legal call. Measured 2026-10-22: the token was minted, and the object
    # check was skipped entirely. Reviewing a plan without an object stays
    # legal, because that is a real use ("does this plan even parse?"); it just
    # yields no token, so the plan cannot be executed on the strength of a
    # review that never looked at the data.
    token_bound_to_object <- !is.null(current_fingerprint) ||
        !is.null(plan$context_fingerprint)
    confirmation_token <- if (valid && requires_confirmation && token_bound_to_object) {
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
        # Bind the token to the CODE that will run, not only to the plan text.
        # ExecuteAIPlan() compares this against the registry it was handed; see
        # sclet_ai_registry_fingerprint() for why the handler body is included.
        registry_fingerprint = sclet_ai_registry_fingerprint(registry),
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
        # The loop keeps the attempt with the fewest errors, so a repair that
        # fixes the reported error while breaking an unrelated step scores WORSE
        # than doing nothing: the candidate is discarded and the round is spent.
        # Measured in this benchmark: one run went 1 -> 4 errors, another 1 -> 1
        # (the same error repeated), and both burned a round to end up where
        # they started. Saying this explicitly is cheaper than a smarter
        # selector, because the failure is in what the model was asked to do.
        "A repair that introduces a NEW error is worse than leaving the reported",
        "one alone: the attempt is scored by its total error count, so the",
        "version with fewer errors is kept and yours may be discarded. Change",
        "only the steps implicated by the errors below; leave every other step",
        "exactly as it was, including its id, action name, and parameters.",
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
        "If a `reference_missing`/`labels_missing` error names a placeholder",
        "such as \"<user-named celldex dataset>\", the parameter needs a real",
        "choice, not a rewrite of the placeholder: use a recorded choice from",
        "context$analysis_story$user_decisions if there is one, otherwise",
        "DELETE that step and name the missing choice in the rationale.",
        "The errors above state what each rejected parameter actually",
        "received; use that to correct the shape."
    )
}

# A recorded user choice (a reference dataset, a communication engine, a
# design column) is settled for the planner. Without this line the model could
# read `analysis_story$user_decisions` and still ask again or substitute its
# own preference, which is what happened before the clarification channel
# existed: the engine retried the same blocked parameter until the repair
# budget ran out. Split into a helper so the text is unit-testable -- the
# prompt itself is assembled inline.
sclet_ai_user_decision_prompt_line <- function() {
    paste(
        "Choices recorded in `context$analysis_story$user_decisions` -- a",
        "reference dataset, a communication engine, a design or cluster column --",
        "are settled. Use them in the matching parameter, do not ask for that",
        "choice again, and do not substitute your own preference for it."
    )
}

# B20 (2026-10-21): `analysis_story$data_support$unexplored` reached the
# context but never the reasoning. Measured on the GSE158055 ambient fixture:
# six planning trials with the block injected, and of the four that did NOT
# plan `run_decontx`, not one mentioned ambient anywhere in its rationale --
# the suggestion was in the context and simply did not enter the thought
# process. Root cause, checked by reading this prompt end to end: NOT ONE LINE
# of it tells the model that `data_support` exists or that `unexplored` is
# something to respond to. The block was data delivered without a covering
# instruction, which is why it behaved like a rumour rather than a checklist.
#
# The instruction has to be a RESPONSE requirement, not a "you may find this
# useful". The spec explicitly warns against the weaker form (M13 / B20):
# "consider the suggestions" is advice a model can ignore silently, and silence
# is exactly the failure we measured. Each suggestion must be either planned or
# declined with a reason -- and a reason must NAME the suggestion, so an
# unmentioned suggestion is detectable rather than merely suspicious.
#
# Kept as a helper (like sclet_ai_user_decision_prompt_line()) so a test can
# pin the exact wording without scraping a 200-line paste().
sclet_ai_unexplored_prompt_line <- function() {
    paste(
        "`context$analysis_story$data_support$unexplored` lists analyses this",
        "object can support that the ledger has no result for -- data value not",
        "yet collected. It is part of your brief, not background reading:",
        "EACH entry must be either PLANNED as a step, or explicitly DECLINED in",
        "your rationale with the entry's action named and the reason it does not",
        "serve the goal. A suggestion that appears in neither the plan nor the",
        "rationale is treated as an oversight. Declining is legitimate -- the",
        "goal is the user's, and an unexplored analysis is not automatically",
        "wanted -- but the choice has to be stated, not silent. `data_support`",
        "also reports `blocked` (what this object cannot answer, with the gate",
        "name) and `shapes` (what the data is); read them before proposing."
    )
}

# Does the goal ask for a contrast between experimental groups, and can the
# plan express it?
#
# Kang IFN-beta real-data run (2026-10-18): the goal asked to "compare
# stimulated versus control cells within donors"; the plan ran
# `run_de_test()` with no `ident.1`/`ident.2`, which is FindAllMarkers -- a
# comparison of CLUSTERS, not conditions. It validated (the action is legal),
# executed, and the session reported success with `ev:findallmarkers`. Nothing
# in the pipeline noticed that the numbers do not answer the question, which is
# worse than refusing: a report about cluster markers reads like a report about
# stimulation.
#
# `run_de_test` CAN express a condition contrast since 2026-10-18 (`group` +
# `ident.1`/`ident.2`, with `covariates` for a paired analysis). What it cannot
# do is infer one: with no idents it runs FindAllMarkers, a comparison of
# CLUSTERS. This check therefore BLOCKS execution when all of the following
# hold: the goal names both a comparison and an experimental group word; the
# object has a non-identity column that looks like such a group; and no DE step
# names both `ident.1` and `ident.2`. It returns a reason list so the caller can
# report what was missing instead of guessing.
sclet_ai_goal_contrast_check <- function(goal, plan, object) {
    empty <- list(blocked = FALSE, reason = NULL, group_columns = character())
    if (is.null(goal) || !length(goal) || !nzchar(trimws(as.character(goal)[[1L]]))) {
        return(empty)
    }
    goal <- as.character(goal)[[1L]]
    compares <- grepl("compar|versus|vs\\.?|between|differ", goal, ignore.case = TRUE)
    group_word <- grepl("stimulat|control|condition|treat|disease|healthy|case",
        goal, ignore.case = TRUE)
    if (!(compares && group_word)) return(empty)

    cd <- SummarizedExperiment::colData(object)
    candidates <- grep("^(condition|stim|stimulation|treatment|treat|group|disease|status|healthy|case)$",
        colnames(cd), ignore.case = TRUE, value = TRUE)
    if (!length(candidates)) return(empty)

    de_steps <- Filter(function(action) {
        identical(as.character(action$action %||% "")[[1L]], "run_de_test")
    }, plan$actions %||% list())
    expressible <- any(vapply(de_steps, function(action) {
        params <- action$params %||% list()
        !is.null(params$ident.1) && !is.null(params$ident.2) &&
            length(params$ident.1) > 0L && length(params$ident.2) > 0L
    }, logical(1L)))
    if (expressible) return(empty)

    list(
        blocked = TRUE,
        reason = paste0(
            "goal_contrast_unavailable: the goal asks for a comparison between ",
            paste(candidates, collapse = "/"),
            " groups, but no step compares those groups. `run_de_test` with no ",
            "ident.1/ident.2 is FindAllMarkers, a comparison of CLUSTERS. ",
            "A condition-aware contrast is expressible instead: set ",
            "group = <the column>, ident.1/ident.2 = <its levels>, and ",
            "covariates = <the subject column> for a paired analysis. Reporting ",
            "cluster markers as the requested contrast would answer a different ",
            "question."
        ),
        group_columns = candidates
    )
}

# Drop a redundant all-cluster marker test when the plan already contains the
# targeted contrast the goal asked for.
#
# Full-Kang plans (2026-10-18) contained `run_de_test()` with no idents (i.e.
# FindAllMarkers) AND `run_de_test(group = condition, ...)`. The second answers
# the question; the first costs a marker pass over 43k cells and can make a
# report look like it produced two independent analyses. It is removed only when
# all of the following hold, so nothing a plan depends on disappears:
#
#   * a targeted DE step exists (group, or both idents);
#   * the goal does not ask for markers;
#   * no other step references the all-cluster step's output
#     (`$\{step.output...\}` or its detest_id);
#   * it is not the only DE step.
#
# Removals are recorded in `plan$metadata$removed_redundant_steps` so the report
# can say what was dropped instead of silently shrinking the plan.
sclet_ai_deduplicate_de_steps <- function(plan, goal = NULL) {
    actions <- plan$actions %||% list()
    if (length(actions) < 2L) return(plan)
    de <- Filter(function(action) {
        identical(as.character(action$action %||% "")[[1L]], "run_de_test")
    }, actions)
    if (length(de) < 2L) return(plan)
    params_of <- function(action) action$params %||% list()
    targeted <- Filter(function(action) {
        params <- params_of(action)
        !is.null(params$group) ||
            (!is.null(params$ident.1) && !is.null(params$ident.2))
    }, de)
    if (!length(targeted)) return(plan)
    goal_text <- if (is.null(goal)) "" else as.character(goal)[[1L]]
    if (nzchar(goal_text) && grepl("marker", goal_text, ignore.case = TRUE)) {
        return(plan)
    }
    all_cluster <- Filter(function(action) {
        params <- params_of(action)
        is.null(params$group) && is.null(params$ident.1) && is.null(params$ident.2)
    }, de)
    if (!length(all_cluster)) return(plan)

    ## Success criteria and stop conditions are part of the plan too: removing a
    ## step they name produced `success_criterion_unresolvable: c9 references
    ## unknown step step_10` on Kang (2026-10-18). A reference there is often
    ## embedded in prose ("the marker table from $\{m1.output.detest_id\} is
    ## non-empty"), and the token collector only recognises a string that IS the
    ## token, so scan every string for the `${<step>.` prefix instead.
    ##
    ## A reference is also spelled as a bare id in a structured field:
    ## `depends_on = "step_10"` on a later action, or a success criterion's
    ## `step_id`. Scanning only params plus token strings missed both, and on
    ## the 2026-10-19 Kang reliability repeats the dedupe removed a step_10
    ## that a criterion and a `depends_on` still named -- 2 of 3 trials died
    ## with "depends on unknown step(s)" / "references unknown step step_10".
    ## A bare id is compared for equality; a token is searched for as a prefix.
    references_id <- function(value, id) {
        if (is.character(value)) {
            return(any(grepl(paste0("\\$\\{", id, "\\."), value)) ||
                any(value == id))
        }
        if (is.list(value)) {
            return(any(vapply(value, references_id, logical(1L), id = id)))
        }
        FALSE
    }
    ## An action must not count as referencing itself: drop its own `id` and
    ## keep params, depends_on and every other field it carries.
    action_candidates <- lapply(actions, function(action)
        action[setdiff(names(action), "id")])
    candidates <- c(action_candidates,
        plan$success_criteria %||% list(), plan$stop_conditions %||% list())
    all_ids <- vapply(actions, function(action)
        as.character(action$id %||% "")[[1L]], character(1))
    referenced_ids <- all_ids[vapply(all_ids, function(id)
        nzchar(id) && any(vapply(candidates, references_id, logical(1L), id = id)),
        logical(1L))]
    removable <- Filter(function(action) {
        id <- as.character(action$id %||% "")[[1L]]
        nzchar(id) && !(id %in% referenced_ids)
    }, all_cluster)
    if (!length(removable)) return(plan)
    remove_ids <- vapply(removable, function(action)
        as.character(action$id)[[1L]], character(1))
    plan$actions <- Filter(function(action)
        !(as.character(action$id %||% "")[[1L]] %in% remove_ids), actions)
    metadata <- plan$metadata %||% list()
    metadata$removed_redundant_steps <- c(
        metadata$removed_redundant_steps %||% character(),
        paste0(remove_ids, ": run_de_test (all-cluster markers) superseded by ",
            paste(vapply(targeted, function(action)
                as.character(action$id)[[1L]], character(1)), collapse = ",")))
    plan$metadata <- metadata
    plan
}

# Which registered action a goal phrase asks for.
#
# A plan can be valid and still not answer the question that was asked: on the
# first real-data pilot (PBMC 10k, 2026-10-18) the goal said "annotate the cell
# types" and one run's plan simply omitted run_annotation, validated, and
# completed. Nothing in the pipeline noticed, because validation only asks
# whether each proposed step is executable. This map is the deterministic part
# of "did the plan cover what the goal named"; it is deliberately keyword-based
# and used for feedback and reporting, never as a validation gate (a loose goal
# phrasing must not block an executable plan).
sclet_ai_goal_action_map <- function() {
    list(
        run_annotation = "annotat|cell[ -]?types?|celltype|label transfer",
        # "differential" alone also matches "differential abundance" (a milo
        # question); keep DE to the expression-test phrasings.
        run_de_test = "marker genes?|markers|differential expression|differentially expressed|DE genes?",
        run_enrichment = "enrich|ontolog|pathway|gene ontology|KEGG|GO term",
        # Bare "correct" was too broad: "correct the labels" is an annotation
        # change, not a batch correction.
        run_integration = "batch effect|batch correct|integrat|harmony|scvi",
        run_velocity = "velocity|splicing|scvelo",
        run_fate_analysis = "fate|cellrank|terminal state",
        run_trajectory = "trajectory|pseudotime|slingshot",
        run_cci = "communication|ligand|receptor|cell[ -]?cell interaction",
        run_doublet_detection = "doublet",
        run_rare_cell_detection = "rare cell|rare population",
        run_milo = "abundance|composition|differential abundance|milo",
        # The two malignant actions are separate because the CNV run is the
        # expensive, reference-dependent half; a goal about copy number alone
        # should not be forced through the call, and a call without it is refused.
        run_cnv_inference = "copy[ -]?number|CNV|inferCNV|genomic instability",
        call_malignant_cells = "malignant|malignancy|tumou?r cells?",
        # Only PLANNER-VISIBLE actions belong here. This map feeds the repair
        # message ("the goal named X, the plan does not do it"), so an entry for
        # an action the planner cannot see would demand a step it cannot
        # propose and spend a billed repair round on it (measured 2026-10-23:
        # "gene set enrichment" demanded the driver-only run_geneset_scoring).
        # The five driver-only mainlines are described in the skill instead.
        run_perturbation_priority = "perturbation priorit|augur|prioriti[sz]e cell types?",
        # Imputation is planner-visible and has no other pattern competing for
        # its goal. The malignant classifier gets NO entry: `call_malignant_cells`
        # already matches "malignant", and a second requirement would make the
        # coverage judge demand both actions for one goal.
        run_imputation = "imput|ALRA"
    )
}

# What the goal named that the plan does not cover.
#
# Returns list(goal, required, missing). `required` is the set of actions the
# goal phrasing asks for; `missing` is required minus what the plan proposes.
# An empty/NULL goal reports nothing rather than guessing.
sclet_ai_goal_coverage <- function(plan, goal = NULL) {
    goal <- goal %||% plan$metadata$goal
    if (is.null(goal) || !length(goal) || !nzchar(trimws(as.character(goal)[[1L]]))) {
        return(list(goal = NULL, required = character(), missing = character()))
    }
    goal <- as.character(goal)[[1L]]
    present <- vapply(plan$actions %||% list(), function(action) {
        as.character(action$action %||% "")[[1L]]
    }, character(1L))
    map <- sclet_ai_goal_action_map()
    required <- names(map)[vapply(map, function(pattern) {
        grepl(pattern, goal, ignore.case = TRUE)
    }, logical(1L))]
    list(
        goal = goal,
        required = required,
        missing = setdiff(required, present)
    )
}

# The instruction sent as a repair round when the goal names an analysis the
# plan does not contain. Phrased as a validation-style defect because that is
# the channel the repair loop already understands.
sclet_ai_goal_coverage_repair_message <- function(coverage) {
    paste0(
        "goal_coverage_missing: the user's goal names ",
        paste(coverage$missing, collapse = ", "),
        ", but the plan has no step using ",
        if (length(coverage$missing) == 1L) "it" else "them",
        ". Add the missing analysis step(s) WITH every required parameter and a ",
        "concrete value: annotation needs ref + labels, integration needs batch, ",
        "and so on. A step added by name alone cannot pass validation and will ",
        "be discarded, leaving the same gap. If the required value is a USER ",
        "CHOICE (for example which celldex reference to use), take it from ",
        "context$analysis_story$user_decisions if it is recorded there; if it is ",
        "not recorded, do NOT invent a reference -- leave that step out and name ",
        "the missing choice in your rationale so the user can answer it. If the ",
        "goal cannot be met on this object, say so instead of silently returning ",
        "a plan that does not answer the question."
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
    # Domain guidelines, injected from inst/skills. The 2026-10-16 sweep
    # (sclet-benchmark-ai) measured the gap this closes: the engine validated
    # and refused correctly, but the model ordered run_rare_cell_detection
    # before its evidence producers 2/2 times and passed an OrgDb where the
    # ontology database belongs -- knowledge that lives in the single-cell
    # guidelines, not in the registry. Selection is keyword-based on the
    # goal; everything matches a NULL goal. The injection costs at most
    # ~17k characters, small against the planning context itself (8,760
    # characters for 32738 x 2700, measured); the option is the escape
    # hatch for cost-sensitive callers.
    skill_block <- local({
        if (!isTRUE(getOption("sclet.ai.use_skills", TRUE))) return("")
        topics <- if (!is.null(goal) && length(goal) == 1L && nzchar(trimws(goal))) {
            sclet_ai_goal_skill_topics(goal)
        } else {
            NULL
        }
        block <- sclet_ai_skill_context(topics)
        if (!nzchar(block)) return("")
        paste0(
            "\n\n## Single-cell analysis guidelines (shipped skills)\n\n",
            "The following guidelines describe how single-cell analyses are ",
            "ordered and what each registered action requires. They are ",
            "guidelines, not overrides: the engine's validation remains the ",
            "authority, and these tell you how to stay on the right side of ",
            "it on the first attempt.\n\n",
            block, "\n"
        )
    })
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
        sclet_ai_user_decision_prompt_line(),
        sclet_ai_unexplored_prompt_line(),
        "Pass each parameter a real value of the declared type: a bare JSON",
        "number for integer/numeric (ncomponents = 50, NOT \"50\", NOT",
        "\"<integer>\", NOT \"e.g. 50\"), a single JSON string for character, and",
        "a JSON boolean for logical. Never pass a type label, a placeholder, or",
        "an example where the value belongs.",
        "Which level set a group-valued parameter takes depends on whether the",
        "step names a `group` colData column, and getting this wrong is a",
        "measured failure (2026-10-21, the only failure in 12 same-tree",
        "trials): a contrast written as levels of a scheme the object does not",
        "have, refused with contrast_invalid.",
        "(a) WITHOUT a `group` parameter (cluster markers, start_cluster), the",
        "values are clustering codes: use only the `code` values from",
        "`context$capabilities$identity_levels$codes`. When",
        "`all_levels_available` is FALSE or no code fits, omit the parameter if",
        "the action supports omission (run_de_test without ident.1 runs the",
        "all-cluster marker path) rather than guessing.",
        "(b) WITH `group = <column>`, ident.1 and ident.2 are two LEVELS OF THAT",
        "colData COLUMN, not clustering codes. Read them from",
        "`context$capabilities$contrast_levels$columns`, find the entry whose",
        "`column` equals your group, and copy two of its `levels` verbatim. If",
        "the entry is absent, or `levels_available` is FALSE, do not guess a",
        "level name: state the missing choice in your rationale (it belongs to",
        "the user) or omit the contrast.",
        "Cite evidence only by the ids that appear in the evidence ledger",
        "(they start with \"ev:\"). Never cite an AI result id",
        "(\"sclet-ai-...\") as evidence: it is a record of a model answer, not",
        "a measured result, and the engine will refuse to record it.",
        "Before run_integration (or any design-aware step), the batch/condition",
        "column must be confirmed by a `confirm_design` action in the same plan",
        "(e.g. {action: \"confirm_design\", params: {batch: \"batch\", condition:",
        "\"condition\"}}). Without it the integration step is rejected with",
        "design_semantics_not_confirmed, and no repair round can clear it.",
        "A value enclosed in angle brackets (e.g. \"<user-named celldex",
        "dataset>\") is a placeholder, never a value: if you need a choice that",
        "belongs to the user (such as which annotation reference), either use a",
        "recorded choice from context$analysis_story$user_decisions, or leave",
        "that step out of the plan and state the missing choice in your",
        "rationale so the user can answer it. A placeholder in a parameter",
        "slot invalidates the whole plan and no repair round can recover it.",
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
    # Repair block goes AFTER the skill block: the repair attempt's most
    # urgent information is what it got wrong, and the skills are background.
    system_prompt <- if (length(repair_block)) {
        paste(c(base_prompt, skill_block, "", repair_block), collapse = "\n")
    } else {
        paste0(base_prompt, skill_block)
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
    # Snapshot the unexplored suggestions the context carried, so the round-
    # progression judge can later check each one was planned or declined by
    # name (B20). Reading them off the context here rather than off the ledger
    # in the judge matters: the judge runs AFTER execution, and by then the
    # ledger already contains this round's evidence -- so "was it unexplored
    # BEFORE this plan" is only answerable from what the planner saw.
    unexplored_snapshot <- tryCatch({
        story <- context$analysis_story %||% list()
        ds <- story$data_support %||% list()
        entries <- ds$unexplored %||% list()
        lapply(entries, function(x) list(
            action = as.character(x$action %||% ""),
            what = as.character(x$what %||% "")
        ))
    }, error = function(e) list())

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
            goal = goal %||% NULL,
            unexplored_suggestions = unexplored_snapshot
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' \dontrun{
#' # needs a provider API key
#' plan <- AIPlanAnalysis(sce, goal = "cluster the cells and find markers")
#' plan$validation$valid
#' }
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

    # Goal coverage: a separate, one-round feedback path. It cannot ride the
    # repair loop above, because that loop only accepts an attempt with FEWER
    # validation errors -- a valid plan missing a step has zero errors already,
    # so a complete candidate would tie and be discarded.
    coverage <- sclet_ai_goal_coverage(best)
    goal_rounds_used <- 0L
    if (length(coverage$missing) && max_repair_rounds >= 1L) {
        candidate <- sclet_ai_plan_once(
            object = object, goal = goal, model = model,
            structured_output = structured_output,
            fallback_on_structure_error = fallback_on_structure_error,
            registry = registry,
            repair_errors = sclet_ai_goal_coverage_repair_message(coverage), ...
        )
        candidate_coverage <- sclet_ai_goal_coverage(candidate)
        # Only a VALID candidate that covers more of the goal replaces the plan;
        # otherwise the executable plan the user already has is kept.
        if (isTRUE(candidate$validation$valid) &&
                length(candidate_coverage$missing) < length(coverage$missing)) {
            best <- candidate
            coverage <- candidate_coverage
            goal_rounds_used <- 1L
        }
    }
    best <- sclet_ai_deduplicate_de_steps(best, goal)
    best$metadata$goal_coverage <- coverage
    best$metadata$goal_repair_rounds_used <- goal_rounds_used
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

# The reduction name `run_integration` registers for a method.
#
# Measured on this machine 2026-10-14 and re-confirmed 2026-10-19 (Harmony on
# the Kang fixture: reducedDimNames() -> c("PCA","HARMONY"), DefaultReduction()
# -> "HARMONY"). fastMNN keeps "corrected" whatever `name` says. This is the one
# place that mapping lives; the descriptor's `allowed_reductions` lists the
# possible names for validation, and this names the one a plan will actually get.
sclet_ai_integration_reduction <- function(method) {
    switch(as.character(method)[[1L]],
        fastMNN = "corrected",
        Harmony = "HARMONY",
        scVI = "scVI",
        "corrected"
    )
}
