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

# Parse an R-style integer range written as one string: "1:10", "c(1:10)",
# "1 : 10", "-3:3".
#
# Why this exists: `dims` (find_neighbors, run_umap, run_velocity) is declared
# `integer_vector`, and a real model writes it as R source -- `dims="1:30"` --
# in a large fraction of runs. The variable is untrusted plan output, so the
# contract has to say which spellings are accepted; the alternative measured in
# the benchmark was one wasted repair round per run (records 25-26 of
# ai_plan_repair, records 1 and 4 of run_annotation_bench, fate_chain_bench).
#
# Only a range is accepted, and only for vector-typed parameters. A scalar
# parameter still refuses "50" (a string where a number belongs), a gene
# selector still refuses "hvg", and the upper bound keeps a typo such as
# "1:10000000" from expanding into an allocation.
sclet_ai_range_vector <- function(value) {
    if (!is.character(value) || length(value) != 1L || is.na(value)) return(NULL)
    text <- gsub("[[:space:]]+", "", value)
    match <- regmatches(text, regexec("^c?\\(?(-?[0-9]+):(-?[0-9]+)\\)?$", text))[[1L]]
    if (!length(match)) return(NULL)
    from <- suppressWarnings(as.integer(match[[2L]]))
    to <- suppressWarnings(as.integer(match[[3L]]))
    if (is.na(from) || is.na(to) || abs(to - from) > 100000L) return(NULL)
    seq.int(from, to)
}

sclet_ai_vector_param <- function(value, type = c("character", "integer")) {
    type <- match.arg(type)
    if (is.list(value)) value <- unlist(value, use.names = FALSE)
    if (type == "character") return(as.character(value))
    ranged <- sclet_ai_range_vector(value)
    if (!is.null(ranged)) return(as.integer(ranged))
    # Untrusted model output can put a non-numeric value (or a type label such
    # as "integer") where an integer vector is expected. A bare as.integer()
    # still returns NA there but emits "NAs introduced by coercion", which
    # duplicates an error the schema check already reports properly
    # ("parameter <x> has invalid type"). Coerce quietly; the invalid value is
    # still rejected downstream by that schema check.
    suppressWarnings(as.integer(value))
}

# Is this selector value a promise that can only be resolved once the
# producing step has run, rather than a value that can be checked now?
#
# `variable_features` is declared in `find_variable_features`' output_schema as
# the marker "hvg": the static form must not carry 2000 gene names into every
# planning context, so the real set only exists after the step executes. A
# downstream step therefore writes
# `features = "${step_2.output.variable_features}"`, and
# `sclet_ai_project_plan_value()` deliberately leaves that reference
# unprojected -- projecting it would substitute the meaningless marker and hand
# validation a literal it could not act on.
#
# So exactly one form reaches here: the well-formed reference. `prerequisites`
# says yes on the strength of the shape, and
# `sclet_ai_resolve_deferred_selector()` turns it into the real gene set in the
# handler, which is what keeps validation and execution in agreement. A literal
# `"hvg"` typed by a model is not a reference and is rejected here, exactly as
# the gene selector rejects it.
sclet_ai_is_deferred_gene_selector <- function(value) {
    !is.null(sclet_ai_deferred_selector_field(value))
}

# The `<field>` of a deferred gene-selector reference, or NULL when the value is
# an ordinary selector. Only the exact `${step.output.field}` shape counts: a
# bare `${step}` or `${step.output}` is a malformed reference that must keep its
# own error message, and a gene vector or a prose description is not a promise.
sclet_ai_deferred_selector_field <- function(value) {
    if (is.null(value) || !is.character(value) || length(value) != 1L) {
        return(NULL)
    }
    if (!grepl("^\\$\\{[^}]+\\}$", value)) {
        return(NULL)
    }
    path <- strsplit(sub("^\\$\\{", "", sub("\\}$", "", value)), "\\.", fixed = FALSE)[[1L]]
    if (length(path) != 3L || !identical(path[[2L]], "output")) {
        return(NULL)
    }
    if (!identical(path[[3L]], "variable_features")) {
        return(NULL)
    }
    path[[3L]]
}

# Resolve a deferred selector to real row names, using the object the step is
# about to run against. This is what makes the validation-time acceptance
# honest: `prerequisites` says yes because the reference is well-formed, and
# here the marker becomes the actual gene set the model was pointing at.
#
# A value that is not a deferred reference is returned untouched, so ordinary
# gene vectors, range notation and row names keep flowing to the selector
# unchanged.
sclet_ai_resolve_deferred_selector <- function(object, value, arg = "features") {
    field <- sclet_ai_deferred_selector_field(value)
    if (is.null(field) || !identical(field, "variable_features")) {
        return(value)
    }
    features <- tryCatch(VariableFeatures(object), error = function(e) e)
    if (inherits(features, "error") || !length(features)) {
        stop(
            "Deferred selector for `", arg, "` could not be resolved: ",
            "the object carries no variable-feature set. Run ",
            "find_variable_features (or FindVariableFeatures) before this step.",
            call. = FALSE
        )
    }
    as.character(features)
}

#' Prerequisite-friendly wrapper around [sclet_ai_gene_selector()]
#'
#' Returns `TRUE` when the selector resolves (or is absent) and the problem
#' message when it does not, so an action's `prerequisites` can validate a
#' selector the schema cannot: `features`/`subset_row` are declared
#' `character_vector`, so a prose description such as "top HVGs from step_2"
#' passes the schema check. Caught here, the model is told at validation time
#' and the repair loop can fix it; caught only by the handler, the plan
#' validates and then dies at execution with no chance to repair. Both real
#' failures in the last benchmark batch were exactly this.
#'
#' @param object A `SingleCellExperiment`.
#' @param value The raw selector value.
#' @param arg Parameter name, used in the message.
#' @return `TRUE`, or a character problem message.
sclet_ai_selector_problem <- function(object, value, arg) {
    if (is.null(value)) {
        return(TRUE)
    }
    resolved <- tryCatch(
        sclet_ai_gene_selector(object, value, arg = arg),
        error = function(e) e
    )
    if (inherits(resolved, "error")) {
        return(conditionMessage(resolved))
    }
    if (!length(resolved)) {
        return(paste0(arg, " selected no features"))
    }
    TRUE
}

#' Resolve a model-supplied gene selector against the object
#'
#' Gene selectors (`features` for `scale_data`/`run_integration`,
#' `subset_row` for `run_pca`) are declared `character_vector`, so a string
#' passes validation. Model output sometimes expresses the selection as R range
#' notation in a string -- `"1:20"`, `"1:20,30:40"` -- which is unambiguous (it
#' means positions) but not a gene name. Handed straight to the action it
#' subsets by a name that does not exist and base R raises an opaque
#' "subscript out of bounds", reproduced on the real pipeline:
#'
#'   ScaleData(object, features = "1:20", assay = "logcounts")
#'   #> Error: subscript out of bounds
#'
#' That error names neither the parameter nor the problem, and it happens after
#' validation has already passed, so the repair loop never sees it. Resolving
#' the notation explicitly fixes the case that has a single sensible reading,
#' and turns every other unresolvable selector into an error that names the
#' parameter, the offending entries, and the accepted forms.
#'
#' @param object A `SingleCellExperiment`.
#' @param value The raw parameter value (`NULL` passes through).
#' @param arg Parameter name, used in error messages.
#' @return `NULL`, or a character vector of row names of `object`.
sclet_ai_gene_selector <- function(object, value, arg = "features") {
    if (is.null(value)) {
        return(NULL)
    }
    if (is.list(value)) value <- unlist(value, use.names = FALSE)
    available <- rownames(object)
    if (is.null(available) || !length(available)) {
        stop("cannot select ", arg, ": the object has no row names", call. = FALSE)
    }
    n <- length(available)

    unsupported <- function(bad) {
        stop(
            arg, " is not usable: ", paste(utils::head(bad, 5L), collapse = ", "),
            " (", length(bad), " entr", if (length(bad) == 1L) "y" else "ies", ") ",
            if (length(bad) == 1L) "is" else "are",
            " neither a row name of the object nor a valid position. Pass row ",
            "names, or positions within 1..", n, ", or range notation such as ",
            "\"1:20\".",
            call. = FALSE
        )
    }

    if (is.numeric(value)) {
        idx <- suppressWarnings(as.integer(value))
        if (anyNA(idx) || any(idx < 1L | idx > n)) {
            unsupported(value)
        }
        return(available[idx])
    }

    value <- as.character(value)
    missing_rows <- setdiff(value, available)
    if (!length(missing_rows)) {
        return(value)
    }

    # Range notation, and only range notation: digits, colons, commas.
    if (length(value) == 1L &&
        grepl("^[0-9]+\\s*:\\s*[0-9]+(\\s*,\\s*[0-9]+\\s*:\\s*[0-9]+)*$", value)) {
        parts <- strsplit(value, ",", fixed = TRUE)[[1L]]
        idx <- unlist(lapply(parts, function(part) {
            ends <- suppressWarnings(as.integer(strsplit(trimws(part), ":", fixed = TRUE)[[1L]]))
            if (length(ends) != 2L || anyNA(ends)) {
                return(integer())
            }
            seq.int(ends[[1L]], ends[[2L]])
        }))
        if (length(idx) && !anyNA(idx) && all(idx >= 1L & idx <= n)) {
            return(available[idx])
        }
    }

    unsupported(missing_rows)
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
                problems <- c(problems, paste0(
                    "parameter ", name, " must be one of: ", paste(spec, collapse = ", "),
                    " (got ", sclet_ai_describe_param_value(value), ")"
                ))
            }
            next
        }
        type <- as.character(spec)[[1L]]
        valid <- sclet_ai_param_type_valid(value, type)
        if (!isTRUE(valid)) {
            problems <- c(problems, paste0(
                "parameter ", name, " has invalid type: expected ", type,
                ", got ", sclet_ai_describe_param_value(value)
            ))
        }
    }
    unique(problems)
}

#' Describe a parameter value for an error message, without echoing it whole
#'
#' "parameter dims has invalid type" states the rule but not the violation, so
#' neither the model reading the repair prompt nor a human reading a benchmark
#' log can tell whether the value was the number as text (`"30"`), a
#' placeholder (`"<integer>"`), or something else entirely. Naming the received
#' shape makes the error actionable for the model and diagnosable afterwards.
#' Long values are truncated: the point is the shape, not the payload.
#'
#' Is a scalar a placeholder rather than a value?
#'
#' A plan that writes `<user-named celldex dataset>` into `ref` has not supplied
#' a reference; it has recorded an open question in the parameter slot. The
#' distinction matters because the two failures need different recovery: an
#' invalid dataset name is a typo the model can fix in a repair round, while a
#' placeholder means the choice belongs to the user and the clarification
#' channel must ask for it. Measured on the PBMC 10k real-data run (2026-10-18):
#' one plan put placeholders in `ref` and `labels`, validation reported one
#' `reference_invalid`, no question was produced, and the whole plan was
#' discarded.
#'
#' @param value A parameter value.
#' @return TRUE for a single string that is enclosed in angle brackets or names
#'   itself as user-named/user-supplied.
#' @noRd
sclet_ai_value_is_placeholder <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
        (grepl("^[[:space:]]*<[^>]*>[[:space:]]*$", value) ||
            grepl("user[- ]?(named|supplied)", value, ignore.case = TRUE))
}

#' @param value Any parameter value.
#' @param max_chars Maximum characters of a scalar value to show.
#' @return A short single-line description, e.g.
#'   `the text "30"` or `a numeric vector of length 3`.
sclet_ai_describe_param_value <- function(value, max_chars = 60L) {
    preview <- function(x) {
        x <- as.character(x)[[1L]]
        if (is.na(x)) {
            return("NA (a missing value)")
        }
        if (nchar(x) > max_chars) {
            x <- paste0(substr(x, 1L, max_chars), "...")
        }
        paste0("\"", x, "\"")
    }
    if (is.null(value)) {
        return("NULL")
    }
    if (is.list(value)) {
        nms <- names(value)
        return(paste0(
            "a list with ", length(value), " element(s)",
            if (length(nms)) paste0(" named ", paste(utils::head(nms, 5L), collapse = ", ")) else ""
        ))
    }
    n <- length(value)
    if (is.character(value)) {
        return(if (n == 1L) paste0("the text ", preview(value)) else paste0("a character vector of length ", n))
    }
    if (is.logical(value)) {
        return(if (n == 1L) paste0("the logical value ", as.character(value)[[1L]]) else paste0("a logical vector of length ", n))
    }
    if (is.numeric(value)) {
        return(if (n == 1L) paste0("the number ", as.character(value)[[1L]]) else paste0("a numeric vector of length ", n))
    }
    if (is.factor(value)) {
        return(if (n == 1L) paste0("the factor level ", preview(as.character(value))) else paste0("a factor of length ", n))
    }
    paste0("an object of class ", paste(class(value), collapse = "/"))
}

#' Is one parameter value acceptable for one declared schema type?
#'
#' Extracted so that the type grammar has exactly one definition. The
#' validator needs to ask a narrower question than
#' [sclet_ai_action_param_problems()] -- "which *present* parameters have a
#' type the schema rejects?" -- and two copies of this switch would drift.
#'
#' @param value The parameter value.
#' @param type A single schema type name.
#' @return A single logical value. Unknown type names are permissive.
sclet_ai_param_type_valid <- function(value, type) {
    isTRUE(switch(
        type,
        character_vector = (is.character(value) && length(value) >= 1L) ||
            (is.list(value) && length(value) >= 1L && all(vapply(value, function(x) is.character(x) && length(x) == 1L, logical(1)))),
        integer_vector = ((is.numeric(value) && length(value) >= 1L && all(value == as.integer(value))) ||
            (is.list(value) && length(value) >= 1L && all(vapply(value, function(x) is.numeric(x) && length(x) == 1L && x == as.integer(x), logical(1)))) ||
            !is.null(sclet_ai_range_vector(value))),
        numeric_vector = (is.numeric(value) && length(value) >= 1L) ||
            !is.null(sclet_ai_range_vector(value)),
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
    ))
}

#' Names of present parameters whose value fails the declared schema type
#'
#' Deliberately narrower than [sclet_ai_action_param_problems()]: it reports
#' only *present* parameters of the wrong type. A parameter that is merely
#' absent is not reported, because an action's prerequisite is frequently the
#' thing that supplies the most useful message about a missing parameter
#' (e.g. `start_cluster_missing`), and skipping the prerequisite there would
#' throw away good diagnostics.
#'
#' @param params The parameter list.
#' @param schema The action's `input_schema`.
#' @return Character vector of parameter names (empty if none).
sclet_ai_action_param_type_problems <- function(params, schema) {
    if (is.null(params) || !is.list(params) || !length(schema)) {
        return(character())
    }
    schema_names <- names(schema)
    if (is.null(schema_names) || any(!nzchar(schema_names))) {
        return(character())
    }
    bad <- character()
    for (name in intersect(names(params), schema_names)) {
        value <- params[[name]]
        if (is.null(value)) {
            next
        }
        spec <- schema[[name]]
        if (is.list(spec) && !is.null(spec$type)) {
            spec <- spec$type
        }
        # A multi-value character spec is an enumeration, not a type.
        if (is.character(spec) && length(spec) > 1L) {
            next
        }
        if (!sclet_ai_param_type_valid(value, as.character(spec)[[1L]])) {
            bad <- c(bad, name)
        }
    }
    bad
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
#'   `"integration"`, `"annotation"`, `"rare_cell"`, `"trajectory"`,
#'   `"velocity"`, `"fate"`, `"spatial"`, `"multimodal"`, `"enrichment"`,
#'   `"malignant"`, `"communication"`, `"milo"`, and `"all"`. Defaults to `"read"`.
#' @return A `sclet_ai_execution_registry` containing the requested actions.
#' @export
AIDefaultExecutionRegistry <- function(object, include = "read") {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    include <- unique(as.character(include))
    allowed_groups <- c("read", "preprocess", "dimred", "graph", "cluster", "integration", "annotation", "malignant", "rare_cell", "trajectory", "velocity", "fate", "spatial", "multimodal", "enrichment", "communication", "milo", "all")
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
                    required_hvg = TRUE,
                    # A referenceable name for the set itself. The static value
                    # is a marker, not the genes: `output_schema` is exposed
                    # statically through `available_actions` before any step has
                    # run, so promising real gene names here would push 2000 of
                    # them into every planning context. The real set is written
                    # into the execution summary by
                    # `sclet_ai_execution_summary()` once it exists, and the
                    # projection resolves ${step.output.variable_features} from
                    # there. Downstream steps (scale_data$features,
                    # run_pca$subset_row) are `character_vector`, so a prose
                    # description used to pass the schema check and die in the
                    # gene selector -- measured as 4 rejections in 4 spellings
                    # on pbmc3k (sclet-benchmark-ai, execute_plan).
                    variable_features = "hvg"
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
                    features = sclet_ai_gene_selector(
                        object,
                        sclet_ai_resolve_deferred_selector(object, params$features, "features"),
                        arg = "features"
                    ),
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
                    # features is a character_vector, so a prose description
                    # passes the schema check. Resolving it here -- the same way
                    # the handler will -- means the model is told at validation
                    # time and the repair loop can fix it, instead of the plan
                    # validating and then dying at execution.
                    if (sclet_ai_is_deferred_gene_selector(params$features)) {
                        return(TRUE)
                    }
                    sclet_ai_selector_problem(object, params$features, "features")
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
    if ("integration" %in% include) {
        actions <- c(actions, list(
            run_integration = AIAction(
                name = "run_integration",
                description = "Run batch-correction integration via fastMNN, Harmony, or scVI and register the result.",
                handler = function(object, params) {
                    method <- match.arg(params$method %||% "fastMNN", c("fastMNN", "Harmony", "scVI"))
                    batch <- params$batch
                    features <- sclet_ai_gene_selector(object, params$features, arg = "features")
                    # RunIntegration() chooses the dimensions internally and
                    # has no `dims` argument. Reject a supplied value during
                    # plan validation rather than silently accepting an inert
                    # parameter that gives the user/model a false impression.
                    if (!is.null(params$dims)) {
                        stop("`run_integration` does not support `dims`; omit it so RunIntegration() can choose internally.", call. = FALSE)
                    }
                    RunIntegration(
                        object,
                        method = method,
                        batch = batch,
                        features = features,
                        layer = params$layer %||% NULL,
                        reduction = params$reduction %||% NULL,
                        name = params$name %||% tolower(method)
                    )
                },
                input_schema = list(
                    batch = list(type = "character", required = TRUE),
                    method = c("fastMNN", "Harmony", "scVI"),
                    name = "character",
                    # Intentionally no `dims`: RunIntegration chooses this
                    # internally, so the plan must not advertise a no-op input.
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
                    # A `confirm_design` step earlier in the SAME plan will
                    # write the confirmation at execution time; validation must
                    # accept that promise (deferred projection), otherwise an
                    # integration plan can never validate on an unconfirmed
                    # object. Execution re-checks against the real state.
                    planned_design <- isTRUE(planned$design_confirmed)
                    if (is.null(confirmation) && !planned_design) {
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
                    # Measured per method on this machine (2026-10-14, offline):
                    # fastMNN registers "corrected", Harmony registers
                    # "HARMONY". scVI (unused here) registers "scVI". fastMNN
                    # ignores the `name` parameter and keeps "corrected".
                    allowed_reductions = c("corrected", "HARMONY", "scVI")
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
                description = "Run a differential-expression or all-cluster marker analysis and register evidence. For an experimental-group contrast (condition/treatment/disease) set group = the colData column and ident.1/ident.2 to its levels; add covariates = a colData column such as donor for a paired pseudobulk analysis.",
                handler = function(object, params) {
                    ident.1 <- if (!is.null(params$ident.1)) sclet_ai_vector_param(params$ident.1, "character") else NULL
                    ident.2 <- if (!is.null(params$ident.2)) sclet_ai_vector_param(params$ident.2, "character") else NULL
                    group <- if (!is.null(params$group)) sclet_ai_vector_param(params$group, "character") else NULL
                    covariates <- if (!is.null(params$covariates)) sclet_ai_vector_param(params$covariates, "character") else NULL
                    if (!is.null(group)) {
                        if (length(group) != 1L) stop("`group` must be one colData column name.")
                        ident.1 <- sclet_ai_resolve_group_value(object, group, ident.1)
                        ident.2 <- sclet_ai_resolve_group_value(object, group, ident.2)
                    } else {
                        ident.1 <- sclet_ai_resolve_ident_value(object, ident.1)
                        ident.2 <- sclet_ai_resolve_ident_value(object, ident.2)
                    }
                    result <- RunDEtest(
                        object,
                        ident.1 = ident.1,
                        ident.2 = ident.2,
                        group = group,
                        covariates = covariates,
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
                            idents <- if (!is.null(group)) {
                                as.character(SummarizedExperiment::colData(result)[[group]])
                            } else {
                                Idents(result)
                            }
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
                                id = if (is.null(ident.1)) {
                                    "ev:findallmarkers"
                                } else if (!is.null(group)) {
                                    paste0("ev:detest_", group, "_",
                                        paste(ident.1, collapse = "_"), "_vs_",
                                        paste(ident.2, collapse = "_"))
                                } else {
                                    paste0("ev:detest_", paste(ident.1, collapse = "_"),
                                        if (length(ident.2)) paste0("_vs_", paste(ident.2, collapse = "_")))
                                },
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
                    group = "character",
                    covariates = "character_vector",
                    name = "character",
                    min.pct = "number",
                    logfc.threshold = "number"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    ident.1 <- if (is.null(params$ident.1)) NULL else sclet_ai_vector_param(params$ident.1, "character")
                    ident.2 <- if (is.null(params$ident.2)) NULL else sclet_ai_vector_param(params$ident.2, "character")
                    group <- if (is.null(params$group)) NULL else sclet_ai_vector_param(params$group, "character")
                    covariates <- if (is.null(params$covariates)) NULL else sclet_ai_vector_param(params$covariates, "character")
                    if (!is.null(group)) {
                        cd <- SummarizedExperiment::colData(object)
                        if (length(group) != 1L || !group %in% colnames(cd)) {
                            return(paste0("group_missing: 'group' must name one colData column; available: ",
                                paste(colnames(cd), collapse = ", ")))
                        }
                        values <- unique(as.character(cd[[group]]))
                        if (is.null(ident.1) || is.null(ident.2) ||
                                length(ident.1) != 1L || length(ident.2) != 1L) {
                            return(paste0("contrast_missing: a '", group,
                                "' contrast needs one ident.1 and one ident.2 (levels: ",
                                paste(values, collapse = ", "), ")"))
                        }
                        if (!ident.1 %in% values || !ident.2 %in% values) {
                            return(paste0("contrast_invalid: ident.1/ident.2 must be levels of '",
                                group, "': ", paste(values, collapse = ", ")))
                        }
                        if (identical(ident.1, ident.2)) {
                            return("contrast_invalid: ident.1 and ident.2 must differ")
                        }
                        if (length(covariates)) {
                            missing <- setdiff(covariates, colnames(cd))
                            if (length(missing)) {
                                return(paste0("covariates_missing: ",
                                    paste(missing, collapse = ", ")))
                            }
                            if (!requireNamespace("limma", quietly = TRUE)) {
                                return("covariates_unavailable: covariate-adjusted DE needs the 'limma' package; install it or omit covariates (unadjusted cell-level test)")
                            }
                        }
                        return(TRUE)
                    }
                    active_ident <- ActiveIdent(object)
                    planned_ident <- (planned %||% list())$active_ident
                    if (!is.null(planned_ident) &&
                        !identical(as.character(planned_ident), as.character(active_ident))) {
                        # An earlier step in this plan will (re)set the identity
                        # -- typically find_clusters writing colLabels -- so the
                        # object's *current* identity is not the one this step
                        # will see, and neither its group sizes nor its label
                        # set can be checked before that step runs. Refusing the
                        # plan here would make the ordinary "cluster, then find
                        # markers" request unplannable from a raw object, which
                        # is the common case. Defer: the producing step's output
                        # contract (required_states/active_ident) and the
                        # execution-time check are the authority.
                        return(TRUE)
                    }
                    ident.1 <- sclet_ai_resolve_ident_value(object, ident.1)
                    ident.2 <- sclet_ai_resolve_ident_value(object, ident.2)
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
                    required_states = list(detest = "any"),
                    # A referenceable name for the detest record this step
                    # registers. Static value is a marker, not the id: the
                    # schema is exposed through `available_actions` before any
                    # step has run, and the real id is written into the
                    # execution summary as an observation (see
                    # `sclet_ai_execution_summary()`), which is also where the
                    # `${step.output.detest_id}` projection resolves it from.
                    # Needed because `run_enrichment` requires `detest_id` and
                    # no registered action used to expose one -- measured as
                    # the reason every model-proposed DE -> enrichment chain
                    # died at validation (sclet-benchmark-ai, 2026-10-16,
                    # three spellings across three legs). Same defect class as
                    # `variable_features`, fixed there in 1bddedc.
                    detest_id = "detest_record_id"
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
                        # Bare, `celldex::Name` and `<celldex::Name>` are one
                        # reference; the literal spelling a model writes must
                        # not decide whether an otherwise-planned run executes.
                        ref_name <- sclet_celldex_reference_name(ref_arg)
                        if (!is.null(ref_name) && !requireNamespace("celldex", quietly = TRUE)) {
                            # Name the missing package rather than calling the
                            # reference unsupported: the latter reads as the
                            # user's mistake and suggests no repair.
                            stop(paste0(
                                "celldex is required to resolve ref='", ref_arg, "'; install it via ",
                                "BiocManager::install('celldex') or pass a SummarizedExperiment ref directly."
                            ), call. = FALSE)
                        }
                        fun <- if (is.null(ref_name)) NULL else sclet_celldex_dataset(ref_name)
                        if (is.null(fun)) {
                            stop(paste0(
                                "reference '", ref_arg, "' is not a supported celldex dataset name. ",
                                "Provide a celldex dataset such as HumanPrimaryCellAtlasData ",
                                "or celldex::MonacoImmuneData, or a SummarizedExperiment ref directly."
                            ), call. = FALSE)
                        }
                        ref <- do.call(fun, list())
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
                            id = sclet_ai_annotation_evidence_id(name),
                            kind = "deterministic_summary",
                            values = aggregate_values,
                            claim_level = "consistent_with"
                        )
                        # The parent is the annotation record id, not an evidence
                        # id, so it never resolves as an ancestor; it is kept
                        # because it also enters the node's dependency_group,
                        # which keeps two different annotation runs from
                        # collapsing into one line of support.
                        result <- RecordAIEvidence(
                            result,
                            evidence,
                            source = name,
                            parents = if (!is.null(ann_rec$id)) as.character(ann_rec$id) else character(),
                            scope = NULL
                        )
                    }
                    # The output schema only requires that *some* annotation
                    # state exists, which an earlier step could have supplied.
                    # Verify this action's own record and its bounded evidence
                    # node, so a run that annotated nothing -- or a ledger that
                    # refused the node -- cannot be reported as completed.
                    if (is.null(ann_rec)) {
                        stop(
                            "annotation_no_result: RunSingleR returned without a completed annotation record.",
                            call. = FALSE
                        )
                    }
                    recorded <- tryCatch(
                        sclet_ai_evidence_get_all(result)[[sclet_ai_annotation_evidence_id(name)]],
                        error = function(e) NULL
                    )
                    if (is.null(recorded)) {
                        stop(
                            "annotation_evidence_missing: SingleR completed but its bounded evidence node was not recorded.",
                            call. = FALSE
                        )
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
                    ref_name <- if (ref_is_name) sclet_celldex_reference_name(ref_arg) else NULL
                    if (ref_is_name && is.null(ref_name)) {
                        # A placeholder is not a reference; refuse it at
                        # validation instead of after six upstream steps have
                        # run. It is reported as MISSING, not invalid, because
                        # the choice belongs to the user: the clarification
                        # catalog maps `reference_missing` to the structured
                        # question, while `reference_invalid` has no question
                        # and used to discard the whole plan (2026-10-18).
                        if (sclet_ai_value_is_placeholder(ref_arg)) {
                            return(paste0(
                                "reference_missing: 'ref' is a placeholder ('",
                                ref_arg, "'), not a dataset. Ask the user which ",
                                "reference to use (a celldex dataset such as ",
                                "celldex::MonacoImmuneData, or a ",
                                "SummarizedExperiment/matrix with dimnames)."
                            ))
                        }
                        # A genuine near-miss or typo is still invalid: the
                        # accepted spellings are named in the message.
                        return(paste0(
                            "reference_invalid: 'ref' is a string but not a celldex dataset name: '",
                            ref_arg, "'. Accepted: HumanPrimaryCellAtlasData, ",
                            "celldex::MonacoImmuneData, <celldex::MonacoImmuneData>, ",
                            "or a SummarizedExperiment/matrix with dimnames."
                        ))
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
                    if (labels_is_name && sclet_ai_value_is_placeholder(labels_arg)) {
                        return(paste0(
                            "labels_missing: 'labels' is a placeholder ('",
                            labels_arg, "'), not a reference colData column. Ask ",
                            "the user which labels column to use (e.g. ",
                            "label.main, label.fine) or a per-sample label vector."
                        ))
                    }
                    if (labels_is_vec && ref_is_object) {
                        # One resolver shared with RunSingleR(): a labels string
                        # naming a colData column of the reference is the
                        # documented spelling, and counting it as a one-element
                        # label vector rejected it at validation (measured on
                        # 2026-10-17 with a SummarizedExperiment reference).
                        resolved <- sclet_resolve_singleR_labels(ref_arg, labels_arg)
                        if (!is.null(resolved$problem)) {
                            return(paste0("labels_invalid: ", resolved$problem))
                        }
                        ref_n <- ncol(ref_arg)
                        if (length(resolved$labels) != ref_n) {
                            return(paste0("labels_mismatch: length(labels)=", length(resolved$labels), " does not match ncol(ref)=", ref_n))
                        }
                    }
                    real_req_ns <- base::requireNamespace
                    if (!real_req_ns("SingleR", quietly = TRUE)) {
                        return("optional_package_missing: SingleR is required for run_annotation; install via BiocManager::install('SingleR')")
                    }
                    if (ref_is_name && !is.null(ref_name)) {
                        if (!real_req_ns("celldex", quietly = TRUE)) {
                            return(paste0("optional_package_missing: celldex is required to resolve ref='", ref_arg, "'; install via BiocManager::install('celldex')"))
                        }
                        # `sclet_celldex_reference_name()` accepts the SHAPE of a
                        # dataset name; only celldex knows which names exist. A
                        # near-miss that ends in "Data" (e.g. BlueprintENCODEData
                        # for BlueprintEncodeData) otherwise passes validation and
                        # fails in the handler -- after the upstream steps have
                        # run, which is the class of failure this benchmark paid
                        # for twice.
                        if (is.null(sclet_celldex_dataset(ref_name))) {
                            return(paste0(
                                "reference_invalid: '", ref_arg, "' is not a dataset exported by ",
                                "celldex. Use a real celldex dataset (e.g. HumanPrimaryCellAtlasData, ",
                                "MonacoImmuneData, BlueprintEncodeData) or pass a ",
                                "SummarizedExperiment/matrix reference directly."
                            ))
                        }
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
                description = paste(
                    "Score doublets with scDblFinder, store per-cell doublet class/score in colData, and",
                    "register one bounded evidence node (aggregate counts and score quantiles) so the call",
                    "is citable from the ledger. The class label is one classifier's estimate, not a",
                    "measurement of doubleness, and no cell is removed or filtered.",
                    "No per-cell score or class vector enters the evidence node."
                ),
                handler = function(object, params) {
                    scored <- RunDoubletFinder(object)
                    # The preprocess record id is fixed by RunDoubletFinder, so
                    # the evidence id is stable and a re-run degrades to a no-op
                    # inside the recorder rather than duplicating or erroring.
                    sclet_ai_record_doublet_evidence(scored, name = "doublet_finder")
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
                    required_states = list(preprocess = "any", ai_evidence = "any"),
                    note = paste(
                        "Per-cell doublet calls are added to colData. No cell is removed or filtered by",
                        "this action. The registered evidence carries only the number of cells, the number",
                        "called doublet, that fraction, and score quantiles; a single classifier's call is",
                        "recorded at claim level 'associated' because nothing independent corroborates it."
                    )
                ),
                mutates_object = TRUE,
                allowed_state_types = c("preprocess", "ai_evidence"),
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
                        min_cluster_size = params$min_cluster_size %||% 3L,
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
                    # Density artefacts below this size are dissolved into
                    # their densest neighbour instead of being reported as
                    # populations (the 2419-for-2700-cells defect).
                    min_cluster_size = "integer",
                    name = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    reduction <- params$reduction %||% "PCA"
                    if (!is.character(reduction) || length(reduction) != 1L || !nzchar(reduction)) {
                        return("reduction_invalid: 'reduction' must be a single non-empty reduction name")
                    }
                    # planned$reductions first (deferred projection): a valid
                    # plan runs run_pca in an earlier step and the reduction
                    # only exists in the planned view at validation time.
                    # Reading reducedDimNames() alone made any plan putting
                    # run_rare_cell_detection after run_pca fail with
                    # "reduction 'PCA' is not available" -- found offline by
                    # the doublet/rare benchmark dry run 2026-10-14 before a
                    # billed call, with nothing wrong in the plan.
                    available <- unique(c(
                        as.character(planned$reductions %||% character()),
                        tryCatch(SingleCellExperiment::reducedDimNames(object),
                            error = function(e) character())))
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
    if ("spatial" %in% include) {
        actions <- c(actions, list(
            check_spatial_readiness = AIAction(
                name = "check_spatial_readiness",
                description = paste(
                    "Read-only readiness report for the spatial analyses recorded on this object.",
                    "It reports which spatial state records exist and what their stored model estimates are,",
                    "and it never runs a spatial algorithm, never ranks or recommends a spatial method, and",
                    "never writes an evidence node or mutates the object."
                ),
                handler = function(object, params) {
                    check_spatial_readiness(object, analysis = params$analysis)
                },
                input_schema = list(
                    analysis = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    records <- sclet_get_state_records(object, "spatial")
                    if (!length(records)) {
                        return(paste0(
                            "no_spatial_record: no spatial state record is available on this object; ",
                            "no spatial analysis has been recorded yet."
                        ))
                    }
                    analysis <- params$analysis
                    if (!is.null(analysis) && !analysis %in% names(records)) {
                        return(paste0(
                            "analysis_not_available: spatial analysis '", analysis,
                            "' is not a recorded spatial analysis. Available spatial records: ",
                            paste(names(records), collapse = ", ")
                        ))
                    }
                    TRUE
                },
                returns = "value",
                requires_confirmation = FALSE,
                output_schema = list(
                    type = "object",
                    required_states = list(spatial = "any"),
                    note = paste(
                        "Read-only: reports spatial model estimates and record counts only. It does not rank",
                        "or recommend a spatial method, and no per-cell or per-spot value is included. No",
                        "evidence node is written and the object is not mutated."
                    )
                ),
                mutates_object = FALSE,
                allowed_state_types = c("spatial"),
                estimated_cost = "low",
                idempotent = TRUE
            ),
            summarize_spatial_evidence = AIAction(
                name = "summarize_spatial_evidence",
                description = paste(
                    "Read-only bounded summary of the stored evidence for one recorded spatial analysis.",
                    "It reports aggregate model estimates and record counts only; it never ranks or recommends",
                    "a spatial method, and it never writes an evidence node or mutates the object."
                ),
                handler = function(object, params) {
                    summarize_spatial_evidence(object, analysis = params$analysis)
                },
                input_schema = list(
                    analysis = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    records <- sclet_get_state_records(object, "spatial")
                    if (!length(records)) {
                        return(paste0(
                            "no_spatial_record: no spatial state record is available on this object; ",
                            "no spatial analysis has been recorded yet."
                        ))
                    }
                    analysis <- params$analysis
                    if (!is.null(analysis) && !analysis %in% names(records)) {
                        return(paste0(
                            "analysis_not_available: spatial analysis '", analysis,
                            "' is not a recorded spatial analysis. Available spatial records: ",
                            paste(names(records), collapse = ", ")
                        ))
                    }
                    TRUE
                },
                returns = "value",
                requires_confirmation = FALSE,
                output_schema = list(
                    type = "object",
                    required_states = list(spatial = "any"),
                    note = paste(
                        "Read-only: reports spatial model estimates and record counts only. It does not rank",
                        "or recommend a spatial method, and no per-cell or per-spot value is included. No",
                        "evidence node is written and the object is not mutated."
                    )
                ),
                mutates_object = FALSE,
                allowed_state_types = c("spatial"),
                estimated_cost = "low",
                idempotent = TRUE
            ),
            compare_spatial_analyses = AIAction(
                name = "compare_spatial_analyses",
                description = paste(
                    "Read-only comparison of the spatial analyses recorded on this object, side by side.",
                    "It reports the stored estimates of each record so a difference can be seen; it does not",
                    "rank or recommend a spatial method, never nominates a \"best\" analysis, and never writes",
                    "an evidence node or mutates the object."
                ),
                handler = function(object, params) {
                    compare_spatial_analyses(object)
                },
                input_schema = list(),
                prerequisites = function(object, params, planned = NULL) {
                    records <- sclet_get_state_records(object, "spatial")
                    if (!length(records)) {
                        return(paste0(
                            "no_spatial_record: no spatial state record is available on this object; ",
                            "no spatial analysis has been recorded yet."
                        ))
                    }
                    TRUE
                },
                returns = "value",
                requires_confirmation = FALSE,
                output_schema = list(
                    type = "object",
                    required_states = list(spatial = "any"),
                    note = paste(
                        "Read-only: reports spatial model estimates and record counts only. It does not rank",
                        "or recommend a spatial method, and no per-cell or per-spot value is included. No",
                        "evidence node is written and the object is not mutated."
                    )
                ),
                mutates_object = FALSE,
                allowed_state_types = c("spatial"),
                estimated_cost = "low",
                idempotent = TRUE
            )
        ))
    }
    if ("multimodal" %in% include) {
        actions <- c(actions, list(
            check_multimodal_readiness = AIAction(
                name = "check_multimodal_readiness",
                description = paste(
                    "Read-only readiness report for the multimodal (assay-platform) records on this object.",
                    "It reports which multimodal state records exist and what their stored estimates are, and",
                    "it never registers a modality, never runs a multimodal workflow, never ranks or recommends",
                    "a modality or method, and never writes an evidence node or mutates the object."
                ),
                handler = function(object, params) {
                    check_multimodal_readiness(object, modalities = params$modalities)
                },
                input_schema = list(
                    modalities = "character_vector"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    records <- sclet_get_state_records(object, "multimodal")
                    if (!length(records)) {
                        return(paste0(
                            "no_multimodal_record: no multimodal state record is available on this object; ",
                            "no multimodal analysis has been recorded yet."
                        ))
                    }
                    TRUE
                },
                returns = "value",
                requires_confirmation = FALSE,
                output_schema = list(
                    type = "object",
                    required_states = list(multimodal = "any"),
                    note = paste(
                        "Read-only: reports modality-platform labels, multimodal model estimates and record",
                        "counts only. It does not rank or recommend a modality or method, and no per-cell value,",
                        "barcode or feature name is included. No evidence node is written and the object is not",
                        "mutated."
                    )
                ),
                mutates_object = FALSE,
                allowed_state_types = c("multimodal"),
                estimated_cost = "low",
                idempotent = TRUE
            ),
            summarize_multimodal_evidence = AIAction(
                name = "summarize_multimodal_evidence",
                description = paste(
                    "Read-only bounded summary of the stored evidence for one recorded multimodal analysis.",
                    "It reports aggregate model estimates and record counts only; it never ranks or recommends",
                    "a modality or method, and it never writes an evidence node or mutates the object."
                ),
                handler = function(object, params) {
                    summarize_multimodal_evidence(object, name = params$name)
                },
                input_schema = list(
                    name = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    records <- sclet_get_state_records(object, "multimodal")
                    if (!length(records)) {
                        return(paste0(
                            "no_multimodal_record: no multimodal state record is available on this object; ",
                            "no multimodal analysis has been recorded yet."
                        ))
                    }
                    name <- params$name
                    if (!is.null(name) && !name %in% names(records)) {
                        return(paste0(
                            "analysis_not_available: multimodal analysis '", name,
                            "' is not a recorded multimodal analysis. Available multimodal records: ",
                            paste(names(records), collapse = ", ")
                        ))
                    }
                    TRUE
                },
                returns = "value",
                requires_confirmation = FALSE,
                output_schema = list(
                    type = "object",
                    required_states = list(multimodal = "any"),
                    note = paste(
                        "Read-only: reports modality-platform labels, multimodal model estimates and record",
                        "counts only. It does not rank or recommend a modality or method, and no per-cell value,",
                        "barcode or feature name is included. No evidence node is written and the object is not",
                        "mutated."
                    )
                ),
                mutates_object = FALSE,
                allowed_state_types = c("multimodal"),
                estimated_cost = "low",
                idempotent = TRUE
            ),
            compare_multimodal_evidence = AIAction(
                name = "compare_multimodal_evidence",
                description = paste(
                    "Read-only comparison of the multimodal analyses recorded on this object, side by side.",
                    "It reports the stored estimates of each record so a difference can be seen; it does not",
                    "rank or recommend a modality or method, never nominates a \"best\" analysis, and never",
                    "writes an evidence node or mutates the object."
                ),
                handler = function(object, params) {
                    compare_multimodal_evidence(object)
                },
                input_schema = list(),
                prerequisites = function(object, params, planned = NULL) {
                    records <- sclet_get_state_records(object, "multimodal")
                    if (!length(records)) {
                        return(paste0(
                            "no_multimodal_record: no multimodal state record is available on this object; ",
                            "no multimodal analysis has been recorded yet."
                        ))
                    }
                    TRUE
                },
                returns = "value",
                requires_confirmation = FALSE,
                output_schema = list(
                    type = "object",
                    required_states = list(multimodal = "any"),
                    note = paste(
                        "Read-only: reports modality-platform labels, multimodal model estimates and record",
                        "counts only. It does not rank or recommend a modality or method, and no per-cell value,",
                        "barcode or feature name is included. No evidence node is written and the object is not",
                        "mutated."
                    )
                ),
                mutates_object = FALSE,
                allowed_state_types = c("multimodal"),
                estimated_cost = "low",
                idempotent = TRUE
            )
        ))
    }
    # Malignant-cell identification as actions (merged 2026-10-18).
    #
    # The package grew RunCNVInference()/CallMalignantCells(), but the AI layer
    # could not reach either -- a planner had no way to name a copy-number run or
    # the two-line call. These two actions mirror the functions, with the
    # invariants the feature documents kept in the actions:
    #
    #   * the CNV reference must be EXOGENOUS and named by the caller; the
    #     action refuses to fall back to clusters the analysis itself produced;
    #   * the call is a TWO-line call; running it without a CNV record would
    #     silently downgrade it to signature-only under the stronger name, so the
    #     prerequisite requires the record.
    if ("malignant" %in% include) {
        actions <- c(actions, list(
            run_cnv_inference = AIAction(
                name = "run_cnv_inference",
                description = paste0(
                    "Infer per-cell copy-number burden with inferCNV against an ",
                    "EXOGENOUS reference: name normal cells with `reference_ids` ",
                    "or a colData column with `reference_col`. Never use a cluster ",
                    "the analysis itself produced as the reference."
                ),
                handler = function(object, params) {
                    name <- params$name %||% "cnv"
                    reference_ids <- if (!is.null(params$reference_ids)) {
                        sclet_ai_vector_param(params$reference_ids, "character")
                    } else {
                        NULL
                    }
                    object <- RunCNVInference(
                        object,
                        reference_cells = reference_cells,
                        reference_col = params$reference_col %||% NULL,
                        species = params$species %||% "human",
                        assay = params$assay %||% "counts",
                        hmm = if (is.null(params$hmm)) TRUE else isTRUE(params$hmm),
                        cutoff = params$cutoff %||% 1,
                        num_threads = params$num_threads %||% 1,
                        name = name
                    )
                    records <- sclet_get_analysis(object, "malignant_cnv")
                    record <- records[[name]]
                    if (is.null(record) || !identical(record$status, "completed")) {
                        # A recorded `unavailable` state is an honest answer, not a
                        # failure; no evidence node is written for it.
                        return(object)
                    }
                    burden <- record$burden$burden
                    values <- list(
                        n_cells = as.integer(record$n_cells_submitted),
                        n_cells_retained = as.integer(record$n_cells_retained),
                        n_cells_dropped = as.integer(length(record$dropped_cells)),
                        n_reference_ids = as.integer(record$n_reference_ids),
                        n_genes = as.integer(record$n_genes),
                        burden_q10 = as.numeric(stats::quantile(burden, 0.1,
                            na.rm = TRUE, names = FALSE)),
                        burden_q50 = as.numeric(stats::quantile(burden, 0.5,
                            na.rm = TRUE, names = FALSE)),
                        burden_q90 = as.numeric(stats::quantile(burden, 0.9,
                            na.rm = TRUE, names = FALSE))
                    )
                    tryCatch(RecordAIEvidence(object, list(
                            id = paste0("ev:cnv_", sclet_ai_evidence_id_token(name)),
                            kind = "deterministic_summary",
                            values = values,
                            claim_level = "associated"
                        ),
                        source = name,
                        parents = character(),
                        scope = NULL), error = function(e) object)
                },
                input_schema = list(
                    reference_ids = "character_vector",
                    reference_col = "character",
                    species = c("human", "mouse"),
                    assay = "character",
                    name = "character",
                    hmm = "logical",
                    cutoff = "number",
                    num_threads = "integer"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    if (!requireNamespace("inferCNV", quietly = TRUE)) {
                        return("optional_package_missing: inferCNV is required for copy-number inference; install Bioconductor 'inferCNV'")
                    }
                    assays <- SummarizedExperiment::assayNames(object)
                    assay <- params$assay %||% "counts"
                    if (!assay %in% assays) {
                        return(paste0("counts_assay_missing: assay '", assay,
                            "' not found; inferCNV needs raw counts (available: ",
                            paste(assays, collapse = ", "), ")"))
                    }
                    if (is.null(params$reference_col) && is.null(params$reference_ids)) {
                        return(paste0(
                            "reference_missing: copy-number inference needs an ",
                            "EXOGENOUS reference; name normal cells with ",
                            "reference_ids, or a colData column with reference_col. ",
                            "Do not use clusters the analysis itself produced."))
                    }
                    if (!is.null(params$reference_col)) {
                        cd <- colnames(SummarizedExperiment::colData(object))
                        if (!is.character(params$reference_col) ||
                                !params$reference_col %in% cd) {
                            return(paste0("reference_col_missing: '",
                                params$reference_col, "' is not a colData column (available: ",
                                paste(cd, collapse = ", "), ")"))
                        }
                    }
                    if (!is.null(params$reference_ids)) {
                        cells <- sclet_ai_vector_param(params$reference_ids, "character")
                        known <- intersect(cells, colnames(object))
                        if (length(known) < 5L) {
                            return(paste0("reference_too_small: ", length(known),
                                " supplied reference cell(s) exist on the object; at ",
                                "least 5 are needed for the CNV cutoff, otherwise the ",
                                "CNV line is unavailable and the call degrades silently"))
                        }
                    }
                    TRUE
                },
                returns = "sce",
                output_schema = list(
                    required_states = list(annotation = "any"),
                    cnv_id = "cnv_record_id"
                ),
                allowed_state_types = "annotation",
                mutates_object = TRUE,
                requires_confirmation = TRUE,
                estimated_cost = "high",
                idempotent = FALSE
            ),
            call_malignant_cells = AIAction(
                name = "call_malignant_cells",
                description = paste0(
                    "Combine copy-number and transcriptional evidence into the ",
                    "four-state malignant call (malignant / likely_malignant / ",
                    "undetermined / non_malignant). Requires a completed CNV run ",
                    "and either a signature score column or signature genes; ",
                    "`identity_markers` names a competing identity programme to ",
                    "demote."
                ),
                handler = function(object, params) {
                    name <- params$name %||% "malignant_call"
                    signature <- NULL
                    if (!is.null(params$signature)) {
                        genes <- sclet_ai_vector_param(params$signature, "character")
                        signature <- stats::setNames(list(genes), "malignant")
                    }
                    identity_markers <- if (!is.null(params$identity_markers)) {
                        sclet_ai_vector_param(params$identity_markers, "character")
                    } else {
                        NULL
                    }
                    object <- CallMalignantCells(
                        object,
                        cnv_id = params$cnv_id %||% "cnv",
                        signature = signature,
                        signature_score = params$signature_score %||% NULL,
                        transcriptomic_cut = params$transcriptomic_cut %||% NULL,
                        cnv_cut = params$cnv_cut %||% NULL,
                        reference_quantile = params$reference_quantile %||% 0.99,
                        identity_genes = identity_genes,
                        identity_demote = params$identity_demote %||% 0.1,
                        assay = params$assay %||% NULL,
                        name = name
                    )
                    state_col <- SummarizedExperiment::colData(object)[[
                        paste0(name, "_state")]]
                    state_counts <- table(as.character(state_col))
                    state_rec <- tryCatch(
                        sclet_get_state_record(object, "annotation", name),
                        error = function(e) NULL)
                    state_summary <- state_rec$summary %||% list()
                    values <- list(
                        n_cells = as.integer(ncol(object)),
                        n_malignant = as.integer(state_counts[["malignant"]] %||% 0L),
                        n_likely_malignant = as.integer(state_counts[["likely_malignant"]] %||% 0L),
                        n_undetermined = as.integer(state_counts[["undetermined"]] %||% 0L),
                        n_non_malignant = as.integer(state_counts[["non_malignant"]] %||% 0L),
                        # Read from the call's own state record, where
                        # CallMalignantCells() computed them, instead of
                        # re-deriving a second answer here.
                        cnv_available = as.integer(isTRUE(state_summary$cnv_available)),
                        identity_conflict = as.integer(isTRUE(state_summary$identity_conflict))
                    )
                    tryCatch(RecordAIEvidence(object, list(
                            id = paste0("ev:malignant_", sclet_ai_evidence_id_token(name)),
                            kind = "deterministic_summary",
                            values = values,
                            claim_level = "associated"
                        ),
                        source = name,
                        parents = character(),
                        scope = NULL), error = function(e) object)
                },
                input_schema = list(
                    cnv_id = "character",
                    signature = "character_vector",
                    signature_score = "character",
                    transcriptomic_cut = "number",
                    cnv_cut = "number",
                    reference_quantile = "number",
                    identity_markers = "character_vector",
                    identity_demote = "number",
                    assay = "character",
                    name = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    cnv_id <- params$cnv_id %||% "cnv"
                    records <- sclet_get_analysis(object, "malignant_cnv")
                    planned_cnv <- planned$cnv_ids %||% character()
                    have_real <- !is.null(records) && !is.null(records[[cnv_id]])
                    have_planned <- length(planned_cnv) &&
                        cnv_id %in% planned_cnv
                    if (!have_real && !have_planned) {
                        return(paste0(
                            "cnv_record_missing: no CNV run '", cnv_id,
                            "' on the object or earlier in this plan; run ",
                            "run_cnv_inference first. The malignant call is a ",
                            "TWO-line call and will not be silently downgraded ",
                            "to signature-only."))
                    }
                    if (is.null(params$signature) && is.null(params$signature_score)) {
                        return(paste0(
                            "signature_missing: the malignant call needs a ",
                            "transcriptional line -- pass signature_score (a colData ",
                            "column) or signature"))
                    }
                    if (!is.null(params$signature_score)) {
                        cd <- colnames(SummarizedExperiment::colData(object))
                        if (!is.character(params$signature_score) ||
                                !params$signature_score %in% cd) {
                            return(paste0("signature_score_missing: '",
                                params$signature_score, "' is not a colData column (available: ",
                                paste(cd, collapse = ", "), ")"))
                        }
                    }
                    TRUE
                },
                returns = "sce",
                output_schema = list(
                    required_states = list(annotation = "any")
                ),
                allowed_state_types = "annotation",
                mutates_object = TRUE,
                requires_confirmation = TRUE,
                estimated_cost = "low",
                idempotent = FALSE
            )
        ))
    }
    if ("enrichment" %in% include) {
        actions <- c(actions, list(
            run_enrichment = AIAction(
                name = "run_enrichment",
                description = paste(
                    "Summarize a completed differential-expression marker table into enrichment categories",
                    "with clusterProfiler, and register bounded aggregate evidence for it. The marker table",
                    "and its cluster column are required: enrichment is an interpretation layer over the DE",
                    "result, not a second measurement of it. The database, ontology, key type and the",
                    "significance cutoffs are the user's biological modelling choices and are never picked",
                    "or recommended here. No gene symbol and no category name or accession is recorded as",
                    "evidence; the full result stays in the state record for get_enrichment()."
                ),
                handler = function(object, params) {
                    name <- params$name %||% "enrichment"
                    # A re-run that produces nothing must not be reported as a
                    # success: RunEnrichment() warns and returns the object
                    # unchanged when no marker passes the cutoffs, when the
                    # marker table has no p-value column (whose comparison
                    # yields logical(0), the same zero-row outcome), or when it
                    # lacks a logFC column. Comparing the record's created_at
                    # before and after distinguishes "this call wrote a result"
                    # from "a stale record from an earlier call is still here".
                    before <- tryCatch(
                        sclet_get_state_record(object, "enrichment", name),
                        error = function(e) NULL
                    )
                    before_created <- before$created_at %||% NA_character_
                    # The marker-table gate may have been skipped at
                    # validation time: in a chained plan the detest record was
                    # only a promise then. Now the upstream step has really run,
                    # so apply the same helper here. Without this, a chain whose
                    # DE table has no cluster column would validate, run, write
                    # a degenerate one-group record and report success, because
                    # RunEnrichment() warns rather than errors in that case and
                    # the record it writes is non-empty.
                    shape_problem <- sclet_ai_enrichment_marker_problem(
                        object, params$detest_id)
                    if (!is.null(shape_problem)) {
                        stop(shape_problem, call. = FALSE)
                    }
                    result <- RunEnrichment(
                        object,
                        id = params$detest_id,
                        p_val_cutoff = params$p_val_cutoff %||% 0.05,
                        logfc_cutoff = params$logfc_cutoff %||% 0.25,
                        db = params$db %||% "GO",
                        orgDb = params$orgDb %||% "org.Hs.eg.db",
                        keyType = params$keyType %||% "SYMBOL",
                        ont = params$ont %||% "BP",
                        name = name
                    )
                    after <- tryCatch(
                        sclet_get_state_record(result, "enrichment", name),
                        error = function(e) NULL
                    )
                    wrote_record <- !is.null(after) &&
                        !identical(as.character(after$created_at %||% NA_character_),
                            as.character(before_created))
                    if (!wrote_record) {
                        stop(paste0(
                            "enrichment_no_result: run_enrichment produced no enrichment record '", name,
                            "'. Either no marker passed p_val_cutoff / logfc_cutoff, or the marker table",
                            " lacks a p-value or logFC column for the cutoffs to be applied. Lower the",
                            " cutoffs, or re-run RunDEtest() so the marker table carries them. Nothing",
                            " was registered as evidence and no claim is supported by this step."
                        ), call. = FALSE)
                    }
                    sclet_ai_record_enrichment_evidence(
                        result,
                        name = name,
                        detest_id = params$detest_id
                    )
                },
                input_schema = list(
                    detest_id = "character",
                    db = "character",
                    orgDb = "character",
                    keyType = "character",
                    ont = "character",
                    p_val_cutoff = "number",
                    logfc_cutoff = "number",
                    name = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    detest <- tryCatch(get_detest(object, params$detest_id),
                        error = function(e) NULL)
                    # TRUE only when the record is on the object right now, so
                    # its contents can be inspected. A chained plan has only a
                    # promise (see below), which is not the same thing and must
                    # not be treated as if it were.
                    record_available <- !is.null(detest)
                    if (!record_available) {
                        # A plan that chains run_de_test -> run_enrichment
                        # cannot have its upstream record at validation time,
                        # because nothing has executed yet. The projection in
                        # sclet_ai_planned_state_records() is the promise that
                        # an earlier step registers it, and the handler
                        # re-checks the marker table for real at execution.
                        planned_detest <- as.character(
                            planned$state_records$detest %||% character())
                        requested <- params$detest_id
                        # A well-formed output reference counts as pointing at
                        # the promise: `detest_id =
                        # "${<step>.output.detest_id}"` resolves at execution
                        # from the producing step's summary, which carries the
                        # real record id (run_de_test now exposes it).
                        # Demanding the LITERAL id here rejected every chained
                        # plan that used the reference -- the exact spelling
                        # the 2026-10-16 benchmark measured as the only
                        # legitimate one once detest_id became referenceable.
                        # The reference shape itself ties the id to a step of
                        # this plan; execution stays the authority on whether
                        # that step really registered a detest record.
                        is_output_reference <- is.character(requested) &&
                            length(requested) == 1L &&
                            grepl("^\\$\\{[A-Za-z0-9_.-]+\\.output\\.[A-Za-z0-9_.-]+\\}$", requested)
                        if (is_output_reference) {
                            requested <- planned_detest
                        }
                        if (length(requested) == 1L && nzchar(requested)) {
                            if (!requested %in% planned_detest) {
                                return(paste0(
                                    "detest_state_missing: run_enrichment is computed from a completed ",
                                    "differential-expression marker table; run RunDEtest() first. ",
                                    "Requested detest id '", as.character(requested),
                                    "' is neither an existing record nor produced by an earlier step of ",
                                    "this plan. Planned detest records: ",
                                    if (length(planned_detest)) paste(planned_detest, collapse = ", ") else "none", "."
                                ))
                            }
                        } else if (!length(planned_detest)) {
                            return(paste0(
                                "detest_state_missing: run_enrichment is computed from a completed ",
                                "differential-expression marker table; run RunDEtest() first."
                            ))
                        }
                    }
                    # The marker table's shape can only be judged once the
                    # record exists. When it is merely planned, this gate is
                    # skipped here and the handler applies the same helper
                    # after the upstream step has really run.
                    if (record_available) {
                        problem <- sclet_ai_enrichment_marker_problem(
                            object, params$detest_id)
                        if (!is.null(problem)) return(problem)
                    }
                    db <- params$db %||% "GO"
                    if (!is.character(db) || length(db) != 1L || !db %in% c("GO", "KEGG")) {
                        return(paste0(
                            "db_invalid: run_enrichment accepts db = \"GO\" or \"KEGG\"; got '",
                            paste(as.character(db), collapse = ", "), "'. The database is a modelling ",
                            "choice and is never selected automatically."
                        ))
                    }
                    org_db <- params$orgDb %||% "org.Hs.eg.db"
                    missing_packages <- c(
                        if (!base::requireNamespace("clusterProfiler", quietly = TRUE)) "clusterProfiler" else NULL,
                        if (is.character(org_db) && length(org_db) == 1L && nzchar(org_db) &&
                            !base::requireNamespace(org_db, quietly = TRUE)) org_db else NULL
                    )
                    if (length(missing_packages)) {
                        return(paste0(
                            "enrichment_package_missing: run_enrichment needs ",
                            paste(missing_packages, collapse = " and "),
                            "; install it via BiocManager::install()."
                        ))
                    }
                    TRUE
                },
                returns = "sce",
                requires_confirmation = TRUE,
                output_schema = list(
                    required_states = list(enrichment = "any", ai_evidence = "any"),
                    note = paste(
                        "Enrichment categories are derived from the DE marker table and are model output,",
                        "not independent measurements: they inherit that table's thresholds and confounders,",
                        "and the recorded evidence declares the DE node as its parent so a claim ceiling",
                        "keeps counting the two as a single line of support. No gene symbol, category name",
                        "or category accession enters the evidence layer, and no cutoff, database or",
                        "ontology is recommended by this action."
                    )
                ),
                mutates_object = TRUE,
                allowed_state_types = c("detest", "enrichment", "ai_evidence"),
                estimated_cost = "medium",
                idempotent = TRUE
            )
        ))
    }
    if ("communication" %in% include) {
        actions <- c(actions, list(
            run_cci = AIAction(
                name = "run_cci",
                description = paste(
                    "Infer cell-cell interactions with the engine the user names, and register bounded",
                    "aggregate evidence for it. The method is required and never defaulted: CellChat,",
                    "CellPhoneDB and NicheNet use different ligand-receptor databases, expression tests and",
                    "sender/receiver framing, so choosing between them is a modelling decision, not a",
                    "parameter to fill in. No ligand, receptor or pathway identity is recorded as",
                    "evidence; the interaction table stays in the analysis layer for get_cci()."
                ),
                handler = function(object, params) {
                    method <- params$method
                    # Record ids are engine-specific and the caller may omit
                    # name, so the newly written record is located by comparing
                    # created_at across communication records rather than
                    # assumed. This also detects a run that wrote no record.
                    before <- sclet_get_state_records(object, "communication")
                    before_created <- vapply(before, function(record) {
                        as.character(record$created_at %||% "")
                    }, character(1L))
                    result <- RunCCI(
                        object,
                        method = method,
                        keep_object = isTRUE(params$keep_object %||% FALSE),
                        group = params$group,
                        assay_name = params$assay_name,
                        layer = params$layer,
                        name = params$name
                    )
                    after <- sclet_get_state_records(result, "communication")
                    created_after <- vapply(after, function(record) {
                        as.character(record$created_at %||% "")
                    }, character(1L))
                    written <- names(created_after)[vapply(names(created_after), function(id) {
                        previous <- if (id %in% names(before_created)) {
                            before_created[[id]]
                        } else {
                            ""
                        }
                        !identical(created_after[[id]], previous)
                    }, logical(1L))]
                    if (!length(written)) {
                        stop(paste0(
                            "communication_no_result: run_cci with method '", method,
                            "' wrote no communication record, so nothing can be reported or claimed from",
                            " it. Check that the group column resolves to at least two populations with",
                            " cells in the expression source, then re-run."
                        ), call. = FALSE)
                    }
                    sclet_ai_record_communication_evidence(
                        result,
                        name = written[[1L]],
                        method = method
                    )
                },
                input_schema = list(
                    method = "character",
                    group = "character",
                    assay_name = "character",
                    layer = "character",
                    name = "character",
                    keep_object = "logical"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    engines <- sclet_ai_clarification_engines()
                    method <- params$method
                    if (is.null(method) || !is.character(method) || length(method) != 1L ||
                            !nzchar(trimws(method))) {
                        return(paste0(
                            "method_missing: run_cci requires an explicit method because the engines are",
                            " not interchangeable (different ligand-receptor databases and statistical tests).",
                            " Pass one of: ", paste(engines, collapse = ", "),
                            ". No engine is chosen for the user."
                        ))
                    }
                    if (!method %in% engines) {
                        return(paste0(
                            "method_invalid: run_cci accepts method = \"CellChat\", \"CellPhoneDB\" or",
                            " \"NicheNet\"; got '", method, "'. No engine is substituted automatically."
                        ))
                    }
                    columns <- colnames(SummarizedExperiment::colData(object))
                    group_key <- params$group %||% tryCatch(ActiveIdent(object), error = function(e) NULL)
                    if (!is.null(group_key) && !as.character(group_key) %in% columns &&
                            !identical(as.character(group_key), "colLabels")) {
                        group_key <- NULL
                    }
                    if (is.null(group_key)) {
                        candidates <- grep("cluster|leiden|louvain|ident", columns,
                            ignore.case = TRUE, value = TRUE)
                        if (length(candidates)) group_key <- candidates[[1L]]
                    }
                    if (is.null(group_key)) {
                        return(paste0(
                            "group_column_missing: run_cci groups cells into sender and receiver",
                            " populations and none is resolvable; set an active identity or pass 'group'.",
                            " Available colData columns: ", paste(columns, collapse = ", "), "."
                        ))
                    }
                    if (!length(SummarizedExperiment::assayNames(object))) {
                        return(
                            "expression_source_missing: run_cci needs an expression assay and this object has none."
                        )
                    }
                    required_packages <- c(CellChat = "CellChat", CellPhoneDB = "basilisk",
                        NicheNet = "nichenetr")
                    if (identical(method, "CellPhoneDB")) {
                        # CellPhoneDB is the one engine that crosses into Python.
                        missing_bridge <- c("reticulate", "basilisk")[
                            !vapply(c("reticulate", "basilisk"),
                                base::requireNamespace, logical(1L), quietly = TRUE)]
                        if (length(missing_bridge)) {
                            return(paste0(
                                "communication_backend_unavailable: CellPhoneDB runs through ",
                                paste(missing_bridge, collapse = " and "),
                                "; install them to use this engine."
                            ))
                        }
                    }
                    needed <- required_packages[[method]]
                    if (!base::requireNamespace(needed, quietly = TRUE)) {
                        return(paste0(
                            "communication_package_missing: run_cci with method '", method, "' needs the '",
                            needed, "' package; install it via BiocManager::install('", needed, "')."
                        ))
                    }
                    TRUE
                },
                returns = "sce",
                requires_confirmation = TRUE,
                output_schema = list(
                    required_states = list(communication = "any", ai_evidence = "any"),
                    note = paste(
                        "Inferred interactions are model estimates computed from expression plus a curated",
                        "ligand-receptor database; a predicted pair is not evidence that cells signalled, and",
                        "an absent pair is uninformative because each engine filters against its own database.",
                        "Only group-flow counts and strength quantiles reach the evidence layer, and every",
                        "communication run on an object shares a dependency group, so agreement between two",
                        "engines over the same expression matrix and labels is not two independent lines of",
                        "support. This action never chooses an engine, and removes or filters no cell."
                    )
                ),
                mutates_object = TRUE,
                allowed_state_types = c("communication", "ai_evidence"),
                estimated_cost = "high",
                idempotent = TRUE
            )
        ))
    }
    if ("milo" %in% include) {
        actions <- c(actions, list(
            run_milo = AIAction(
                name = "run_milo",
                description = paste(
                    "Test differential abundance of neighbourhoods with Milo/miloR and register one",
                    "bounded evidence node. The sample and condition design roles are mandatory: this",
                    "action never guesses which colData column plays which role, and they must be",
                    "declared through ConfirmAIDesignSemantics() naming the columns yourself. Only",
                    "neighbourhood counts and FDR quantiles reach the evidence layer; the design frame,",
                    "the model formula and every sample or condition label stay in the state record."
                ),
                handler = function(object, params) {
                    name <- params$name %||% "milo"
                    # Same silent-path discipline as enrichment and CCI: Milo
                    # can reach testNhoods and write a record with nothing
                    # usable in it, so created_at is compared before and after
                    # and an absent record is a typed failure, not a success.
                    before <- sclet_get_state_records(object, "milo")
                    before_created <- vapply(before, function(record) {
                        as.character(record$created_at %||% "")
                    }, character(1L))
                    result <- RunMilo(
                        object,
                        sample_col = params$sample_col,
                        condition_col = params$condition_col,
                        covariates = params$covariates,
                        contrasts = params$contrasts,
                        prop = params$prop %||% 0.1,
                        name = name
                    )
                    after <- sclet_get_state_records(result, "milo")
                    created_after <- vapply(after, function(record) {
                        as.character(record$created_at %||% "")
                    }, character(1L))
                    written <- names(created_after)[vapply(names(created_after), function(id) {
                        previous <- if (id %in% names(before_created)) before_created[[id]] else ""
                        !identical(created_after[[id]], previous)
                    }, logical(1L))]
                    if (!length(written)) {
                        stop(paste0(
                            "milo_no_result: run_milo wrote no milo record for name '", name,
                            "', so nothing can be reported or claimed from it. Check that the design ",
                            "columns resolve and that the neighbourhood graph was built, then re-run."
                        ), call. = FALSE)
                    }
                    sclet_ai_record_milo_evidence(result, name = written[[1L]])
                },
                input_schema = list(
                    sample_col = list(type = "character", required = TRUE),
                    condition_col = list(type = "character", required = TRUE),
                    covariates = "character_vector",
                    contrasts = "character_vector",
                    prop = "number",
                    name = "character"
                ),
                prerequisites = function(object, params, planned = NULL) {
                    if (!base::requireNamespace("miloR", quietly = TRUE)) {
                        return(paste0(
                            "milo_package_missing: run_milo needs miloR; ",
                            "install it via BiocManager::install('miloR')."
                        ))
                    }
                    columns <- colnames(SummarizedExperiment::colData(object))
                    # Deferred projection, same as every other consuming gate:
                    # a plan that runs run_pca earlier in the SAME plan has not
                    # executed yet, so `reducedDimNames(object)` is empty at
                    # validation time and this gate refused every chained milo
                    # plan (measured 2026-10-16, three records). The planned
                    # view carries the reductions earlier steps declare, and
                    # execution stays the authority through RunMilo()'s own
                    # reducedDim check.
                    reductions <- planned$reductions %||%
                        tryCatch(SingleCellExperiment::reducedDimNames(object),
                            error = function(e) character())
                    if (!"PCA" %in% reductions) {
                        return(paste0(
                            "reduction_missing: run_milo builds its neighbourhood graph on a PCA ",
                            "embedding; run RunPCA() first. Available reductions: ",
                            paste(reductions, collapse = ", "), "."
                        ))
                    }
                    for (role in c("sample_col", "condition_col")) {
                        value <- params[[role]]
                        plain <- role
                        if (is.null(value) || !is.character(value) || length(value) != 1L ||
                                !nzchar(trimws(value))) {
                            return(paste0(
                                plain, "_missing: run_milo requires the '", plain,
                                "' column name to be supplied explicitly; the design role is never ",
                                "guessed from column names. Available colData columns: ",
                                paste(columns, collapse = ", "), "."
                            ))
                        }
                        if (!value %in% columns) {
                            return(paste0(
                                plain, "_not_in_colData: '", value, "' is not a colData column. ",
                                "Available colData columns: ", paste(columns, collapse = ", "), "."
                            ))
                        }
                    }
                    if (is.null(params$contrasts) || !length(params$contrasts)) {
                        return(paste0(
                            "contrasts_missing: run_milo needs explicit contrasts, e.g. ",
                            "\"trtA - ctrl\". The condition contrast is a statement about the design, ",
                            "not something to infer. Distinct condition values: ",
                            paste(utils::head(unique(as.character(
                                SummarizedExperiment::colData(object)[[params$condition_col]]
                            )), 8L), collapse = ", "), "."
                        ))
                    }
                    # Design semantics gate, shared with run_integration rather
                    # than forked: the user's own declaration must cover both
                    # roles, mapped onto the role names Milo uses.
                    confirmation <- tryCatch(sclet_ai_find_design_confirmation(object, design = list(
                        sample = params$sample_col,
                        condition = params$condition_col
                    )), error = function(e) NULL)
                    if (is.null(confirmation)) {
                        return(paste0(
                            "design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, ",
                            "design = list(sample = '", params$sample_col,
                            "', condition = '", params$condition_col,
                            "')) before running Milo. Differential abundance attributes neighbourhoods ",
                            "to samples and contrasts between conditions, so both roles must be the ",
                            "user's own declaration; no column is inferred here."
                        ))
                    }
                    TRUE
                },
                returns = "sce",
                requires_confirmation = TRUE,
                output_schema = list(
                    required_states = list(milo = "any", ai_evidence = "any"),
                    note = paste(
                        "Differential abundance results are associations between neighbourhoods and",
                        "conditions under a chosen graph and test; they are not causal and not",
                        "measured effects. Counts and FDR quantiles only reach the evidence layer;",
                        "the design frame, model formula and all sample or condition labels stay in",
                        "the state record. Every Milo node on an object shares a dependency group, so",
                        "a second contrast over the same graph is not a second line of support."
                    )
                ),
                mutates_object = TRUE,
                allowed_state_types = c("reduction", "milo", "ai_evidence"),
                estimated_cost = "high",
                idempotent = TRUE
            )
        ))
    }
    AIExecutionRegistry(actions)
}

# Resolve the CellRank trajectory record a fate action should act on.
#
# Returns list(id = <record id>, record = <record>, message = NULL) on success,
# or list(id = NULL, record = NULL, message = "<typed reason>") when no usable
# record exists. Only completed records whose method is "CellRank" (or whose
# artifacts$analysis_key is "cellrank", matching
# sclet_ai_record_fate_evidence()) count. A record id supplied by the user must
# name one of those exactly.
sclet_ai_resolve_cellrank_record <- function(object, cellrank_id = NULL) {
    records <- tryCatch(sclet_get_state_records(object, "trajectory"), error = function(e) list())
    is_cellrank <- vapply(records, function(record) {
        identical(as.character(record$status)[[1L]], "completed") &&
            (identical(as.character(record$method %||% ""), "CellRank") ||
                identical(as.character(record$artifacts$analysis_key %||% ""), "cellrank") ||
                identical(tolower(as.character(record$artifacts$analysis_key %||% "")), "cellrank"))
    }, logical(1L))
    ids <- names(records)[is_cellrank]
    if (!length(ids)) {
        return(list(id = NULL, record = NULL, message = paste0(
            "cellrank_record_missing: no completed CellRank trajectory record is available; ",
            "run run_fate_analysis first. Available trajectory records: ",
            if (length(records)) paste(names(records), collapse = ", ") else "(none)"
        )))
    }
    if (!is.null(cellrank_id) && (!is.character(cellrank_id) || length(cellrank_id) != 1L ||
        is.na(cellrank_id) || !nzchar(trimws(cellrank_id)))) {
        return(list(id = NULL, record = NULL, message = paste0(
            "cellrank_id_unknown: cellrank_id must be a single non-empty record id. ",
            "Available CellRank records: ", paste(ids, collapse = ", ")
        )))
    }
    if (!is.null(cellrank_id) && !cellrank_id %in% ids) {
        return(list(id = NULL, record = NULL, message = paste0(
            "cellrank_id_unknown: trajectory record '", cellrank_id,
            "' is not a completed CellRank record. Available CellRank records: ",
            paste(ids, collapse = ", ")
        )))
    }
    id <- if (is.null(cellrank_id)) ids[[1L]] else cellrank_id
    list(id = id, record = records[[id]], message = NULL)
}

# The terminal-state codes actually present in the record's stored
# terminal-state column.
#
# Codes are derived exactly as sclet_ai_record_fate_evidence() derives them --
# paste0("cluster_", seq_along(unique(<column values>))) -- so a code the user
# supplies can only be accepted when it names a state the stored result really
# carries. The set is returned verbatim; nothing here orders or scores states.
sclet_ai_terminal_state_codes <- function(object, record) {
    terminal_state_column <- as.character(record$artifacts$terminal_state_col %||% "")
    if (length(terminal_state_column) != 1L || !nzchar(terminal_state_column)) return(character())
    coldata <- tryCatch(SummarizedExperiment::colData(object), error = function(e) NULL)
    if (is.null(coldata) || !terminal_state_column %in% colnames(coldata)) return(character())
    labels <- as.character(coldata[[terminal_state_column]])
    labels <- unique(labels[!is.na(labels) & nzchar(labels)])
    if (!length(labels)) return(character())
    paste0("cluster_", seq_along(labels))
}

# Validate user-supplied terminal-state codes, returning them unchanged.
#
# This is the handler-side twin of the terminal_state_unknown prerequisite: it
# re-checks rather than repairing, and preserves the user's own ordering and
# selection. Codes that the stored result does not carry are dropped, so a
# confirmation can never smuggle in a state the model output does not contain.
sclet_ai_validated_terminal_state_codes <- function(object, record, codes) {
    codes <- sclet_ai_vector_param(codes, "character")
    valid <- sclet_ai_terminal_state_codes(object, record)
    unique(codes[codes %in% valid])
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

sclet_ai_record_velocity_evidence <- function(object, name) {
    record <- sclet_get_state_record(object, "velocity", name)
    if (is.null(record)) return(object)
    velocity_coldata <- intersect(
        c("velocity_pseudotime", "velocity_confidence", "root_cells", "end_points"),
        colnames(SummarizedExperiment::colData(object))
    )
    column_summaries <- list()
    for (column in velocity_coldata) {
        values <- as.numeric(SummarizedExperiment::colData(object)[[column]])
        values <- values[is.finite(values)]
        if (!length(values)) next
        quantiles <- round(unname(stats::quantile(values, c(0.25, 0.5, 0.75))), 4L)
        column_summaries[[column]] <- list(
            n_observed = as.integer(length(values)),
            mean = round(mean(values), 4L),
            quantile_25 = as.numeric(quantiles[[1L]]),
            median = as.numeric(quantiles[[2L]]),
            quantile_75 = as.numeric(quantiles[[3L]]),
            max = round(max(values), 4L)
        )
    }
    mode <- record$params$mode %||% "deterministic"
    mode_code <- match(mode, c("deterministic", "stochastic", "dynamical"))
    if (is.na(mode_code)) mode_code <- 0L
    reductions_written <- record$reductions %||% character()
    values <- list(
        mode_code = as.integer(mode_code),
        n_reductions_written = as.integer(length(reductions_written)),
        n_fields_written = as.integer(length(velocity_coldata)),
        velocity_is_model_estimate = TRUE,
        raw_values_included = FALSE
    )
    if (length(column_summaries)) {
        values$n_summarized_fields <- as.integer(length(column_summaries))
        for (column in names(column_summaries)) {
            values[[column]] <- column_summaries[[column]]
        }
    }
    evidence <- list(
        id = paste0("ev:velocity_", name),
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

#' Field names a plan may reference in an earlier step's output
#'
#' One definition, shared by the runtime output summary and by validation.
#' They must agree. When they did not, validation resolved a reference against
#' the action's declared `output_schema` (which has keys like `required_hvg`)
#' while execution resolved it against the runtime summary (which does not), so
#' a plan containing `${step_2.output.required_hvg}` *passed validation and
#' then failed at execution* with "Plan output field is not available":
#'
#'   step_3: Plan output field is not available: ${step_2.output.required_hvg}
#'
#' Reproduced in two of four real benchmark runs. Validation that accepts what
#' execution cannot do is worse than stricter validation, because the failure
#' surfaces only after the plan has been paid for.
#'
#' @param descriptor An action descriptor (`NULL` for the base keys only).
#' @return Character vector of referenceable field names.
#' Contract keys an action's output summary may expose
#'
#' One definition, shared by the summary builder and by the validator's view of
#' which field names a plan reference or a success criterion may name. They must
#' agree: when they did not, validation accepted references execution could not
#' resolve.
#'
#' `note` is excluded deliberately. Those entries hold long prose and `note` is
#' hard-denied by the privacy classifier in `R/ai-privacy.R`, so carrying it
#' onto the ledger path would be wrong on both counts. Keys whose declared value
#' is `NULL` are excluded too, because assigning `NULL` into a list removes the
#' element rather than creating one, which would silently desynchronise the two
#' functions.
#'
#' @param descriptor An action descriptor (`NULL` for no contract keys).
#' @return Character vector of exposable contract key names.
sclet_ai_exposable_contract_keys <- function(descriptor = NULL) {
    contract <- (descriptor %||% list())$output_schema %||% list()
    keys <- names(contract)
    if (!length(keys)) {
        return(character())
    }
    keep <- vapply(keys, function(k) {
        !identical(k, "note") && !is.null(contract[[k]])
    }, logical(1))
    unique(keys[keep])
}

sclet_ai_referenceable_output_fields <- function(descriptor = NULL) {
    if (!is.null(descriptor) && !identical(descriptor$returns, "sce")) {
        # A `returns = "value"` action's execution summary is the value itself
        # (or a class/length pair when it is long or non-atomic), so no named
        # field can be promised before the value exists. Offering the SCE field
        # list regardless made a criterion such as `field = "n_cells"` on a
        # check_qc step pass validation and then report `not_met` at execution
        # with "field n_cells not found in output of step" -- the same
        # validation-accepts-what-execution-cannot-do defect this function
        # exists to prevent.
        return(character())
    }

    unique(c(
        "class", "n_features", "n_cells", "assays", "reductions",
        "active_assay", "active_reduction", "reduction", "fingerprint",
        "active_ident", "active_graph",
        sclet_ai_exposable_contract_keys(descriptor)
    ))
}

sclet_ai_execution_summary <- function(value, descriptor = NULL) {
    if (inherits(value, "SingleCellExperiment")) {
        active_reduction <- tryCatch(DefaultReduction(value), error = function(e) NULL)
        summary <- list(
            class = class(value),
            n_features = nrow(value),
            n_cells = ncol(value),
            assays = SummarizedExperiment::assayNames(value),
            reductions = SingleCellExperiment::reducedDimNames(value),
            active_assay = sclet_get_active_assay(value),
            fingerprint = tryCatch(GetAnalysisLedger(value)$fingerprint, error = function(e) NULL)
        )
        # Only when the object really has one: a NULL element still counts as
        # "present" for the contract-fill loop below, so declaring `reduction`
        # unconditionally suppressed the schema's own value on a plain object
        # (measured in test-ai-phase3's produce/consume chain).
        if (!is.null(active_reduction)) {
            summary$active_reduction <- active_reduction
            # Alias plans naturally chain after integration:
            # `find_neighbors(reduction = "${step.output.reduction}")`.
            summary$reduction <- active_reduction
        }
        # Observed state, not declared state: the object really has these. They
        # are skipped when absent so the key means "the object has one", not
        # "the object has one and it is NULL".
        active_ident <- tryCatch(sclet_get_active_ident(value), error = function(e) NULL)
        if (!is.null(active_ident)) summary$active_ident <- active_ident
        active_graph <- tryCatch(DefaultGraph(value), error = function(e) NULL)
        if (!is.null(active_graph)) summary$active_graph <- active_graph

        # The variable-feature set, as an observation rather than a declaration.
        # Written before the contract-key fill so the static marker in
        # `output_schema` ("hvg") is replaced by the real gene names and never
        # copied into the summary -- the fill in the loop below only adds keys
        # that are still absent. Only this action promises the field, so the
        # selector is a name check rather than a descriptor lookup: the summary
        # is built from the returned object, and `find_variable_features` is
        # the action that produces an HVG set on it.
        hvg_set <- if (isTRUE((descriptor %||% list())$name == "find_variable_features")) {
            tryCatch(VariableFeatures(value), error = function(e) NULL)
        } else {
            NULL
        }
        if (!is.null(hvg_set) && length(hvg_set)) summary$variable_features <- as.character(hvg_set)

        # The detest record id this step registered, as an observation rather
        # than a declaration -- same pattern as `variable_features` above. The
        # static marker in `output_schema` ("detest_record_id") is replaced by
        # the real id, read back from the returned object, so a downstream
        # `run_enrichment(detest_id = "${step.output.detest_id}")` resolves.
        # Only run_de_test promises the field, hence the name check.
        detest_ids <- if (isTRUE((descriptor %||% list())$name == "run_de_test")) {
            tryCatch(names(sclet_get_state_records(value, "detest")),
                error = function(e) character())
        } else {
            NULL
        }
        if (length(detest_ids)) summary$detest_id <- detest_ids[[1L]]

        # Chainable ids for the malignant workflow, observed from the returned
        # object (the step records them under their analysis maps).
        cnv_ids <- if (isTRUE((descriptor %||% list())$name == "run_cnv_inference")) {
            tryCatch(names(sclet_get_analysis(value, "malignant_cnv")),
                error = function(e) character())
        } else {
            NULL
        }
        if (length(cnv_ids)) summary$cnv_id <- cnv_ids[[length(cnv_ids)]]


        # Declared contract keys fill whatever the observations did not already
        # establish. An observation must never be overwritten by a declaration.
        contract <- (descriptor %||% list())$output_schema %||% list()
        for (name in sclet_ai_exposable_contract_keys(descriptor)) {
            if (name %in% names(summary)) {
                next
            }
            declared <- contract[[name]]
            if (is.list(declared)) {
                # Keep the names, not the nested structure. A contract value such
                # as list(clustering = "louvain_clusters") holds internal
                # state-type labels whose KEY matches the privacy classifier's
                # group pattern, so carrying it verbatim made the sanitizer
                # redact `required_states.clustering` as group labels -- a false
                # positive on a state name, but a real warning on a path that
                # leads to the ledger. The summary's job is to record what the
                # step promised, which the names do.
                declared <- utils::head(names(declared), 10L)
            }
            summary[[name]] <- declared
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
#'   `status`. Dry-run results additionally include `success_criteria_preview`,
#'   `stop_conditions_preview`, and `human_confirmations` (threaded through from
#'   `validation$human_confirmations`), all purely descriptive and unevaluated.
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
        current_plan_fingerprint <- sclet_ai_plan_fingerprint(plan)
        if (!identical(as.character(validation$plan_id), as.character(plan$plan_id)) ||
            !identical(as.character(validation$plan_fingerprint), current_plan_fingerprint) ||
            (!is.null(validation$context_fingerprint) &&
                !identical(as.character(validation$context_fingerprint), as.character(current_fingerprint)))) {
            sclet_ai_error(
                "sclet_ai_stale_plan",
                "The supplied plan validation does not match the current object or plan. Revalidate before execution."
            )
        }
    }
    # Backstop to the plan_empty validation rule, for a caller that bypasses
    # ValidateAIPlan() with a hand-built validation object: executing zero
    # steps must fail loudly rather than report status="completed" with no
    # results, which reads as a finished analysis. Measured 2026-10-16; see
    # the plan_empty comment in ValidateAIPlan().
    if (!length(plan$actions %||% list())) {
        sclet_ai_error(
            "sclet_ai_empty_plan",
            "The plan has no actions; nothing to execute. Re-propose the plan."
        )
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
                success_assessment = list(),
                success_criteria_preview = plan$success_criteria %||% list(),
                stop_conditions_preview = plan$stop_conditions %||% list(),
                human_confirmations = validation$human_confirmations %||% list()
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
