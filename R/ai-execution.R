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

# Does a named assay hold raw counts, as inferCNV requires?
#
# `run_cnv_inference` defaults to the `counts` assay, but the model may name
# another one. The handler refuses a transformed assay with a clear message;
# validation did not, so the refusal surfaced only at execution, after a plan
# had been paid for (measured 2026-10-19: the first real-data malignant
# benchmark planned `assay = "logcounts"`, validated, and failed the step).
# Returning a typed problem here gives the repair loop the same sentence
# before anything runs.
#
# Deliberately a sample, not a full scan: a sparse assay's nonzeros are exact,
# and proving the negative case does not require reading a dense matrix twice.
# The handler keeps its full check as the authority.
sclet_ai_assay_count_problem <- function(object, assay) {
    if (!assay %in% SummarizedExperiment::assayNames(object)) return(NULL)
    x <- SummarizedExperiment::assay(object, assay)
    values <- if (methods::is(x, "sparseMatrix")) {
        x@x
    } else {
        rows <- seq_len(min(nrow(x), 100L))
        cols <- seq_len(min(ncol(x), 100L))
        as.numeric(as.matrix(x[rows, cols, drop = FALSE]))
    }
    values <- values[!is.na(values)]
    if (!length(values)) return(NULL)
    if (any(values < 0) || any(values != round(values))) {
        return(paste0(
            "counts_assay_not_count: assay '", assay, "' does not hold raw ",
            "counts (it has negative or non-integer values); inferCNV needs ",
            "raw counts. Use assay = \"counts\" (or the name of the object's ",
            "count assay)."))
    }
    NULL
}


# Is this selector value a promise that can only be resolved once the
# producing step has run, rather than a value that can be checked now?
#
# `variable_features` is declared in `find_variable_features`' output_schema as
# the marker "hvg": the static form must not carry 2000 gene names into every
# planning context, so the real set only exists after the step executes. A
# downstream step therefore writes
# `features = "$\{step_2.output.variable_features\}"`, and
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
# an ordinary selector. Only the exact `$\{step.output.field\}` shape counts: a
# bare `$\{step\}` or `$\{step.output\}` is a malformed reference that must keep its
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
#' Report which action parameters do not match the registry schema
#'
#' Compares a plan step's parameters against the input schema of the action it
#' names and returns a human-readable description of what does not line up:
#' unknown parameters, missing required ones, and values of the wrong type.
#'
#' @param params The `params` list of a plan step.
#' @param schema The action's `input_schema`, as held by the registry.
#' @return A character string describing the problems, empty when the
#'   parameters match the schema.
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
    # The catalog is assembled one block per registry position instead of as a
    # single 2,800-line literal. Each block keeps its own `%in% include` guard
    # and is appended in the order below, so group selection and action order
    # are exactly what they were before the split.
    actions <- list()
    actions <- c(actions, sclet_ai_registry_read(object, include))
    actions <- c(actions, sclet_ai_registry_preprocess(object, include))
    actions <- c(actions, sclet_ai_registry_dimred_graph(object, include))
    actions <- c(actions, sclet_ai_registry_cluster_design(object, include))
    actions <- c(actions, sclet_ai_registry_integration(object, include))
    actions <- c(actions, sclet_ai_registry_annotation(object, include))
    actions <- c(actions, sclet_ai_registry_rare_cell(object, include))
    actions <- c(actions, sclet_ai_registry_trajectory_velocity_fate(object, include))
    actions <- c(actions, sclet_ai_registry_spatial_multimodal(object, include))
    actions <- c(actions, sclet_ai_registry_malignant(object, include))
    actions <- c(actions, sclet_ai_registry_enrichment_communication_milo(object, include))
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
    sclet_ai_record_evidence_for_run(object, evidence, source = name,
        parents = character(), scope = NULL)
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
    sclet_ai_record_evidence_for_run(object, evidence, source = name,
        parents = character(), scope = NULL)
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
#' a plan containing `$\{step_2.output.required_hvg\}` *passed validation and
#' then failed at execution* with "Plan output field is not available":
#'
#'   step_3: Plan output field is not available: $\{step_2.output.required_hvg\}
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
            # `find_neighbors(reduction = "$\{step.output.reduction\}")`.
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
        # `run_enrichment(detest_id = "$\{step.output.detest_id\}")` resolves.
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

# Bound one resolved step's parameters for the ledger.
#
# The execution record is stored ON the object, so a parameter's value cannot
# be copied wholesale (a 2000-gene `features` vector is not a ledger entry),
# yet the provenance skill promises the resolved parameters are readable so a
# reader can re-run the same analysis. Scalars are kept verbatim; longer
# vectors keep their type and length only -- enough to see that `dims = 1:30`
# or `features = <2000 names>` was what ran, without copying either, and
# without putting free-form values of unknown origin into the record.
sclet_ai_bounded_step_params <- function(params) {
    if (is.null(params) || !length(params)) {
        return(list())
    }
    bounded <- lapply(params, function(value) {
        if (is.null(value)) {
            return(NULL)
        }
        if (is.list(value)) {
            return(sclet_ai_bounded_step_params(value))
        }
        if (is.atomic(value) && length(value) == 1L) {
            return(value)
        }
        if (is.atomic(value)) {
            return(list(type = typeof(value), length = length(value)))
        }
        list(type = class(value)[[1L]])
    })
    bounded <- bounded[!vapply(bounded, is.null, logical(1L))]
    bounded
}

sclet_ai_record_execution <- function(object, plan, results, status, dry_run) {
    execution_id <- sclet_ai_execution_id(plan$plan_id)
    # Per-step parameters and timing, not just counts. `results` carries the
    # duration already; the resolved parameters were computed inside the
    # execution loop and used to be discarded, so the ledger could say a step
    # completed but not what it was told to do (2026-10-19: neither
    # `nfeatures` nor `logfc.threshold` was recoverable from any record).
    steps <- lapply(results, function(r) {
        list(
            id = r$id,
            action = r$action,
            status = r$status,
            duration_sec = r$duration_sec %||% NA_real_,
            params = r$params %||% list()
        )
    })
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
            requires_confirmation = isTRUE(plan$requires_confirmation),
            steps = steps,
            # B20: the unexplored suggestions the planner SAW, snapshotted at
            # plan time. The report has to be able to say "this run was told
            # about X and did not do it" -- without the snapshot that sentence
            # is unwritable, because after execution the ledger already contains
            # this round's evidence and "was it unexplored before" is no longer
            # answerable. Empty when the context carried no data_support block.
            unexplored_suggestions = (plan$metadata %||% list())$
                unexplored_suggestions %||% list()
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
    # An unreadable ledger must not be reported as "no evidence node found":
    # that sentence blames the analysis for a failure of the read, and a caller
    # chasing it will look for a step that never ran. Keep the difference.
    evidence_read_error <- NULL
    after_evidence <- tryCatch(
        sclet_ai_evidence_get_all(object_after),
        error = function(e) {
            evidence_read_error <<- conditionMessage(e)
            list()
        }
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
                    reason = if (!is.null(evidence_read_error)) {
                        paste0("could not read the evidence ledger (", evidence_read_error,
                            "), so this criterion was not evaluated")
                    } else {
                        paste("no evidence node found for", crit$evidence_id)
                    },
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
        step_params <- sclet_ai_bounded_step_params(resolved_params)
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
                duration_sec = 0,
                params = step_params
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
                duration_sec = elapsed,
                params = step_params
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
                    duration_sec = elapsed,
                    params = step_params
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
                    duration_sec = elapsed,
                    params = step_params
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
            duration_sec = elapsed,
            params = step_params
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
