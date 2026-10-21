sclet_ai_require <- function() {
    if (!requireNamespace("aisdk", quietly = TRUE)) {
        sclet_ai_error(
            "sclet_ai_missing_dependency",
            "The optional package 'aisdk' is required for AI calls. Install it before using sclet AI functions."
        )
    }
    invisible(TRUE)
}

sclet_ai_resolve_model <- function(model = NULL) {
    if (!is.null(model)) {
        return(model)
    }
    sclet_ai_require()
    model <- tryCatch({
        ns <- asNamespace("aisdk")
        getter <- if (exists("get_model", envir = ns, inherits = FALSE)) {
            get("get_model", envir = ns)
        } else if (exists("model", envir = ns, inherits = FALSE)) {
            get("model", envir = ns)
        } else {
            NULL
        }
        if (is.function(getter)) getter() else NULL
    }, error = function(e) NULL)
    if (!is.null(model)) {
        return(model)
    }
    model_name <- Sys.getenv("OPENAI_MODEL", unset = "")
    if (!nzchar(model_name)) {
        sclet_ai_error(
            "sclet_ai_missing_model",
            "No AI model was supplied and aisdk has no configured default model. Set an aisdk default model or OPENAI_MODEL."
        )
    }
    # Current aisdk accepts provider:model identifiers directly. Older aisdk
    # releases may expose only provider factories, so retain that fallback.
    if (grepl(":", model_name, fixed = TRUE)) {
        return(model_name)
    }
    tryCatch({
        provider <- aisdk::create_openai()
        provider$language_model(model_name)
    }, error = function(e) sclet_ai_error(
        "sclet_ai_provider_error",
        paste("Could not create AI model:", conditionMessage(e)),
        parent = e
    ))
}

sclet_ai_tool_description <- function(tools) {
    if (!is.list(tools) || !length(tools)) {
        return("")
    }
    descriptions <- vapply(tools, function(tool) {
        if (is.function(tool)) {
            return("- unnamed R tool")
        }
        if (is.list(tool)) {
            name <- tool$name %||% "unnamed"
            description <- tool$description %||% "No description"
            return(paste0("- ", name, ": ", description))
        }
        "- unnamed tool"
    }, character(1))
    paste(c("Available read-only tools:", descriptions), collapse = "\n")
}

sclet_ai_native_parameter <- function(spec) {
    sclet_ai_require()
    if (is.list(spec) && !is.null(spec$schema)) {
        return(spec$schema)
    }
    if (is.character(spec) && length(spec) > 1L &&
        exists("z_enum", asNamespace("aisdk"), inherits = FALSE)) {
        return(aisdk::z_enum(spec))
    }
    if (identical(spec, "logical") || identical(spec, "boolean")) {
        return(aisdk::z_boolean())
    }
    if (identical(spec, "numeric") || identical(spec, "number")) {
        return(aisdk::z_number())
    }
    if (identical(spec, "integer")) {
        return(aisdk::z_number())
    }
    if (identical(spec, "list") || identical(spec, "object")) {
        return(aisdk::z_object(.additional_properties = TRUE))
    }
    aisdk::z_string(nullable = TRUE)
}

# Convert sclet's provider-neutral descriptors to aisdk native Tool objects.
# This is deliberately kept internal so the stable sclet registry does not
# depend on aisdk R6 internals.
sclet_ai_native_tools <- function(tools) {
    sclet_ai_require()
    if (!is.list(tools) || !length(tools) ||
        !exists("tool", asNamespace("aisdk"), inherits = FALSE)) {
        return(list())
    }
    lapply(tools, function(spec) {
        if (!is.list(spec) || !is.function(spec$handler)) {
            sclet_ai_error("sclet_ai_tool_error", "Invalid sclet AI tool descriptor.")
        }
        parameters <- spec$input_schema %||% list()
        schema <- if (length(parameters) &&
            exists("z_object", asNamespace("aisdk"), inherits = FALSE)) {
            fields <- lapply(parameters, sclet_ai_native_parameter)
            do.call(aisdk::z_object, fields)
        } else {
            NULL
        }
        tryCatch(
            aisdk::tool(
                name = spec$name,
                description = spec$description,
                parameters = schema,
                execute = spec$handler
            ),
            error = function(e) sclet_ai_error(
                "sclet_ai_tool_error",
                paste("Could not create native aisdk tool:", conditionMessage(e)),
                parent = e
            )
        )
    })
}

# Claim levels a finding may carry. The single source of truth: the
# structured-output schema below and the local validator in ai-result.R both
# read this list, so they cannot drift apart.
sclet_ai_finding_claim_levels <- function() {
    c(
        "observed", "measured", "estimated", "associated",
        "consistent_with", "suggestive", "hypothesis", "causal"
    )
}

# Severities a finding may carry, likewise shared with the local validator.
sclet_ai_finding_severities <- function() {
    c("info", "low", "warning", "high", "error")
}

# Where a plan-level criterion may derive its value from.
sclet_ai_criterion_sources <- function() {
    c("action_output", "evidence")
}

# The comparisons a criterion may ask for.
sclet_ai_criterion_checks <- function() {
    c("exists", "equals", "gte", "lte", "in")
}

sclet_ai_evidence_citation_instruction <- function() {
    paste(
        "When the supplied ledger records ai_evidence nodes (records with",
        "type 'ai_evidence', ids such as 'ev:...'), every finding that rests on",
        "one of them must list that node's id in the finding's evidence_refs",
        "field, verbatim. Cite only ids that appear in the supplied ledger;",
        "never invent an evidence id.",
        sep = " "
    )
}

sclet_ai_criterion_instruction <- function() {
    paste(
        "When a proposed step produces something that can be verified after",
        "execution -- a reduction, a graph, a clustering with named groups, a",
        "marker table, an evidence node -- include a success_criterion naming",
        "it, and include a stop_condition when a wrong result should halt the",
        "plan. Use source = \"action_output\" with the step's id for anything",
        "this plan produces, or source = \"evidence\" with an evidence_id it is",
        "expected to write. Use check = \"exists\" when presence is enough;",
        "\"equals\", \"gte\", \"lte\" and \"in\" additionally require both field",
        "and value. Never require a criterion: if the plan produces nothing",
        "verifiable, omit both lists rather than inventing criteria.",
        sep = " "
    )
}

# A reusable native structured-output schema for future callers. Tool/JSON mode
# is selected by aisdk itself; the local result validator remains authoritative.
sclet_ai_native_result_schema <- function() {
    sclet_ai_require()
    if (!exists("z_object", asNamespace("aisdk"), inherits = FALSE)) {
        return(NULL)
    }
    do.call(
        aisdk::z_object,
        list(
            answer = aisdk::z_string(nullable = TRUE),
            findings = aisdk::z_array(aisdk::z_object(
                id = aisdk::z_string(nullable = TRUE),
                statement = aisdk::z_string(nullable = TRUE),
                description = aisdk::z_string(nullable = TRUE),
                severity = aisdk::z_enum(sclet_ai_finding_severities(), nullable = TRUE),
                claim_level = aisdk::z_enum(sclet_ai_finding_claim_levels(), nullable = TRUE),
                evidence_refs = aisdk::z_array(aisdk::z_string(), nullable = TRUE),
                .additional_properties = TRUE
            )),
            evidence = aisdk::z_array(aisdk::z_any()),
            warnings = aisdk::z_array(aisdk::z_any()),
            recommendations = aisdk::z_array(aisdk::z_any()),
            proposed_actions = aisdk::z_array(aisdk::z_any()),
            success_criteria = aisdk::z_array(sclet_ai_native_criterion_schema()),
            stop_conditions = aisdk::z_array(sclet_ai_native_criterion_schema())
        )
    )
}

# The declared shape of one plan-level criterion (used by both
# `success_criteria` and `stop_conditions`, which normalize identically).
sclet_ai_native_criterion_schema <- function() {
    aisdk::z_object(
        id = aisdk::z_string(nullable = TRUE),
        description = aisdk::z_string(nullable = TRUE),
        source = aisdk::z_enum(sclet_ai_criterion_sources(), nullable = TRUE),
        step_id = aisdk::z_string(nullable = TRUE),
        evidence_id = aisdk::z_string(nullable = TRUE),
        check = aisdk::z_enum(sclet_ai_criterion_checks(), nullable = TRUE),
        field = aisdk::z_string(nullable = TRUE),
        value = aisdk::z_any(),
        .additional_properties = TRUE
    )
}

sclet_ai_schema_error <- function(path, expected, actual) {
    sclet_ai_error(
        "sclet_ai_invalid_output",
        paste0("AI output schema mismatch at `", path, "`: expected ", expected,
            ", got ", actual, ".")
    )
}

sclet_ai_schema_type_name <- function(value) {
    if (is.null(value)) return("null")
    if (is.factor(value)) return("string")
    if (is.logical(value)) return("boolean")
    if (is.integer(value)) return("integer")
    if (is.numeric(value)) return("number")
    if (is.character(value)) return("string")
    if (is.list(value)) return("object")
    "unsupported value"
}

sclet_ai_schema_types <- function(schema) {
    type <- schema$type %||% NULL
    if (is.null(type)) return(character())
    if (is.character(type)) return(type)
    if (is.list(type) && all(vapply(type, is.character, logical(1)))) {
        return(unlist(type, use.names = FALSE))
    }
    character()
}

sclet_ai_validate_schema_node <- function(value, schema, path = "$") {
    if (is.null(schema) || !is.list(schema)) return(invisible(TRUE))
    types <- sclet_ai_schema_types(schema)
    actual <- sclet_ai_schema_type_name(value)
    non_null_types <- setdiff(types, "null")
    if (identical(actual, "null")) {
        if (!"null" %in% types && length(types)) {
            sclet_ai_schema_error(path, paste(types, collapse = " or "), actual)
        }
        return(invisible(TRUE))
    }
    if (length(non_null_types) && !"any" %in% non_null_types) {
        compatible <- vapply(non_null_types, function(type) {
            switch(type,
                object = identical(actual, "object"),
                array = is.list(value) && is.null(names(value)),
                string = identical(actual, "string"),
                boolean = identical(actual, "boolean"),
                number = is.numeric(value) && length(value) == 1L && is.finite(value),
                integer = is.numeric(value) && length(value) == 1L && is.finite(value) && value == floor(value),
                FALSE
            )
        }, logical(1))
        if (!any(compatible)) {
            sclet_ai_schema_error(path, paste(non_null_types, collapse = " or "), actual)
        }
    }

    if (!is.null(schema$enum)) {
        choices <- unlist(schema$enum, recursive = TRUE, use.names = FALSE)
        if (!any(vapply(choices, function(choice) {
            identical(value, choice) ||
                (is.numeric(value) && length(value) == 1L && is.numeric(choice) &&
                    length(choice) == 1L && isTRUE(value == choice))
        }, logical(1)))) {
            sclet_ai_schema_error(path, paste0("one of {", paste(choices, collapse = ", "), "}"),
                paste0("value `", paste(as.character(value), collapse = ","), "`"))
        }
    }
    if ("object" %in% non_null_types && !is.null(schema$properties)) {
        if (is.null(names(value)) && length(value)) {
            sclet_ai_schema_error(path, "a named object", "an unnamed list")
        }
        properties <- schema$properties
        required <- schema$required %||% list()
        required <- unlist(required, recursive = TRUE, use.names = FALSE)
        # aisdk's z_object marks every property required, including fields
        # declared `nullable = TRUE`. The local result contract treats those
        # as omittable ("a criterion that needs no value still validates"),
        # so a nullable field may be absent locally; when present it may be
        # null or must match its type.
        nullable_fields <- names(properties)[vapply(properties, function(prop) {
            "null" %in% sclet_ai_schema_types(prop)
        }, logical(1))]
        required <- setdiff(required, nullable_fields)
        missing <- setdiff(required, names(value))
        if (length(missing)) {
            sclet_ai_error("sclet_ai_invalid_output",
                paste0("AI output is missing required fields at `", path, "`: ", paste(missing, collapse = ", "), "."))
        }
        extra <- setdiff(names(value), names(properties))
        if (length(extra) && identical(schema$additionalProperties, FALSE)) {
            sclet_ai_error("sclet_ai_invalid_output",
                paste0("AI output has unsupported fields at `", path, "`: ", paste(extra, collapse = ", "), "."))
        }
        for (field in intersect(names(value), names(properties))) {
            sclet_ai_validate_schema_node(value[[field]], properties[[field]], paste0(path, ".", field))
        }
    }
    if ("array" %in% non_null_types && !is.null(schema$items)) {
        for (i in seq_along(value)) {
            sclet_ai_validate_schema_node(value[[i]], schema$items, paste0(path, "[", i, "]"))
        }
    }
    invisible(TRUE)
}

sclet_ai_validate_schema <- function(result, schema) {
    if (is.null(schema)) {
        return(result)
    }
    if (is.function(schema)) {
        checked <- tryCatch(schema(result), error = function(e) e)
        if (inherits(checked, "error")) {
            sclet_ai_error(
                "sclet_ai_invalid_output",
                paste("AI output failed schema validation:", conditionMessage(checked)),
                parent = checked
            )
        }
        if (isFALSE(checked)) {
            sclet_ai_error("sclet_ai_invalid_output", "AI output failed schema validation.")
        }
        return(result)
    }
    if (inherits(schema, "z_schema")) {
        # A normalized sclet_ai_result always carries the wrapper fields
        # task/context/metadata on top of the structured payload the schema
        # describes. Validate the schema-relevant payload only, so the local
        # check does not reject fields the schema never meant to constrain.
        value <- result
        if (inherits(result, "sclet_ai_result")) {
            declared <- names(schema$properties %||% list())
            undeclared_wrapper <- setdiff(c("task", "context", "metadata"), declared)
            value <- result[setdiff(names(result), undeclared_wrapper)]
        }
        sclet_ai_validate_schema_node(value, schema)
        return(result)
    }
    if (is.list(schema) && !is.null(schema$required)) {
        missing <- setdiff(as.character(schema$required), names(result))
        if (length(missing)) {
            sclet_ai_error(
                "sclet_ai_invalid_output",
                paste("AI output is missing required fields:", paste(missing, collapse = ", "))
            )
        }
    }
    result
}

sclet_ai_create_agent <- function(...) {
    factory <- getOption("sclet.ai.create_agent", NULL)
    if (is.function(factory)) {
        return(do.call(factory, list(...)))
    }
    aisdk::create_agent(...)
}

sclet_ai_generate_object <- function(...) {
    factory <- getOption("sclet.ai.generate_object", NULL)
    if (is.function(factory)) {
        return(do.call(factory, list(...)))
    }
    aisdk::generate_object(...)
}

#' Call an aisdk-backed sclet AI task
#'
#' @param task Task identifier.
#' @param context Structured context returned by `GetAnalysisLedger()`.
#' @param tools Optional sclet read-only tool descriptors.
#' @param schema Optional local validator or aisdk structured-output schema.
#' @param model Optional aisdk model identifier or model object.
#' @param system_prompt Optional additional system instructions.
#' @param structured_output Logical. Ask the provider for a structured
#'   response instead of free text.
#' @param fallback_on_structure_error Logical. Fall back to a text answer
#'   when the provider's structured output does not validate.
#' @param privacy Privacy level for the outbound payload: `"standard"`,
#'   `"strict"` or `"off"`.
#' @param privacy_consent Logical. Record that the user consented to the
#'   privacy level used for this call.
#' @param enforce_privacy Logical. Refuse a call whose payload the privacy
#'   layer cannot bound.
#' @param ... Additional adapter options, including `max_steps`.
#' @return A `sclet_ai_result` object.
sclet_ai_call <- function(
    task,
    context,
    tools = list(),
    schema = NULL,
    model = NULL,
    system_prompt = NULL,
    structured_output = FALSE,
    fallback_on_structure_error = TRUE,
    privacy = getOption("sclet.ai.privacy", "standard"),
    privacy_consent = FALSE,
    enforce_privacy = getOption("sclet.ai.enforce_privacy", TRUE),
    ...
) {
    if (!is.list(context)) {
        sclet_ai_error("sclet_ai_context_error", "`context` must be a structured list.")
    }
    dots <- list(...)
    mock <- getOption("sclet.ai.call", default = NULL)
    if (is.null(mock)) mock <- getOption("sclet.ai.adapter", NULL)
    privacy_policy <- NULL
    if (!length(enforce_privacy) || !is.logical(enforce_privacy) || is.na(enforce_privacy)) {
        enforce_privacy <- TRUE
    }
    if (isTRUE(enforce_privacy)) {
        sanitized <- sclet_ai_sanitize_context(
            object = NULL,
            context = context,
            privacy = privacy,
            user_consent = privacy_consent
        )
        context <- sanitized$payload
        privacy_policy <- sanitized$policy
        if (!is.null(privacy_policy)) {
            if (!is.function(mock)) {
                if (isTRUE(privacy_policy$requires_user_consent)) {
                    warning(
                        call. = FALSE,
                        paste(
                            "sclet AI outbound payload requires user consent before use.",
                            "Fields requiring consent were redacted:",
                            paste(sort(unique(privacy_policy$redacted_fields)), collapse = ", "),
                            "Re-run with privacy_consent = TRUE to include them under explicit user approval."
                        )
                    )
                }
                if (!is.null(privacy_policy$warnings) && length(privacy_policy$warnings)) {
                    for (w in unique(privacy_policy$warnings)) {
                        warning(call. = FALSE, paste("sclet AI payload sanitizer:", w))
                    }
                }
            }
        }
    }
    if (is.function(mock)) {
        response <- tryCatch(
            mock(task = task, context = context, tools = tools, schema = schema, model = model, ...),
            error = function(e) sclet_ai_error(
                "sclet_ai_provider_error",
                paste("Mock AI adapter failed:", conditionMessage(e)),
                parent = e
            )
        )
        result <- sclet_ai_normalize_response(
            response, task = task, context = context,
            metadata = list(
                mock = TRUE,
                native_tools = FALSE,
                payload_policy = privacy_policy
            )
        )
        return(sclet_ai_validate_schema(result, schema))
    }

    model <- sclet_ai_resolve_model(model)
    sclet_ai_require()
    prompt <- paste(
        "You are sclet's R-native single-cell analysis assistant.",
        "Use only the supplied ledger facts and clearly label uncertainty.",
        "Do not claim causality unless the evidence explicitly supports it.",
        if (!is.null(system_prompt)) system_prompt else "",
        sclet_ai_evidence_citation_instruction(),
        sclet_ai_criterion_instruction(),
        "Structured analysis ledger:",
        paste(utils::capture.output(dput(context)), collapse = "\n"),
        sep = "\n\n"
    )
    structured_output <- isTRUE(structured_output)
    if (structured_output && !length(tools) &&
        exists("generate_object", asNamespace("aisdk"), inherits = FALSE)) {
        native_schema <- if (inherits(schema, "z_schema")) {
            schema
        } else {
            sclet_ai_native_result_schema()
        }
        started <- Sys.time()
        response <- tryCatch(
            sclet_ai_generate_object(
                model = model,
                prompt = task,
                schema = native_schema,
                schema_name = task,
                system = prompt,
                mode = "tool",
                max_retries = getOption("sclet.ai.max_retries", 1L)
            ),
            error = function(e) sclet_ai_error(
                "sclet_ai_provider_error",
                paste("aisdk structured call failed:", conditionMessage(e)),
                parent = e
            )
        )
        if (!is.list(response) || is.null(response$object)) {
            if (inherits(response, "GenerateObjectResult")) {
                msg <- sprintf(
                    "aisdk structured call produced no object (valid=%s, attempts=%s, finish_reason=%s).",
                    if (is.logical(response$valid)) response$valid else "NA",
                    if (!is.null(response$attempts)) response$attempts else "NA",
                    if (!is.null(response$finish_reason)) response$finish_reason else "unknown"
                )
                rt <- response$raw_text %||% ""
                if (is.character(rt) && nzchar(rt)) {
                    msg <- paste0(
                        msg,
                        "\n\nraw_text (truncated to 800 chars):\n",
                        substr(rt, 1, 800)
                    )
                }
                sclet_ai_error("sclet_ai_provider_error", msg)
            }
            sclet_ai_error("sclet_ai_provider_error",
                "aisdk structured call produced no object.")
        }
        raw <- response$object
        result <- sclet_ai_normalize_response(
            raw,
            task = task,
            context = context,
            metadata = list(
                model = if (is.character(model)) model else NULL,
                duration_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
                provider = "aisdk",
                native_tools = FALSE,
                structured_output = TRUE,
                payload_policy = privacy_policy
            )
        )
        return(sclet_ai_validate_schema(result, schema))
    }
    native_tools <- sclet_ai_native_tools(tools)
    max_steps <- dots$max_steps %||% getOption("sclet.ai.max_steps", 6L)
    agent <- tryCatch(
        sclet_ai_create_agent(
            name = "sclet_copilot",
            description = "A single-cell analysis auditor, planner, and interpreter",
            system_prompt = prompt,
            tools = native_tools,
            model = model
        ),
        error = function(e) sclet_ai_error(
            "sclet_ai_provider_error",
            paste("Could not create aisdk agent:", conditionMessage(e)),
            parent = e
        )
    )
    started <- Sys.time()
    response <- tryCatch(
        agent$run(task, max_steps = max_steps),
        error = function(e) sclet_ai_error(
            "sclet_ai_provider_error",
            paste("aisdk AI call failed:", conditionMessage(e)),
            parent = e
        )
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    text <- if (is.list(response) && !is.null(response$text)) response$text else response
    if (!is.character(text) && !is.list(text)) {
        extracted <- tryCatch(response$text, error = function(e) NULL)
        if (!is.null(extracted)) text <- extracted
    }
    if (!is.character(text) && !is.list(text)) text <- as.character(text)
    result <- sclet_ai_normalize_response(
        text,
        task = task,
        context = context,
        metadata = list(
            model = if (is.character(model)) model else NULL,
            duration_sec = elapsed,
            provider = "aisdk",
            native_tools = length(native_tools) > 0L,
            max_steps = max_steps,
            payload_policy = privacy_policy
        )
    )
    if (!isTRUE(structured_output)) {
        return(sclet_ai_validate_schema(result, schema))
    }

    native_schema <- if (inherits(schema, "z_schema")) {
        schema
    } else {
        sclet_ai_native_result_schema()
    }
    structured_started <- Sys.time()
    structured_response <- tryCatch(
        sclet_ai_generate_object(
            model = model,
            prompt = paste(
                "Convert the following sclet analysis assistant response into the",
                "requested structured result. Preserve evidence and uncertainty.",
                "Do not add facts that are absent from the response or ledger.",
                "Assistant response:",
                text,
                sep = "\n\n"
            ),
            schema = native_schema,
            schema_name = task,
            system = paste(
                "Return only an object matching the supplied schema.",
                "Use claim_level values conservatively; causal claims require evidence_refs.",
                sclet_ai_evidence_citation_instruction(),
                sclet_ai_criterion_instruction(),
                sep = "\n"
            ),
            mode = "tool",
            max_retries = getOption("sclet.ai.max_retries", 1L)
        ),
        error = function(e) e
    )
    if (inherits(structured_response, "error")) {
        if (!isTRUE(fallback_on_structure_error)) {
            sclet_ai_error(
                "sclet_ai_provider_error",
                paste("aisdk structured summary failed:", conditionMessage(structured_response)),
                parent = structured_response
            )
        }
        result$warnings <- c(
            result$warnings,
            list(list(
                code = "structured_output_fallback",
                message = conditionMessage(structured_response)
            ))
        )
        result$metadata$structured_output <- FALSE
        result$metadata$structured_output_error <- conditionMessage(structured_response)
        return(sclet_ai_validate_schema(result, schema))
    }
    if (is.list(structured_response) && !is.null(structured_response$object)) {
        structured_object <- structured_response$object
        structured_result <- sclet_ai_normalize_response(
            structured_object,
            task = task,
            context = context,
            metadata = list(
                model = if (is.character(model)) model else NULL,
                duration_sec = elapsed + as.numeric(difftime(Sys.time(), structured_started, units = "secs")),
                provider = "aisdk",
                native_tools = length(native_tools) > 0L,
                structured_output = TRUE,
                tool_loop_steps = max_steps,
                payload_policy = privacy_policy
            )
        )
        return(sclet_ai_validate_schema(structured_result, schema))
    }
    if (inherits(structured_response, "GenerateObjectResult")) {
        msg <- sprintf(
            "aisdk structured summary produced no object (valid=%s, attempts=%s, finish_reason=%s).",
            if (is.logical(structured_response$valid)) structured_response$valid else "NA",
            if (!is.null(structured_response$attempts)) structured_response$attempts else "NA",
            if (!is.null(structured_response$finish_reason)) structured_response$finish_reason else "unknown"
        )
        rt <- structured_response$raw_text %||% ""
        if (is.character(rt) && nzchar(rt)) {
            msg <- paste0(
                msg,
                "\n\nraw_text (truncated to 800 chars):\n",
                substr(rt, 1, 800)
            )
        }
        if (!isTRUE(fallback_on_structure_error)) {
            sclet_ai_error("sclet_ai_provider_error", msg, parent = structured_response)
        }
        result$warnings <- c(
            result$warnings,
            list(list(
                code = "structured_output_fallback_no_object",
                message = msg
            ))
        )
        result$metadata$structured_output <- FALSE
        result$metadata$structured_output_error <- msg
        return(sclet_ai_validate_schema(result, schema))
    }
    structured_result <- sclet_ai_normalize_response(
        structured_response,
        task = task,
        context = context,
        metadata = list(
            model = if (is.character(model)) model else NULL,
            duration_sec = elapsed + as.numeric(difftime(Sys.time(), structured_started, units = "secs")),
            provider = "aisdk",
            native_tools = length(native_tools) > 0L,
            structured_output = TRUE,
            tool_loop_steps = max_steps,
            payload_policy = privacy_policy
        )
    )
    sclet_ai_validate_schema(structured_result, schema)
}
