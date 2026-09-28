sclet_ai_tool <- function(name, description, handler, input_schema = list()) {
    list(
        name = name,
        description = description,
        input_schema = input_schema,
        output_schema = list(type = "object"),
        read_only = TRUE,
        handler = handler
    )
}

#' Build the read-only tool registry for an SCE
#'
#' @param object A `SingleCellExperiment` object.
#' @param mode Registry mode. The first implementation exposes only read-only
#'   tools; `"execute"` returns an empty registry until execution safeguards are
#'   implemented.
#' @return A named list of tool descriptors.
sclet_ai_tool_registry <- function(object, mode = c("read", "execute")) {
    mode <- match.arg(mode)
    if (identical(mode, "execute")) {
        return(list())
    }
    list(
        get_status = sclet_ai_tool(
            "get_status",
            "Return deterministic Status() and health information.",
            function() Status(object)
        ),
        get_ledger = sclet_ai_tool(
            "get_ledger",
            "Return the bounded AI-facing analysis ledger.",
            function(detail = "summary", target = NULL) {
                GetAnalysisLedger(object, detail = detail, target = target)
            },
            list(detail = c("summary", "full"), target = "character")
        ),
        get_capabilities = sclet_ai_tool(
            "get_capabilities",
            "Return available assays, reductions, analyses and read-only capabilities.",
            function() GetAnalysisLedger(object)$capabilities
        ),
        get_analysis_record = sclet_ai_tool(
            "get_analysis_record",
            "Return one analysis record by type or id without exposing full matrices.",
            function(type = NULL, id = NULL) {
                ledger <- GetAnalysisLedger(object, detail = "full")
                records <- ledger$analyses
                if (!length(records)) return(NULL)
                hits <- records[vapply(records, function(x) {
                    identical(x$type, type) || identical(x$id, id) ||
                        (!is.null(type) && grepl(type, x$type %||% "", fixed = TRUE))
                }, logical(1))]
                if (!length(hits)) NULL else if (length(hits) == 1L) hits[[1L]] else hits
            },
            list(type = "character", id = "character")
        ),
        get_analysis_lineage = sclet_ai_tool(
            "get_analysis_lineage",
            "Return conservative parent references for recorded analyses.",
            function(target = NULL) {
                lineage <- GetAnalysisLedger(object)$lineage
                if (is.null(target)) return(lineage)
                lineage[grepl(target, names(lineage), fixed = TRUE)]
            },
            list(target = "character")
        ),
        get_quality_checks = sclet_ai_tool(
            "get_quality_checks",
            "Return deterministic health and missing-prerequisite checks.",
            function() GetAnalysisLedger(object)$quality_checks
        )
    )
}

sclet_ai_call_tools <- function(tools, name, args = list()) {
    if (!is.list(tools) || is.null(tools[[name]])) {
        sclet_ai_error("sclet_ai_tool_error", paste0("Unknown AI tool: ", name))
    }
    tool <- tools[[name]]
    if (!is.function(tool$handler)) {
        sclet_ai_error("sclet_ai_tool_error", paste0("AI tool has no handler: ", name))
    }
    do.call(tool$handler, args)
}

#' Ask the AI to summarize current sclet status
#'
#' @param object A `SingleCellExperiment` object.
#' @param model Optional aisdk model.
#' @param structured_output Logical. Use the tool-loop plus structured-output
#'   summary stage. Defaults to `TRUE`.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return A structured `sclet_ai_result`.
#' @export
AIStatus <- function(object, model = NULL, structured_output = TRUE, ...) {
    context <- sclet_ai_context(object, target = "status", detail = "summary")
    sclet_ai_call(
        task = "status_review",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Summarize the current analysis status for a scientist.",
            "Separate completed capabilities, missing prerequisites, and warnings.",
            "Do not invent biological conclusions."
        ),
        ...
    )
}

#' Ask the AI to review deterministic QC and readiness information
#'
#' @param object A `SingleCellExperiment` object.
#' @param model Optional aisdk model.
#' @param structured_output Logical. Use the tool-loop plus structured-output
#'   summary stage. Defaults to `TRUE`.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return A structured `sclet_ai_result`.
#' @export
AIReviewQC <- function(object, model = NULL, structured_output = TRUE, ...) {
    context <- sclet_ai_context(object, target = "qc", detail = "summary")
    sclet_ai_call(
        task = "qc_review",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Review only the supplied deterministic QC and readiness facts.",
            "Rank risks, identify missing checks, and propose conservative next steps.",
            "Do not claim that a biological population is real from QC alone."
        ),
        ...
    )
}

#' Ask the AI to explain a recorded analysis
#'
#' @param object A `SingleCellExperiment` object.
#' @param type Optional analysis type filter.
#' @param id Optional analysis id filter.
#' @param model Optional aisdk model.
#' @param structured_output Logical. Use the tool-loop plus structured-output
#'   summary stage. Defaults to `TRUE`.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return A structured `sclet_ai_result`.
#' @export
AIExplainAnalysis <- function(object, type = NULL, id = NULL, model = NULL, structured_output = TRUE, ...) {
    target <- type %||% id
    context <- sclet_ai_context(object, target = target, detail = "full")
    sclet_ai_call(
        task = "analysis_explanation",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Explain the selected analysis in terms of inputs, method, outputs,",
            "evidence and limitations. Distinguish measured results from hypotheses.",
            "Do not upgrade association to causality."
        ),
        ...
    )
}

#' Ask the AI to recommend a conservative next analysis step
#'
#' @param object A `SingleCellExperiment` object.
#' @param model Optional aisdk model.
#' @param structured_output Logical. Use the tool-loop plus structured-output
#'   summary stage. Defaults to `TRUE`.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return A structured `sclet_ai_result` containing recommendations and proposed actions.
#' @export
AIRecommendNextStep <- function(object, model = NULL, structured_output = TRUE, ...) {
    context <- sclet_ai_context(object, target = "planning", detail = "summary")
    sclet_ai_call(
        task = "next_step_recommendation",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Recommend the next analysis only when its prerequisites are visible.",
            "List missing inputs and risks. Return proposed actions for review,",
            "but never imply that they were executed."
        ),
        ...
    )
}
