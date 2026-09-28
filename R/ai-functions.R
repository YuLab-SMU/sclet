#' Audit a feature analysis chain with the R-native AI adapter
#'
#' @param sce A `SingleCellExperiment` object.
#' @param features Character vector of features to audit.
#' @param raw_layer Raw assay name.
#' @param integrated_layer Corrected/integrated assay name. If `NULL`, infer it
#'   from the integration record when possible.
#' @param model Optional aisdk model.
#' @param structured_output Logical. Use the tool-loop plus structured-output
#'   summary stage. Defaults to `TRUE`.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return A structured `sclet_ai_result`.
#' @export
AIAuditAnalysisChain <- function(
    sce,
    features,
    raw_layer = "counts",
    integrated_layer = NULL,
    model = NULL,
    structured_output = TRUE,
    ...
) {
    features <- intersect(as.character(features), rownames(sce))
    if (!length(features)) {
        stop("None of the specified features are found in the SCE object.")
    }
    if (is.null(integrated_layer)) {
        integration <- get_integration(sce)
        if (!is.null(integration)) {
            integrated_layer <- integration$artifacts$assay %||%
                integration$artifacts$layer
        }
    }
    assays <- SummarizedExperiment::assayNames(sce)
    if (!raw_layer %in% assays) {
        stop(sprintf("Raw layer '%s' not found.", raw_layer))
    }
    if (is.null(integrated_layer) || !integrated_layer %in% assays) {
        stop("An integrated assay is required for AIAuditAnalysisChain().")
    }
    context <- sclet_ai_context(sce, target = "report", detail = "full")
    context$requested_audit <- list(
        features = features,
        raw_layer = raw_layer,
        integrated_layer = integrated_layer
    )
    sclet_ai_call(
        task = "analysis_chain_audit",
        context = context,
        tools = sclet_ai_tool_registry(sce),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Audit whether the requested features are robust across the named",
            "expression layers. Discuss design limitations and evidence scope.",
            "Do not calculate statistics that are not present in the ledger.",
            sep = "\n"
        ),
        ...
    )
}

# Keep the original character-returning facade while routing calls through the
# structured adapter. This definition is collated after the legacy copilot.R.
sclet_copilot <- function(sce, question, model = NULL) {
    result <- sclet_ai_call(
        task = "copilot",
        context = GetAnalysisLedger(sce, detail = "full"),
        tools = sclet_ai_tool_registry(sce),
        model = model,
        system_prompt = paste(
            "Answer the user's question using only the ledger and read-only",
            "tools. Explain evidence and limitations explicitly.",
            sep = "\n"
        )
    )
    result$answer %||% ""
}

#' Record a structured AI result in the analysis ledger
#'
#' @param object A `SingleCellExperiment` object.
#' @param result A `sclet_ai_result` returned by an AI function.
#' @param id Optional stable record id. Defaults to a timestamped task id.
#' @param active Logical. Reserved for future AI state activation; defaults to
#'   `FALSE` so recording a review does not change the active analysis view.
#' @return The updated `SingleCellExperiment` object.
#' @export
RecordAIResult <- function(object, result, id = NULL, active = FALSE) {
    validate_sclet_ai_result(result)
    if (is.null(id)) {
        stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
        id <- paste(result$task, stamp, sep = "_")
    }
    if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
        stop("`id` must be a single non-empty character string.")
    }
    context <- result$context
    record <- list(
        id = id,
        type = paste0("ai_", result$task),
        method = result$metadata$provider %||% "aisdk",
        inputs = list(
            task = result$task,
            context_schema_version = context$schema_version %||% NULL,
            context_fingerprint = context$fingerprint %||% NULL
        ),
        params = result$metadata,
        artifacts = list(
            answer = result$answer,
            findings = result$findings,
            evidence = result$evidence,
            recommendations = result$recommendations,
            proposed_actions = result$proposed_actions
        ),
        summary = list(
            n_findings = length(result$findings),
            n_recommendations = length(result$recommendations),
            n_proposed_actions = length(result$proposed_actions)
        ),
        created_at = Sys.time()
    )
    # AI review records live in the general analysis registry for compatibility;
    # they do not alter the active computational analysis state.
    sclet_set_analysis(object, id, record)
}
