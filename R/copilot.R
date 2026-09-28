#' Summarize the AI-facing analysis ledger for an LLM
#'
#' This compatibility formatter uses the same bounded ledger as the newer
#' R-native AI functions. New code should prefer `GetAnalysisLedger()`.
#'
#' @param sce A `SingleCellExperiment` object.
#' @return A character string containing the structured ledger.
#' @export
SummarizeContextForLLM <- function(sce) {
    ledger <- GetAnalysisLedger(sce, detail = "full", include_artifacts = FALSE)
    paste(
        paste(
            "You are sclet_copilot, an expert single-cell analysis assistant.",
            "Use only ledger facts; label uncertainty and do not invent causality.",
            sep = "\n"
        ),
        "=== SCLET AI LEDGER ===",
        paste(utils::capture.output(dput(ledger)), collapse = "\n"),
        sep = "\n\n"
    )
}

#' sclet AI Copilot compatibility facade
#'
#' @param sce A `SingleCellExperiment` object.
#' @param question Character string containing the user task.
#' @param model Optional aisdk model identifier or model object.
#' @param structured_output Logical. Use the tool-loop plus structured-output
#'   summary stage. Defaults to `TRUE`.
#' @param ... Additional arguments passed to the AI adapter.
#' @return A character response from the structured sclet AI adapter.
#' @export
sclet_copilot <- function(sce, question, model = NULL, structured_output = TRUE, ...) {
    result <- sclet_ai_call(
        task = "copilot",
        context = GetAnalysisLedger(sce, detail = "full"),
        tools = sclet_ai_tool_registry(sce),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Answer the user's question using only the ledger and read-only",
            "tools. Explain evidence and limitations explicitly.",
            sep = "\n"
        ),
        ...
    )
    result$answer %||% ""
}

#' Audit an analysis chain through the R-native AI layer
#'
#' @param sce A `SingleCellExperiment` object.
#' @param features Candidate features to audit.
#' @param raw_layer Raw assay name.
#' @param integrated_layer Corrected/integrated assay name. If `NULL`, infer it
#'   from the integration state when possible.
#' @param ... Additional arguments passed to `AIAuditAnalysisChain()`.
#' @return A character report for backward compatibility.
#' @export
AuditAnalysisChain <- function(
    sce,
    features,
    raw_layer = "counts",
    integrated_layer = NULL,
    ...
) {
    result <- AIAuditAnalysisChain(
        sce = sce,
        features = features,
        raw_layer = raw_layer,
        integrated_layer = integrated_layer,
        ...
    )
    result$answer %||% ""
}
