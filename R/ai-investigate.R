# Read-only AI investigation facade for advanced analysis questions

sclet_ai_question_requires_design <- function(question) {
    grepl(
        "batch|integration|condition|sample|subject|patient|donor|trajectory|root|reference|\u6279\u6b21|\u6574\u5408|\u6761\u4ef6|\u6837\u672c|\u75c5\u4eba|\u4f9b\u4f53|\u8f68\u8ff9|\u6839\u8282\u70b9|\u53c2\u8003",
        question,
        ignore.case = TRUE
    )
}

sclet_ai_investigate_design <- function(object, design = NULL, question = NULL) {
    if (!is.null(design)) {
        return(list(
            status = "user_supplied",
            design = design,
            clarification_required = FALSE
        ))
    }
    if (is.null(design) && !is.null(question) && !sclet_ai_question_requires_design(question)) {
        return(list(
            status = "not_required",
            design = NULL,
            clarification_required = FALSE
        ))
    }
    readiness <- tryCatch(check_integration_readiness(object), error = function(e) NULL)
    if (!is.null(readiness) && identical(readiness$status, "clarification_required")) {
        return(list(
            status = "clarification_required",
            design = readiness$design,
            questions = readiness$questions,
            blocked_actions = readiness$blocked_actions,
            clarification_required = TRUE
        ))
    }
    list(
        status = "not_supplied",
        design = NULL,
        clarification_required = FALSE
    )
}

sclet_ai_investigate_context <- function(object, question, design, scope) {
    profile <- if (exists("GetAIProfile", mode = "function", inherits = TRUE)) {
        tryCatch(GetAIProfile(object), error = function(e) NULL)
    } else NULL
    ledger <- GetAnalysisLedger(
        object,
        detail = "summary",
        include_artifacts = FALSE,
        include_data = FALSE
    )
    diagnostics <- list(
        integration_readiness = tryCatch(
            check_integration_readiness(object, design = design),
            error = function(e) list(status = "not_available", reason = conditionMessage(e))
        )
    )
    list(
        schema_version = "advanced-investigation-1.0",
        question = question,
        scope = scope,
        design = design,
        profile = profile,
        ledger = ledger,
        diagnostics = diagnostics,
        execution = list(
            allowed = FALSE,
            mutation = FALSE,
            actions_executed = FALSE
        ),
        privacy = list(
            complete_matrix_included = FALSE,
            assay_values_included = FALSE,
            metadata_values_included = FALSE,
            raw_identifiers_included = FALSE
        )
    )
}

#' Investigate an advanced single-cell analysis question without executing actions
#'
#' `AIInvestigate()` is a read-only facade for dataset-specific questions. It
#' combines a bounded dataset profile, the analysis ledger, and deterministic
#' readiness diagnostics. It never changes the `SingleCellExperiment` and never
#' executes a proposed action.
#'
#' @param object A `SingleCellExperiment` object.
#' @param question A single non-empty scientific question.
#' @param design Optional named list describing sample, batch, condition or
#'   subject semantics supplied by the user.
#' @param scope Investigation scope.
#' @param model Optional aisdk model; omitted uses the configured default.
#' @param structured_output Logical; request structured provider output.
#' @param fallback_on_structure_error Logical; retain a normalized answer if
#'   structured output is unavailable.
#' @param ... Additional adapter arguments.
#' @return A read-only `sclet_ai_result`.
#' @export
AIInvestigate <- function(
    object,
    question,
    design = NULL,
    scope = c("diagnosis", "advanced", "report"),
    model = NULL,
    structured_output = TRUE,
    fallback_on_structure_error = TRUE,
    ...
) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    if (!is.character(question) || length(question) != 1L || is.na(question) || !nzchar(trimws(question))) {
        stop("question must be one non-empty character string", call. = FALSE)
    }
    scope <- match.arg(scope)
    design_state <- sclet_ai_investigate_design(object, design, question = question)
    context <- sclet_ai_investigate_context(
        object,
        question = question,
        design = design_state$design,
        scope = scope
    )
    system_prompt <- paste(
        "You are investigating an advanced single-cell analysis question.",
        "Use only deterministic facts in the supplied profile, ledger, and diagnostics.",
        "Do not execute actions, invent metadata semantics, or claim causal effects.",
        "Return findings with evidence references, missing evidence, alternative explanations,",
        "and candidate routes. If design semantics are missing, return clarification_required.",
        sep = " "
    )
    result <- sclet_ai_call(
        task = question,
        context = context,
        tools = list(),
        model = model,
        system_prompt = system_prompt,
        structured_output = structured_output,
        fallback_on_structure_error = fallback_on_structure_error,
        ...
    )
    result$metadata <- utils::modifyList(
        result$metadata %||% list(),
        list(
            interface = "AIInvestigate",
            scope = scope,
            read_only = TRUE,
            actions_executed = FALSE,
            clarification_required = isTRUE(design_state$clarification_required)
        )
    )
    if (isTRUE(design_state$clarification_required)) {
        result$warnings <- c(
            result$warnings,
            list(list(
                code = "clarification_required",
                questions = design_state$questions,
                blocked_actions = design_state$blocked_actions
            ))
        )
    }
    validate_sclet_ai_result(result)
    result
}
