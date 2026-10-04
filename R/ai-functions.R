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

#' Ask the AI a bounded question about a sclet analysis
#'
#' `AskAI()` is the beginner-facing read-only conversational interface. It sends
#' the question together with the bounded analysis ledger, never the complete
#' expression matrix, and returns a structured `sclet_ai_result`.
#'
#' @param object A `SingleCellExperiment` object.
#' @param question A single natural-language question.
#' @param model Optional aisdk model.
#' @param structured_output Logical. Use structured output when available.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return A structured `sclet_ai_result`.
#' @export
AskAI <- function(object, question, model = NULL, structured_output = TRUE, ...) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    if (!is.character(question) || length(question) != 1L ||
        is.na(question) || !nzchar(trimws(question))) {
        stop("`question` must be one non-empty character value.")
    }
    context <- GetAnalysisLedger(object, detail = "full")
    sclet_ai_call(
        task = "copilot",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Answer the beginner's question using only the supplied sclet ledger.",
            "Explain what is known, what is uncertain, and what the beginner should do next.",
            "Do not execute actions, invent matrix-level statistics, or upgrade hypotheses",
            "to biological conclusions.",
            "Beginner question:", question,
            sep = "\n"
        ),
        ...
    )
}

#' Run a beginner-friendly AI-led sclet analysis
#'
#' `RunAIAnalysis()` hides registry, validation, dry-run, and confirmation details
#' behind one guided entry point. It asks the AI for a plan, validates the plan
#' against a built-in allowlist, previews the actions, and executes only after
#' explicit approval. With `confirm = "ask"`, interactive sessions prompt the
#' user; non-interactive sessions stop after the dry run.
#'
#' @param object A `SingleCellExperiment` object.
#' @param goal A single natural-language analysis goal.
#' @param model Optional aisdk model.
#' @param confirm One of `"ask"`, `"yes"`, or `"no"`. Logical `TRUE`/`FALSE`
#'   are also accepted. `"ask"` prompts only in interactive sessions.
#' @param registry Optional `AIExecutionRegistry`. If `NULL`, the beginner
#'   registry includes the native preprocessing, dimensional-reduction, graph,
#'   and clustering actions.
#' @param include Action groups used to build the default registry.
#' @param structured_output Logical. Use structured output for the AI plan.
#' @param fallback_on_structure_error Logical. Retain the tool-loop plan if
#'   structured output fails.
#' @param record Logical. Record the execution and AI plan in the ledger.
#' @param ... Additional arguments passed to `AIPlanAnalysis()`.
#' @return An object of class `sclet_ai_analysis` containing the original or
#'   updated object, plan, validation, dry-run preview, execution, and report.
#' @export
RunAIAnalysis <- function(
    object,
    goal,
    model = NULL,
    confirm = "ask",
    registry = NULL,
    include = c("read", "preprocess", "dimred", "graph", "cluster"),
    structured_output = TRUE,
    fallback_on_structure_error = TRUE,
    record = TRUE,
    ...
) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    if (!is.character(goal) || length(goal) != 1L || is.na(goal) || !nzchar(trimws(goal))) {
        stop("`goal` must be one non-empty character value.")
    }
    if (is.logical(confirm)) {
        if (length(confirm) != 1L || is.na(confirm)) stop("`confirm` must be one value.")
        confirm <- if (confirm) "yes" else "no"
    } else {
        confirm <- match.arg(confirm, c("ask", "yes", "no"))
    }
    if (is.null(registry)) {
        registry <- AIDefaultExecutionRegistry(object, include = include)
    } else if (!inherits(registry, "sclet_ai_execution_registry")) {
        registry <- AIExecutionRegistry(registry)
    }

    plan <- AIPlanAnalysis(
        object,
        goal = goal,
        model = model,
        structured_output = structured_output,
        fallback_on_structure_error = fallback_on_structure_error,
        ...
    )
    validation <- ValidateAIPlan(plan, object = object, registry = registry, strict = FALSE)
    if (!isTRUE(validation$valid)) {
        return(structure(
            list(
                object = object,
                plan = plan,
                registry = registry,
                validation = validation,
                preview = NULL,
                execution = NULL,
                status = "invalid_plan",
                report = list(
                    status = "invalid_plan",
                    goal = goal,
                    errors = validation$errors,
                    warnings = validation$warnings,
                    clarification = sclet_ai_format_clarification(validation$errors)
                )
            ),
            class = c("sclet_ai_analysis", "list")
        ))
    }

    preview <- ExecuteAIPlan(
        object,
        plan,
        registry,
        validation = validation,
        dry_run = TRUE,
        record = FALSE
    )
    approved <- identical(confirm, "yes")
    if (identical(confirm, "ask")) {
        if (interactive()) {
            cat("\nAI proposes", length(plan$actions), "registered actions for:\n", goal, "\n")
            cat("Execute this plan? [y/N] ")
            answer <- tolower(trimws(readLines(n = 1L, warn = FALSE)))
            approved <- answer %in% c("y", "yes")
        } else {
            message("Non-interactive session: returning the dry-run preview without execution.")
        }
    }
    if (!approved) {
        return(structure(
            list(
                object = object,
                plan = plan,
                registry = registry,
                validation = validation,
                preview = preview,
                execution = NULL,
                status = "dry_run",
                report = list(
                    status = "dry_run",
                    goal = goal,
                    plan_id = plan$plan_id,
                    actions = vapply(plan$actions, function(x) x$action, character(1)),
                    errors = validation$errors,
                    warnings = validation$warnings
                )
            ),
            class = c("sclet_ai_analysis", "list")
        ))
    }

    execution <- ExecuteAIPlan(
        object,
        plan,
        registry,
        validation = validation,
        dry_run = FALSE,
        confirmation = validation$confirmation_token,
        record = record
    )
    updated <- execution$object
    if (isTRUE(record) && inherits(plan$ai_result, "sclet_ai_result")) {
        updated <- RecordAIResult(
            updated,
            plan$ai_result,
            id = paste0("ai_plan_", plan$plan_id)
        )
    }
    structure(
        list(
            object = updated,
            plan = plan,
            registry = registry,
            validation = validation,
            preview = preview,
            execution = execution,
            status = execution$status,
            report = list(
                status = execution$status,
                goal = goal,
                plan_id = plan$plan_id,
                actions = vapply(plan$actions, function(x) x$action, character(1)),
                errors = validation$errors,
                warnings = validation$warnings,
                execution_id = execution$execution_id
            )
        ),
        class = c("sclet_ai_analysis", "list")
    )
}

#' @export
print.sclet_ai_analysis <- function(x, ...) {
    cat("sclet beginner AI analysis [", x$status, "]\n", sep = "")
    if (!is.null(x$report$goal)) cat("Goal: ", x$report$goal, "\n", sep = "")
    if (!is.null(x$report$actions)) {
        cat("Actions: ", paste(x$report$actions, collapse = " -> "), "\n", sep = "")
    }
    invisible(x)
}

sclet_ai_clarification_catalog <- function() {
    data.frame(
        question_id = c(
            "design_batch", "annotation_reference", "annotation_labels",
            "trajectory_root", "trajectory_root_unknown", "trajectory_group",
            "reduction", "counts_assay", "unclassified"
        ),
        pattern = c(
            "design_semantics_not_confirmed", "reference_missing", "labels_missing",
            "start_cluster_missing", "start_cluster_unknown", "group_column_missing",
            "reduction_missing", "counts_assay_missing", ""
        ),
        blocked_action = c(
            "run_integration", "run_annotation", "run_annotation",
            "run_trajectory", "run_trajectory", "run_trajectory",
            "run_trajectory", "run_doublet_detection", "unknown"
        ),
        related_function = c(
            "ConfirmAIDesignSemantics", "run_annotation", "run_annotation",
            "run_trajectory", "run_trajectory", "run_trajectory",
            "run_trajectory", "NormalizeData", NA_character_
        ),
        stringsAsFactors = FALSE
    )
}

sclet_ai_clarification_question_ids <- function() {
    as.character(sclet_ai_clarification_catalog()$question_id)
}

sclet_ai_clarification_text <- function(question_id) {
    switch(question_id,
        design_batch = paste0(
            "An integration step was blocked because the technical batch column has not been confirmed. ",
            "Which colData column represents the technical batch? Confirm it with ",
            "ConfirmAIDesignSemantics(object, design = list(batch = '<column name>')). ",
            "This assistant will not choose the column for you."
        ),
        annotation_reference = paste0(
            "An annotation step was blocked because no reference dataset was supplied. ",
            "Which reference should be used (for example a celldex dataset name or a ",
            "SummarizedExperiment)? Pass it explicitly as the 'ref' parameter of run_annotation."
        ),
        annotation_labels = paste0(
            "An annotation step was blocked because no labels column was supplied for the reference. ",
            "Which labels column of the reference should be used? Pass it explicitly as the ",
            "'labels' parameter of run_annotation."
        ),
        trajectory_root = paste0(
            "A trajectory step was blocked because the trajectory root was not specified. ",
            "Which cluster represents the origin of the trajectory? Pass it explicitly as the ",
            "'start_cluster' parameter of run_trajectory; automatic root selection is not used ",
            "because the origin is a biological assumption."
        ),
        trajectory_root_unknown = paste0(
            "A trajectory step was blocked because the requested start cluster is not one of the ",
            "current cluster identities. Which existing cluster represents the origin?"
        ),
        trajectory_group = paste0(
            "A trajectory step was blocked because the requested cluster column is not present in ",
            "colData. Which colData column holds the cluster labels?"
        ),
        reduction = paste0(
            "A step was blocked because the requested dimensionality reduction is not available. ",
            "Run the corresponding reduction first, or name a reduction that already exists."
        ),
        counts_assay = paste0(
            "A step was blocked because the object has no 'counts' assay. ",
            "Supply a counts assay, for example by reading raw data with Read10X()."
        ),
        "The plan was rejected by an existing validation rule and needs a human decision before it can run."
    )
}


#' Turn plan validation errors into structured clarification questions
#'
#' This is a presentation helper only. It parses the error strings already
#' produced by `ValidateAIPlan()` and maps recognized patterns onto structured
#' questions. It performs no validation of its own, proposes no answer and no
#' default column, and it never discards an unrecognized error: anything that
#' does not match a known pattern is preserved verbatim as a generic question.
#'
#' @param validation_errors Character vector of errors from `ValidateAIPlan()`.
#' @param readiness_results Optional named list of already-computed readiness
#'   results. Names are question ids; values may supply `blocked_action` or
#'   `text` overrides. Nothing is recomputed here.
#' @return A list with `status`, `questions` and the untouched `raw_errors`.
#' @keywords internal
#' @noRd
sclet_ai_format_clarification <- function(validation_errors, readiness_results = list()) {
    errors <- as.character(validation_errors %||% character())
    catalog <- sclet_ai_clarification_catalog()
    patterns <- catalog$pattern[nzchar(catalog$pattern)]
    questions <- list()
    seen <- character()
    for (error in errors) {
        hit <- patterns[vapply(patterns, function(p) grepl(p, error, fixed = TRUE), logical(1L))]
        matched <- if (length(hit)) {
            catalog$question_id[catalog$pattern %in% hit][[1L]]
        } else {
            "unclassified"
        }
        override <- readiness_results[[matched]]
        blocked_action <- catalog$blocked_action[catalog$question_id == matched][[1L]]
        text <- sclet_ai_clarification_text(matched)
        if (is.list(override)) {
            if (!is.null(override$blocked_action)) blocked_action <- as.character(override$blocked_action)[[1L]]
            if (!is.null(override$text)) text <- as.character(override$text)[[1L]]
        }
        # keep one question per distinct id, but never drop an unrecognized error
        key <- if (identical(matched, "unclassified")) paste0("unclassified:", error) else matched
        if (key %in% seen) next
        seen <- c(seen, key)
        related <- catalog$related_function[catalog$question_id == matched][[1L]]
        questions[[length(questions) + 1L]] <- list(
            id = matched,
            text = text,
            blocked_action = blocked_action,
            related_function = if (is.na(related)) NULL else related,
            recognized = !identical(matched, "unclassified"),
            source_error = error
        )
    }
    list(
        status = if (length(questions)) "clarification_required" else "no_clarification_needed",
        questions = questions,
        n_questions = length(questions),
        raw_errors = errors
    )
}

#' Resolve AI clarification questions interactively
#'
#' Presents the structured questions produced by
#' `sclet_ai_format_clarification()` to the user, collects an answer for each,
#' applies the answer through the existing confirmation mechanism where one
#' exists, and records each accepted answer as a `user_decision` evidence node.
#' This closes the loop described in the clarification contract: the AI was
#' blocked, the user was asked, and the user answered.
#'
#' The function never answers on the user's behalf. It proposes no default
#' column, no default cluster and no default reference, and an empty answer is
#' treated as "not answered" rather than being filled in. Answers are only
#' applied when the question type has a real confirmation function
#' (`design_batch` goes through `ConfirmAIDesignSemantics()`); other question
#' types describe plan parameters, so they are recorded for the audit trail
#' and returned for the caller to feed back into its plan.
#'
#' In a non-interactive session without an `ask` function the function refuses
#' to prompt and returns `status = "needs_interactive"` without modifying the
#' object or recording anything.
#'
#' @title ResolveAIClarifications
#' @param object A `SingleCellExperiment` object.
#' @param clarification A clarification list as returned by
#'   `sclet_ai_format_clarification()`, or `result$report$clarification` from
#'   `RunAIAnalysis()`.
#' @param ask Optional function taking a question list and returning the user's
#'   answer as a character scalar. Defaults to a terminal prompt. Supply this
#'   to embed a different UI or to script the interaction.
#' @param apply_answer Optional function `(object, question, answer)` returning
#'   the updated object, used for question types that have no built-in
#'   confirmation function.
#' @param record Logical. Record accepted answers as `user_decision` evidence.
#' @return A list with the updated `object`, a `status`, and a `transcript` of
#'   per-question `resolved` / `skipped` / `invalid` outcomes.
#' @export
ResolveAIClarifications <- function(
    object,
    clarification,
    ask = NULL,
    apply_answer = NULL,
    record = TRUE) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    questions <- clarification$questions
    if (is.null(questions) || !length(questions)) {
        return(list(
            object = object,
            status = "nothing_to_resolve",
            resolved = character(),
            skipped = character(),
            invalid = character(),
            transcript = list(),
            raw_errors = clarification$raw_errors %||% character()
        ))
    }
    if (!is.null(ask) && !is.function(ask)) stop("ask must be a function", call. = FALSE)
    if (!is.null(apply_answer) && !is.function(apply_answer)) {
        stop("apply_answer must be a function", call. = FALSE)
    }
    if (is.null(ask) && !interactive()) {
        return(list(
            object = object,
            status = "needs_interactive",
            resolved = character(),
            skipped = vapply(questions, function(q) as.character(q$id), character(1L)),
            invalid = character(),
            transcript = list(),
            questions = questions,
            raw_errors = clarification$raw_errors %||% character()
        ))
    }

    resolved <- character()
    skipped <- character()
    invalid <- character()
    transcript <- list()

    for (question in questions) {
        question_id <- as.character(question$id)
        blocked <- question$blocked_action %||% "unknown"
        if (is.null(ask)) {
            cat("\n[sclet] ", question$text, "\n", sep = "")
            if (!is.null(question$related_function)) {
                cat("       Related function: ", question$related_function, "\n", sep = "")
            }
            cat("       Your answer (empty to skip): ")
            answer <- tryCatch(readLines(n = 1L, warn = FALSE), error = function(e) character())
            answer <- if (length(answer)) trimws(as.character(answer)[[1L]]) else ""
        } else {
            answer <- tryCatch(ask(question), error = function(e) NULL)
            answer <- if (is.null(answer)) "" else trimws(as.character(answer)[[1L]])
        }
        if (!nzchar(answer)) {
            skipped <- c(skipped, question_id)
            transcript[[length(transcript) + 1L]] <- list(
                id = question_id, outcome = "skipped", applied = FALSE
            )
            next
        }
        # a design confirmation answer must name a real column before anything
        # is written; an unusable answer is reported, never confirmed
        if (identical(question_id, "design_batch")) {
            columns <- colnames(SummarizedExperiment::colData(object))
            if (!answer %in% columns) {
                invalid <- c(invalid, question_id)
                transcript[[length(transcript) + 1L]] <- list(
                    id = question_id, outcome = "invalid", applied = FALSE,
                    reason = "answer is not an existing colData column"
                )
                next
            }
            object <- ConfirmAIDesignSemantics(object, design = list(batch = answer))
        } else if (is.function(apply_answer)) {
            object <- apply_answer(object, question, answer)
        }
        if (isTRUE(record)) {
            object <- tryCatch(
                sclet_ai_record_clarification_response(
                    object,
                    question_id = question_id,
                    answer = answer,
                    blocked_action = blocked
                ),
                error = function(e) object
            )
        }
        resolved <- c(resolved, question_id)
        transcript[[length(transcript) + 1L]] <- list(
            id = question_id, outcome = "resolved", applied = TRUE,
            blocked_action = blocked, recorded = isTRUE(record)
        )
    }

    status <- if (length(skipped) || length(invalid)) {
        if (length(resolved)) "partially_resolved" else "unresolved"
    } else {
        "resolved"
    }
    list(
        object = object,
        status = status,
        resolved = resolved,
        skipped = skipped,
        invalid = invalid,
        transcript = transcript,
        questions = questions,
        raw_errors = clarification$raw_errors %||% character()
    )
}

#' Record a user's answer to a clarification question as evidence
#'
#' Records that a clarification question was asked and answered, so an audit can
#' reconstruct "the AI was blocked, the user was asked, the user answered". This
#' never performs the confirmation itself: call `ConfirmAIDesignSemantics()` (or
#' whatever the question's `related_function` is) first, then record the response.
#'
#' Free-text answers are deliberately not stored. Only a numeric code identifying
#' the question, a numeric code identifying which `colData` column the user named,
#' and boolean flags are written, so the payload stays within the existing
#' `sclet_ai_evidence_value_ok()` rules without adding a special case to bypass
#' them.
#'
#' @param object A `SingleCellExperiment` object.
#' @param question_id Question id from `sclet_ai_format_clarification()`.
#' @param answer The user's answer. Only used to check whether it names an
#'   existing `colData` column; the text itself is never recorded.
#' @param blocked_action Optional action name that was blocked.
#' @return The SCE with a `user_decision` evidence node.
#' @keywords internal
#' @noRd
sclet_ai_record_clarification_response <- function(object, question_id, answer, blocked_action = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    if (!is.character(question_id) || length(question_id) != 1L || !nzchar(trimws(question_id))) {
        stop("question_id must be a single non-empty string", call. = FALSE)
    }
    if (!is.character(answer) || length(answer) != 1L) {
        stop("answer must be a single character string", call. = FALSE)
    }
    known <- sclet_ai_clarification_question_ids()
    question_code <- match(question_id, known)
    if (is.na(question_code)) question_code <- length(known)
    columns <- colnames(SummarizedExperiment::colData(object))
    answer_index <- match(answer, columns)
    actions <- sort(unique(names(AIDefaultExecutionRegistry(object, include = "all"))))
    action_code <- if (is.null(blocked_action)) 0L else match(as.character(blocked_action)[[1L]], actions)
    if (is.na(action_code)) action_code <- 0L
    values <- list(
        question_code = as.integer(question_code),
        question_is_known = question_id %in% known,
        blocked_action_code = as.integer(action_code),
        answer_column_index = if (is.na(answer_index)) NA_integer_ else as.integer(answer_index),
        answer_column_present = !is.na(answer_index),
        free_text_recorded = FALSE,
        raw_values_included = FALSE
    )
    evidence <- list(
        id = paste0("ev:clarification_q", question_code, "_a", action_code),
        kind = "user_decision",
        values = values,
        claim_level = "observed"
    )
    RecordAIEvidence(object, evidence, source = NULL, parents = character(), scope = NULL)
}
