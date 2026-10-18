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

# The evidence ids a set of findings cites, read by EXACT name.
#
# `finding$evidence_refs` is not safe here: R's `$` does partial matching on
# lists, and the recorded finding also carries `evidence_refs_attached`
# (logical). A finding without refs therefore returned the attached FLAG, and
# `as.character(FALSE)` became the literal id "FALSE", which
# `ValidateAIEvidenceRefs()` then rejected -- so every result whose findings
# cited nothing failed to record. `[[` does not partial-match.
sclet_ai_finding_refs <- function(findings) {
    refs <- unlist(lapply(findings, function(finding) {
        if (!is.list(finding)) return(character())
        if (!"evidence_refs" %in% names(finding)) return(character())
        as.character(finding[["evidence_refs"]] %||% character())
    }), use.names = FALSE)
    unique(refs[nzchar(refs)])
}

# Build the findings as they will be recorded. Applies the claim audit's
# derived levels and, when the engine attached an evidence reference on a
# finding's behalf, carries that reference onto the recorded finding. The
# attachment is re-run here rather than read back from the audit because the
# audit reports only THAT an attachment happened, not the resulting id -- and a
# raised `claim_level` whose stated basis is absent from the ledger would be
# exactly the kind of unsupported claim the audit exists to prevent. See
# `.dev/spec-ai-evidence-authorship.md`.
sclet_ai_recorded_findings <- function(object, result, audit, audit_claims) {
    findings <- result$findings
    if (!isTRUE(audit_claims) || is.null(audit) || !length(audit$findings)) {
        return(findings)
    }
    lapply(seq_along(findings), function(i) {
        finding <- findings[[i]]
        entry <- audit$findings[[min(i, length(audit$findings))]]
        # Exact-name read via the helper even though a raw finding cannot carry
        # `evidence_refs_attached` yet: this loop WRITES that flag below, so a
        # second pass over recorded findings would partial-match, and the guard
        # costs nothing here (see `sclet_ai_finding_refs()` for the trap).
        if (is.list(finding) && isTRUE(entry$evidence_refs_attached) &&
            !length(sclet_ai_finding_refs(list(finding)))) {
            attachment <- sclet_ai_attach_finding_evidence(object, finding)
            if (isTRUE(attachment$attached)) {
                finding$evidence_refs <- attachment$finding$evidence_refs
            }
        }
        if (is.list(finding)) {
            finding$claim_level <- entry$claim_level
            finding$requested_claim_level <- entry$requested_claim_level
            finding$derived_from_ceiling <- isTRUE(entry$derived_from_ceiling)
            finding$ceiling_reasons <- entry$ceiling_reasons %||% character()
            finding$evidence_refs_attached <- isTRUE(entry$evidence_refs_attached)
            # The rule name travels with the flag, so a ledger reader can tell
            # WHY: the sole node was attached, the model cited for itself, or
            # the engine declined to guess among several.
            finding$evidence_attachment_rule <- as.character(entry$attachment_rule %||% "none")
        }
        finding
    })
}

#' Record a structured AI result in the analysis ledger
#'
#' @param object A `SingleCellExperiment` object.
#' @param result A `sclet_ai_result` returned by an AI function.
#' @param id Optional stable record id. Defaults to a timestamped task id.
#' @param active Logical. Reserved for future AI state activation; defaults to
#'   `FALSE` so recording a review does not change the active analysis view.
#' @param audit_claims Logical. Audit each finding against the claim ceiling of
#'   the recorded evidence and refuse to record an over-claimed result. When
#'   `TRUE`, also refuse to record findings whose `evidence_refs` do not resolve
#'   to evidence nodes on `object`. Defaults to `TRUE` so recorded claims are
#'   supported by resolvable evidence. Set to `FALSE` to disable both checks.
#' @return The updated `SingleCellExperiment` object. Each record summary
#'   includes `cited_user_decisions`, the evidence ids of cited user-decision
#'   nodes (an empty character vector when none are cited).
#' @export
RecordAIResult <- function(object, result, id = NULL, active = FALSE, audit_claims = TRUE) {
    validate_sclet_ai_result(result)
    audit <- NULL
    if (isTRUE(audit_claims)) {
        audit <- sclet_ai_audit_result_claims(object, result)
        if (identical(audit$status, "overclaimed")) {
            stop("refusing to record an over-claimed AI result: ",
                paste(audit$problems, collapse = "; "),
                call. = FALSE)
        }
    }
    # The recorded findings first, because when claim auditing is on the engine
    # may have ATTACHED an evidence reference on a finding's behalf
    # (see `.dev/spec-ai-evidence-authorship.md`). That reference has to travel
    # with the recorded finding -- otherwise the ledger would show a raised
    # claim whose stated basis is missing -- and it has to be included in the
    # resolvability check below, so an attached id is held to exactly the same
    # bar as one the model cited.
    recorded_findings <- sclet_ai_recorded_findings(object, result, audit, audit_claims)
    all_refs <- sclet_ai_finding_refs(recorded_findings)
    if (isTRUE(audit_claims) && length(all_refs)) {
        # Deliberately no `fingerprint` argument: this check only asks whether
        # each cited evidence id resolves to a real, completed evidence node on
        # the object, not whether that node was recorded against the exact same
        # ledger fingerprint as `result$context`. Fingerprint staleness is a
        # separate concern, and `result$context$fingerprint` is frequently a
        # hand-built placeholder rather than a real ledger fingerprint; treating
        # a mismatch here as an error would wrongly reject results built from
        # synthetic/test contexts whose evidence is otherwise perfectly valid.
        ref_check <- ValidateAIEvidenceRefs(object, all_refs)
        if (length(ref_check$errors)) {
            stop("refusing to record an AI result citing unresolvable evidence: ",
                paste(ref_check$errors, collapse = "; "),
                call. = FALSE)
        }
    }
    cited_user_decisions <- character()
    if (length(all_refs)) {
        all_nodes <- tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())
        cited_user_decisions <- all_refs[vapply(all_refs, function(id) {
            node <- all_nodes[[id]]
            !is.null(node) && sclet_ai_evidence_kind_is_user_decision(node$kind %||% "")
        }, logical(1L))]
    }
    if (is.null(id)) {
        stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
        id <- paste(result$task, stamp, sep = "_")
    }
    if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
        stop("`id` must be a single non-empty character string.")
    }
    # When claim auditing is on, the recorded `claim_level` of each finding is
    # the engine-derived value, never the model's word. The original request is
    # kept on the finding as `requested_claim_level` so a downgrade is visible
    # rather than silent, and the ceiling's `reasons` travel with the finding
    # so a level is never quoted without its basis. This is the change that
    # makes `claim_level` on a recorded finding a function of the evidence,
    # not of the model's self-report.
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
            findings = recorded_findings,
            evidence = result$evidence,
            recommendations = result$recommendations,
            proposed_actions = result$proposed_actions
        ),
        summary = list(
            n_findings = length(recorded_findings),
            n_recommendations = length(result$recommendations),
            n_proposed_actions = length(result$proposed_actions),
            cited_user_decisions = cited_user_decisions
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
#' @param record Logical. Record the execution and AI plan in the ledger. When
#'   `FALSE`, post-execution interpretation is not recorded in the ledger.
#' @param interpret Logical. When `TRUE` and the execution completes, create an
#'   evidence-linked report by calling `AIExplainAnalysis()` on the updated
#'   object. The result is stored in `report$interpretation` only when it is a
#'   valid `sclet_ai_result`; if `record = TRUE`, it is also recorded in the
#'   ledger. Defaults to `FALSE` so no additional LLM call is made unless
#'   explicitly requested. Interpretation errors or an invalid result are
#'   captured in `report$interpretation_error`, and the completed execution is
#'   retained.
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
    interpret = FALSE,
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
                    warnings = validation$warnings,
                    goal_coverage = plan$metadata$goal_coverage %||% NULL
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
    plan_record_error <- NULL
    if (isTRUE(record) && inherits(plan$ai_result, "sclet_ai_result")) {
        # Recording the plan narrative is bookkeeping, not part of the analysis:
        # a citation problem in it must not discard an execution that already
        # ran. Measured on real data (PBMC 10k, 2026-10-18): the plan narrative
        # cited `sclet-ai-*` (AI result ids) instead of `ev:*`, RecordAIResult
        # correctly refused, and the refusal propagated out of RunAIAnalysis --
        # a completed 9/9 analysis was reported as a failed run. Degrade to an
        # unrecorded plan narrative plus a reported error.
        recorded_plan <- tryCatch(
            RecordAIResult(
                updated,
                plan$ai_result,
                id = paste0("ai_plan_", plan$plan_id)
            ),
            error = function(e) conditionMessage(e)
        )
        if (is.character(recorded_plan)) {
            plan_record_error <- recorded_plan
        } else {
            updated <- recorded_plan
        }
    }
    report <- list(
        status = execution$status,
        goal = goal,
        plan_id = plan$plan_id,
        actions = vapply(plan$actions, function(x) x$action, character(1)),
        errors = validation$errors,
        warnings = validation$warnings,
        execution_id = execution$execution_id,
        success_assessment = execution$success_assessment,
        # Did the plan cover what the goal named? A plan can validate and run
        # while silently omitting a step the goal asked for (measured on the
        # first PBMC 10k pilot: annotation was dropped and nothing surfaced it).
        # Reported here, and fed back as one planning round inside
        # AIPlanAnalysis when the gap is detected.
        goal_coverage = plan$metadata$goal_coverage %||% NULL,
        plan_record_error = plan_record_error
    )
    if (isTRUE(interpret) && execution$status %in% c("completed", "completed_with_errors")) {
        interpretation_state <- new.env(parent = emptyenv())
        interpretation_state$error <- NULL
        interpretation_state$result <- NULL
        tryCatch({
            interpretation_state$result <- AIExplainAnalysis(updated, id = execution$execution_id, model = model)
        }, error = function(e) {
            interpretation_state$error <- conditionMessage(e)
        })
        if (!is.null(interpretation_state$error)) {
            report$interpretation_error <- interpretation_state$error
        } else if (!inherits(interpretation_state$result, "sclet_ai_result") ||
                   !isTRUE(validate_sclet_ai_result(interpretation_state$result, error = FALSE))) {
            report$interpretation_error <- paste0(
                "AIExplainAnalysis() returned an invalid interpretation result (class: ",
                paste(class(interpretation_state$result), collapse = "/"),
                "); expected a valid sclet_ai_result."
            )
            if ("interpretation" %in% names(report)) {
                report$interpretation <- NULL
            }
            interpretation_state$result <- NULL
        } else {
            if (isTRUE(record)) {
                recorded_interpretation <- tryCatch(
                    RecordAIResult(
                        updated,
                        interpretation_state$result,
                        id = paste0("ai_interpretation_", plan$plan_id)
                    ),
                    error = function(e) conditionMessage(e)
                )
                if (is.character(recorded_interpretation)) {
                    report$interpretation_error <- recorded_interpretation
                } else {
                    updated <- recorded_interpretation
                    report$interpretation <- interpretation_state$result
                }
            } else {
                report$interpretation <- interpretation_state$result
            }
        }
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
            report = report
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
            "reduction", "counts_assay", "cci_method", "cci_method",
            "unclassified"
        ),
        pattern = c(
            "design_semantics_not_confirmed", "reference_missing", "labels_missing",
            "start_cluster_missing", "start_cluster_unknown", "group_column_missing",
            "reduction_missing", "counts_assay_missing", "method_missing",
            "method_invalid", ""
        ),
        blocked_action = c(
            "run_integration", "run_annotation", "run_annotation",
            "run_trajectory", "run_trajectory", "run_trajectory",
            "run_trajectory", "run_doublet_detection", "run_cci", "run_cci",
            "unknown"
        ),
        related_function = c(
            "ConfirmAIDesignSemantics", "run_annotation", "run_annotation",
            "run_trajectory", "run_trajectory", "run_trajectory",
            "run_trajectory", "NormalizeData", "run_cci", "run_cci",
            NA_character_
        ),
        stringsAsFactors = FALSE
    )
}

# Engine vocabulary `run_cci` validates against. One definition, so a recorded
# clarification answer can be checked against the same list the action enforces.
sclet_ai_clarification_engines <- function() {
    c("CellChat", "CellPhoneDB", "NicheNet")
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
        cci_method = paste0(
            "A cell-cell communication step was blocked because the engine was not chosen. ",
            "Which engine should be used: CellChat, CellPhoneDB or NicheNet? ",
            "Pass it explicitly as the 'method' parameter of run_cci; the engines use ",
            "different ligand-receptor databases, so this assistant will not choose one for you."
        ),
        "The plan was rejected by an existing validation rule and needs a human decision before it can run."
    )
}

# Can a clarification answer be written into the ledger and into the planning
# context?
#
# The clarification layer records "the user was asked and answered" as a
# user_decision evidence node, and that node's summary is what the planner sees
# in `analysis_story$user_decisions` on the next round. Writing the answer
# itself is therefore only allowed when it belongs to a bounded public
# vocabulary: an engine name from the fixed engine list, a celldex dataset
# name, a colData / assay / reduction name that already exists on the object, or
# an already-anonymized `cluster_N` code. Anything else -- free prose, a sample
# name, a path -- is recorded as "answered, value withheld", so the audit can
# see that a choice was made without the ledger absorbing the text.
sclet_ai_clarification_public_answer <- function(question_id, answer, object) {
    if (!is.character(answer) || length(answer) != 1L || is.na(answer) || !nzchar(answer)) {
        return(list(public = FALSE, code = 0L, text = "answered, value withheld"))
    }
    question_id <- as.character(question_id)[[1L]]
    public <- FALSE
    code <- 0L
    if (identical(question_id, "cci_method")) {
        engines <- sclet_ai_clarification_engines()
        code <- match(answer, engines)
        public <- !is.na(code)
        if (is.na(code)) code <- 0L
    } else if (identical(question_id, "annotation_reference")) {
        public <- !is.null(sclet_celldex_reference_name(answer))
    } else if (identical(question_id, "annotation_labels")) {
        public <- grepl("^[A-Za-z][A-Za-z0-9._]{0,63}$", answer)
    } else if (identical(question_id, "trajectory_root")) {
        public <- grepl("^(cluster|group|sample|condition|route)_[0-9]+$", answer)
    } else if (identical(question_id, "design_batch") ||
            identical(question_id, "trajectory_group") ||
            identical(question_id, "group_column")) {
        public <- answer %in% colnames(SummarizedExperiment::colData(object))
    } else if (identical(question_id, "reduction")) {
        public <- answer %in% SingleCellExperiment::reducedDimNames(object)
    } else if (identical(question_id, "counts_assay")) {
        public <- answer %in% SummarizedExperiment::assayNames(object)
    }
    text <- if (isTRUE(public)) {
        paste0("user chose ", question_id, " = ", answer)
    } else {
        paste0("user answered ", question_id, " (value withheld)")
    }
    list(public = isTRUE(public), code = as.integer(code), text = text)
}

# Why an answer cannot be accepted for a question, or NULL when it can.
#
# This checks the answer against the vocabulary the blocked action itself would
# accept; it does not decide anything the action does not. Used by
# ResolveAIClarifications() so a typo is reported as `invalid` instead of being
# recorded as if the user had made a usable choice.
sclet_ai_clarification_answer_problem <- function(question_id, answer, object) {
    question_id <- as.character(question_id)[[1L]]
    if (identical(question_id, "design_batch")) {
        if (!answer %in% colnames(SummarizedExperiment::colData(object))) {
            return("answer is not an existing colData column")
        }
        return(NULL)
    }
    if (identical(question_id, "cci_method")) {
        if (!answer %in% sclet_ai_clarification_engines()) {
            return(paste0(
                "answer is not one of the registered engines: ",
                paste(sclet_ai_clarification_engines(), collapse = ", ")
            ))
        }
        return(NULL)
    }
    if (identical(question_id, "annotation_reference")) {
        name <- sclet_celldex_reference_name(answer)
        if (is.null(name)) {
            return("answer is not a celldex dataset name")
        }
        if (requireNamespace("celldex", quietly = TRUE) &&
                is.null(sclet_celldex_dataset(name))) {
            return(paste0("answer names '", name, "', which celldex does not export"))
        }
        return(NULL)
    }
    if (identical(question_id, "annotation_labels")) {
        if (!grepl("^[A-Za-z][A-Za-z0-9._]{0,63}$", answer)) {
            return("answer is not a single identifier-shaped labels column name")
        }
        return(NULL)
    }
    NULL
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
#' treated as "not answered" rather than being filled in. An answer that the
#' blocked action would not accept is reported as `invalid` and is neither
#' applied nor recorded (for example an engine name that `run_cci` does not
#' register, or a labels column that is not an identifier).
#'
#' Answers are applied through an existing confirmation mechanism where one
#' exists (`design_batch` goes through `ConfirmAIDesignSemantics()`). Every
#' accepted answer is recorded as a `user_decision` evidence node; when the
#' answer is a public vocabulary term (an engine, a celldex dataset, an existing
#' column, an anonymized cluster code) the node's summary names it, and that
#' summary is what the planner reads on the next round as
#' `analysis_story$user_decisions`. That is the channel by which "the engine
#' asked, the user answered" stops being a transcript the caller has to thread
#' back into its own goal text.
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
#'   per-question `resolved` / `skipped` / `invalid` outcomes. Each resolved
#'   entry distinguishes `applied` (something was written to the object through
#'   a confirmation or an applier) from `recorded` (the answer was written to
#'   the ledger as a `user_decision` node); a plan-parameter answer is recorded
#'   without being applied.
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
        # An unusable answer is reported as invalid, never recorded as if the
        # user had made a usable choice. The design branch also performs the
        # confirmation itself; the other branches only need the answer checked
        # against the vocabulary their action accepts.
        problem <- sclet_ai_clarification_answer_problem(question_id, answer, object)
        if (!is.null(problem)) {
            invalid <- c(invalid, question_id)
            transcript[[length(transcript) + 1L]] <- list(
                id = question_id, outcome = "invalid", applied = FALSE,
                reason = problem
            )
            next
        }
        applied <- FALSE
        if (identical(question_id, "design_batch")) {
            object <- ConfirmAIDesignSemantics(object, design = list(batch = answer))
            applied <- TRUE
        } else if (is.function(apply_answer)) {
            object <- apply_answer(object, question, answer)
            applied <- TRUE
        }
        # `applied` says whether something was written to the OBJECT through a
        # confirmation or an applier. A recorded-only answer (a plan parameter
        # such as the communication engine) is `recorded = TRUE` but
        # `applied = FALSE`: nothing on the object changed, and reporting it as
        # applied made the transcript read like a modification that never
        # happened.
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
            id = question_id, outcome = "resolved", applied = applied,
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
#' The answer is written into the node's summary -- which is what the planner
#' reads back as `analysis_story$user_decisions` -- only when it is a public
#' vocabulary term: an engine name, a celldex dataset name, an existing
#' `colData`/assay/reduction name, or an already-anonymized `cluster_N` code.
#' Any other answer is recorded as "answered, value withheld", so the audit sees
#' that a choice was made while the ledger does not absorb free text.
#'
#' The numeric `values` never carry the text: only codes identifying the
#' question, the blocked action, the named `colData` column and the answer's
#' position in the public vocabulary, so the payload stays within the existing
#' `sclet_ai_evidence_value_ok()` rules without a special case.
#'
#' @param object A `SingleCellExperiment` object.
#' @param question_id Question id from `sclet_ai_format_clarification()`.
#' @param answer The user's answer. Recorded verbatim only when it is a public
#'   vocabulary term (see Details); otherwise the value is withheld.
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
    # The answer is only written when it is a public vocabulary term; see
    # sclet_ai_clarification_public_answer(). The summary is what the planner
    # sees as `analysis_story$user_decisions` on the next round, so this is the
    # channel that turns "the user was asked" into "the user's choice is known".
    public <- sclet_ai_clarification_public_answer(question_id, answer, object)
    values <- list(
        question_code = as.integer(question_code),
        question_is_known = question_id %in% known,
        blocked_action_code = as.integer(action_code),
        answer_column_index = if (is.na(answer_index)) NA_integer_ else as.integer(answer_index),
        answer_column_present = !is.na(answer_index),
        answer_is_public_vocabulary = isTRUE(public$public),
        answer_code = as.integer(public$code),
        free_text_recorded = FALSE,
        raw_values_included = FALSE
    )
    evidence <- list(
        id = sclet_ai_clarification_evidence_id(object, question_code, action_code),
        kind = "user_decision",
        values = values,
        claim_level = "observed",
        summary = public$text
    )
    RecordAIEvidence(object, evidence, source = NULL, parents = character(), scope = NULL)
}

# The evidence id for one clarification answer, made unique by observation.
#
# The obvious id -- ev:clarification_q<question>_a<action> -- is fixed per
# question and blocked action, so a SECOND answer to the same question (the
# user changed their mind, or the engine hit the same gate in a later round)
# collided with the first. `RecordAIEvidence()` refuses a duplicate id, and
# `ResolveAIClarifications()` records inside a tryCatch, so the second answer
# was silently dropped: the transcript said "recorded" while the ledger kept
# only the first, stale choice. The first answer keeps the plain id (existing
# records and tests depend on it); later answers gain a numeric suffix.
sclet_ai_clarification_evidence_id <- function(object, question_code, action_code) {
    base <- paste0("ev:clarification_q", question_code, "_a", action_code)
    existing <- tryCatch(
        names(sclet_ai_evidence_get_all(object)),
        error = function(e) character()
    )
    if (!base %in% existing) return(base)
    index <- 1L
    repeat {
        index <- index + 1L
        candidate <- paste0(base, "_", index)
        if (!candidate %in% existing) return(candidate)
    }
}
