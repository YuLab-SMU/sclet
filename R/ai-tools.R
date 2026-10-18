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

# Annotate a read-only AI result with the evidence claim ceiling
#
# The read-only entry points (`AIInvestigate()`, `AIExplainAnalysis()`,
# `AIRecommendNextStep()`) hand their findings to a user who may act on them,
# so those findings must carry the same ceiling audit that `RecordAIResult()`
# applies at the write gate -- the gap recorded in
# `.dev/ai-advanced-analysis-spec.md` ("known still-unresolved gaps"). The
# audit itself lives in `R/ai-claims.R`; this helper is only the wiring shared
# by the three call sites, so they behave identically.
#
# A failure in the annotation can never break a read-only call. When the audit
# helper is missing or raises, the result is returned un-annotated with
# `metadata$claim_ceiling_applied = FALSE`, so a caller can tell "annotated"
# from "annotation unavailable" rather than have the difference swallowed. The
# audited value is accepted only when it is a superset of the result passed in,
# so a mis-shaped audit return cannot corrupt the result, and the result class
# is restored when the audit drops it.
sclet_ai_annotate_result_ceiling <- function(object, result) {
    applied <- FALSE
    if (exists("sclet_ai_annotate_findings_ceiling", mode = "function", inherits = TRUE)) {
        annotated <- tryCatch(
            sclet_ai_annotate_findings_ceiling(object, result),
            error = function(e) NULL
        )
        if (is.list(annotated) && is.list(annotated$findings) &&
            all(names(result) %in% names(annotated))) {
            if (inherits(result, "sclet_ai_result") && !inherits(annotated, "sclet_ai_result")) {
                class(annotated) <- class(result)
            }
            result <- annotated
            applied <- TRUE
        }
    }
    metadata <- result$metadata %||% list()
    metadata$claim_ceiling_applied <- applied
    result$metadata <- metadata
    result
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
    evidence_lines <- sclet_ai_evidence_prompt_lines(object)
    result <- sclet_ai_call(
        task = "analysis_explanation",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        system_prompt = sclet_ai_explanation_prompt(evidence_lines),
        ...
    )
    if (!is.list(result)) {
        # A provider that returns something other than a structured result is
        # reported by RunAIAnalysis as an invalid interpretation; no annotation
        # or citation work can apply to it, and touching it would replace the
        # caller's error with a `$` crash.
        return(result)
    }
    result <- sclet_ai_attach_ref_coverage(
        sclet_ai_sanitize_finding_refs(object, result))

    # One revision round, only when it is justified: there is evidence to cite,
    # enough findings for a ratio to mean something, and fewer than half of them
    # are cited. Measured on Kang (2026-10-18): adding the id list to the prompt
    # still produced 0.22-0.3 coverage, because the model kept writing findings
    # with no reference. Asking it once to attach an id or remove the finding is
    # cheaper than shipping a report that is mostly assertions.
    coverage <- result$metadata$evidence_ref_coverage
    if (length(evidence_lines) && isTRUE(
            result$metadata$evidence_ref_coverage$n_findings >= 3L) &&
            isTRUE(coverage$fraction < 0.5)) {
        uncited <- vapply(Filter(function(finding)
            !length(finding$evidence_refs %||% character()), result$findings %||% list()),
            function(finding) as.character(finding$statement %||%
                finding$description %||% "")[[1L]], character(1))
        revision_prompt <- paste0(
            sclet_ai_explanation_prompt(evidence_lines),
            "\n\nYour previous explanation had findings WITHOUT evidence_refs: ",
            paste0("- ", uncited, collapse = "\n"),
            "\nRevise the explanation: attach at least one exact id from the ",
            "list to each finding, or remove that finding and state it in ",
            "`answer` as a hypothesis. Keep supported findings unchanged."
        )
        revised <- sclet_ai_call(
            task = "analysis_explanation",
            context = context,
            tools = sclet_ai_tool_registry(object),
            model = model,
            structured_output = structured_output,
            system_prompt = revision_prompt,
            ...
        )
        revised <- sclet_ai_attach_ref_coverage(
            sclet_ai_sanitize_finding_refs(object, revised))
        revised_coverage <- revised$metadata$evidence_ref_coverage
        if (isTRUE(revised_coverage$fraction > coverage$fraction)) {
            revised$metadata$explanation_revision_rounds <- 1L
            result <- revised
        }
    }
    sclet_ai_annotate_result_ceiling(object, result)
}

# The explanation prompt exists as a function so the citation rules are
# testable. On real data only 3 of 13 findings cited an evidence id (Kang
# contrast run, 2026-10-18); the generic "explain evidence" prompt did not ask
# for a per-finding reference, and an un-referenced finding is indistinguishable
# from one the model made up.
sclet_ai_explanation_prompt <- function(evidence_lines = character()) {
    base <- paste(
        "Explain the selected analysis in terms of inputs, method, outputs,",
        "evidence and limitations. Distinguish measured results from hypotheses.",
        "Do not upgrade association to causality.",
        "Every finding must list the exact evidence ids it rests on in",
        "`evidence_refs`: they start with `ev:` and appear in",
        "context$evidence. Do not cite an AI result id (`sclet-ai-...`) as",
        "evidence. If you cannot support a finding with at least one listed id,",
        "do NOT emit it as a finding: state it in `answer` as a hypothesis",
        "instead. An unsupported finding is worse than a shorter answer."
    )
    if (length(evidence_lines)) {
        # Measured on Kang (2026-10-18): with the ids only inside the JSON
        # context, coverage was 0.73; the model still invented/omitted refs. A
        # short explicit list names what it may cite, so it does not have to
        # search the context for ids.
        base <- paste0(base, "\n\n",
            "The evidence ids available on this object are listed below. Use ",
            "these exact strings; do not add, shorten or translate them:\n",
            paste(evidence_lines, collapse = "\n"))
    }
    base
}

# The evidence ids an explanation may cite, as prompt lines.
#
# Bounded so a large ledger does not dominate the prompt: the first `max_ids`
# nodes in ledger order, plus a count of the rest. User-decision nodes are
# included (they are real ids) but marked, because they support provenance, not
# measured results.
sclet_ai_evidence_prompt_lines <- function(object, max_ids = 20L) {
    nodes <- tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())
    if (!length(nodes)) return(character())
    ids <- names(nodes)
    shown <- utils::head(ids, max_ids)
    lines <- vapply(shown, function(id) {
        node <- nodes[[id]]
        paste0("- ", id, " (", node$kind %||% "?", ", claim ",
            node$claim_level %||% "?",
            if (sclet_ai_evidence_kind_is_user_decision(node$kind %||% "")) {
                ", user decision -- provenance only"
            } else {
                ""
            }, ")")
    }, character(1))
    more <- if (length(ids) > length(shown)) {
        paste0("- ... ", length(ids) - length(shown), " more ids exist but are not listed")
    } else {
        character()
    }
    c(lines, more)
}

# Count how many findings carry an evidence id, so the caller (report,
# benchmark, or a claim audit) can tell "cited" from "asserted".
# A report-level warning when most findings lack an evidence id.
#
# Counting alone was not enough: a run could return 9 findings with 2 cited
# (0.22) and nothing in the report said the rest were assertions. The threshold
# is explicit and only fires when there are enough findings for a ratio to mean
# something.
sclet_ai_coverage_warning <- function(coverage, min_fraction = 0.5, min_findings = 3L) {
    if (is.null(coverage)) return(NULL)
    n <- as.integer(coverage$n_findings %||% 0L)
    fraction <- coverage$fraction
    if (is.na(fraction) || n < min_findings || fraction >= min_fraction) return(NULL)
    paste0(
        "only ", coverage$n_with_refs %||% 0L, " of ", n,
        " findings cite an evidence id; the rest are assertions and should be ",
        "read as hypotheses"
    )
}

sclet_ai_attach_ref_coverage <- function(result) {
    findings <- result$findings %||% list()
    n_findings <- length(findings)
    with_refs <- sum(vapply(findings, function(finding) {
        length(finding$evidence_refs %||% character()) > 0L
    }, logical(1L)))
    metadata <- result$metadata %||% list()
    metadata$evidence_ref_coverage <- list(
        n_findings = as.integer(n_findings),
        n_with_refs = as.integer(with_refs),
        fraction = if (n_findings > 0L) with_refs / n_findings else NA_real_
    )
    result$metadata <- metadata
    result
}

#' Ask the AI to recommend a conservative next analysis step
#'
#' @param object A `SingleCellExperiment` object.
#' @param model Optional aisdk model.
#' @param structured_output Logical. Use the tool-loop plus structured-output
#'   summary stage. Defaults to `TRUE`.
#' @param registry Optional `AIExecutionRegistry`. Its action names,
#'   descriptions, and input schemas are exposed to the model via
#'   `context$available_actions` so recommended actions can reference real,
#'   registered action names rather than invented ones. Defaults to
#'   `AIDefaultExecutionRegistry(object, include = "all")`.
#' @param ... Additional arguments passed to `sclet_ai_call()`.
#' @return A structured `sclet_ai_result` containing recommendations and proposed actions.
#' @export
AIRecommendNextStep <- function(object, model = NULL, structured_output = TRUE, registry = NULL, ...) {
    if (is.null(registry)) {
        registry <- AIDefaultExecutionRegistry(object, include = "all")
    }
    context <- sclet_ai_context(
        object,
        target = "planning",
        detail = "summary",
        registry = registry
    )
    result <- sclet_ai_call(
        task = "next_step_recommendation",
        context = context,
        tools = sclet_ai_tool_registry(object),
        model = model,
        structured_output = structured_output,
        system_prompt = paste(
            "Recommend the next analysis only when its prerequisites are visible.",
            "List missing inputs and risks. Return proposed actions for review,",
            "but never imply that they were executed.",
            "Any proposed action's `action` field must be an exact, verbatim key",
            "from `context$available_actions`, with a `params` list using only the",
            "parameter names listed in that action's input_schema. Do not invent",
            "action names. If no listed action fits, say so in your answer instead.",
            "For group-valued parameters (ident.1, ident.2, start_cluster), use",
            "only the `code` values from `context$capabilities$identity_levels$codes`.",
            "When `all_levels_available` is FALSE or no code fits, omit the",
            "parameter if the action supports omission rather than guessing.",
            "When a parameter must carry a value that an earlier step produces,",
            "take it from the producing action's `output_schema`, or reference it",
            "with the plan's ${step_id.output.field} syntax. Never put a prose",
            "description where a concrete value or reference is required."
        ),
        ...
    )
    sclet_ai_annotate_result_ceiling(object, result)
}
