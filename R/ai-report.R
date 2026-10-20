# Render a recorded AI result as a Markdown report.
#
# Every report written so far (M1, M2) was assembled by the driver script that
# happened to produce the analysis, so each driver wrote its own renderer and
# the report could drift from the ledger it claimed to describe. This renders
# from the LEDGER ONLY and re-checks it: the citations are resolved, each
# recorded claim level is re-derived from the evidence, and any disagreement is
# printed instead of smoothed over. A report that cannot be checked is a
# narrative, not a record.

sclet_ai_report_record <- function(ledger, id = NULL) {
    records <- ledger$analyses %||% list()
    keys <- names(records)
    candidate <- vapply(keys, function(key) {
        grepl("^ai_", key) && !identical(as.character(records[[key]]$type %||% ""), "ai_execution")
    }, logical(1L))
    candidates <- keys[candidate]
    if (is.null(id)) {
        if (!length(candidates)) {
            return(NULL)
        }
        if (length(candidates) > 1L) {
            # Latest wins: a report is a record of one moment, and the caller
            # can always name an older one explicitly.
            stamps <- vapply(candidates, function(key) {
                as.character(records[[key]]$created_at %||% "")
            }, character(1L))
            candidates <- candidates[order(stamps, decreasing = TRUE)]
        }
        return(records[[candidates[[1L]]]])
    }
    record <- records[[id]] %||% NULL
    if (is.null(record)) {
        stop("no recorded AI result with id '", id, "' on this object.", call. = FALSE)
    }
    record
}

sclet_ai_report_refs <- function(findings) {
    refs <- unlist(lapply(findings, function(f) {
        as.character(f$evidence_refs %||% character())
    }), use.names = FALSE)
    unique(refs[nzchar(refs)])
}

# One table row per finding, with the engine's verdict beside it. The point of
# the "audit" column is that a reader never sees a claim level without seeing
# whether the engine moved it.
sclet_ai_report_finding_rows <- function(findings) {
    vapply(findings, function(f) {
        statement <- gsub("\\|", "/", as.character(f$statement %||% ""))
        refs <- as.character(f$evidence_refs %||% character())
        refs <- if (length(refs)) paste0("`", refs, "`", collapse = ", ") else "-"
        recorded <- as.character(f$claim_level %||% "hypothesis")
        requested <- as.character(f$requested_claim_level %||% recorded)
        audit <- if (identical(requested, recorded)) {
            "ok"
        } else {
            paste0(requested, " -> ", recorded)
        }
        sprintf("| %s | %s | %s | %s |", statement, refs, recorded, audit)
    }, character(1L))
}

sclet_ai_report_evidence_rows <- function(ledger, refs) {
    nodes <- ledger$state_records$ai_evidence %||% list()
    if (!length(refs)) {
        return(character())
    }
    vapply(refs, function(ref) {
        node <- nodes[[ref]]$summary %||% NULL
        if (is.null(node)) {
            return(sprintf("| `%s` | **unresolved** | - | - |", ref))
        }
        parents <- as.character(node$parents %||% character())
        group <- as.character(node$dependency_group %||% "")
        sprintf("| `%s` | %s | %s | %s |", ref,
            as.character(node$claim_level %||% "-"),
            if (length(parents)) paste0("`", parents, "`", collapse = ", ") else "-",
            if (nzchar(group)) substr(group, 1L, 12L) else "-")
    }, character(1L))
}

# The self-check: refs must resolve, and each recorded level must equal the
# level the engine derives today. Both are printed, pass or fail.
sclet_ai_report_checks <- function(object, ledger, record) {
    findings <- record$artifacts$findings %||% list()
    refs <- sclet_ai_report_refs(findings)
    records <- ledger$state_records$ai_evidence %||% list()
    unresolved <- setdiff(refs, names(records))
    mismatch <- character()
    for (finding in findings) {
        derived <- tryCatch(
            sclet_ai_derive_finding_claim_level(object, finding)$claim_level,
            error = function(e) NA_character_
        )
        recorded <- as.character(finding$claim_level %||% "hypothesis")
        if (!identical(as.character(derived), recorded)) {
            mismatch <- c(mismatch, paste0(recorded, " vs ", derived))
        }
    }
    list(
        refs = refs,
        unresolved = unresolved,
        level_mismatch = mismatch,
        n_findings = length(findings),
        ok = !length(unresolved) && !length(mismatch)
    )
}

#' Render a recorded AI result as a Markdown report
#'
#' Reads a result recorded by [RecordAIResult()] from the analysis ledger and
#' renders it as Markdown: the findings with their recorded claim levels, the
#' evidence each one cites, the claim ceiling for that evidence, the warnings
#' and recommendations, and the provenance of the run that produced it.
#'
#' The report is rendered from the ledger and nothing else, and it re-checks
#' itself: every cited evidence id must resolve, and every recorded claim level
#' is re-derived from the evidence. A disagreement is printed in the report's
#' verification section rather than smoothed over, so a report that cannot be
#' checked says so.
#'
#' No model is called and no raw values are included: evidence nodes carry
#' bounded aggregates, and the renderer prints what the ledger holds.
#'
#' @param object A `SingleCellExperiment`.
#' @param id Record id of the AI result to render. Defaults to the most
#'   recently recorded AI result on the object.
#' @param path Optional file path. When supplied the report is written there;
#'   the Markdown text is returned either way.
#' @return The report as a single string of Markdown.
#' @examples
#' \dontrun{
#' rendered <- RenderAIReport(sce)
#' RenderAIReport(sce, id = "ai_report_pbmc10k_m1", path = "report.md")
#' }
#' @export
RenderAIReport <- function(object, id = NULL, path = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.", call. = FALSE)
    }
    ledger <- GetAnalysisLedger(object, detail = "full", include_artifacts = TRUE)
    record <- sclet_ai_report_record(ledger, id = id)
    if (is.null(record)) {
        stop("no recorded AI result on this object; record one with RecordAIResult() first.",
            call. = FALSE)
    }
    findings <- record$artifacts$findings %||% list()
    checks <- sclet_ai_report_checks(object, ledger, record)
    ceiling <- NULL
    if (length(checks$refs)) {
        ceiling <- tryCatch(
            check_claim_ceiling(object, evidence_ids = checks$refs),
            error = function(e) NULL
        )
    }
    dims <- dim(object)
    lines <- c(
        paste0("# ", sub("^ai_", "", as.character(record$type %||% "ai_result")),
            " - ", as.character(record$id)),
        "",
        sprintf("- Recorded: %s", as.character(record$created_at %||% "-")),
        sprintf("- Data: %s genes x %s cells", dims[1L], dims[2L]),
        sprintf("- Ledger fingerprint at record time: `%s`; now: `%s`",
            as.character(record$inputs$context_fingerprint %||% "-"),
            as.character(ledger$fingerprint %||% "-")),
        sprintf("- Findings: %s", checks$n_findings),
        "",
        if (nzchar(as.character(record$artifacts$answer %||% ""))) {
            c("## Summary", "", as.character(record$artifacts$answer), "")
        },
        "## Findings (claim levels derived by the engine)",
        "",
        "| Finding | Evidence cited | Recorded level | Engine audit |",
        "|---|---|---|---|",
        if (length(findings)) sclet_ai_report_finding_rows(findings) else "| (no findings) | - | - | - |",
        "",
        "## Evidence",
        "",
        "| Evidence | claim_level | parents | dependency_group |",
        "|---|---|---|---|",
        if (length(checks$refs)) sclet_ai_report_evidence_rows(ledger, checks$refs) else "| (none cited) | - | - | - |",
        ""
    )
    if (!is.null(ceiling)) {
        lines <- c(lines,
            "## Claim ceiling",
            "",
            sprintf("- Ceiling: `%s` (%s independent line(s) of support)",
                as.character(ceiling$ceiling_claim_level %||% "-"),
                as.character(ceiling$n_independent_groups %||% "-")),
            sprintf("- Reason: %s", as.character(ceiling$ceiling_reasons %||% "-")),
            ""
        )
    }
    raw_warnings <- record$artifacts$warnings %||% list()
    warnings <- if (is.list(raw_warnings)) {
        vapply(raw_warnings, function(x) {
            if (is.list(x)) as.character(x$text %||% x$description %||% "") else as.character(x)
        }, character(1L))
    } else {
        as.character(raw_warnings)
    }
    warnings <- warnings[nzchar(warnings)]
    recommendations <- record$artifacts$recommendations %||% list()
    if (length(warnings) || length(recommendations)) {
        lines <- c(lines, "## Warnings and recommendations", "")
        if (length(warnings)) {
            lines <- c(lines, paste0("- Warning: ", warnings))
        }
        if (length(recommendations)) {
            lines <- c(lines, paste0("- Recommendation: ", vapply(recommendations, function(x) {
                if (is.list(x)) as.character(x$text %||% x$description %||% "") else as.character(x)
            }, character(1L))))
        }
        lines <- c(lines, "")
    }
    steps <- sclet_ai_report_steps(ledger)
    if (length(steps)) {
        lines <- c(lines, "## Executed steps (from the execution records)", "", steps, "")
    }
    lines <- c(lines,
        "## Verification",
        "",
        if (length(checks$refs)) {
            sprintf("- Citations resolve: %s of %s",
                if (length(checks$unresolved)) "FAIL" else "PASS",
                length(checks$refs))
        } else {
            "- Citations resolve: nothing cited"
        },
        if (length(checks$unresolved)) {
            paste0("  - unresolved: ", paste0("`", checks$unresolved, "`", collapse = ", "))
        },
        sprintf("- Claim levels re-derived: %s", if (length(checks$level_mismatch)) {
            paste0("FAIL (", paste(checks$level_mismatch, collapse = "; "), ")")
        } else {
            "every recorded level matches the engine"
        }),
        sprintf("- Result: %s", if (checks$ok) "all checks passed" else "**problems found, see above**"),
        "",
        "> Rendered and re-checked from the ledger by `RenderAIReport()`; numbers and citations come from the records, not from a model."
    )
    text <- paste(lines, collapse = "\n")
    if (!is.null(path)) {
        writeLines(lines, path)
    }
    text
}

# The executed steps of every AI run on the object, newest run last. This is
# what makes "did a QC step run" answerable: a read-only action leaves no
# analysis record of its own, but the execution record keeps one row per step.
sclet_ai_report_steps <- function(ledger) {
    records <- ledger$analyses %||% list()
    executions <- Filter(function(record) {
        identical(as.character(record$type %||% ""), "ai_execution")
    }, records)
    if (!length(executions)) {
        return(character())
    }
    stamps <- vapply(executions, function(record) {
        as.character(record$created_at %||% "")
    }, character(1L))
    executions <- executions[order(stamps)]
    rows <- character()
    for (record in executions) {
        steps <- record$params$steps %||% list()
        if (!length(steps)) next
        rows <- c(rows, sprintf("### %s", as.character(record$id %||% "execution")))
        rows <- c(rows, "", "| step | action | status | seconds |", "|---|---|---|---|")
        for (step in steps) {
            rows <- c(rows, sprintf("| %s | %s | %s | %s |",
                as.character(step$id %||% ""),
                as.character(step$action %||% ""),
                as.character(step$status %||% ""),
                format(round(as.numeric(step$duration_sec %||% NA_real_), 2L))))
        }
        rows <- c(rows, "")
    }
    rows
}
