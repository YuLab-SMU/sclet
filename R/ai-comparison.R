sclet_ai_comparison_records <- function(object, ids = NULL) {
    ledger <- GetAnalysisLedger(object, detail = "summary", include_artifacts = FALSE, include_data = FALSE)
    records <- ledger$analyses %||% list()
    state_records <- ledger$state_records %||% list()
    if (length(state_records)) {
        flat <- unlist(state_records, recursive = FALSE, use.names = FALSE)
        keep <- flat[vapply(flat, function(r) {
            is.list(r) && grepl("integration", tolower(r$type %||% r$summary[["type"]] %||% ""), ignore.case = TRUE)
        }, logical(1L))]
        records <- c(records, keep)
    }
    if (!length(records)) return(list(ledger = ledger, records = list()))
    records <- records[vapply(records, function(record) {
        if (!is.list(record)) return(FALSE)
        type <- tolower(as.character(record$type %||% ""))
        id <- tolower(as.character(record$id %||% ""))
        grepl("integration|harmony|fastmnn|mnn|scvi|batch|correct", paste(type, id))
    }, logical(1L))]
    if (length(records)) {
        record_ids <- vapply(records, function(record) as.character(record$id %||% ""), character(1L))
        records <- records[!duplicated(record_ids)]
    }
    if (!is.null(ids)) {
        ids <- unique(as.character(ids))
        records <- records[vapply(records, function(record) {
            as.character(record$id %||% "") %in% ids
        }, logical(1L))]
    }
    list(ledger = ledger, records = records)
}

sclet_ai_metric <- function(name, value, method = "recorded_summary", baseline_ref = NULL, evidence_id = NULL) {
    list(
        name = name,
        value = value,
        direction = "not_specified",
        scope = list(),
        method = method,
        parameters = list(),
        baseline_ref = baseline_ref,
        uncertainty = list(status = "not_reported"),
        evidence_id = evidence_id
    )
}

sclet_ai_record_metrics <- function(record, baseline_ref = NULL) {
    summary <- record$summary %||% record$summary %||% list()
    metrics <- summary$metrics %||% list()
    if (!is.list(metrics)) return(list())
    result <- list()
    add_structured <- function(m) {
        if (!is.list(m) || is.null(m$name)) return()
        m$baseline_ref <- m$baseline_ref %||% baseline_ref
        if (!"direction" %in% names(m)) m$direction <- "not_specified"
        if (!"uncertainty" %in% names(m)) m$uncertainty <- list(status = "not_reported")
        if (is.null(m$value) && identical(m$status %||% m$uncertainty$status, "not_available")) {
            m$value <- NA_real_
            if (is.null(m$uncertainty$status)) m$uncertainty$status <- "not_available"
            if (!is.null(m$reason) && is.null(m$uncertainty$reason)) m$uncertainty$reason <- m$reason
        }
        result[[length(result) + 1L]] <<- m
    }
    for (i in seq_along(metrics)) {
        m <- metrics[[i]]
        if (is.list(m) && !is.null(m$name)) {
            add_structured(m)
        } else {
            name <- names(metrics)[[i]]
            value <- m
            if (is.null(value) || is.list(value) || !is.numeric(value) && !is.logical(value)) next
            result[[length(result) + 1L]] <- sclet_ai_metric(name, unname(value), "recorded_summary", baseline_ref)
        }
    }
    extras <- intersect(names(summary), c("batch_mixing", "biological_preservation", "cluster_stability", "runtime_sec", "memory_mb"))
    for (name in extras) {
        value <- summary[[name]]
        if (is.null(value) || is.list(value) || !is.numeric(value) && !is.logical(value)) next
        result[[length(result) + 1L]] <- sclet_ai_metric(name, unname(value), "legacy_summary", baseline_ref)
    }
    result
}

sclet_ai_comparison_route_metrics <- function(record, baseline_ref) {
    mets <- sclet_ai_record_metrics(record, baseline_ref = baseline_ref)
    stats::setNames(mets, vapply(mets, function(m) m$name, character(1L)))
}

sclet_ai_comparison_metric_value <- function(route_metrics, name) {
    m <- route_metrics[[name]]
    if (is.null(m)) return(NA_real_)
    val <- m$value
    if (length(val) != 1L || !is.numeric(val) || !is.finite(val)) return(NA_real_)
    unname(val)
}

sclet_ai_comparison_sort_routes <- function(ids, metrics_by_route, criterion, baseline) {
    if (!length(ids)) return(ids)
    score_name <- switch(criterion,
        biological_preservation = "biological_preservation",
        stability = "cluster_stability",
        cost = "runtime_sec",
        evidence = NULL)
    if (is.null(score_name)) return(ids)
    values <- vapply(ids, function(id) {
        sclet_ai_comparison_metric_value(metrics_by_route[[id]] %||% list(), score_name)
    }, numeric(1L))
    direction_is_higher_better <- switch(criterion,
        biological_preservation = TRUE,
        stability = TRUE,
        cost = FALSE,
        TRUE)
    na_last <- rep(NA_real_, length(values))
    if (direction_is_higher_better) {
        order_vec <- order(-values, na.last = TRUE)
    } else {
        order_vec <- order(values, na.last = TRUE)
    }
    ids[order_vec]
}

sclet_ai_comparison_tradeoffs <- function(ids, metrics_by_route, baseline, threshold = 0.1) {
    tradeoffs <- list()
    if (!length(ids) || is.na(baseline) || !(baseline %in% ids)) return(tradeoffs)
    baseline_batch <- sclet_ai_comparison_metric_value(metrics_by_route[[baseline]], "batch_mixing")
    baseline_bio <- sclet_ai_comparison_metric_value(metrics_by_route[[baseline]], "biological_preservation")
    if (!is.finite(baseline_batch) || !is.finite(baseline_bio)) return(tradeoffs)
    for (id in setdiff(ids, baseline)) {
        batch_val <- sclet_ai_comparison_metric_value(metrics_by_route[[id]], "batch_mixing")
        bio_val <- sclet_ai_comparison_metric_value(metrics_by_route[[id]], "biological_preservation")
        if (!is.finite(batch_val) || !is.finite(bio_val)) next
        batch_delta <- batch_val - baseline_batch
        bio_delta <- bio_val - baseline_bio
        if (batch_delta > 0 && bio_delta < (-threshold)) {
            tradeoffs[[length(tradeoffs) + 1L]] <- list(
                type = "batch_vs_bio",
                routes = c(id, baseline),
                batch_delta = unname(batch_delta),
                biological_delta = unname(bio_delta),
                threshold = threshold,
                message = "batch mixing improved but biological preservation declined; no automatic selection."
            )
        }
    }
    tradeoffs
}

sclet_ai_comparison_normalize_missing <- function(metrics, known_names = c("batch_mixing", "biological_preservation", "cluster_stability", "runtime_sec")) {
    for (name in known_names) {
        if (!(name %in% names(metrics))) {
            metrics[[name]] <- list(
                name = name,
                value = NA_real_,
                direction = "not_specified",
                scope = list(),
                method = "not_recorded",
                parameters = list(),
                baseline_ref = NULL,
                uncertainty = list(status = "not_available", reason = "metric_not_recorded"),
                evidence_id = NULL
            )
        } else {
            m <- metrics[[name]]
            if (is.null(m$value) || (!is.finite(m$value) && !identical(m$uncertainty$status, "not_available"))) {
                m$value <- NA_real_
                if (is.null(m$uncertainty$status)) m$uncertainty$status <- "not_available"
                if (is.null(m$uncertainty$reason)) m$uncertainty$reason <- "value_missing_or_not_finite"
            }
            metrics[[name]] <- m
        }
    }
    metrics
}

#' Compare existing integration analysis routes without executing a route
#'
#' @param object A `SingleCellExperiment` object.
#' @param ids Optional analysis ids to compare.
#' @param criterion Comparison criterion.
#' @return A read-only route comparison or a typed not-available result.
#' @export
CompareAIAnalyses <- function(
    object,
    ids = NULL,
    criterion = c("evidence", "stability", "biological_preservation", "cost")
) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    criterion <- match.arg(criterion)
    collected <- sclet_ai_comparison_records(object, ids)
    records <- collected$records
    if (!length(records)) {
        return(list(
            status = "not_available",
            reason = "no_comparable_integration_records",
            criterion = criterion,
            routes = list(),
            baseline = NULL,
            metrics = list(),
            tradeoffs = list(),
            recommendation = NULL,
            fingerprints = list(ledger = collected$ledger$fingerprint),
            warnings = list(),
            execution = list(allowed = FALSE, performed = FALSE)
        ))
    }
    ids_found <- vapply(records, function(record) as.character(record$id %||% "unknown"), character(1L))
    baseline <- ids_found[grepl("raw|baseline|original", tolower(ids_found))][1L]
    if (is.na(baseline) || !length(baseline)) baseline <- if (length(ids_found)) ids_found[[1L]] else NULL
    baseline <- unname(baseline)
    routes_raw <- lapply(records, function(record) {
        list(
            id = record$id %||% NULL,
            type = record$type %||% NULL,
            method = record$method %||% NULL,
            status = record$status %||% record$summary$status %||% "unknown",
            inputs = record$inputs %||% list(),
            params = record$params %||% list(),
            summary = record$summary %||% list(),
            artifacts = record$artifacts %||% list(),
            evidence_refs = record$evidence_refs %||% character(),
            raw_values_included = FALSE
        )
    })
    names(routes_raw) <- ids_found
    metrics_by_route <- lapply(records, function(r) {
        id <- as.character(r$id %||% "unknown")
        structured <- sclet_ai_comparison_route_metrics(r, baseline_ref = baseline)
        sclet_ai_comparison_normalize_missing(structured)
    })
    names(metrics_by_route) <- ids_found
    sorted_ids <- sclet_ai_comparison_sort_routes(ids_found, metrics_by_route, criterion, baseline)
    routes <- routes_raw[sorted_ids]
    all_metrics <- list()
    for (id in sorted_ids) {
        rmets <- metrics_by_route[[id]]
        for (mname in names(rmets)) {
            m <- rmets[[mname]]
            label <- paste(id, mname, sep = "::")
            all_metrics[[label]] <- m
        }
    }
    tradeoffs <- sclet_ai_comparison_tradeoffs(sorted_ids, metrics_by_route, baseline)
    warnings <- list("This comparison reuses recorded summaries and does not execute integration.")
    if (identical(criterion, "stability")) {
        stab_values <- vapply(sorted_ids, function(id) {
            sclet_ai_comparison_metric_value(metrics_by_route[[id]], "cluster_stability")
        }, numeric(1L))
        if (all(is.na(stab_values))) {
            warnings <- c(warnings, list("cluster_stability metric not available for any route; ordering falls back to evidence order."))
        }
    }
    if (!length(tradeoffs)) {
        tradeoffs <- list(
            status = "requires_metric_alignment",
            message = "Compare metrics across routes; no route is automatically selected.",
            criterion = criterion
        )
    }
    list(
        status = "available",
        criterion = criterion,
        routes = routes,
        baseline = baseline,
        metrics = all_metrics,
        tradeoffs = tradeoffs,
        recommendation = NULL,
        fingerprints = list(ledger = collected$ledger$fingerprint),
        warnings = warnings,
        execution = list(allowed = FALSE, performed = FALSE)
    )
}
