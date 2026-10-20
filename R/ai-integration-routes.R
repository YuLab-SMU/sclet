sclet_ai_route_metric <- function(name, value, direction, method = "recorded_summary",
                                     baseline_ref = NULL, evidence_id = NULL,
                                     scope = list(), parameters = list(),
                                     status = NULL, reason = NULL) {
    out <- list(
        name = name,
        value = value,
        direction = direction,
        scope = scope,
        method = method,
        parameters = parameters,
        baseline_ref = baseline_ref,
        uncertainty = list(status = status %||% "not_reported"),
        evidence_id = evidence_id
    )
    if (!is.null(status)) out$status <- status
    if (!is.null(reason)) out$uncertainty$reason <- reason
    out
}

sclet_ai_compute_pca_eta_squared <- function(object, group_col, max_pcs = 10L,
                                              reduction = NULL) {
    if (is.null(group_col)) return(NA_real_)
    reduction_names <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) NULL)
    if (!length(reduction_names)) return(NA_real_)
    if (!is.null(reduction)) {
        if (!reduction %in% reduction_names) return(NA_real_)
        pca_name <- reduction
    } else {
        pca_idx <- grep("pca", tolower(reduction_names), ignore.case = TRUE)[1L]
        if (is.na(pca_idx)) return(NA_real_)
        pca_name <- reduction_names[[pca_idx]]
    }
    pca <- tryCatch(SingleCellExperiment::reducedDim(object, pca_name), error = function(e) NULL)
    if (is.null(pca)) return(NA_real_)
    if (!is.matrix(pca) && !inherits(pca, "Matrix")) return(NA_real_)
    n_use <- min(max_pcs, ncol(pca))
    if (n_use < 1L) return(NA_real_)
    cd <- SummarizedExperiment::colData(object)
    if (!(group_col %in% colnames(cd))) return(NA_real_)
    groups <- as.character(cd[[group_col]])
    if (length(groups) != nrow(pca)) return(NA_real_)
    eta_values <- vapply(seq_len(n_use), function(i) {
        sclet_ai_diag_eta_squared(as.numeric(pca[, i]), groups)
    }, numeric(1))
    mean(eta_values, na.rm = TRUE)
}

# The embedding a route is judged on.
#
# A route's own output is the thing that must be measured: fastMNN writes
# "corrected", Harmony "HARMONY", scVI "scVI", and the raw baseline has no
# integration record at all, so it is judged on the PCA. Before this existed,
# the metric always read the reduction whose name matched "pca" -- which on a
# corrected object is still the UNCORRECTED PCA. Measured on Kang IFN-beta
# (2026-10-20, 5,000 cells): the registered `batch_mixing` was 0.8778 for the
# raw baseline, for fastMNN and for Harmony alike, while the batch eta^2 on
# each route's own embedding was 0.1222 / 0.0064 / 0.0427 -- a 19-fold
# difference the comparison reported as "no route changed anything".
sclet_ai_route_reduction <- function(object, route_id) {
    available <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    if (!length(available)) return(NULL)
    record <- tryCatch(
        sclet_get_state_record(object, "integration", tolower(route_id)),
        error = function(e) NULL
    )
    candidate <- record$artifacts$reduction %||% NULL
    if (is.character(candidate) && length(candidate) == 1L &&
        candidate %in% available) {
        return(candidate)
    }
    # The baseline route has no integration record: it is judged on the PCA.
    pca_idx <- grep("pca", tolower(available), ignore.case = TRUE)[1L]
    if (!is.na(pca_idx)) return(available[[pca_idx]])
    NULL
}

sclet_ai_compute_route_metrics <- function(object, route_id, baseline_ref = "raw", design = NULL, runtime_sec = NA_real_) {
    metrics <- list()
    batch <- design$batch
    condition <- design$condition
    if (is.null(design)) {
        return(list(
            sclet_ai_route_metric("batch_mixing", NA_real_, "higher_is_better",
                status = "not_available", reason = "missing_design_batch",
                baseline_ref = baseline_ref),
            sclet_ai_route_metric("biological_preservation", NA_real_, "higher_is_better",
                status = "not_available", reason = "missing_design",
                baseline_ref = baseline_ref),
            sclet_ai_route_metric("cluster_stability", NA_real_, "higher_is_better",
                status = "not_available",
                reason = "stability_resampling_skipped_in_phase_c",
                baseline_ref = baseline_ref),
            sclet_ai_route_metric("runtime_sec", runtime_sec, "lower_is_better",
                method = "wall_clock", baseline_ref = baseline_ref)
        ))
    }
    # The embedding this route is judged on. Recorded in every metric's
    # parameters, because a number without the embedding it came from is not
    # reproducible: "batch_mixing = 0.99" is a different measurement when it
    # is read off the corrected embedding than off the raw PCA.
    reduction <- sclet_ai_route_reduction(object, route_id)
    eta_batch <- if (!is.null(batch)) {
        sclet_ai_compute_pca_eta_squared(object, batch, reduction = reduction)
    } else NA_real_
    if (is.finite(eta_batch)) {
        batch_metric <- sclet_ai_route_metric(
            "batch_mixing",
            unname(1 - eta_batch),
            "higher_is_better",
            method = "embedding_eta_squared_complement",
            parameters = list(max_pcs = 10L, batch_column = "anonymized",
                reduction = reduction %||% NA_character_),
            baseline_ref = baseline_ref
        )
    } else {
        batch_metric <- sclet_ai_route_metric(
            "batch_mixing", NA_real_, "higher_is_better",
            status = "not_available", reason = "pca_or_batch_unavailable",
            parameters = list(max_pcs = 10L, reduction = reduction %||% NA_character_),
            baseline_ref = baseline_ref
        )
    }
    metrics[[length(metrics) + 1L]] <- batch_metric
    eta_condition <- if (!is.null(condition)) {
        sclet_ai_compute_pca_eta_squared(object, condition, reduction = reduction)
    } else NA_real_
    if (is.finite(eta_condition)) {
        bio_metric <- sclet_ai_route_metric(
            "biological_preservation",
            unname(eta_condition),
            "higher_is_better",
            method = "embedding_eta_squared",
            parameters = list(max_pcs = 10L, condition_column = "anonymized",
                reduction = reduction %||% NA_character_),
            baseline_ref = baseline_ref
        )
    } else {
        bio_metric <- sclet_ai_route_metric(
            "biological_preservation", NA_real_, "higher_is_better",
            status = "not_available", reason = "condition_unavailable",
            parameters = list(max_pcs = 10L, reduction = reduction %||% NA_character_),
            baseline_ref = baseline_ref
        )
    }
    metrics[[length(metrics) + 1L]] <- bio_metric
    stability_metric <- sclet_ai_route_metric(
        "cluster_stability", NA_real_, "higher_is_better",
        status = "not_available",
        reason = "stability_resampling_skipped_in_phase_c",
        baseline_ref = baseline_ref
    )
    metrics[[length(metrics) + 1L]] <- stability_metric
    runtime_metric <- sclet_ai_route_metric(
        "runtime_sec", runtime_sec, "lower_is_better",
        method = "wall_clock", baseline_ref = baseline_ref
    )
    metrics[[length(metrics) + 1L]] <- runtime_metric
    metrics
}

sclet_ai_register_baseline_route <- function(object, design) {
    metrics <- sclet_ai_compute_route_metrics(object, "raw", baseline_ref = "raw",
        design = design, runtime_sec = 0)
    id <- "raw"
    object <- sclet_set_analysis_state(
        object = object,
        type = "integration",
        id = id,
        method = "baseline_route",
        inputs = list(route = "raw", design = names(design)),
        summary = list(
            route = "raw",
            baseline = TRUE,
            design_semantics = "confirmed",
            metrics = metrics
        ),
        active = FALSE
    )
    sclet_log_command(
        object,
        "RunIntegrationRoutes_baseline",
        params = list(route = "raw"),
        outputs = list(analysis_id = id)
    )
}

# Attach route metrics to an integration state record.
#
# Two callers need this and they must not diverge: the plan-driven action
# (`run_integration`) and the route runner (`sclet_ai_run_one_route`). The
# metrics are computed on the route's own embedding and stored in the record's
# summary, so a later comparison reads numbers that were measured on the thing
# they describe.
sclet_ai_record_route_metrics <- function(object, route_name, batch = NULL,
                                           runtime_sec = NA_real_,
                                           design = NULL) {
    if (is.null(design)) {
        design <- sclet_ai_design_for_metrics(object, batch)
    }
    record <- tryCatch(
        sclet_get_state_record(object, "integration", tolower(route_name)),
        error = function(e) NULL
    )
    if (is.null(record)) {
        return(object)
    }
    metrics <- sclet_ai_compute_route_metrics(object, tolower(route_name),
        baseline_ref = "raw", design = design, runtime_sec = runtime_sec)
    summary <- record$summary %||% list()
    summary$design_semantics <- "confirmed"
    summary$metrics <- metrics
    sclet_set_analysis_state(
        object,
        type = "integration",
        id = tolower(route_name),
        method = record$method %||% "RunIntegration",
        inputs = record$inputs %||% list(),
        summary = summary,
        artifacts = record$artifacts %||% list(),
        active = FALSE
    )
}

# The design a route's metrics are computed under: the batch column the
# action was given, plus the condition column the user CONFIRMED for the same
# object. The condition role is never guessed from a column name -- it comes
# from the recorded design confirmation, or it is absent and the biological
# preservation metric reports itself as not available.
sclet_ai_design_for_metrics <- function(object, batch = NULL) {
    if (is.null(batch)) {
        # The batch column is not a guess: it is the one the recorded routes
        # were run with (RunIntegration stores it in the route's inputs).
        records <- tryCatch(
            sclet_get_state_records(object, "integration"),
            error = function(e) list())
        for (id in names(records)) {
            recorded_batch <- records[[id]]$inputs$batch %||% NULL
            if (is.character(recorded_batch) && length(recorded_batch) == 1L &&
                nzchar(recorded_batch)) {
                batch <- recorded_batch
                break
            }
        }
    }
    confirmation <- if (is.null(batch)) {
        NULL
    } else {
        sclet_ai_find_design_confirmation(object, batch = batch)
    }
    recorded <- confirmation$summary$design %||% list()
    design <- list()
    if (!is.null(batch) && nzchar(batch)) {
        design$batch <- batch
    } else if (!is.null(recorded$batch)) {
        design$batch <- as.character(recorded$batch)
    }
    if (!is.null(recorded$condition)) {
        design$condition <- as.character(recorded$condition)
    }
    design
}

sclet_ai_run_one_route <- function(object, route, design) {
    started <- Sys.time()
    batch <- design$batch
    route_name <- tolower(route)
    method <- route
    if (identical(route_name, "fastmnn")) method <- "fastMNN"
    if (identical(route_name, "harmony")) method <- "Harmony"
    if (identical(route_name, "scvi")) method <- "scVI"
    res <- RunIntegration(
        object,
        method = method,
        batch = batch,
        name = route_name
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    metrics <- sclet_ai_compute_route_metrics(
        res, route_name, baseline_ref = "raw", design = design,
        runtime_sec = elapsed
    )
    existing <- tryCatch(sclet_get_state_record(res, "integration", route_name), error = function(e) NULL)
    if (!is.null(existing) && is.list(existing)) {
        existing_summary <- existing$summary %||% list()
        existing_summary$design_semantics <- "confirmed"
        existing_summary$metrics <- metrics
        res <- sclet_set_analysis_state(
            res,
            type = "integration",
            id = route_name,
            method = existing$method %||% paste0("RunIntegration_", method),
            inputs = utils::modifyList(existing$inputs %||% list(), list(route = route, design = names(design))),
            summary = existing_summary,
            artifacts = existing$artifacts %||% list(),
            active = FALSE
        )
    }
    res
}

#' Run multiple integration routes under confirmed design and register comparison metrics
#'
#' Run raw baseline plus one or more integration methods and register each route
#' as a state-recorded analysis with bounded metrics. Each non-interactive sessions
#' default to dry-run unless confirm is not "ask".
#'
#' @param object A SingleCellExperiment.
#' @param design Named list with at least batch. Semantics are taken as confirmed by
#'   caller.
#' @param routes Character vector of routes to execute. Default includes
#'   c("raw", "fastMNN").
#' @param confirm Confirmation policy. "ask" in non-interactive sessions returns a
#'   dry-run preview.
#' @param baseline Baseline route id; defaults to "raw".
#' @param ... Reserved for future arguments.
#' @return A list describing object, plan, validation, preview, execution, status, and
#'   report.
#' @export
RunIntegrationRoutes <- function(
    object,
    design,
    routes = c("raw", "fastMNN"),
    confirm = c("ask", "always", "never"),
    baseline = "raw",
    ...
) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    confirm <- match.arg(confirm)
    routes <- unique(as.character(routes))
    valid_routes <- c("raw", "fastMNN", "Harmony", "scVI")
    if (!length(routes) || any(!routes %in% valid_routes)) {
        stop("routes must be a subset of: ", paste(valid_routes, collapse = ", "), call. = FALSE)
    }
    if (missing(design) || !is.list(design)) {
        return(list(
            status = "clarification_required",
            questions = list(
                list(
                    id = "design_batch",
                    text = "design must be a named list identifying at least a $batch column.",
                    required = TRUE
                ),
                list(
                    id = "design_condition",
                    text = "Optionally provide $condition for biological preservation checks.",
                    required = FALSE
                )
            ),
            assumptions = list(assumed_batch = NULL, assumed_condition = NULL),
            blocked_actions = routes,
            object = object,
            execution = list(allowed = FALSE, performed = FALSE)
        ))
    }
    if (is.null(design$batch)) {
        return(list(
            status = "clarification_required",
            questions = list(list(
                id = "design_batch",
                text = "design$batch is required before running integration.",
                required = TRUE
            )),
            assumptions = list(),
            blocked_actions = routes,
            object = object,
            execution = list(allowed = FALSE, performed = FALSE)
        ))
    }
    if (!(design$batch %in% colnames(SummarizedExperiment::colData(object)))) {
        return(list(
            status = "clarification_required",
            questions = list(list(
                id = "design_batch_missing",
                text = paste0("design$batch = '", design$batch, "' is not a colData column."),
                required = TRUE
            )),
            blocked_actions = routes,
            object = object,
            execution = list(allowed = FALSE, performed = FALSE)
        ))
    }
    plan <- list(
        plan_id = paste0("integration_routes_", format(Sys.time(), "%Y%m%d%H%M%S")),
        task = "multi_route_integration",
        routes = routes,
        baseline = baseline,
        design = design,
        requires_confirmation = TRUE
    )
    preview <- list(
        n_routes = length(routes),
        routes = routes,
        baseline = baseline,
        estimated_costs = stats::setNames(lapply(routes, function(r) {
            switch(tolower(r),
                raw = "low",
                fastmnn = "medium",
                harmony = "medium",
                scvi = "high",
                "medium")
        }), routes)
    )
    interactive_confirmed <- FALSE
    if (identical(confirm, "never")) {
        return(list(
            object = object,
            plan = plan,
            validation = list(valid = TRUE, design_confirmed = TRUE),
            preview = preview,
            execution = list(allowed = FALSE, performed = FALSE, n_routes = 0L,
                reason = "confirm_never"),
            status = "cancelled",
            report = list()
        ))
    }
    if (identical(confirm, "ask")) {
        if (!interactive()) {
            return(list(
                object = object,
                plan = plan,
                validation = list(valid = TRUE, design_confirmed = TRUE),
                preview = preview,
                execution = list(allowed = FALSE, performed = FALSE, n_routes = 0L,
                    reason = "non_interactive_dry_run_confirm_ask"),
                status = "dry_run",
                report = list(dry_run_note = "Re-run with confirm=\"always\" to execute.")
            ))
        }
        interactive_confirmed <- TRUE
    }
    if (identical(confirm, "always")) {
        interactive_confirmed <- TRUE
    }
    if (!isTRUE(interactive_confirmed)) {
        return(list(
            object = object,
            plan = plan,
            validation = list(valid = TRUE, design_confirmed = TRUE),
            preview = preview,
            execution = list(allowed = FALSE, performed = FALSE),
            status = "dry_run",
            report = list()
        ))
    }
    current <- object
    n_exec <- 0L
    route_summaries <- list()
    for (route in routes) {
        if (identical(route, "raw")) {
            current <- sclet_ai_register_baseline_route(current, design)
        } else {
            current <- sclet_ai_run_one_route(current, route, design)
        }
        n_exec <- n_exec + 1L
        state_rec <- tryCatch({
            id <- if (identical(route, "raw")) "raw" else tolower(route)
            sclet_get_state_record(current, "integration", id)
        }, error = function(e) NULL)
        route_summaries[[length(route_summaries) + 1L]] <- list(
            route = route,
            id = if (identical(route, "raw")) "raw" else tolower(route),
            summary = state_rec$summary %||% list()
        )
    }
    report <- list(
        routes_executed = vapply(route_summaries, function(x) x$route, character(1)),
        metrics_summary = lapply(route_summaries, function(x) {
            mets <- x$summary$metrics %||% list()
            out <- lapply(mets, function(m) {
                list(name = m$name, value = m$value %||% NA_real_, status = m$status %||% m$uncertainty$status %||% "reported")
            })
            stats::setNames(out, vapply(mets, function(m) m$name, character(1)))
        }),
        tradeoffs_detected = list(),
        recommendation = NULL
    )
    list(
        object = current,
        plan = plan,
        validation = list(valid = TRUE, design_confirmed = TRUE, baseline = baseline),
        preview = preview,
        execution = list(allowed = TRUE, performed = TRUE, n_routes = n_exec),
        status = "completed",
        report = report
    )
}

# The recorded route metrics, in a shape the comparison can read.
sclet_ai_integration_route_metrics <- function(object) {
    records <- tryCatch(
        sclet_get_state_records(object, "integration"),
        error = function(e) list())
    if (!length(records)) {
        return(list())
    }
    out <- list()
    for (id in names(records)) {
        rec <- records[[id]]
        # The comparison record is not a route: it holds the comparison table,
        # not a route's own metrics.
        if (identical(as.character(rec$method %||% ""), "sclet_ai_route_comparison")) next
        metrics <- rec$summary$metrics %||% list()
        if (!length(metrics)) next
        value_of <- function(name) {
            for (m in metrics) {
                if (identical(as.character(m$name %||% ""), name)) {
                    return(m$value %||% NA_real_)
                }
            }
            NA_real_
        }
        out[[id]] <- list(
            id = id,
            route = as.character(rec$summary$route %||% id),
            baseline = isTRUE(rec$summary$baseline),
            batch_mixing = value_of("batch_mixing"),
            biological_preservation = value_of("biological_preservation"),
            runtime_sec = value_of("runtime_sec"),
            reduction = rec$artifacts$reduction %||% NA_character_
        )
    }
    out
}

# Register ONE evidence node for the route comparison.
#
# One node, not one per route: the routes are alternative methods applied to
# the same expression matrix under the same design, so they are one line of
# support, not corroboration. Route names stay in the state records; the node
# carries anonymized `route_N` codes and numbers, because an evidence node is
# allowed to leave the process.
sclet_ai_record_integration_comparison <- function(object, name = "integration_comparison") {
    routes <- sclet_ai_integration_route_metrics(object)
    if (!length(routes)) {
        stop("no_integration_route_recorded: there is no recorded route to compare.",
            call. = FALSE)
    }
    # The baseline is the object's own uncorrected embedding. A recorded
    # "raw" route (RunIntegrationRoutes writes one) is used when it exists;
    # otherwise the baseline is measured here, on the same cells and the same
    # PCA the routes started from, so the comparison has something to compare
    # against even when the plan only ran the integration steps.
    if (!any(vapply(routes, function(r) isTRUE(r$baseline), logical(1L)))) {
        design <- sclet_ai_design_for_metrics(object)
        base_metrics <- sclet_ai_compute_route_metrics(object, "raw",
            baseline_ref = "raw", design = design, runtime_sec = 0)
        value_of <- function(name) {
            for (m in base_metrics) {
                if (identical(as.character(m$name %||% ""), name)) {
                    return(m$value %||% NA_real_)
                }
            }
            NA_real_
        }
        routes[["raw"]] <- list(
            id = "raw", route = "raw", baseline = TRUE,
            batch_mixing = value_of("batch_mixing"),
            biological_preservation = value_of("biological_preservation"),
            runtime_sec = value_of("runtime_sec"),
            reduction = sclet_ai_route_reduction(object, "raw") %||% NA_character_
        )
    }
    # The baseline comes first, so route_1 is always "no correction" and the
    # comparison reads in the order the question is asked.
    baseline_ids <- names(routes)[vapply(routes, function(r) isTRUE(r$baseline), logical(1L))]
    other_ids <- setdiff(names(routes), baseline_ids)
    routes <- routes[c(baseline_ids, other_ids)]
    ids <- names(routes)
    codes <- paste0("route_", seq_along(ids))
    names(codes) <- ids
    pick <- function(field) {
        vapply(ids, function(id) as.numeric(routes[[id]][[field]] %||% NA_real_), numeric(1L))
    }
    batch_mixing <- pick("batch_mixing")
    preservation <- pick("biological_preservation")
    runtime <- pick("runtime_sec")
    baseline_idx <- which(vapply(ids, function(id) isTRUE(routes[[id]]$baseline), logical(1L)))
    if (!length(baseline_idx)) baseline_idx <- 1L
    base_batch <- batch_mixing[baseline_idx[1L]]
    base_preservation <- preservation[baseline_idx[1L]]
    finite_batch <- is.finite(batch_mixing)
    best_batch_idx <- if (any(finite_batch)) which.max(replace(batch_mixing, !finite_batch, -Inf)) else NA_integer_
    # A route is a win only if it improves batch mixing without giving up
    # biological signal. Anything else is a trade-off, and the comparison says
    # so instead of recommending the largest number.
    improved <- if (is.na(best_batch_idx)) {
        FALSE
    } else {
        is.finite(base_batch) && batch_mixing[best_batch_idx] > base_batch
    }
    retained <- if (is.na(best_batch_idx) || !is.finite(base_preservation)) {
        NA
    } else {
        is.finite(preservation[best_batch_idx]) &&
            preservation[best_batch_idx] >= base_preservation
    }
    tradeoff <- if (is.na(best_batch_idx) || !is.finite(base_preservation)) {
        NA
    } else {
        !isTRUE(retained)
    }
    values <- list(
        n_routes = as.integer(length(ids)),
        route_codes = as.character(codes),
        # Aligned with route_codes by position: the code hides which method
        # produced which number, the state record keeps the mapping.
        batch_mixing_per_route = as.numeric(batch_mixing),
        biological_preservation_per_route = as.numeric(preservation),
        runtime_sec_per_route = as.numeric(runtime),
        # A numeric index, not a route code: a code would have to be a string,
        # and an NA string is not a value an evidence node may carry.
        best_batch_mixing_route_index = as.integer(best_batch_idx),
        batch_mixing_improved = improved,
        biological_preservation_retained = retained,
        tradeoff_detected = tradeoff,
        raw_values_included = FALSE
    )
    # The comparison itself is a recorded analysis: its own state record holds
    # the readable table (route names included) and is the evidence's source,
    # so the evidence node can stay anonymous.
    object <- sclet_set_analysis_state(
        object,
        type = "integration",
        id = name,
        method = "sclet_ai_route_comparison",
        inputs = list(routes = as.character(ids)),
        summary = list(
            n_routes = as.integer(length(ids)),
            routes = lapply(ids, function(id) routes[[id]]),
            batch_mixing_improved = improved,
            biological_preservation_retained = retained,
            tradeoff_detected = tradeoff
        ),
        active = FALSE
    )
    evidence <- list(
        id = "ev:integration_comparison",
        kind = "deterministic_summary",
        values = values,
        claim_level = "associated"
    )
    sclet_ai_record_evidence_for_run(object, evidence, source = name,
        parents = as.character(ids), scope = list(route_family = "route_1"))
}
