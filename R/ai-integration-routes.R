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

sclet_ai_compute_pca_eta_squared <- function(object, group_col, max_pcs = 10L) {
    if (is.null(group_col)) return(NA_real_)
    reduction_names <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) NULL)
    if (!length(reduction_names)) return(NA_real_)
    pca_idx <- grep("pca", tolower(reduction_names), ignore.case = TRUE)[1L]
    if (is.na(pca_idx)) return(NA_real_)
    pca_name <- reduction_names[[pca_idx]]
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
    eta_batch <- if (!is.null(batch)) {
        sclet_ai_compute_pca_eta_squared(object, batch)
    } else NA_real_
    if (is.finite(eta_batch)) {
        batch_metric <- sclet_ai_route_metric(
            "batch_mixing",
            unname(1 - eta_batch),
            "higher_is_better",
            method = "pca_eta_squared_complement",
            parameters = list(max_pcs = 10L, batch_column = "anonymized"),
            baseline_ref = baseline_ref
        )
    } else {
        batch_metric <- sclet_ai_route_metric(
            "batch_mixing", NA_real_, "higher_is_better",
            status = "not_available", reason = "pca_or_batch_unavailable",
            baseline_ref = baseline_ref
        )
    }
    metrics[[length(metrics) + 1L]] <- batch_metric
    eta_condition <- if (!is.null(condition)) {
        sclet_ai_compute_pca_eta_squared(object, condition)
    } else NA_real_
    if (is.finite(eta_condition)) {
        bio_metric <- sclet_ai_route_metric(
            "biological_preservation",
            unname(eta_condition),
            "higher_is_better",
            method = "pca_eta_squared",
            parameters = list(max_pcs = 10L, condition_column = "anonymized"),
            baseline_ref = baseline_ref
        )
    } else {
        bio_metric <- sclet_ai_route_metric(
            "biological_preservation", NA_real_, "higher_is_better",
            status = "not_available", reason = "condition_unavailable",
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
