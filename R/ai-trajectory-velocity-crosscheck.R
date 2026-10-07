#' Compare trajectory and velocity pseudotime evidence summaries
#'
#' This function reads only the bounded evidence nodes emitted by completed
#' trajectory and velocity runs (`ev:trajectory_<name>` and
#' `ev:velocity_<name>`). It compares their pseudotime spread summaries
#' (relative IQR / max ratios) and reports whether the two estimates are of
#' similar magnitude. It never reads raw per-cell `colData` columns, never
#' claims either pseudotime is correct, and never writes a new evidence node.
#'
#' @param object A `SingleCellExperiment` object.
#' @param trajectory_id Optional analysis-run name for the trajectory evidence
#'   node (matching the `name` parameter of `run_trajectory`). When `NULL`, the
#'   single available trajectory run is used; if more than one exists, the
#'   function returns `status = "ambiguous_run"`.
#' @param velocity_id Optional analysis-run name for the velocity evidence node
#'   (matching the `name` parameter of `run_velocity`). When `NULL`, the single
#'   available velocity run is used; if more than one exists, the function
#'   returns `status = "ambiguous_run"`.
#' @return A bounded list with `status`, spread comparisons, and explicit
#'   caveats. When either evidence node is missing or incomplete, returns
#'   `status = "not_available"` with a typed `reason`. When multiple runs of
#'   either type exist and no id was given, returns `status = "ambiguous_run"`
#'   listing the available ids. The result contains no per-cell values and
#'   always sets `raw_values_included = FALSE`.
#' @export
compare_trajectory_velocity_evidence <- function(object, trajectory_id = NULL, velocity_id = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    evidence <- sclet_ai_evidence_get_all(object)
    traj_ids <- grep("^ev:trajectory_", names(evidence), value = TRUE)
    vel_ids <- grep("^ev:velocity_", names(evidence), value = TRUE)

    if (!length(traj_ids) && !length(vel_ids)) {
        return(list(
            status = "not_available",
            reason = "no_trajectory_or_velocity_evidence_records",
            raw_values_included = FALSE
        ))
    }
    if (!length(traj_ids)) {
        return(list(
            status = "not_available",
            reason = "trajectory_evidence_missing",
            available_velocity_runs = sub("^ev:velocity_", "", vel_ids),
            raw_values_included = FALSE
        ))
    }
    if (!length(vel_ids)) {
        return(list(
            status = "not_available",
            reason = "velocity_evidence_missing",
            available_trajectory_runs = sub("^ev:trajectory_", "", traj_ids),
            raw_values_included = FALSE
        ))
    }

    if (is.null(trajectory_id)) {
        if (length(traj_ids) > 1L) {
            return(list(
                status = "ambiguous_run",
                available_trajectory_runs = sub("^ev:trajectory_", "", traj_ids),
                available_velocity_runs = sub("^ev:velocity_", "", vel_ids),
                raw_values_included = FALSE
            ))
        }
        trajectory_id <- sub("^ev:trajectory_", "", traj_ids)
    }
    if (is.null(velocity_id)) {
        if (length(vel_ids) > 1L) {
            return(list(
                status = "ambiguous_run",
                available_trajectory_runs = sub("^ev:trajectory_", "", traj_ids),
                available_velocity_runs = sub("^ev:velocity_", "", vel_ids),
                raw_values_included = FALSE
            ))
        }
        velocity_id <- sub("^ev:velocity_", "", vel_ids)
    }

    traj_node <- evidence[[paste0("ev:trajectory_", trajectory_id)]]
    vel_node <- evidence[[paste0("ev:velocity_", velocity_id)]]

    if (is.null(traj_node)) {
        return(list(
            status = "not_available",
            reason = "trajectory_evidence_missing",
            available_trajectory_runs = sub("^ev:trajectory_", "", traj_ids),
            raw_values_included = FALSE
        ))
    }
    if (is.null(vel_node)) {
        return(list(
            status = "not_available",
            reason = "velocity_evidence_missing",
            available_velocity_runs = sub("^ev:velocity_", "", vel_ids),
            raw_values_included = FALSE
        ))
    }

    traj_values <- traj_node$values
    vel_values <- vel_node$values

    if (is.null(traj_values$first_lineage_median)) {
        return(list(
            status = "not_available",
            reason = "trajectory_lineage_not_available",
            trajectory_run = trajectory_id,
            raw_values_included = FALSE
        ))
    }
    if (is.null(vel_values$velocity_pseudotime)) {
        return(list(
            status = "not_available",
            reason = "velocity_pseudotime_column_not_available",
            velocity_run = velocity_id,
            raw_values_included = FALSE
        ))
    }

    traj_spread <- sclet_ai_tv_spread(
        traj_values$first_lineage_iqr,
        traj_values$first_lineage_max
    )
    vel_spread <- sclet_ai_tv_spread(
        vel_values$velocity_pseudotime$quantile_75 - vel_values$velocity_pseudotime$quantile_25,
        vel_values$velocity_pseudotime$max
    )
    spread_ratio <- sclet_ai_tv_ratio(vel_spread, traj_spread)
    comparable <- is.finite(spread_ratio)

    if (!comparable) {
        label <- "relative spread could not be compared because one or both bounded summaries were degenerate"
    } else if (spread_ratio >= 0.5 && spread_ratio <= 2.0) {
        label <- "relative spread of the two pseudotime estimates is of similar order of magnitude"
    } else {
        label <- "relative spread of the two pseudotime estimates differs by more than 2x; this may reflect genuinely different dynamics captured by each method, not necessarily an error"
    }

    list(
        status = "available",
        trajectory_run = trajectory_id,
        velocity_run = velocity_id,
        trajectory_relative_spread = traj_spread,
        velocity_relative_spread = vel_spread,
        spread_ratio = spread_ratio,
        comparable_spread = comparable,
        spread_consistency_label = label,
        claim_level = "consistent_with",
        caveats = list(
            "Both pseudotime values are relative orderings, not absolute time, and are model- or algorithm-dependent (slingshot lineage ordering vs scVelo-estimated dynamics).",
            "This comparison only checks whether the bounded spread summaries are of similar magnitude; it does not establish per-cell agreement, directional (root-to-tip) agreement, or that either method is biologically correct."
        ),
        raw_values_included = FALSE
    )
}

sclet_ai_tv_spread <- function(iqr, max_val) {
    if (!is.numeric(iqr) || !is.numeric(max_val) || !length(iqr) || !length(max_val)) {
        return(NA_real_)
    }
    if (!is.finite(max_val) || max_val == 0) {
        return(NA_real_)
    }
    as.numeric(iqr) / as.numeric(max_val)
}

sclet_ai_tv_ratio <- function(num, denom) {
    if (!is.finite(num) || !is.finite(denom) || denom == 0) {
        return(NA_real_)
    }
    num / denom
}
