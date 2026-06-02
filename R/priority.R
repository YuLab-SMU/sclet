#' Run perturbation priority analysis with Augur
#'
#' Prioritize cell types by their responsiveness to experimental perturbation
#' using the Augur framework. Augur trains a machine-learning model for each
#' cell type to predict experimental condition from molecular measurements, then
#' ranks cell types by cross-validated AUC.
#'
#' @param sce A SingleCellExperiment object.
#' @param condition_col Character. Column in `colData(sce)` that contains
#'   condition labels (e.g. treatment vs control).
#' @param group_col Character. Column in `colData(sce)` that contains cell
#'   type or cluster labels. If `NULL`, defaults to `ActiveIdent(sce)`.
#'   The special value `"colLabels"` uses `SingleCellExperiment::colLabels(sce)`.
#' @param layer Character. Expression layer to use for the Augur analysis.
#'   If `NULL`, uses the active layer and avoids scaled layers.
#' @param assay Character. Assay to use if `layer` is not specified.
#' @param name Character. Analysis record id. Defaults to `"augur"`.
#' @param ... Additional arguments passed to `Augur::calculate_auc()`.
#'
#' @return A `SingleCellExperiment` object with perturbation priority results
#'   registered in the analysis-state layer. The raw Augur result list is also
#'   stored under `metadata(sce)$sclet$analyses$augur` for backwards
#'   compatibility.
#' @export
RunPerturbationPriority <- function(sce, condition_col, group_col = NULL,
    layer = NULL, assay = NULL, name = "augur", ...) {
    if (!requireNamespace("Augur", quietly = TRUE)) {
        stop(
            "Package 'Augur' is needed for this function to work. ",
            "Please install it with: ",
            'remotes::install_github("neurorestore/Augur")',
            call. = FALSE
        )
    }
    if (!is(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment object.")
    }
    if (is.null(condition_col) || !nzchar(condition_col) ||
        !condition_col %in% colnames(SummarizedExperiment::colData(sce))) {
        stop(
            "A valid `condition_col` must be provided that exists in colData(sce)."
        )
    }
    cd <- SummarizedExperiment::colData(sce)

    if (is.null(group_col)) {
        group_col <- ActiveIdent(sce)
    }
    if (!is.null(group_col) && identical(group_col, "colLabels")) {
        labels <- SingleCellExperiment::colLabels(sce)
        sce$sclet_priority_group_temp <- as.character(labels)
        group_col <- "sclet_priority_group_temp"
    }
    if (is.null(group_col) || !nzchar(group_col) ||
        !group_col %in% colnames(cd)) {
        stop(
            "A valid `group_col` must be provided that exists in colData(sce)."
        )
    }
    condition_values <- unique(as.character(cd[[condition_col]]))
    if (length(condition_values) < 2L) {
        stop("`condition_col` must contain at least two distinct condition levels.")
    }

    src <- sclet_resolve_expression_source(
        sce, layer = layer, assay = assay,
        prefer_nonscaled = TRUE, fallback_layers = c("logcounts"),
        context = "perturbation priority"
    )

    expr_matrix <- SummarizedExperiment::assay(sce, src$assay)
    meta <- as.data.frame(cd)

    cli::cli_alert_info("Running Augur perturbation priority on {.val {group_col}} by {.val {condition_col}}")
    augur_result <- Augur::calculate_auc(
        input = expr_matrix,
        meta = meta,
        cell_type_col = group_col,
        label_col = condition_col,
        ...
    )

    sce <- sclet_set_analysis(
        sce,
        "augur",
        list(
            id = name,
            condition_col = condition_col,
            group_col = group_col,
            assay = src$assay,
            AUC = augur_result$AUC,
            result = augur_result,
            timestamp = Sys.time()
        )
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "priority",
        id = name,
        method = "Augur",
        inputs = list(
            condition_col = condition_col,
            group_col = group_col,
            assay = src$assay,
            layer = src$layer
        ),
        artifacts = list(
            analysis_key = "augur",
            AUC = augur_result$AUC
        ),
        params = list(...)
    )
    sce <- sclet_log_command(
        sce,
        "RunPerturbationPriority",
        params = list(
            condition_col = condition_col,
            group_col = group_col,
            assay = src$assay,
            name = name
        ),
        outputs = list(
            analysis = "augur",
            state = "priority"
        )
    )

    sce
}

#' Run rare cell detection
#'
#' Reserve interface for rare-cell identification workflows. Currently a
#' placeholder that records the analysis intent so the book documentation
#' and API surface can be built around this mainline before the backend is
#' finalized.
#'
#' @param sce A SingleCellExperiment object.
#' @param method Reserved for future backend selection.
#' @param name Analysis record id. Defaults to `"rareq"`.
#' @param ... Reserved for future parameters.
#'
#' @return Updated `SingleCellExperiment` with a placeholder priority record.
#' @export
RunRareCellDetection <- function(sce, method = NULL, name = "rareq", ...) {
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    message(
        "RunRareCellDetection() is currently a reserved interface. ",
        "Rare cell detection backends will be added in a future release."
    )

    placeholder <- list(
        id = name,
        method = if (is.null(method)) "reserved" else method,
        status = "reserved",
        timestamp = Sys.time()
    )

    sce <- sclet_set_analysis(sce, "rare_cells", placeholder)
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "rare_cells",
        id = name,
        method = if (is.null(method)) "reserved" else method,
        inputs = list(),
        artifacts = list(analysis_key = "rare_cells"),
        summary = list(status = "reserved")
    )
    sce <- sclet_log_command(
        sce,
        "RunRareCellDetection",
        params = list(name = name),
        outputs = list(analysis = "rare_cells", state = "rare_cells")
    )

    sce
}

#' Run a state priority and perturbation sensitivity workflow
#'
#' This is the semantic shell over the state-priority mainline. It currently
#' dispatches perturbation priority analysis through Augur, and other
#' prioritization modules will be added as backends mature.
#'
#' @param sce A SingleCellExperiment object.
#' @param condition_col Character. Column in colData containing condition labels.
#' @param group_col Character. Column in colData containing cell type labels.
#' @param run_perturbation_priority Logical. Whether to run Augur priority.
#' @param run_rare_cell_detection Logical. Whether to run rare cell detection.
#' @param name Workflow analysis id. Defaults to `"state_priority"`.
#' @param ... Additional arguments passed to the selected backends.
#'
#' @return Updated `SingleCellExperiment` with state priority workflow record.
#' @export
RunStatePriorityWorkflow <- function(sce, condition_col = NULL,
    group_col = NULL, run_perturbation_priority = TRUE,
    run_rare_cell_detection = FALSE, name = "state_priority", ...) {

    component_ids <- list()

    if (isTRUE(run_perturbation_priority)) {
        priority_id <- paste0(name, "_priority")
        sce <- RunPerturbationPriority(
            sce,
            condition_col = condition_col,
            group_col = group_col,
            name = priority_id,
            ...
        )
        component_ids$perturbation_priority <- priority_id
    }

    if (isTRUE(run_rare_cell_detection)) {
        rare_id <- paste0(name, "_rare")
        sce <- RunRareCellDetection(sce, name = rare_id)
        component_ids$rare_cells <- rare_id
    }

    if (!length(component_ids)) {
        stop(
            "RunStatePriorityWorkflow() needs at least one ",
            "prioritization task enabled."
        )
    }

    workflow_analysis <- list(
        id = name,
        method = "state_priority_mainline",
        inputs = list(
            condition_col = condition_col,
            group_col = group_col
        ),
        artifacts = component_ids,
        summary = list(
            has_perturbation_priority = !is.null(component_ids$perturbation_priority),
            has_rare_cells = !is.null(component_ids$rare_cells)
        ),
        created_at = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "state_priority", workflow_analysis)
    sce <- sclet_log_command(
        sce,
        "RunStatePriorityWorkflow",
        params = list(
            condition_col = condition_col,
            group_col = group_col,
            run_perturbation_priority = run_perturbation_priority,
            run_rare_cell_detection = run_rare_cell_detection,
            name = name
        ),
        outputs = list(
            workflow = name,
            perturbation_priority = component_ids$perturbation_priority,
            rare_cells = component_ids$rare_cells
        )
    )

    sce
}
