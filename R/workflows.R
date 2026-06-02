#' Run a trajectory-oriented analysis workflow
#'
#' This high-level entry point organizes the trajectory mainline around the
#' user's analysis goal instead of exposing only individual algorithms. It can
#' orchestrate diffusion map, Slingshot, RNA velocity, and CellRank/CellFate
#' using the existing atomic interfaces in `sclet`.
#'
#' @param sce A SingleCellExperiment object.
#' @param group Grouping variable used for trajectory inference. If `NULL`, uses
#'   `ActiveIdent(sce)`. The special value `"colLabels"` is materialized to a
#'   stable `colData` column for downstream methods that require a column name.
#' @param reduction Reduction used by Slingshot/velocity/fate. If `NULL`, uses
#'   `DefaultReduction(sce)` and falls back to `"PCA"` when available.
#' @param run_diffusion_map Logical. Whether to compute `RunDiffusionMap()`.
#' @param run_slingshot Logical. Whether to run `RunSlingshot()`.
#' @param run_velocity One of `"auto"`, `"never"`, or `"always"`.
#' @param run_fate One of `"auto"`, `"never"`, or `"always"`.
#' @param diffusion_components Number of diffusion components.
#' @param velocity_mode Velocity mode passed to `RunVelocity()`.
#' @param start_cluster,end_cluster Optional lineage root/end clusters for Slingshot.
#' @param name Workflow analysis id. Defaults to `"trajectory_workflow"`.
#' @return Updated SingleCellExperiment object.
#' @export
RunTrajectoryWorkflow <- function(
    sce,
    group = NULL,
    reduction = NULL,
    run_diffusion_map = TRUE,
    run_slingshot = TRUE,
    run_velocity = c("auto", "never", "always"),
    run_fate = c("auto", "never", "always"),
    diffusion_components = 50,
    velocity_mode = c("deterministic", "stochastic", "dynamical"),
    start_cluster = NULL,
    end_cluster = NULL,
    name = "trajectory_workflow"
) {
    run_velocity <- match.arg(run_velocity)
    run_fate <- match.arg(run_fate)
    velocity_mode <- match.arg(velocity_mode)

    reduction <- sclet_workflow_resolve_reduction(sce, reduction)
    group_info <- sclet_workflow_resolve_group(sce, group = group, prefix = name)
    sce <- group_info$object
    group_col <- group_info$group

    component_ids <- list()

    if (isTRUE(run_diffusion_map)) {
        sce <- RunDiffusionMap(sce, n_components = diffusion_components, reduction = reduction)
        component_ids$diffusion_map <- "DM"
    }

    if (isTRUE(run_slingshot)) {
        slingshot_id <- paste0(name, "_slingshot")
        sce <- RunSlingshot(
            sce,
            group = group_col,
            reduction = reduction,
            start_cluster = start_cluster,
            end_cluster = end_cluster,
            name = slingshot_id
        )
        component_ids$trajectory <- slingshot_id
    }

    should_run_velocity <- switch(
        run_velocity,
        auto = all(c("spliced", "unspliced") %in% SummarizedExperiment::assayNames(sce)),
        never = FALSE,
        always = TRUE
    )
    if (isTRUE(should_run_velocity)) {
        velocity_id <- paste0(name, "_velocity")
        sce <- RunVelocity(
            sce,
            mode = velocity_mode,
            use.dimred = reduction,
            name = velocity_id
        )
        component_ids$velocity <- velocity_id
    } else if (identical(run_velocity, "always")) {
        stop("Assays 'spliced' and 'unspliced' are required when `run_velocity = \"always\"`.")
    }

    should_run_fate <- switch(
        run_fate,
        auto = has_velocity(sce),
        never = FALSE,
        always = TRUE
    )
    if (isTRUE(should_run_fate)) {
        if (!has_velocity(sce)) {
            stop("Velocity results are required to run cell fate inference.")
        }
        fate_id <- paste0(name, "_fate")
        sce <- RunCellFate(
            sce,
            reduction = reduction,
            cluster_key = group_col,
            name = fate_id
        )
        component_ids$fate <- fate_id
    }

    workflow_analysis <- list(
        id = name,
        method = "trajectory_mainline",
        inputs = list(
            group = group_col,
            reduction = reduction
        ),
        artifacts = component_ids,
        summary = list(
            has_diffusion_map = isTRUE(run_diffusion_map),
            has_slingshot = isTRUE(run_slingshot),
            has_velocity = has_velocity(sce),
            has_fate = !is.null(component_ids$fate)
        ),
        created_at = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "trajectory_workflow", workflow_analysis)
    sce <- sclet_log_command(
        sce,
        "RunTrajectoryWorkflow",
        params = list(
            group = group_col,
            reduction = reduction,
            run_diffusion_map = run_diffusion_map,
            run_slingshot = run_slingshot,
            run_velocity = run_velocity,
            run_fate = run_fate,
            diffusion_components = diffusion_components,
            velocity_mode = velocity_mode,
            start_cluster = start_cluster,
            end_cluster = end_cluster,
            name = name
        ),
        outputs = list(
            workflow = name,
            trajectory = component_ids$trajectory,
            velocity = component_ids$velocity,
            fate = component_ids$fate
        )
    )

    sce
}

#' Run a program/regulon interpretation workflow
#'
#' This high-level entry point organizes mechanism-oriented analyses around the
#' user's question, combining pathway/signature scoring and regulon activity in
#' a single workflow shell.
#'
#' @param sce A SingleCellExperiment object.
#' @param gene_sets Optional named gene-set list for `RunGeneSetScoring()`.
#' @param scoring_method Scoring backend for `RunGeneSetScoring()`.
#' @param assay_use Assay used for gene-set scoring and SCENIC.
#' @param scoring_as_altExp Logical. Whether to store gene-set scores in an altExp.
#' @param run_scenic One of `"auto"`, `"never"`, or `"always"`.
#' @param tfs_path,motif_annotations_path,database_paths Inputs for `RunSCENIC()`.
#' @param name Workflow analysis id.
#' @return Updated SingleCellExperiment object.
#' @export
RunProgramWorkflow <- function(
    sce,
    gene_sets = NULL,
    scoring_method = c("UCell", "AUCell", "GSVA"),
    assay_use = "counts",
    scoring_as_altExp = FALSE,
    run_scenic = c("auto", "never", "always"),
    tfs_path = NULL,
    motif_annotations_path = NULL,
    database_paths = NULL,
    name = "program_workflow"
) {
    scoring_method <- match.arg(scoring_method)
    run_scenic <- match.arg(run_scenic)

    component_ids <- list()

    if (!is.null(gene_sets)) {
        scoring_id <- paste0(name, "_scores")
        sce <- RunGeneSetScoring(
            sce,
            gene_sets = gene_sets,
            method = scoring_method,
            assay_use = assay_use,
            name = scoring_id,
            as_altExp = scoring_as_altExp
        )
        component_ids$geneset_scoring <- scoring_id
    }

    scenic_inputs_ready <- !is.null(tfs_path) &&
        !is.null(motif_annotations_path) &&
        !is.null(database_paths)
    should_run_scenic <- switch(
        run_scenic,
        auto = scenic_inputs_ready,
        never = FALSE,
        always = TRUE
    )
    if (isTRUE(should_run_scenic)) {
        if (!isTRUE(scenic_inputs_ready)) {
            stop("`tfs_path`, `motif_annotations_path`, and `database_paths` are required when SCENIC is enabled.")
        }
        scenic_id <- paste0(name, "_scenic")
        sce <- RunSCENIC(
            sce,
            tfs_path = tfs_path,
            motif_annotations_path = motif_annotations_path,
            database_paths = database_paths,
            assay_use = assay_use,
            name = scenic_id
        )
        component_ids$scenic <- scenic_id
    }

    if (!length(component_ids)) {
        stop("RunProgramWorkflow() needs at least one program-level task: gene-set scoring and/or SCENIC.")
    }

    workflow_analysis <- list(
        id = name,
        method = "program_mainline",
        inputs = list(
            assay = assay_use,
            gene_sets = if (is.null(gene_sets)) NULL else names(gene_sets)
        ),
        artifacts = component_ids,
        summary = list(
            has_geneset_scoring = !is.null(component_ids$geneset_scoring),
            has_scenic = !is.null(component_ids$scenic)
        ),
        created_at = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "program_workflow", workflow_analysis)
    sce <- sclet_log_command(
        sce,
        "RunProgramWorkflow",
        params = list(
            scoring_method = scoring_method,
            assay_use = assay_use,
            scoring_as_altExp = scoring_as_altExp,
            run_scenic = run_scenic,
            has_gene_sets = !is.null(gene_sets),
            name = name
        ),
        outputs = list(
            workflow = name,
            geneset_scoring = component_ids$geneset_scoring,
            scenic = component_ids$scenic
        )
    )

    sce
}

#' Run a reference mapping workflow
#'
#' This high-level entry point organizes reference annotation and label transfer
#' into a single stable workflow shell.
#'
#' @param object A SingleCellExperiment object.
#' @param ref Reference dataset.
#' @param labels Labels in the reference dataset.
#' @param method Mapping backend. One of `"SingleR"` or `"KNN"`.
#' @param assay.type Assay used for `SingleR`.
#' @param layer Layer used for mapping.
#' @param features Features used for `KNN`.
#' @param k Number of neighbors for `KNN`.
#' @param name Workflow analysis id.
#' @return Updated SingleCellExperiment object.
#' @export
RunReferenceWorkflow <- function(
    object,
    ref,
    labels = NULL,
    method = c("SingleR", "KNN"),
    assay.type = NULL,
    layer = NULL,
    features = NULL,
    k = 5,
    name = "reference_workflow"
) {
    method <- match.arg(method)
    mapping_id <- paste0(name, "_mapping")

    object <- RunReferenceMapping(
        object = object,
        ref = ref,
        labels = labels,
        method = method,
        assay.type = assay.type,
        layer = layer,
        features = features,
        k = k,
        name = mapping_id
    )

    mapping_record <- get_mapping(object, id = mapping_id)
    workflow_analysis <- list(
        id = name,
        method = "reference_mainline",
        inputs = list(
            backend = method,
            assay = assay.type,
            layer = layer
        ),
        artifacts = list(
            mapping = mapping_id,
            labels_col = mapping_record$artifacts$labels_col,
            score_col = mapping_record$artifacts$score_col
        ),
        summary = mapping_record$summary,
        created_at = Sys.time()
    )
    object <- sclet_set_analysis(object, "reference_workflow", workflow_analysis)
    object <- sclet_log_command(
        object,
        "RunReferenceWorkflow",
        params = list(
            method = method,
            assay.type = assay.type,
            layer = layer,
            features = if (is.null(features)) NULL else length(features),
            k = k,
            name = name
        ),
        outputs = list(
            workflow = name,
            mapping = mapping_id
        )
    )

    object
}

#' Plot a trajectory workflow overview
#'
#' Compose the key panels from a trajectory workflow into a single overview so
#' users can follow the mainline directly instead of manually assembling
#' lineage, velocity, and fate plots one by one.
#'
#' @param object A SingleCellExperiment object.
#' @param id Workflow analysis id. Defaults to the stored `"trajectory_workflow"` analysis.
#' @param reduction Reduction used for plotting. Defaults to the reduction
#'   recorded in the workflow.
#' @param group Grouping column used for lineage plotting. Defaults to the group
#'   recorded in the workflow.
#' @param panels Character vector selecting panels from `"lineage"`,
#'   `"velocity"`, `"terminal_states"`, and `"fate_probability"`.
#' @return A combined plot object from `aplot::plot_list()`.
#' @export
plot_trajectory_overview <- function(
    object,
    id = NULL,
    reduction = NULL,
    group = NULL,
    panels = c("lineage", "velocity", "terminal_states", "fate_probability")
) {
    if (!requireNamespace("aplot", quietly = TRUE)) {
        stop("Package 'aplot' is needed for this function to work. Please install it.")
    }

    workflow <- sclet_get_analysis(object, "trajectory_workflow")
    if (is.null(workflow)) {
        stop("No trajectory workflow record found. Please run RunTrajectoryWorkflow() first.")
    }
    if (!is.null(id) && !identical(workflow$id, id)) {
        stop("Trajectory workflow id '", id, "' not found in object.")
    }

    reduction <- reduction %||% workflow$inputs$reduction
    group <- group %||% workflow$inputs$group
    panels <- match.arg(
        panels,
        choices = c("lineage", "velocity", "terminal_states", "fate_probability"),
        several.ok = TRUE
    )

    plot_list <- list()
    trajectory_id <- workflow$artifacts$trajectory
    velocity_id <- workflow$artifacts$velocity
    fate_id <- workflow$artifacts$fate

    if ("lineage" %in% panels && !is.null(trajectory_id) && has_trajectory(object, id = trajectory_id)) {
        if (is.null(group)) {
            stop("No grouping column recorded for the trajectory workflow. Please provide `group`.")
        }
        plot_list$Lineage <- lineage_plot(
            object,
            group = group,
            reduction = reduction,
            mapping = ggplot2::aes(color = .data[[group]])
        )
    }

    if ("velocity" %in% panels && !is.null(velocity_id) && has_velocity(object, id = velocity_id)) {
        plot_list$Velocity <- VelocityPlot(
            object,
            reduction = reduction,
            group.by = group,
            id = velocity_id
        )
    }

    if ("terminal_states" %in% panels && !is.null(fate_id) && has_cellrank(object, id = fate_id)) {
        plot_list$`Terminal States` <- plot_fate_terminal_states(
            object,
            reduction = reduction,
            id = fate_id
        )
    }

    if ("fate_probability" %in% panels && !is.null(fate_id) && has_cellrank(object, id = fate_id)) {
        cellrank <- get_cellrank(object, id = fate_id)
        fate_name <- NULL
        if (!is.null(cellrank$artifacts$fate_probability_names)) {
            fate_name <- cellrank$artifacts$fate_probability_names[[1]]
        } else if (!is.null(cellrank$fate_probability_names)) {
            fate_name <- cellrank$fate_probability_names[[1]]
        }
        if (!is.null(fate_name)) {
            plot_list$`Fate Probability` <- plot_fate_probability(
                object,
                fate = fate_name,
                reduction = reduction,
                id = fate_id
            )
        }
    }

    if (!length(plot_list)) {
        stop("No trajectory workflow panels are available to plot.")
    }

    aplot::plot_list(gglist = plot_list, ncol = 2)
}

sclet_workflow_resolve_reduction <- function(object, reduction = NULL) {
    if (!is.null(reduction)) {
        return(reduction)
    }
    reduction <- DefaultReduction(object)
    if (!is.null(reduction)) {
        return(reduction)
    }
    if ("PCA" %in% SingleCellExperiment::reducedDimNames(object)) {
        return("PCA")
    }
    stop("No active reduction found. Please provide `reduction` or run RunPCA()/RunUMAP() first.")
}

sclet_workflow_resolve_group <- function(object, group = NULL, prefix = "workflow") {
    if (is.null(group)) {
        group <- ActiveIdent(object)
    }
    if (is.null(group)) {
        stop("No active identity found. Please provide `group`.")
    }
    if (identical(group, "colLabels")) {
        groups <- SingleCellExperiment::colLabels(object)
        if (is.null(groups)) {
            stop("Active identity is 'colLabels' but `colLabels(object)` is empty.")
        }
        group_col <- paste0(prefix, "_group")
        SummarizedExperiment::colData(object)[[group_col]] <- groups
        return(list(object = object, group = group_col))
    }
    if (!group %in% colnames(SummarizedExperiment::colData(object))) {
        stop(sprintf("Grouping column '%s' not found in colData(object).", group))
    }
    list(object = object, group = group)
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
