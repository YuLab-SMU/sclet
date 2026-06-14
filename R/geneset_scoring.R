#' Run Gene Set Scoring
#'
#' This function calculates gene set or pathway activity scores for individual cells using 
#' methods like AUCell, UCell, or GSVA. The resulting scores are added to the `colData` of the 
#' SingleCellExperiment object and registered in the analysis-state contract.
#'
#' @param sce A SingleCellExperiment object.
#' @param gene_sets A named list of character vectors containing gene symbols.
#' @param features Deprecated alias of `gene_sets`, kept for backward compatibility.
#' @param method Character. The scoring method to use. Options are "AUCell", "UCell", or "GSVA". Defaults to "UCell".
#' @param assay_use Character. The assay to use for scoring. Defaults to "counts" (often preferred for AUCell/UCell) or "logcounts".
#' @param ncores Integer. Number of cores for parallel processing. Defaults to 1.
#' @param name Character. Name of the scoring run, used as prefix for colData columns. Defaults to "Score".
#' @param as_altExp Logical. If TRUE, stores scores as an alternative experiment (altExp) instead of colData. Useful for differential pathway testing. Defaults to FALSE.
#' @param ... Additional arguments passed to the underlying scoring function.
#' 
#' @return A SingleCellExperiment object with gene set scores added to `colData(sce)` and state registered.
#' 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay colData colData<-
#' @export
RunGeneSetScoring <- function(sce, gene_sets = NULL, features = NULL,
                              method = c("UCell", "AUCell", "GSVA"), 
                              assay_use = "counts", ncores = 1, name = "Score", 
                              as_altExp = FALSE, ...) {
    
    method <- match.arg(method)

    if (is.null(gene_sets) && !is.null(features)) {
        gene_sets <- features
        warning("`features` is deprecated; please use `gene_sets` instead.", call. = FALSE)
    }

    if (!is.null(gene_sets) && !is.null(features)) {
        stop("Please supply only one of `gene_sets` or `features`.")
    }

    if (is.null(gene_sets)) {
        stop("`gene_sets` must be provided as a named list of character vectors.")
    }
    
    if (!is.list(gene_sets) || is.null(names(gene_sets))) {
        stop("gene_sets must be a named list of character vectors.")
    }
    
    if (!assay_use %in% assayNames(sce)) {
        stop(sprintf("Assay '%s' not found in the SingleCellExperiment object.", assay_use))
    }
    
    expr_mat <- assay(sce, assay_use)
    expr_mat_for_rank <- NULL
    
    scores <- NULL
    
    if (method == "UCell") {
        if (!requireNamespace("UCell", quietly = TRUE)) {
            stop("Please install 'UCell' package to use method='UCell'.")
        }
        message("Running UCell scoring...")
        expr_mat_for_rank <- sclet_as_dgCMatrix(expr_mat)
        scores <- UCell::ScoreSignatures_UCell(expr_mat_for_rank, features = gene_sets, ncores = ncores, name = "", ...)
        
    } else if (method == "AUCell") {
        if (!requireNamespace("AUCell", quietly = TRUE)) {
            stop("Please install 'AUCell' package to use method='AUCell'.")
        }
        message("Running AUCell scoring...")
        expr_mat_for_rank <- sclet_as_dgCMatrix(expr_mat)
        cells_rankings <- AUCell::AUCell_buildRankings(expr_mat_for_rank, nCores = ncores, plotStats = FALSE, verbose = FALSE)
        cells_AUC <- AUCell::AUCell_calcAUC(gene_sets, cells_rankings, nCores = ncores, ...)
        scores <- t(AUCell::getAUC(cells_AUC))
        
    } else if (method == "GSVA") {
        if (!requireNamespace("GSVA", quietly = TRUE)) {
            stop("Please install 'GSVA' package to use method='GSVA'.")
        }
        message("Running GSVA scoring...")
        
        # In GSVA 1.50+, parameter specification is via gsvaParam
        if (utils::packageVersion("GSVA") >= "1.50.0") {
            param <- GSVA::gsvaParam(expr_mat, gene_sets, ...)
            gsva_res <- GSVA::gsva(param)
        } else {
            gsva_res <- GSVA::gsva(expr_mat, gene_sets, parallel.sz = ncores, ...)
        }
        scores <- t(gsva_res)
    }
    
    artifact_info <- list(
        storage = if (isTRUE(as_altExp)) "altExp" else "colData"
    )

    # Format scores
    if (!as_altExp) {
        # Prefix colnames with name
        colnames(scores) <- paste0(name, "_", colnames(scores))
        
        # Add to colData
        for (cn in colnames(scores)) {
            colData(sce)[[cn]] <- scores[, cn]
        }
        artifact_info$score_columns <- colnames(scores)
        msg <- sprintf("Gene set scoring completed using %s. Added %d columns to colData.", method, ncol(scores))
    } else {
        # Create altExp
        if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
            stop("SingleCellExperiment is required for altExp.")
        }
        alt_sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(scores = t(scores))
        )
        SingleCellExperiment::altExp(sce, name) <- alt_sce
        artifact_info$altExp <- name
        artifact_info$assay <- "scores"
        msg <- sprintf("Gene set scoring completed using %s. Added altExp '%s' with %d features.", method, name, ncol(scores))
    }
    
    # Register state
    geneset_analysis <- list(
        id = name,
        method = method,
        assay_use = assay_use,
        gene_sets = names(gene_sets),
        artifacts = artifact_info,
        as_altExp = as_altExp,
        timestamp = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "geneset_scoring", geneset_analysis)
    sce <- sclet_set_state_record(
        object = sce,
        type = "geneset_scoring",
        id = name,
        active = TRUE,
        value = list(
            method = method,
            inputs = list(
                assay = assay_use,
                gene_sets = names(gene_sets)
            ),
            artifacts = artifact_info,
            params = list(
                ncores = ncores,
                as_altExp = as_altExp
            ),
            summary = list(
                n_programs = ncol(scores),
                n_cells = nrow(scores)
            ),
            created_at = Sys.time()
        )
    )
    sce <- sclet_log_command(
        sce,
        "RunGeneSetScoring",
        params = list(
            method = method,
            assay_use = assay_use,
            ncores = ncores,
            name = name,
            as_altExp = as_altExp,
            gene_sets = names(gene_sets)
        ),
        outputs = artifact_info
    )
    
    message(msg)
    return(sce)
}

#' Plot program or pathway activity on an embedding
#'
#' This function provides a unified visualization entry for program-level
#' activity results produced by either `RunGeneSetScoring()` or `RunSCENIC()`.
#'
#' @param object A SingleCellExperiment object.
#' @param program Character. Program, pathway, signature, or regulon name.
#' @param source One of `"auto"`, `"geneset_scoring"`, or `"scenic"`.
#' @param reduction Character. Reduction to use. If `NULL`, uses `DefaultReduction(object)`.
#' @param id Optional analysis record id for the selected source.
#' @param assay Character. Assay name inside the selected altExp. If `NULL`, uses
#'   the recorded assay for the source.
#' @param point_size Numeric. Point size.
#' @param low,high Colors for the activity gradient.
#' @return A ggplot object.
#' @export
plot_program_activity <- function(object, program, source = c("auto", "geneset_scoring", "scenic"),
                                  reduction = NULL, id = NULL, assay = NULL,
                                  point_size = 0.6, low = "grey90", high = "firebrick") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    source <- match.arg(source)
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("A valid dimensional reduction is required to plot program activity.")
    }

    candidate_sources <- if (identical(source, "auto")) c("geneset_scoring", "scenic") else source
    resolved <- NULL
    errors <- character(0)
    for (candidate in candidate_sources) {
        attempt <- tryCatch(
            sclet_resolve_program_activity(object, program = program, source = candidate, id = id, assay = assay),
            error = function(e) e
        )
        if (!inherits(attempt, "error")) {
            resolved <- attempt
            break
        }
        errors <- c(errors, sprintf("%s: %s", candidate, conditionMessage(attempt)))
    }
    if (is.null(resolved)) {
        stop(paste(c(
            sprintf("Program '%s' could not be resolved.", program),
            errors
        ), collapse = "\n"))
    }

    emb <- SingleCellExperiment::reducedDim(object, reduction)
    df <- data.frame(
        dim1 = emb[, 1],
        dim2 = emb[, 2],
        activity = resolved$activity,
        stringsAsFactors = FALSE
    )
    df <- df[order(df$activity), , drop = FALSE]

    ggplot2::ggplot(df, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data$activity)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_gradient(low = low, high = high) +
        ggplot2::labs(
            x = paste0(reduction, "_1"),
            y = paste0(reduction, "_2"),
            color = resolved$legend,
            title = resolved$title
        ) +
        ggplot2::theme_classic()
}

#' Plot aggregated program activity across groups
#'
#' This function summarizes one or more programs across cell groups and renders
#' the result as a heatmap.
#'
#' @param object A SingleCellExperiment object.
#' @param programs Character vector of program, pathway, signature, or regulon names.
#' @param source One of `"auto"`, `"geneset_scoring"`, or `"scenic"`.
#' @param group.by Grouping variable. If `NULL`, uses `ActiveIdent(object)`. The
#'   special value `"colLabels"` uses `SingleCellExperiment::colLabels(object)`.
#' @param id Optional analysis record id for the selected source.
#' @param assay Character. Assay name inside the selected altExp. If `NULL`, uses
#'   the recorded assay for the source.
#' @param scale_rows Logical. If `TRUE`, z-score scales each program across groups.
#' @param low,mid,high Colors for the fill gradient.
#' @return A ggplot object.
#' @export
plot_program_heatmap <- function(object, programs, source = c("auto", "geneset_scoring", "scenic"),
                                 group.by = NULL, id = NULL, assay = NULL, scale_rows = TRUE,
                                 low = "#2166AC", mid = "white", high = "#B2182B") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    source <- match.arg(source)
    if (length(programs) == 0) {
        stop("`programs` must contain at least one entry.")
    }
    if (is.null(group.by)) {
        group.by <- ActiveIdent(object)
    }
    if (is.null(group.by)) {
        stop("No active identity found. Please provide `group.by`.")
    }

    if (identical(group.by, "colLabels")) {
        group_values <- SingleCellExperiment::colLabels(object)
        if (is.null(group_values)) {
            stop("group.by is 'colLabels' but `colLabels(object)` is empty.")
        }
    } else if (group.by %in% colnames(SummarizedExperiment::colData(object))) {
        group_values <- SummarizedExperiment::colData(object)[[group.by]]
    } else {
        stop(sprintf("Grouping column '%s' not found in colData(object).", group.by))
    }

    activity_list <- lapply(programs, function(program) {
        resolved <- sclet_resolve_program_activity(object, program = program, source = source, id = id, assay = assay)
        resolved$activity
    })
    names(activity_list) <- programs

    keep <- !is.na(group_values)
    if (!any(keep)) {
        stop("No non-missing group labels available to plot.")
    }
    group_values <- as.character(group_values[keep])
    activity_list <- lapply(activity_list, function(x) as.numeric(x[keep]))

    summary_df <- do.call(
        rbind,
        lapply(programs, function(program) {
            values <- activity_list[[program]]
            agg <- stats::aggregate(
                values,
                by = list(group = group_values),
                FUN = mean,
                na.rm = TRUE
            )
            transform(agg, program = program)
        })
    )
    colnames(summary_df)[colnames(summary_df) == "x"] <- "value"
    summary_df$group <- factor(summary_df$group, levels = unique(group_values))
    summary_df$program <- factor(summary_df$program, levels = programs)

    if (isTRUE(scale_rows)) {
        scaled_values <- ave(summary_df$value, summary_df$program, FUN = function(x) {
            s <- stats::sd(x, na.rm = TRUE)
            if (is.na(s) || identical(s, 0)) {
                rep(0, length(x))
            } else {
                as.numeric(scale(x))
            }
        })
        summary_df$value <- scaled_values
        fill_name <- "Scaled activity"
    } else {
        fill_name <- "Mean activity"
    }

    ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$group, y = .data$program, fill = .data$value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.3) +
        ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high, name = fill_name) +
        ggplot2::labs(
            x = if (identical(group.by, "colLabels")) "Group (colLabels)" else group.by,
            y = "Program"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        )
}

sclet_resolve_program_activity <- function(object, program, source, id = NULL, assay = NULL) {
    if (identical(source, "geneset_scoring")) {
        record <- get_geneset_scoring(object, id = id)
        if (is.null(record)) {
            stop("No gene set scoring results found.")
        }
        artifacts <- record$artifacts
        if (is.null(artifacts$storage) && !is.null(record$score_columns)) {
            artifacts$storage <- if (isTRUE(record$as_altExp)) "altExp" else "colData"
        }

        if (identical(artifacts$storage, "colData")) {
            score_columns <- artifacts$score_columns
            if (is.null(score_columns)) {
                score_columns <- record$score_columns
            }
            if (is.null(score_columns)) {
                stop("No score columns recorded for gene set scoring results.")
            }

            candidate_columns <- unique(c(
                program,
                paste0(record$id, "_", program),
                score_columns[basename(score_columns) == program],
                score_columns[endsWith(score_columns, paste0("_", program))]
            ))
            candidate_columns <- candidate_columns[!is.na(candidate_columns) & nzchar(candidate_columns)]
            score_col <- candidate_columns[candidate_columns %in% colnames(SummarizedExperiment::colData(object))][1]
            if (is.na(score_col) || is.null(score_col)) {
                stop(sprintf("Program '%s' not found in recorded gene set score columns.", program))
            }

            return(list(
                activity = as.numeric(SummarizedExperiment::colData(object)[[score_col]]),
                title = program,
                legend = "Score"
            ))
        }

        alt_exp_name <- artifacts$altExp
        assay_name <- if (is.null(assay)) artifacts$assay else assay
        if (is.null(alt_exp_name) || !alt_exp_name %in% SingleCellExperiment::altExpNames(object)) {
            stop("Recorded gene set scoring altExp not found in object.")
        }
        if (is.null(assay_name)) {
            assay_name <- "scores"
        }
        score_mat <- SummarizedExperiment::assay(SingleCellExperiment::altExp(object, alt_exp_name), assay_name)
        if (!program %in% rownames(score_mat)) {
            stop(sprintf("Program '%s' not found in altExp '%s'.", program, alt_exp_name))
        }

        return(list(
            activity = as.numeric(score_mat[program, ]),
            title = program,
            legend = "Score"
        ))
    }

    if (identical(source, "scenic")) {
        scenic_record <- get_scenic(object, id = id)
        if (is.null(scenic_record)) {
            stop("No SCENIC results found.")
        }
        alt_exp_name <- scenic_record$artifacts$altExp
        if (is.null(alt_exp_name)) {
            alt_exp_name <- "SCENIC_AUC"
        }
        if (!alt_exp_name %in% SingleCellExperiment::altExpNames(object)) {
            stop(sprintf("SCENIC altExp '%s' not found in object.", alt_exp_name))
        }
        assay_name <- if (is.null(assay)) scenic_record$artifacts$assay else assay
        if (is.null(assay_name)) {
            assay_name <- "AUC"
        }
        auc_mat <- SummarizedExperiment::assay(SingleCellExperiment::altExp(object, alt_exp_name), assay_name)
        if (!program %in% rownames(auc_mat)) {
            stop(sprintf("Regulon '%s' not found in SCENIC assay '%s'.", program, assay_name))
        }

        return(list(
            activity = as.numeric(auc_mat[program, ]),
            title = program,
            legend = "AUC"
        ))
    }

    stop(sprintf("Unsupported program activity source '%s'.", source))
}

#' Dotplot of program activity across groups
#'
#' Summarizes program activity per group with dot size reflecting the fraction of
#' cells expressing detectable activity and dot color reflecting mean activity.
#'
#' @param object A SingleCellExperiment object.
#' @param programs Character vector of program names.
#' @param source One of `"auto"`, `"geneset_scoring"`, or `"scenic"`.
#' @param group.by Grouping variable. If `NULL`, uses `ActiveIdent(object)`.
#' @param id Optional analysis record id.
#' @param low,high Colors for the fill gradient. Defaults to blue-red.
#'
#' @return A ggplot object.
#' @export
plot_program_dotplot <- function(object, programs, source = c("auto", "geneset_scoring", "scenic"),
                                 group.by = NULL, id = NULL,
                                 low = "#2166AC", high = "#B2182B") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }
    if (length(programs) == 0) {
        stop("`programs` must contain at least one entry.")
    }

    source <- match.arg(source)
    if (is.null(group.by)) {
        group.by <- ActiveIdent(object)
    }
    if (is.null(group.by)) {
        stop("No active identity found. Please provide `group.by`.")
    }
    if (identical(group.by, "colLabels")) {
        group_values <- SingleCellExperiment::colLabels(object)
        if (is.null(group_values)) {
            stop("group.by is 'colLabels' but `colLabels(object)` is empty.")
        }
    } else if (group.by %in% colnames(SummarizedExperiment::colData(object))) {
        group_values <- SummarizedExperiment::colData(object)[[group.by]]
    } else {
        stop(sprintf("Grouping column '%s' not found in colData(object).", group.by))
    }
    group_values <- as.character(group_values)

    dot_data <- lapply(programs, function(prog) {
        activity <- tryCatch(
            get_program(object, program = prog, source = source, id = id),
            error = function(e) NULL
        )
        if (is.null(activity)) return(NULL)
        groups <- unique(group_values)
        do.call(rbind, lapply(groups, function(g) {
            idx <- which(group_values == g)
            act <- activity[idx]
            data.frame(
                program = prog,
                group = g,
                mean_activity = mean(act, na.rm = TRUE),
                pct_expressed = sum(act > 0, na.rm = TRUE) / length(act) * 100,
                stringsAsFactors = FALSE
            )
        }))
    })
    dot_data <- do.call(rbind, dot_data)
    if (is.null(dot_data) || nrow(dot_data) == 0) {
        stop("No program activity data could be resolved.")
    }
    dot_data$program <- factor(dot_data$program, levels = programs)
    dot_data$group <- factor(dot_data$group)

    ggplot2::ggplot(dot_data, ggplot2::aes(
        x = .data$program, y = .data$group,
        size = .data$pct_expressed, fill = .data$mean_activity
    )) +
        ggplot2::geom_point(shape = 21, color = "grey70") +
        ggplot2::scale_size_continuous(name = "Percent\nexpressed", range = c(1, 6)) +
        ggplot2::scale_fill_gradient2(low = low, mid = "white", high = high, name = "Mean\nactivity") +
        ggplot2::labs(x = NULL, y = NULL, title = "Program Activity Dotplot") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
            plot.title = ggplot2::element_text(face = "bold", size = 13),
            panel.grid.major = ggplot2::element_line(linewidth = 0.3)
        )
}
