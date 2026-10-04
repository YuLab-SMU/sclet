# Deterministic, bounded diagnostics for AI-led advanced analysis

sclet_ai_diag_not_available <- function(reason, requested = NULL) {
    list(
        status = "not_available",
        reason = as.character(reason),
        requested = requested
    )
}

sclet_ai_diag_coldata <- function(object, column) {
    if (!inherits(object, "SingleCellExperiment")) {
        return(NULL)
    }
    if (!is.character(column) || length(column) != 1L || !nzchar(column)) {
        return(NULL)
    }
    data <- SummarizedExperiment::colData(object)
    if (!(column %in% colnames(data))) {
        return(NULL)
    }
    data[[column]]
}

sclet_ai_diag_group <- function(object, group, role = "group") {
    values <- if (is.character(group) && length(group) == 1L) {
        sclet_ai_diag_coldata(object, group)
    } else if (length(group) == ncol(object)) {
        group
    } else {
        NULL
    }
    if (is.null(values) || length(values) != ncol(object)) {
        return(sclet_ai_diag_not_available(
            paste0("missing_", role, "_column"),
            requested = group
        ))
    }
    missing <- is.na(values) | !nzchar(as.character(values))
    values <- as.character(values)
    values[missing] <- "__missing__"
    levels <- unique(values)
    encoded <- match(values, levels)
    list(
        status = "available",
        requested = if (is.character(group) && length(group) == 1L) group else NULL,
        n_groups = length(levels),
        group_index = encoded,
        labels = paste0("group_", seq_along(levels)),
        missing_cells = sum(missing),
        raw_values_included = FALSE
    )
}

sclet_ai_diag_numeric_columns <- function(object, exclude = character()) {
    data <- SummarizedExperiment::colData(object)
    columns <- setdiff(colnames(data), exclude)
    columns[vapply(columns, function(column) {
        value <- data[[column]]
        is.numeric(value) || is.integer(value)
    }, logical(1L))]
}

sclet_ai_diag_summary <- function(value) {
    value <- as.numeric(value)
    observed <- is.finite(value)
    value <- value[observed]
    if (!length(value)) {
        return(list(n = 0L, missing_fraction = 1, mean = NULL, median = NULL, sd = NULL))
    }
    list(
        n = as.integer(length(value)),
        missing_fraction = NULL,
        mean = unname(mean(value)),
        median = unname(stats::median(value)),
        sd = if (length(value) > 1L) unname(stats::sd(value)) else 0
    )
}

#' Summarize numeric QC columns by an anonymized metadata grouping
#'
#' @param object A `SingleCellExperiment` object.
#' @param group A `colData` column name or a vector with one value per cell.
#' @return A bounded aggregate summary; raw group values are never returned.
#' @export
summarize_qc_by_group <- function(object, group) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    grouping <- sclet_ai_diag_group(object, group, "group")
    if (!identical(grouping$status, "available")) {
        return(grouping)
    }
    columns <- sclet_ai_diag_numeric_columns(object, exclude = if (is.character(group)) group else character())
    if (!length(columns)) {
        return(sclet_ai_diag_not_available("no_numeric_qc_columns"))
    }
    data <- SummarizedExperiment::colData(object)
    summaries <- lapply(columns, function(column) {
        value <- as.numeric(data[[column]])
        by_group <- lapply(seq_len(grouping$n_groups), function(index) {
            selected <- grouping$group_index == index
            result <- sclet_ai_diag_summary(value[selected])
            result$missing_fraction <- mean(is.na(value[selected]) | !is.finite(value[selected]))
            result
        })
        names(by_group) <- grouping$labels
        list(
            column = column,
            class = class(data[[column]]),
            by_group = by_group,
            raw_values_included = FALSE
        )
    })
    names(summaries) <- columns
    list(
        status = "available",
        group = grouping[c("requested", "n_groups", "labels", "missing_cells", "raw_values_included")],
        metrics = summaries,
        n_cells = ncol(object),
        raw_values_included = FALSE
    )
}

sclet_ai_diag_reduction_name <- function(object, requested = NULL) {
    available <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    if (!length(available)) {
        return(NULL)
    }
    if (!is.null(requested) && requested %in% available) {
        return(requested)
    }
    index <- match(tolower(available), "pca")
    if (any(!is.na(index))) return(available[which(!is.na(index))[1L]])
    index <- grep("pca", tolower(available))
    if (length(index)) available[[index[[1L]]]] else NULL
}

sclet_ai_diag_eta_squared <- function(values, groups) {
    keep <- is.finite(values) & !is.na(groups)
    values <- as.numeric(values[keep])
    groups <- groups[keep]
    if (length(values) < 2L || length(unique(groups)) < 2L) return(NA_real_)
    total <- sum((values - mean(values))^2)
    if (!is.finite(total) || total == 0) return(0)
    between <- sum(vapply(split(values, groups), function(x) {
        length(x) * (mean(x) - mean(values))^2
    }, numeric(1L)))
    min(1, max(0, between / total))
}

#' Summarize PCA association with a metadata column
#'
#' @param object A `SingleCellExperiment` object.
#' @param metadata A `colData` column name or a vector with one value per cell.
#' @param reduction Optional reduction name; defaults to PCA if available.
#' @return A bounded association summary without metadata values.
#' @export
summarize_pca_metadata_association <- function(object, metadata, reduction = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    reduction <- sclet_ai_diag_reduction_name(object, reduction)
    if (is.null(reduction)) return(sclet_ai_diag_not_available("pca_not_available"))
    embedding <- tryCatch(SingleCellExperiment::reducedDim(object, reduction), error = function(e) NULL)
    if (is.null(embedding) || nrow(embedding) != ncol(object)) {
        return(sclet_ai_diag_not_available("pca_dimensions_not_aligned"))
    }
    value <- if (is.character(metadata) && length(metadata) == 1L) {
        sclet_ai_diag_coldata(object, metadata)
    } else if (length(metadata) == ncol(object)) metadata else NULL
    if (is.null(value)) return(sclet_ai_diag_not_available("metadata_not_available", metadata))
    missing <- is.na(value) | !nzchar(as.character(value))
    values <- as.character(value)
    values[missing] <- "__missing__"
    numeric_metadata <- is.numeric(value) || is.integer(value)
    associations <- lapply(seq_len(min(ncol(embedding), 20L)), function(index) {
        pc <- as.numeric(embedding[, index])
        score <- if (numeric_metadata) {
            keep <- is.finite(pc) & is.finite(as.numeric(value))
            if (sum(keep) < 3L) NA_real_ else abs(stats::cor(pc[keep], as.numeric(value)[keep]))
        } else {
            sclet_ai_diag_eta_squared(pc, values)
        }
        list(component = index, association = unname(score), method = if (numeric_metadata) "absolute_pearson" else "eta_squared")
    })
    scores <- vapply(associations, function(x) x$association, numeric(1L))
    max_association <- if (any(is.finite(scores))) max(scores[is.finite(scores)]) else NULL
    list(
        status = "available",
        reduction = reduction,
        metadata = if (is.character(metadata) && length(metadata) == 1L) metadata else NULL,
        metadata_type = if (numeric_metadata) "numeric" else "categorical",
        n_components = ncol(embedding),
        associations = associations,
        max_association = max_association,
        missing_cells = sum(missing),
        raw_values_included = FALSE
    )
}

#' Summarize cluster composition across anonymized samples
#'
#' @param object A `SingleCellExperiment` object.
#' @param sample A `colData` sample column name or vector.
#' @param cluster A `colData` cluster column name or vector.
#' @return An aggregate cluster-by-sample composition summary.
#' @export
summarize_cluster_sample_composition <- function(object, sample, cluster = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    data <- SummarizedExperiment::colData(object)
    if (is.null(cluster)) {
        candidates <- grep("cluster|leiden|louvain|ident", colnames(data), ignore.case = TRUE, value = TRUE)
        cluster <- if (length(candidates)) candidates[[1L]] else NULL
    }
    samples <- sclet_ai_diag_group(object, sample, "sample")
    clusters <- sclet_ai_diag_group(object, cluster, "cluster")
    if (!identical(samples$status, "available")) return(samples)
    if (!identical(clusters$status, "available")) return(clusters)
    counts <- table(clusters$group_index, samples$group_index)
    proportions <- counts / pmax(rowSums(counts), 1)
    list(
        status = "available",
        sample = samples[c("requested", "n_groups", "labels", "raw_values_included")],
        cluster = clusters[c("requested", "n_groups", "labels", "raw_values_included")],
        counts = unname(as.matrix(counts)),
        proportions = unname(as.matrix(proportions)),
        cluster_sizes = as.integer(rowSums(counts)),
        raw_values_included = FALSE
    )
}

#' Check whether an object is ready for rare-cell and doublet diagnosis
#'
#' Reports whether the deterministic inputs required by
#' `RunRareCellDetection()` and `RunDoubletFinder()` are present. A missing
#' doublet run is not blocking, but it is reported because doublet evidence is
#' one of the independent signals required before a small population may be
#' described with anything stronger than a low-confidence guess.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param cluster Optional \code{colData} cluster column name.
#' @return A read-only readiness summary.
#' @export
check_rare_cell_readiness <- function(object, cluster = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    columns <- colnames(SummarizedExperiment::colData(object))
    active_ident <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    idents <- if (!is.null(active_ident)) tryCatch(Idents(object), error = function(e) NULL) else NULL
    has_idents <- !is.null(idents) && length(idents) == ncol(object) &&
        (is.factor(idents) || is.character(idents)) &&
        length(unique(as.character(stats::na.omit(as.character(idents))))) >= 1L
    cluster_col <- sclet_ai_diag_cluster_column(object, cluster)
    has_cluster_column <- !is.null(cluster_col)
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    has_pca <- "PCA" %in% reductions
    doublet_col <- "scDblFinder.class" %in% columns
    questions <- character()
    if (!has_idents && !has_cluster_column) {
        questions <- c(questions, paste0(
            "Which colData column holds the cluster / ident assignment used to define small populations? ",
            "Available colData columns: ", paste(columns, collapse = ", ")
        ))
    }
    if (!has_pca) {
        questions <- c(questions, paste0(
            "No PCA reduction is available; run RunPCA() first because rare-cell detection is ",
            "density based on a dimensional reduction. Available reductions: ",
            paste(reductions, collapse = ", ")
        ))
    }
    status <- if (!has_idents && !has_cluster_column) "not_ready" else if (!has_pca) "not_ready" else "ready_for_diagnostic"
    list(
        status = status,
        checks = list(
            has_idents = has_idents,
            has_cluster_column = has_cluster_column,
            cluster_column_resolved = cluster_col,
            has_pca = has_pca,
            reductions = reductions,
            doublet_evidence_available = doublet_col,
            raw_values_included = FALSE
        ),
        blocked_actions = if (identical(status, "not_ready")) {
            c("run_rare_cell_detection", "run_doublet_detection")
        } else {
            character()
        },
        questions = questions,
        notes = if (doublet_col) {
            character()
        } else {
            "doublet evidence missing: without scDblFinder.class the number of independent signals for any small population stays below the required minimum"
        }
    )
}

sclet_ai_diag_cluster_column <- function(object, requested = NULL) {
    columns <- colnames(SummarizedExperiment::colData(object))
    if (is.character(requested) && length(requested) == 1L && nzchar(requested)) {
        return(if (requested %in% columns) requested else NULL)
    }
    if (length(requested) == ncol(object) && length(requested)) return("__vector__")
    if ("rare_cluster" %in% columns) return("rare_cluster")
    candidates <- grep("cluster|leiden|louvain|ident", columns, ignore.case = TRUE, value = TRUE)
    if (length(candidates)) candidates[[1L]] else NULL
}

sclet_ai_diag_sample_column <- function(object, exclude = character()) {
    columns <- setdiff(colnames(SummarizedExperiment::colData(object)), exclude)
    candidates <- grep("sample|batch|donor|patient|subject|individual|library", columns, ignore.case = TRUE, value = TRUE)
    if (length(candidates)) candidates[[1L]] else NULL
}

sclet_ai_diag_qc_columns <- function(object, exclude = character()) {
    all_columns <- colnames(SummarizedExperiment::colData(object))
    columns <- setdiff(all_columns, exclude)
    columns <- columns[!grepl("^scDblFinder|^decontx", columns)]
    sclet_ai_diag_numeric_columns(object, exclude = setdiff(all_columns, columns))
}

sclet_ai_diag_marker_signal <- function(object, cluster_label, cluster_index, n_clusters) {
    detest <- tryCatch(sclet_get_state_records(object, "detest"), error = function(e) list())
    annotation <- tryCatch(sclet_get_state_records(object, "annotation"), error = function(e) list())
    detest <- if (is.list(detest)) detest else list()
    annotation <- if (is.list(annotation)) annotation else list()
    covered <- FALSE
    for (record in detest) {
        if (!is.list(record)) next
        # FindAllMarkers runs over every cluster by construction, so it covers
        # whichever cluster is being summarized. FindMarkers records the two
        # compared labels explicitly, so match those directly.
        if (identical(as.character(record$method %||% ""), "FindAllMarkers")) {
            covered <- TRUE
            break
        }
        parts <- c(
            as.character(record$id %||% ""),
            as.character(record$inputs$ident.1 %||% ""),
            as.character(record$inputs$ident.2 %||% ""),
            as.character(record$params$ident.1 %||% ""),
            as.character(record$params$ident.2 %||% "")
        )
        if (cluster_label %in% parts) {
            covered <- TRUE
            break
        }
        n_groups <- suppressWarnings(as.integer(record$summary$n_groups %||% NA_integer_))
        if (length(n_groups) == 1L && !is.na(n_groups) && n_groups >= n_clusters) {
            covered <- TRUE
            break
        }
    }
    annotation_cols <- grep("_labels$|pruned[.]labels$|^annotation",
        colnames(SummarizedExperiment::colData(object)), value = TRUE)
    n_labels <- 0L
    for (column in annotation_cols) {
        value <- as.character(SummarizedExperiment::colData(object)[[column]])
        labels <- unique(value[cluster_index])
        labels <- labels[!is.na(labels) & nzchar(labels) & labels != "NA"]
        n_labels <- max(n_labels, length(labels))
    }
    n_annotation_records <- length(annotation)
    available <- covered || n_annotation_records > 0L || n_labels > 0L
    list(
        available = available,
        n_marker_records = as.integer(length(detest)),
        n_annotation_records = as.integer(n_annotation_records),
        cluster_covered_by_marker_analysis = isTRUE(covered),
        n_distinct_annotation_labels = as.integer(n_labels),
        raw_values_included = FALSE
    )
}

#' Aggregate the already computed signals bearing on a small population
#'
#' `summarize_small_cluster_evidence()` only aggregates signals that have
#' already been computed on the object: QC columns, doublet calls, recorded
#' marker/annotation analyses, and sample or batch columns. It deliberately
#' returns no verdict about whether a small population is real or noise; it
#' reports how many independent signal classes are available so that a caller
#' can grade the strength of any claim made about that population.
#'
#' @param object A `SingleCellExperiment` object.
#' @param cluster A `colData` cluster column name, or a vector with one value
#'   per cell. Defaults to a detected cluster-like column.
#' @param size_threshold Maximum number of cells for a population to be
#'   considered small.
#' @param label Optional single cluster label to summarize. When `NULL` and
#'   exactly one small cluster exists, that cluster is summarized; when several
#'   small clusters exist, each is listed under `clusters`.
#' @return A bounded list with `independent_signals` for each of `qc`,
#'   `doublet`, `marker` and `sample_replication`, plus
#'   `n_independent_signals_available`. No population-level verdict is
#'   returned and raw per-cell values are never included.
#' @export
summarize_small_cluster_evidence <- function(object, cluster = NULL, size_threshold = 10L, label = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    size_threshold <- as.numeric(size_threshold)
    if (length(size_threshold) != 1L || !is.finite(size_threshold) || size_threshold < 1) {
        stop("size_threshold must be one positive number", call. = FALSE)
    }
    cluster_col <- sclet_ai_diag_cluster_column(object, cluster)
    if (is.null(cluster_col)) return(sclet_ai_diag_not_available("missing_cluster_column", requested = cluster))
    grouping <- if (identical(cluster_col, "__vector__")) {
        sclet_ai_diag_group(object, cluster, "cluster")
    } else {
        sclet_ai_diag_group(object, cluster_col, "cluster")
    }
    if (!identical(grouping$status, "available")) return(grouping)
    sizes <- tabulate(grouping$group_index, nbins = grouping$n_groups)
    small <- which(sizes <= size_threshold)
    if (!is.null(label)) {
        label <- as.character(label)
        if (length(label) != 1L || !nzchar(label)) stop("label must be a single non-empty string", call. = FALSE)
        index <- match(label, grouping$labels)
        if (is.na(index)) stop("label is not one of the detected clusters", call. = FALSE)
        single <- sclet_ai_diag_small_cluster_signals(object, grouping, index, sizes, cluster_col)
        single$size_threshold <- size_threshold
        single$size_is_below_threshold <- as.integer(sizes[[index]]) <= size_threshold
        single$n_small_clusters <- length(small)
        return(single)
    }
    if (!length(small)) {
        return(list(
            status = "available",
            n_small_clusters = 0L,
            size_threshold = size_threshold,
            clusters = list(),
            raw_values_included = FALSE
        ))
    }
    entries <- lapply(small, function(index) {
        sclet_ai_diag_small_cluster_signals(object, grouping, index, sizes, cluster_col)
    })
    names(entries) <- grouping$labels[small]
    base <- list(
        status = "available",
        size_threshold = size_threshold,
        n_small_clusters = length(small),
        clusters = entries,
        raw_values_included = FALSE
    )
    if (length(small) != 1L) return(base)
    single <- entries[[1L]]
    single$size_threshold <- size_threshold
    single$size_is_below_threshold <- TRUE
    single$n_small_clusters <- 1L
    single$clusters <- entries
    single
}

sclet_ai_diag_small_cluster_signals <- function(object, grouping, index, sizes, cluster_col) {
    selected <- grouping$group_index == index
    size <- as.integer(sizes[[index]])
    n_cells <- ncol(object)
    other <- !selected
    columns <- colnames(SummarizedExperiment::colData(object))
    sample_col <- sclet_ai_diag_sample_column(object, exclude = c(cluster_col, "rare_cluster"))

    qc_columns <- if (identical(cluster_col, "__vector__")) {
        character()
    } else {
        sclet_ai_diag_qc_columns(object, exclude = unique(c(cluster_col, sample_col)))
    }
    qc_etas <- numeric()
    if (length(qc_columns) && any(other) && any(selected)) {
        for (column in qc_columns) {
            value <- as.numeric(SummarizedExperiment::colData(object)[[column]])
            groups <- rep("other_cells", length(value))
            groups[selected] <- "cluster_cells"
            eta <- sclet_ai_diag_eta_squared(value, groups)
            if (is.finite(eta)) qc_etas <- c(qc_etas, eta)
        }
    }
    qc <- if (!length(qc_etas)) {
        list(available = FALSE, reason = "no_numeric_qc_columns", raw_values_included = FALSE)
    } else {
        list(
            available = TRUE,
            n_qc_metrics = as.integer(length(qc_etas)),
            max_eta_squared = round(max(qc_etas), 4L),
            deviates_from_other_cells = max(qc_etas) >= 0.25,
            raw_values_included = FALSE
        )
    }

    doublet <- if ("scDblFinder.class" %in% columns) {
        cls <- as.character(SummarizedExperiment::colData(object)[["scDblFinder.class"]])
        inside <- cls[selected]
        outside <- cls[other]
        list(
            available = TRUE,
            n_doublet = as.integer(sum(inside == "doublet", na.rm = TRUE)),
            doublet_fraction = round(sum(inside == "doublet", na.rm = TRUE) / max(length(inside), 1L), 4L),
            doublet_fraction_other_cells = if (length(outside)) {
                round(sum(outside == "doublet", na.rm = TRUE) / length(outside), 4L)
            } else {
                NA_real_
            },
            raw_values_included = FALSE
        )
    } else {
        list(available = FALSE, reason = "scDblFinder_class_missing", raw_values_included = FALSE)
    }

    marker <- sclet_ai_diag_marker_signal(
        object,
        cluster_label = grouping$labels[[index]],
        cluster_index = selected,
        n_clusters = grouping$n_groups
    )

    sample_replication <- if (is.null(sample_col)) {
        list(available = FALSE, reason = "no_sample_or_batch_column", raw_values_included = FALSE)
    } else {
        composition <- summarize_cluster_sample_composition(object, sample = sample_col,
            cluster = if (identical(cluster_col, "__vector__")) selected else cluster_col)
        if (!identical(composition$status, "available")) {
            list(available = FALSE, reason = "sample_composition_not_available", raw_values_included = FALSE)
        } else {
            row <- composition$proportions[index, , drop = TRUE]
            present <- sum(row > 0)
            list(
                available = TRUE,
                n_samples_present = as.integer(present),
                n_samples_total = as.integer(composition$sample$n_groups),
                present_in_multiple_samples = present >= 2L,
                max_sample_fraction = round(max(row, na.rm = TRUE), 4L),
                raw_values_included = FALSE
            )
        }
    }

    signals <- list(
        qc = qc,
        doublet = doublet,
        marker = marker,
        sample_replication = sample_replication
    )
    n_available <- sum(vapply(signals, function(x) isTRUE(x$available), logical(1L)))
    list(
        status = "available",
        cluster = grouping$labels[[index]],
        size = size,
        fraction_of_total = round(size / max(n_cells, 1L), 6L),
        independent_signals = signals,
        n_independent_signals_available = as.integer(n_available),
        raw_values_included = FALSE
    )
}

sclet_ai_diag_embedding_name <- function(object, requested = NULL) {
    available <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    if (!length(available)) return(NULL)
    if (is.character(requested) && length(requested) == 1L && requested %in% available) {
        return(requested)
    }
    preferred <- c("umap", "tsne", "pca")
    for (target in preferred) {
        index <- match(tolower(available), target)
        if (any(!is.na(index))) return(available[which(!is.na(index))[[1L]]])
    }
    NULL
}

#' Check whether an object is ready for trajectory inference
#'
#' Reports whether a cluster assignment and a usable embedding are present for
#' `RunSlingshot()`. This diagnostic deliberately never proposes a root: which
#' cluster represents the origin of a trajectory is a biological assumption that
#' only the user can make, so the result contains no recommendation, no
#' candidate root, and no ordering of clusters.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param reduction Optional preferred reduction name. Defaults to `"UMAP"`
#'   when available, then any usable embedding.
#' @return A read-only readiness summary.
#' @export
check_trajectory_readiness <- function(object, reduction = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    columns <- colnames(SummarizedExperiment::colData(object))
    active_ident <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    idents <- if (!is.null(active_ident)) tryCatch(Idents(object), error = function(e) NULL) else NULL
    has_idents <- !is.null(idents) && length(idents) == ncol(object) &&
        (is.factor(idents) || is.character(idents)) &&
        length(unique(as.character(stats::na.omit(as.character(idents))))) >= 2L
    cluster_col <- NULL
    candidates <- grep("cluster|leiden|louvain|ident", columns, ignore.case = TRUE, value = TRUE)
    if (length(candidates)) cluster_col <- candidates[[1L]]
    has_cluster_column <- !is.null(cluster_col)
    n_clusters <- 0L
    if (has_idents) {
        n_clusters <- length(unique(as.character(stats::na.omit(as.character(idents)))))
    } else if (has_cluster_column) {
        values <- as.character(SummarizedExperiment::colData(object)[[cluster_col]])
        n_clusters <- length(unique(stats::na.omit(values)))
    }
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    resolved <- sclet_ai_diag_embedding_name(object, reduction)
    has_reduction <- !is.null(resolved)
    questions <- character()
    if (!has_idents && !has_cluster_column) {
        questions <- c(questions, paste0(
            "No cluster assignment is available; trajectory inference needs an explicit cluster column. ",
            "Available colData columns: ", paste(columns, collapse = ", ")
        ))
    }
    if (!has_reduction) {
        questions <- c(questions, paste0(
            "No usable embedding is available; run the dimensionality reduction first. ",
            "Available reductions: ", paste(reductions, collapse = ", ")
        ))
    }
    status <- if (!has_idents && !has_cluster_column) {
        "not_ready"
    } else if (!has_reduction) {
        "not_ready"
    } else {
        "ready_for_diagnostic"
    }
    list(
        status = status,
        checks = list(
            has_idents = has_idents,
            has_cluster_column = has_cluster_column,
            cluster_column_resolved = cluster_col,
            n_clusters = as.integer(n_clusters),
            reduction_resolved = resolved,
            available_reductions = reductions,
            raw_values_included = FALSE
        ),
        blocked_actions = if (identical(status, "not_ready")) "trajectory" else character(),
        questions = questions,
        notes = paste(
            "No root or start cluster is suggested here: choosing the trajectory origin is a biological",
            "assumption reserved for the user, and run_trajectory requires an explicit start_cluster."
        )
    )
}

#' Summarize cluster position within a trajectory embedding
#'
#' Reports how clusters are distributed across the leading embedding
#' dimensions. This is a description of the data only: it never ranks clusters,
#' never suggests a root, and never interprets the cluster order as time.
#'
#' @param object A `SingleCellExperiment` object.
#' @param cluster Optional `colData` cluster column name, or a vector with one
#'   value per cell.
#' @param reduction Optional reduction name; defaults to UMAP, then any usable
#'   embedding.
#' @param n_dims Number of leading embedding dimensions to summarize.
#' @return A bounded list of per-cluster dimension summaries.
#' @export
summarize_trajectory_cluster_order <- function(object, cluster = NULL, reduction = NULL, n_dims = 2L) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    n_dims <- as.integer(n_dims)
    if (length(n_dims) != 1L || is.na(n_dims) || n_dims < 1L) {
        stop("n_dims must be one positive integer", call. = FALSE)
    }
    reduction <- sclet_ai_diag_embedding_name(object, reduction)
    if (is.null(reduction)) {
        return(sclet_ai_diag_not_available("reduction_not_available", requested = reduction))
    }
    embedding <- tryCatch(SingleCellExperiment::reducedDim(object, reduction), error = function(e) NULL)
    if (is.null(embedding) || nrow(embedding) != ncol(object)) {
        return(sclet_ai_diag_not_available("reduction_dimensions_not_aligned", requested = reduction))
    }
    grouping <- if (is.null(cluster)) {
        candidates <- grep("cluster|leiden|louvain|ident",
            colnames(SummarizedExperiment::colData(object)), ignore.case = TRUE, value = TRUE)
        if (!length(candidates)) return(sclet_ai_diag_not_available("missing_cluster_column"))
        sclet_ai_diag_group(object, candidates[[1L]], "cluster")
    } else {
        sclet_ai_diag_group(object, cluster, "cluster")
    }
    if (!identical(grouping$status, "available")) return(grouping)
    n_dims <- min(n_dims, ncol(embedding))
    dims <- lapply(seq_len(n_dims), function(index) {
        values <- as.numeric(embedding[, index])
        by_group <- lapply(seq_len(grouping$n_groups), function(group_index) {
            group_values <- values[grouping$group_index == group_index]
            group_values <- group_values[is.finite(group_values)]
            list(
                mean = if (length(group_values)) unname(mean(group_values)) else NA_real_,
                sd = if (length(group_values) > 1L) unname(stats::sd(group_values)) else 0,
                n = as.integer(length(group_values))
            )
        })
        names(by_group) <- grouping$labels
        by_group
    })
    names(dims) <- paste0("dim_", seq_len(n_dims))
    list(
        status = "available",
        reduction = reduction,
        n_dims = as.integer(n_dims),
        n_groups = grouping$n_groups,
        groups = grouping[c("requested", "labels", "missing_cells", "raw_values_included")],
        by_dimension = dims,
        root_suggested = FALSE,
        raw_values_included = FALSE
    )
}

#' Summarize very small clusters without exposing cluster labels
#'
#' @param object A `SingleCellExperiment` object.
#' @param threshold Maximum number of cells for a small cluster.
#' @param cluster A `colData` cluster column name or vector.
#' @return A bounded small-cluster summary.
#' @export
summarize_small_clusters <- function(object, threshold = 10L, cluster = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    threshold <- as.numeric(threshold)
    if (length(threshold) != 1L || !is.finite(threshold) || threshold < 1) stop("threshold must be one positive number", call. = FALSE)
    data <- SummarizedExperiment::colData(object)
    if (is.null(cluster)) {
        candidates <- grep("cluster|leiden|louvain|ident", colnames(data), ignore.case = TRUE, value = TRUE)
        cluster <- if (length(candidates)) candidates[[1L]] else NULL
    }
    grouping <- sclet_ai_diag_group(object, cluster, "cluster")
    if (!identical(grouping$status, "available")) return(grouping)
    sizes <- tabulate(grouping$group_index, nbins = grouping$n_groups)
    small <- which(sizes <= threshold)
    list(
        status = "available",
        cluster = grouping[c("requested", "n_groups", "labels", "raw_values_included")],
        threshold = threshold,
        small_clusters = paste0("cluster_", small),
        small_sizes = as.integer(sizes[small]),
        small_fractions = unname(sizes[small] / ncol(object)),
        n_small_clusters = length(small),
        raw_values_included = FALSE
    )
}

#' Check whether an object is ready for annotation and differential-expression actions
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param design Optional named list describing an annotation design. A named
#'   element \code{reference} or \code{ref} identifies the reference dataset
#'   (e.g. a celldex function name such as \code{"HumanPrimaryCellAtlasData"});
#'   \code{labels} optionally identifies the labels column;
#'   \code{groups} optionally identifies the cluster column used for DE.
#'   Passing a reference here avoids the \code{clarification_required} status.
#' @return A read-only readiness and clarification summary.
#' @export
check_annotation_readiness <- function(object, design = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    design <- if (is.null(design)) list() else design
    columns <- colnames(SummarizedExperiment::colData(object))
    assay_names <- SummarizedExperiment::assayNames(object)
    has_assay <- length(assay_names) > 0L
    has_counts <- "counts" %in% assay_names || "logcounts" %in% assay_names
    active_ident <- ActiveIdent(object)
    idents <- if (!is.null(active_ident)) Idents(object) else NULL
    has_idents <- !is.null(idents) && length(idents) == ncol(object) &&
        (is.factor(idents) || is.character(idents)) &&
        length(unique(as.character(stats::na.omit(as.character(idents))))) >= 1L
    idents_table <- if (has_idents) table(as.character(idents), useNA = "no") else table(character())
    smallest_group_n <- if (length(idents_table)) min(as.integer(idents_table)) else 0L
    de_min_n <- 2L
    groups_col <- NULL
    if (!is.null(design$groups)) {
        groups_col <- as.character(design$groups[[1L]])
        if (!groups_col %in% columns) groups_col <- NULL
    }
    has_groups_col <- !is.null(groups_col)
    ref_resolved <- NULL
    ref_candidate <- design$reference %||% design$ref
    if (length(ref_candidate) == 1L && is.character(ref_candidate) && nzchar(ref_candidate)) {
        ref_resolved <- ref_candidate
    }
    labels_candidate <- design$labels
    labels_resolved <- NULL
    if (length(labels_candidate) == 1L && is.character(labels_candidate) && nzchar(labels_candidate)) {
        labels_resolved <- labels_candidate
    }
    questions <- character()
    if (is.null(ref_resolved)) {
        questions <- c(questions, "Which reference dataset should be used? (e.g. HumanPrimaryCellAtlasData, MouseRNAseqData, BlueprintENCODEData, DatabaseImmuneCellExpressionData, NovershternHematopoieticData, MonacoImmuneData, or a SummarizedExperiment object identifier)")
        questions <- c(questions, "What is the species and tissue of the dataset?")
    }
    if (is.null(labels_resolved) && !is.null(ref_resolved)) {
        questions <- c(questions, "Which labels column in the reference should be used (e.g. label.main, label.fine)?")
    }
    if (!has_idents && !has_groups_col) {
        questions <- c(questions, "Which colData column holds cluster / group / cell-identity assignments for DE and per-cluster annotation summaries?")
    }
    status <- if (!has_assay || !has_counts) {
        "not_ready"
    } else if (length(questions) > 0L) {
        "clarification_required"
    } else if (!has_idents && !has_groups_col) {
        "clarification_required"
    } else if (has_idents && smallest_group_n < de_min_n) {
        "not_ready"
    } else {
        "ready_for_diagnostic"
    }
    blocked_actions <- if (identical(status, "not_ready")) c("run_annotation", "run_de_test")
        else if (identical(status, "clarification_required")) c("run_annotation", "run_de_test")
        else character()
    list(
        status = status,
        checks = list(
            assay_available = has_assay,
            counts_or_logcounts_present = has_counts,
            active_ident = active_ident,
            has_idents = has_idents,
            smallest_group_n = as.integer(smallest_group_n),
            groups_column_resolved = groups_col,
            reference_resolved = ref_resolved,
            labels_resolved = labels_resolved,
            raw_values_included = FALSE
        ),
        blocked_actions = blocked_actions,
        questions = questions
    )
}

#' Check whether an object is ready for a protected integration analysis
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param design Optional named list with \code{sample}, \code{batch},
#'   \code{condition} and \code{subject} column names.
#' @return A read-only readiness and clarification summary.
#' @export
check_integration_readiness <- function(object, design = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    data <- SummarizedExperiment::colData(object)
    columns <- colnames(data)
    design <- if (is.null(design)) list() else design
    resolve <- function(role) {
        value <- design[[role]]
        if (length(value) == 1L && is.character(value) && value %in% columns) return(list(status = "confirmed", column = value))
        list(status = "missing", column = NULL)
    }
    resolved <- lapply(c("sample", "batch", "condition", "subject"), resolve)
    names(resolved) <- c("sample", "batch", "condition", "subject")
    has_assay <- length(SummarizedExperiment::assayNames(object)) > 0L
    batch <- resolved$batch
    condition <- resolved$condition
    status <- if (!has_assay) {
        "not_ready"
    } else if (is.null(design) || !identical(batch$status, "confirmed")) {
        "clarification_required"
    } else if (!is.null(condition$column) && identical(batch$column, condition$column)) {
        "not_ready"
    } else {
        "ready_for_diagnostic"
    }
    list(
        status = status,
        design = resolved,
        protected_variables = if (!is.null(condition$column)) condition$column else character(),
        checks = list(
            assay_available = has_assay,
            batch_confirmed = identical(batch$status, "confirmed"),
            condition_distinct_from_batch = is.null(condition$column) || !identical(batch$column, condition$column),
            raw_values_included = FALSE
        ),
        blocked_actions = if (identical(status, "clarification_required")) "integration" else character(),
        questions = if (identical(status, "clarification_required")) "Which colData column represents technical batch?" else character()
    )
}

