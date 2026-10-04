#' Build a bounded, read-only dataset profile for AI planning
#'
#' `GetAIProfile()` describes the structure and analysis readiness of a
#' `SingleCellExperiment` without returning assay values, cell identifiers, or
#' metadata values. Metadata design fields are name-based candidates only and
#' require user confirmation before they can be treated as semantic variables.
#'
#' @param object A `SingleCellExperiment` object.
#' @param max_sparse_elements Maximum number of sparse assay elements for the
#'   bounded `nnzero` summary. Larger sparse assays are reported as skipped.
#' @return A value-free, versioned dataset profile.
#' @export
GetAIProfile <- function(object, max_sparse_elements = 1000000L) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    max_sparse_elements <- as.numeric(max_sparse_elements)
    if (length(max_sparse_elements) != 1L ||
        !is.finite(max_sparse_elements) || max_sparse_elements < 0) {
        stop("max_sparse_elements must be one non-negative finite number", call. = FALSE)
    }

    col_data <- SummarizedExperiment::colData(object)
    row_data <- SummarizedExperiment::rowData(object)
    assay_names <- sclet_ai_profile_names(
        SummarizedExperiment::assayNames(object)
    )
    col_data_columns <- sclet_ai_profile_names(colnames(col_data))
    row_data_columns <- sclet_ai_profile_names(colnames(row_data))
    status <- sclet_ai_profile_status(object)
    ledger <- sclet_ai_profile_ledger(object)
    assays <- sclet_ai_profile_assays(object, assay_names, max_sparse_elements)

    design <- sclet_ai_profile_design(col_data, col_data_columns)
    health <- sclet_ai_profile_health(status$health)
    qc_candidates <- sclet_ai_profile_candidates(
        col_data,
        col_data_columns,
        c(
            "ncount", "ncounts", "totalcount", "totalcounts", "ngene",
            "ngenes", "nfeature", "nfeatures", "percentmt", "pctmt",
            "mito", "doublet", "cellcycle"
        )
    )
    structure <- sclet_ai_profile_structure(
        object = object,
        assay_names = assay_names,
        assays = assays,
        col_data_columns = col_data_columns,
        row_data_columns = row_data_columns,
        status = status
    )

    profile <- list(
        schema_version = "1.0",
        dataset = list(
            class = class(object),
            n_cells = as.integer(ncol(object)),
            n_genes = as.integer(nrow(object)),
            assays = assays,
            col_data = list(
                n_columns = length(col_data_columns),
                columns = col_data_columns
            ),
            row_data = list(
                n_columns = length(row_data_columns),
                columns = row_data_columns
            ),
            metadata_keys = sclet_ai_profile_names(names(S4Vectors::metadata(object)))
        ),
        design = design,
        qc = list(
            status = if (length(health)) "observed" else "not_available",
            health = health,
            candidate_metrics = qc_candidates,
            by_group = list(
                status = "not_available",
                reason = "design_candidates_require_user_confirmation"
            )
        ),
        structure = structure,
        analysis_state = sclet_ai_profile_analysis_state(status, ledger),
        diagnostics = sclet_ai_profile_diagnostics(
            design = design,
            col_data = col_data,
            col_data_columns = col_data_columns
        ),
        capabilities = sclet_ai_profile_capabilities(
            assay_names = assay_names,
            assays = assays,
            structure = structure
        ),
        cost_estimates = sclet_ai_profile_costs(
            assays,
            max_sparse_elements = max_sparse_elements
        ),
        privacy = list(
            assay_values_included = FALSE,
            complete_matrix_included = FALSE,
            metadata_values_included = FALSE,
            cell_identifiers_included = FALSE,
            gene_identifiers_included = FALSE,
            bounded_assay_summary = TRUE,
            user_confirmation_required_for_design = TRUE
        )
    )
    profile$fingerprint <- sclet_ai_profile_fingerprint(profile)
    profile
}

sclet_ai_profile_names <- function(value) {
    if (is.null(value)) {
        return(character())
    }
    value <- unique(as.character(value))
    value <- value[!is.na(value) & nzchar(value)]
    sort(value)
}

sclet_ai_profile_status <- function(object) {
    getter <- get0("Status", mode = "function", inherits = TRUE)
    if (is.null(getter)) {
        return(list())
    }
    tryCatch(getter(object, details = FALSE), error = function(e) list())
}

sclet_ai_profile_ledger <- function(object) {
    getter <- get0("GetAnalysisLedger", mode = "function", inherits = TRUE)
    if (is.null(getter)) {
        return(list())
    }
    tryCatch(
        getter(
            object,
            detail = "summary",
            include_artifacts = FALSE,
            include_data = FALSE
        ),
        error = function(e) list()
    )
}

sclet_ai_profile_assays <- function(object, assay_names, max_sparse_elements) {
    if (!length(assay_names)) {
        return(list())
    }
    result <- lapply(assay_names, function(name) {
        value <- tryCatch(
            SummarizedExperiment::assay(object, name, withDimnames = FALSE),
            error = function(e) NULL
        )
        dimensions <- if (is.null(value)) NULL else dim(value)
        if (!is.null(dimensions)) {
            dimensions <- as.integer(dimensions)
        }
        sparse <- !is.null(value) && inherits(value, "sparseMatrix")
        total_elements <- if (is.null(dimensions)) {
            NA_real_
        } else {
            as.double(dimensions[[1L]]) * as.double(dimensions[[2L]])
        }
        bounded_summary <- list(status = "not_requested")
        if (sparse) {
            if (is.finite(total_elements) && total_elements <= max_sparse_elements) {
                nnzero <- tryCatch(Matrix::nnzero(value), error = function(e) NA_real_)
                bounded_summary <- list(
                    status = if (is.na(nnzero)) "not_available" else "available",
                    nnzero = if (is.na(nnzero)) NULL else as.numeric(nnzero),
                    density = if (is.na(nnzero) || total_elements == 0) {
                        NULL
                    } else {
                        as.numeric(nnzero / total_elements)
                    },
                    max_elements = max_sparse_elements
                )
            } else {
                bounded_summary <- list(
                    status = "skipped",
                    reason = "size_limit",
                    max_elements = max_sparse_elements
                )
            }
        }
        list(
            name = name,
            class = if (is.null(value)) character() else class(value),
            dimensions = dimensions,
            sparse = sparse,
            total_elements = total_elements,
            bounded_summary = bounded_summary,
            values_included = FALSE
        )
    })
    names(result) <- assay_names
    result
}

sclet_ai_profile_column_info <- function(col_data, column) {
    value <- col_data[[column]]
    missing_fraction <- tryCatch({
        missing <- is.na(value)
        as.numeric(mean(as.logical(missing)))
    }, error = function(e) NA_real_)
    n_unique <- tryCatch(length(unique(value)), error = function(e) NA_integer_)
    list(
        column = column,
        class = class(value),
        missing_fraction = missing_fraction,
        n_unique = as.integer(n_unique)
    )
}

sclet_ai_profile_candidates <- function(col_data, columns, aliases) {
    if (!length(columns)) {
        return(list())
    }
    normalized <- gsub("[^a-z0-9]", "", tolower(columns))
    matched <- columns[normalized %in% aliases]
    if (!length(matched)) {
        return(list())
    }
    result <- lapply(matched, function(column) {
        info <- sclet_ai_profile_column_info(col_data, column)
        c(
            info,
            list(
                status = "candidate",
                semantic = "unknown",
                evidence = "column_name_pattern"
            )
        )
    })
    names(result) <- matched
    result
}

sclet_ai_profile_design <- function(col_data, columns) {
    roles <- list(
        sample = c("sample", "sampleid", "library", "libraryid", "origident"),
        batch = c("batch", "batchid", "batchname", "run", "runid"),
        condition = c(
            "condition", "conditionid", "group", "groupid", "treatment",
            "treatmentid", "casecontrol", "disease", "status"
        ),
        subject = c(
            "subject", "subjectid", "donor", "donorid", "patient",
            "patientid", "individual", "individualid"
        )
    )
    result <- lapply(roles, function(aliases) {
        candidates <- sclet_ai_profile_candidates(col_data, columns, aliases)
        candidate_columns <- sclet_ai_profile_names(names(candidates))
        status <- if (length(candidate_columns)) {
            "candidate"
        } else if (!length(columns)) {
            "missing"
        } else {
            "unknown"
        }
        list(
            status = status,
            semantic = if (length(candidate_columns)) "unknown" else NULL,
            column = if (length(candidate_columns) == 1L) candidate_columns else NULL,
            columns = candidate_columns,
            candidates = candidates,
            confirmation_required = TRUE
        )
    })
    result$metadata <- list(
        status = if (length(columns)) "observed" else "missing",
        n_columns = length(columns),
        columns = columns,
        values_included = FALSE
    )
    result
}

sclet_ai_profile_health <- function(value) {
    if (!is.list(value)) {
        return(list())
    }
    known <- c(
        "has_spliced_assay", "has_unspliced_assay", "has_velocity",
        "has_trajectory", "has_cellrank", "has_fate", "has_perturbation",
        "has_perturbation_priority", "has_rare_cells"
    )
    result <- lapply(known, function(name) {
        item <- value[[name]]
        if (is.null(item)) NULL else isTRUE(item)
    })
    names(result) <- known
    result[!vapply(result, is.null, logical(1L))]
}

sclet_ai_profile_structure <- function(
    object,
    assay_names,
    assays,
    col_data_columns,
    row_data_columns,
    status
) {
    reductions <- sclet_ai_profile_names(
        tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    )
    reduction_summary <- lapply(reductions, function(name) {
        value <- tryCatch(SingleCellExperiment::reducedDim(object, name), error = function(e) NULL)
        list(
            name = name,
            class = if (is.null(value)) character() else class(value),
            dimensions = if (is.null(value)) NULL else as.integer(dim(value)),
            values_included = FALSE
        )
    })
    names(reduction_summary) <- reductions

    alt_names <- sclet_ai_profile_names(
        tryCatch(SingleCellExperiment::altExpNames(object), error = function(e) character())
    )
    alt_summary <- lapply(alt_names, function(name) {
        value <- tryCatch(SingleCellExperiment::altExp(object, name), error = function(e) NULL)
        list(
            name = name,
            class = if (is.null(value)) character() else class(value),
            n_cells = if (is.null(value)) NULL else as.integer(ncol(value)),
            n_features = if (is.null(value)) NULL else as.integer(nrow(value)),
            assays = if (is.null(value)) character() else sclet_ai_profile_names(
                SummarizedExperiment::assayNames(value)
            ),
            values_included = FALSE
        )
    })
    names(alt_summary) <- alt_names

    available <- if (is.list(status$available)) status$available else list()
    layers <- sclet_ai_profile_names(available$layers)
    graphs <- sclet_ai_profile_names(available$graphs)
    list(
        assays = assays,
        reductions = reduction_summary,
        alternative_experiments = alt_summary,
        layers = layers,
        graphs = graphs,
        col_data_columns = col_data_columns,
        row_data_columns = row_data_columns,
        complete_matrix_included = FALSE
    )
}

sclet_ai_profile_analysis_state <- function(status, ledger) {
    available <- if (is.list(status$available)) status$available else list()
    active <- if (is.list(ledger$active_view)) ledger$active_view else list()
    analyses <- if (is.list(ledger$analyses)) names(ledger$analyses) else character()
    state_records <- if (is.list(ledger$state_records)) names(ledger$state_records) else character()
    workflows <- if (is.list(ledger$workflows)) names(ledger$workflows) else character()
    list(
        ledger = list(
            available = isTRUE(length(ledger)),
            schema_version = if (is.null(ledger$schema_version)) NULL else as.character(ledger$schema_version)
        ),
        active = list(
            present = length(active) > 0L,
            fields = sclet_ai_profile_names(names(active))
        ),
        recorded = list(
            analyses = sclet_ai_profile_names(analyses),
            state_records = sclet_ai_profile_names(state_records),
            workflows = sclet_ai_profile_names(workflows),
            available_analyses = sclet_ai_profile_names(available$analyses)
        ),
        n_commands = if (is.null(status$n_commands)) NULL else as.integer(status$n_commands),
        values_included = FALSE
    )
}

sclet_ai_profile_diagnostics <- function(design, col_data, col_data_columns) {
    cluster <- sclet_ai_profile_candidates(
        col_data,
        col_data_columns,
        c("cluster", "clusters", "leiden", "louvain", "seuratclusters", "celltype", "celllabel")
    )
    annotation <- sclet_ai_profile_candidates(
        col_data,
        col_data_columns,
        c("annotation", "annot", "celltype", "celllabel", "label", "labels")
    )
    diagnostic <- function(role, candidates, reason) {
        columns <- sclet_ai_profile_names(names(candidates))
        list(
            status = if (length(columns)) "candidate" else if (
                identical(design$metadata$status, "missing")
            ) "missing" else "unknown",
            semantic = if (length(columns)) "unknown" else NULL,
            columns = columns,
            evidence = if (length(columns)) "column_name_pattern" else reason,
            confirmation_required = TRUE
        )
    }
    list(
        batch = list(
            status = design$batch$status,
            candidate_columns = design$batch$columns,
            semantic = design$batch$semantic,
            evidence = if (length(design$batch$columns)) {
                "design_batch_candidate"
            } else {
                "no_confirmed_batch_design"
            }
        ),
        rare_cluster = diagnostic(
            "cluster",
            cluster,
            "no_cluster_or_cell_group_candidate"
        ),
        annotation_readiness = diagnostic(
            "annotation",
            annotation,
            "no_annotation_candidate"
        ),
        inference_policy = "name_based_candidates_only"
    )
}

sclet_ai_profile_capabilities <- function(assay_names, assays, structure) {
    list(
        read_only = TRUE,
        execution = FALSE,
        profile_generation = TRUE,
        assays = assay_names,
        reductions = sclet_ai_profile_names(names(structure$reductions)),
        layers = structure$layers,
        graphs = structure$graphs,
        sparse_summary = TRUE,
        complete_matrix_transfer = FALSE,
        metadata_value_transfer = FALSE,
        values_included = FALSE
    )
}

sclet_ai_profile_costs <- function(assays, max_sparse_elements) {
    total_elements <- vapply(assays, function(assay) {
        value <- assay$total_elements
        if (length(value) && is.finite(value)) value else 0
    }, numeric(1L))
    sparse <- vapply(assays, function(assay) isTRUE(assay$sparse), logical(1L))
    list(
        operation = "bounded_structural_profile",
        assay_count = length(assays),
        total_assay_elements = sum(total_elements),
        sparse_assay_count = sum(sparse),
        sparse_summary_limit = max_sparse_elements,
        matrix_transfer_elements = 0,
        complete_matrix_read = FALSE
    )
}

sclet_ai_profile_fingerprint <- function(value) {
    raw <- serialize(value, connection = NULL, version = 2)
    if (!length(raw)) {
        return("sclet-ai-profile-00000000")
    }
    weights <- seq_along(raw) %% 100003
    checksum <- sum((as.double(raw) + 1) * (weights + 1)) %% 2147483647
    sprintf("sclet-ai-profile-%08x", as.integer(checksum))
}
