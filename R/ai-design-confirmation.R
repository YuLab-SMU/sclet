ConfirmAIDesignSemantics <- function(object, design) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    if (!is.list(design) || is.null(names(design))) {
        stop("`design` must be a named list mapping roles (batch/condition/subject/sample) to colData column names.")
    }
    known_roles <- c("batch", "condition", "subject", "sample")
    roles <- names(design)
    if (is.null(roles) || !all(nzchar(roles))) {
        stop("`design` must be a named list with non-empty names.")
    }
    unknown_roles <- setdiff(roles, known_roles)
    if (length(unknown_roles)) {
        stop(
            "Unsupported design role(s): ",
            paste(unknown_roles, collapse = ", "),
            ". Supported roles are: batch, condition, subject, sample."
        )
    }
    if (!("batch" %in% roles)) {
        stop("`design` must contain at minimum a 'batch' role.")
    }
    cd <- SummarizedExperiment::colData(object)
    cd_names <- colnames(cd)
    for (role in roles) {
        col <- design[[role]]
        if (!is.character(col) || length(col) != 1L || is.na(col) || !nzchar(col)) {
            stop("Design role '", role, "' must be a single non-empty column name (character scalar).")
        }
        if (!col %in% cd_names) {
            stop("Design role '", role, "' points to column '", col, "' which is not present in colData(object). Available columns: ", paste(cd_names, collapse = ", "))
        }
    }
    design_sorted <- design[order(roles)]
    id <- paste(
        c("design", paste(sprintf("%s=%s", names(design_sorted), design_sorted), collapse = ",")),
        collapse = ";"
    )
    value_key <- sclet_ai_design_value_key(object, design_sorted)
    design_fp <- as.character(sclet_ai_design_fingerprint(object))
    ledger_fp_pre <- as.character(GetAnalysisLedger(object, detail = "summary", include_artifacts = FALSE, include_data = FALSE)$fingerprint)
    confirmed_at <- as.character(Sys.time())
    sclet_set_analysis_state(
        object = object,
        type = "ai_design_confirmation",
        id = id,
        method = "ConfirmAIDesignSemantics",
        inputs = design_sorted,
        summary = list(
            design = design_sorted,
            design_value_key = value_key,
            object_fingerprint = design_fp,
            ledger_fingerprint_before = ledger_fp_pre,
            confirmed_at = confirmed_at
        ),
        active = FALSE
    )
}

sclet_ai_design_value_key <- function(object, design_role_map) {
    cols <- unique(as.character(unlist(design_role_map)))
    cols <- cols[nzchar(cols)]
    if (!length(cols)) {
        return(NULL)
    }
    cd <- tryCatch(SummarizedExperiment::colData(object), error = function(e) NULL)
    if (is.null(cd) || !all(cols %in% colnames(cd))) {
        return(NULL)
    }
    hashes <- vapply(cols, function(col) {
        v <- cd[[col]]
        if (is.factor(v)) {
            paste(levels(v), as.integer(v), collapse = "|", sep = ":")
        } else if (is.atomic(v)) {
            if (length(v) <= 200L) {
                paste(as.character(v), collapse = "\u0001")
            } else {
                val <- as.character(v)
                n <- length(val)
                paste(
                    c("n=", as.character(n), "|head=", paste(val[seq_len(10L)], collapse = ","), "|tail=", paste(val[seq.int(to = n, length.out = 10L)], collapse = ","), "|uniq=", paste(unique(val), collapse = ",")),
                    collapse = ""
                )
            }
        } else {
            paste(class(v), collapse = ",")
        }
    }, character(1L), USE.NAMES = TRUE)
    out <- as.list(hashes)
    names(out) <- cols
    out
}

sclet_ai_design_fingerprint <- function(object) {
    tryCatch({
        cd <- SummarizedExperiment::colData(object)
        sclet_ai_fingerprint(list(
            n_cells = as.integer(ncol(object)),
            n_features = as.integer(nrow(object)),
            cd_names = as.character(colnames(cd)),
            has_rowdata = !is.null(tryCatch(
                SummarizedExperiment::rowData(object),
                error = function(e) NULL
            ))
        ))
    }, error = function(e) NULL)
}

sclet_ai_find_design_confirmation <- function(object, batch) {
    if (!inherits(object, "SingleCellExperiment")) {
        return(NULL)
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        return(NULL)
    }
    records <- tryCatch(
        sclet_get_state_records(object, "ai_design_confirmation"),
        error = function(e) list()
    )
    if (!length(records)) {
        return(NULL)
    }
    design_fp <- as.character(sclet_ai_design_fingerprint(object))
    for (id in names(records)) {
        rec <- records[[id]]
        summary <- rec$summary %||% list()
        if (!identical(as.character(summary$object_fingerprint %||% ""), design_fp)) {
            next
        }
        design_role_map <- summary$design %||% list()
        if (!identical(as.character(design_role_map$batch %||% ""), as.character(batch))) {
            next
        }
        if (!is.null(summary$design_value_key)) {
            cur_key <- sclet_ai_design_value_key(object, design_role_map)
            if (is.null(cur_key) || !identical(cur_key, summary$design_value_key)) {
                next
            }
        }
        return(rec)
    }
    NULL
}
