## End-to-end RunAIAnalysis() coverage for opt-in domains outside the
## beginner default registry.
##
## Proves that communication and malignant-cell workflows complete a real
## context -> plan -> execute -> report loop through RunAIAnalysis(), not just
## ValidateAIPlan()/ExecuteAIPlan() on hand-built plans.

sclet_ai_test_e2e_cci_object <- function(n_cells = 40L) {
    genes <- c("CXCL12", "ACKR1", "HLA-A", "CD8A", "APP", "CD74", "ACTB", "GAPDH")
    counts <- matrix(stats::rpois(length(genes) * n_cells, lambda = 8L),
        nrow = length(genes), ncol = n_cells,
        dimnames = list(genes, paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    sce$cluster <- factor(rep(c("beta cells", "alpha cells"), length.out = n_cells),
        levels = c("beta cells", "alpha cells"))
    ActiveIdent(sce) <- "cluster"
    sce
}

sclet_ai_test_e2e_cci_table <- function(method = "CellChat") {
    data.frame(
        source_group = c("beta cells", "beta cells", "alpha cells", "alpha cells", "beta cells"),
        target_group = c("alpha cells", "alpha cells", "beta cells", "beta cells", "beta cells"),
        ligand = c("CXCL12", "APP", "HLA-A", "CXCL12", "ACKR1"),
        receptor = c("ACKR1", "CD74", "CD8A", "ACKR1", "CXCL12"),
        pathway = c("CXCL12_ACKR1", "APP_CD74", "HLA-A_CD8A", "CXCL12_ACKR1", "ACKR1_CXCL12"),
        probability = c(0.91, 0.62, 0.44, 0.31, 0.08),
        p_value = c(0.001, 0.02, 0.34, 0.48, 0.99),
        method = method,
        stringsAsFactors = FALSE
    )
}

sclet_ai_test_e2e_cci_store <- function(object, id = "cci1", method = "CellChat",
                                        analysis_key = "cellchat",
                                        table = sclet_ai_test_e2e_cci_table(method),
                                        group = "cluster") {
    object <- sclet:::sclet_set_analysis(object, analysis_key, list(
        id = id,
        method = method,
        object = NULL,
        interactions = table,
        group = group
    ))
    sclet:::sclet_set_analysis_state(
        object = object,
        type = "communication",
        id = id,
        method = method,
        inputs = list(assay = "counts", layer = NULL, group = group),
        artifacts = list(analysis_key = analysis_key),
        params = list(species = "human"),
        summary = list(n_groups = length(unique(table$source_group))),
        active = TRUE
    )
}

sclet_ai_test_e2e_cci_engine_stub <- function(id = "cci1", method = "CellChat",
                                              analysis_key = "cellchat") {
    forced_method <- method
    function(object, method = NULL, ...) {
        sclet_ai_test_e2e_cci_store(object, id = id,
            method = method %||% forced_method,
            analysis_key = analysis_key)
    }
}

sclet_ai_test_e2e_pin_communication_engines <- function(env = parent.frame()) {
    real_require_namespace <- base::requireNamespace
    testthat::local_mocked_bindings(
        requireNamespace = function(package, ...) {
            if (identical(package, "CellChat")) return(TRUE)
            if (identical(package, "basilisk")) return(TRUE)
            if (identical(package, "nichenetr")) return(FALSE)
            real_require_namespace(package, ...)
        },
        .package = "base",
        .env = env
    )
}

sclet_ai_test_e2e_malignant_object <- function() {
    set.seed(1L)
    counts <- matrix(stats::rnbinom(50L * 20L, mu = 5, size = 2),
        nrow = 50L, ncol = 20L,
        dimnames = list(paste0("g", seq_len(50L)), paste0("c", seq_len(20L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    SummarizedExperiment::colData(sce)$reference <- rep(c(TRUE, FALSE), each = 10L)
    SummarizedExperiment::colData(sce)$malignant_score <- c(
        seq(2, 0.5, length.out = 6L),
        seq(0.4, -0.4, length.out = 4L),
        seq(-0.5, -2, length.out = 10L)
    )
    sce
}

sclet_ai_test_e2e_cnv_stub <- function(default_name = "cnv") {
    function(object,
             reference_cells = NULL,
             reference_col = NULL,
             species = "human",
             assay = "counts",
             hmm = TRUE,
             cutoff = 1,
             num_threads = 1,
             name = default_name,
             ...) {
        actual_name <- name %||% default_name
        cdata <- SummarizedExperiment::colData(object)
        ref_col <- reference_col %||% "reference"
        burden <- c(seq(0.12, 0.09, length.out = 6L),
            seq(0.085, 0.07, length.out = 4L),
            seq(0.05, 0.02, length.out = ncol(object) - 10L))
        cdata[[paste0(actual_name, "_burden")]] <- burden
        cdata[[paste0(actual_name, "_is_reference")]] <- as.logical(cdata[[ref_col]])
        SummarizedExperiment::colData(object) <- cdata

        record <- list(
            status = "completed",
            burden = data.frame(cell = colnames(object), burden = burden),
            n_cells_submitted = as.integer(ncol(object)),
            n_cells_retained = as.integer(ncol(object)),
            dropped_cells = character(),
            n_reference_cells = as.integer(sum(as.logical(cdata[[ref_col]]))),
            n_genes = as.integer(nrow(object)),
            params = list(hmm = isTRUE(hmm), cutoff = cutoff),
            created_at = Sys.time()
        )
        existing <- sclet:::sclet_get_analysis(object, "malignant_cnv")
        if (is.null(existing)) existing <- list()
        existing[[actual_name]] <- record
        object <- sclet:::sclet_set_analysis(object, "malignant_cnv", existing)
        object <- sclet:::sclet_set_state_record(
            object = object,
            type = "cnv",
            id = actual_name,
            value = list(
                method = "inferCNV per-cell copy-number burden",
                inputs = list(assay = assay),
                artifacts = list(columns = paste0(actual_name, c("_burden", "_is_reference"))),
                params = record$params,
                summary = list(
                    status = "completed",
                    n_cells = as.integer(ncol(object)),
                    n_cells_retained = as.integer(ncol(object)),
                    n_cells_dropped = 0L,
                    n_reference_cells = record$n_reference_cells,
                    burden_q = unname(stats::quantile(burden,
                        c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
                ),
                created_at = Sys.time()
            )
        )
        object
    }
}

sclet_ai_test_e2e_call_malignant_stub <- function(default_name = "malignant_call") {
    function(object,
             cnv_id = "cnv",
             signature = NULL,
             signature_score = NULL,
             transcriptomic_cut = NULL,
             cnv_cut = NULL,
             reference_quantile = 0.99,
             identity_genes = NULL,
             identity_demote = 0.1,
             assay = NULL,
             name = default_name,
             ...) {
        actual_name <- name %||% default_name
        cdata <- SummarizedExperiment::colData(object)
        score_col <- signature_score %||% "malignant_score"
        scores <- cdata[[score_col]]
        state <- c(
            rep("malignant", 6L),
            rep("likely_malignant", 4L),
            rep("non_malignant", ncol(object) - 10L)
        )
        cdata[[paste0(actual_name, "_state")]] <- factor(state,
            levels = c("malignant", "likely_malignant", "undetermined", "non_malignant"))
        cdata[[paste0(actual_name, "_signature_score")]] <- scores
        SummarizedExperiment::colData(object) <- cdata

        state_counts <- table(state)
        object <- sclet:::sclet_set_state_record(
            object = object,
            type = "malignancy",
            id = actual_name,
            value = list(
                method = "four-state malignant-cell call from CNV + signature evidence",
                inputs = list(cnv_id = cnv_id),
                artifacts = list(columns = paste0(actual_name, c("_state", "_signature_score"))),
                summary = list(
                    n_malignant = as.integer(unname(state_counts[["malignant"]])),
                    n_likely_malignant = as.integer(unname(state_counts[["likely_malignant"]])),
                    n_undetermined = 0L,
                    n_non_malignant = as.integer(unname(state_counts[["non_malignant"]])),
                    cnv_available = TRUE,
                    identity_conflict = FALSE,
                    signature_cut = 0,
                    signature_cut_is_median = 0L,
                    n_signature_positive = as.integer(sum(scores >= 0, na.rm = TRUE)),
                    n_cnv_positive = 10L,
                    cnv_reference_threshold = 0.08,
                    n_pre_demote_supported = 10L
                ),
                created_at = Sys.time()
            )
        )
        object
    }
}

test_that("run_cci survives RunAIAnalysis end to end with bounded evidence", {
    sclet_ai_test_e2e_pin_communication_engines()
    object <- sclet_ai_test_e2e_cci_object()
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Infer cell-cell interactions with CellChat.",
            proposed_actions = list(list(
                id = "cc1",
                action = "run_cci",
                params = list(method = "CellChat", name = "e2e_cci")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(object, include = c("read", "communication"))
    result <- with_mocked_bindings(
        RunAIAnalysis(object, goal = "Infer cell-cell interactions between the clusters.",
            confirm = "yes", registry = registry),
        RunCCI = sclet_ai_test_e2e_cci_engine_stub(id = "e2e_cci", method = "CellChat"),
        .package = "sclet"
    )

    expect_equal(result$status, "completed")
    expect_equal(result$execution$status, "completed")
    expect_equal(vapply(result$plan$actions, function(x) x$action, character(1L)), "run_cci")
    expect_true(has_cci(result$object, "e2e_cci"))

    evidence <- sclet:::sclet_ai_evidence_get_all(result$object)
    node <- evidence[["ev:communication_e2e_cci"]]
    expect_false(is.null(node))
    expect_false(isTRUE(node$values$raw_values_included))
    expect_true(isTRUE(node$values$engine_was_named_by_caller))
    expect_equal(get_cci(result$object, id = "e2e_cci")$method[[1L]], "CellChat")
})

test_that("malignant two-step plan survives RunAIAnalysis end to end with bounded evidence", {
    real_require_namespace <- base::requireNamespace
    testthat::local_mocked_bindings(
        requireNamespace = function(package, ...) {
            if (identical(package, "infercnv")) return(TRUE)
            real_require_namespace(package, ...)
        },
        .package = "base"
    )
    object <- sclet_ai_test_e2e_malignant_object()
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Infer copy-number burden, then call malignant cells.",
            proposed_actions = list(
                list(id = "cnv1", action = "run_cnv_inference",
                    params = list(reference_col = "reference", assay = "counts", name = "e2e_cnv")),
                list(id = "call1", action = "call_malignant_cells",
                    params = list(cnv_id = "e2e_cnv", signature_score = "malignant_score",
                        name = "e2e_malignant"))
            )
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(object, include = c("read", "malignant"))
    result <- with_mocked_bindings(
        RunAIAnalysis(object, goal = "Identify malignant cells in this tumour sample.",
            confirm = "yes", registry = registry),
        RunCNVInference = sclet_ai_test_e2e_cnv_stub(default_name = "e2e_cnv"),
        CallMalignantCells = sclet_ai_test_e2e_call_malignant_stub(default_name = "e2e_malignant"),
        .package = "sclet"
    )

    expect_equal(result$status, "completed")
    expect_equal(result$execution$status, "completed")
    expect_equal(vapply(result$plan$actions, function(x) x$action, character(1L)),
        c("run_cnv_inference", "call_malignant_cells"))

    evidence <- sclet:::sclet_ai_evidence_get_all(result$object)
    cnv_node <- evidence[["ev:cnv_e2e_cnv"]]
    malignant_node <- evidence[["ev:malignant_e2e_malignant"]]
    expect_false(is.null(cnv_node))
    expect_false(is.null(malignant_node))
    expect_equal(malignant_node$claim_level, "associated")

    per_cell <- any(vapply(malignant_node$values, function(v) {
        is.atomic(v) && length(v) >= ncol(result$object)
    }, logical(1L)))
    expect_false(isTRUE(per_cell))

    summary <- sclet:::sclet_get_state_record(result$object, "malignancy", "e2e_malignant")$summary
    expect_equal(summary$n_malignant, 6L)
    expect_equal(summary$n_likely_malignant, 4L)
    expect_true(isTRUE(summary$cnv_available))
})
