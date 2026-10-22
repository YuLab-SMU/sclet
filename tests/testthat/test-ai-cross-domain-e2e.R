## End-to-end RunAIAnalysis() coverage for opt-in domains outside the
## beginner default registry.
##
## Proves that communication and malignant-cell workflows complete a real
## context -> plan -> execute -> report loop through RunAIAnalysis(), not just
## ValidateAIPlan()/ExecuteAIPlan() on hand-built plans.
##
## The communication fixtures live in helper-ai-fixtures.R, which testthat
## loads for this file whether it is run alone or as part of the suite.

## The malignant fixture's design, declared ONCE here and then used to build
## the per-cell values. Assertions compare against these declared numbers, not
## against a helper the stub also calls: a check whose two sides share one
## function cannot fail, and the earlier version of this file had exactly that
## (setting the state mapping to a single label left every assertion green).
##   * 6 malignant, 4 likely_malignant, 1 undetermined, 9 non_malignant
##   * 5 normal cells are the reference population and carry the LOWEST burden
sclet_ai_test_e2e_malignant_design <- function() {
    list(
        n_malignant = 6L,
        n_likely_malignant = 4L,
        n_undetermined = 1L,
        n_non_malignant = 9L,
        n_reference = 5L
    )
}

## Expression scores for the object, in cell order: malignant plateau, a band
## that straddles the cutoff, one NA (the cell the call cannot classify), then
## the non-malignant tail.
##
## These counts are deliberately NOT read from sclet_ai_test_e2e_malignant_design():
## if the scores were derived from the same constants the assertions compare
## against, editing a declared count would silently reshape the fixture to agree
## with itself and the assertion could never fail (that was the defect this file
## had). They are independent literals instead, and the stopifnot below turns
## "someone edited only one of the two places" into a loud fixture error at
## setup time rather than a confusing count mismatch later.
sclet_ai_test_e2e_malignant_scores <- function() {
    n_plateau <- 6L
    n_straddle <- 4L
    n_tail <- 9L
    design <- sclet_ai_test_e2e_malignant_design()
    stopifnot(
        "malignant plateau count disagrees with the declared design" =
            n_plateau == design$n_malignant,
        "likely_malignant band count disagrees with the declared design" =
            n_straddle == design$n_likely_malignant,
        "non_malignant tail count disagrees with the declared design" =
            n_tail == design$n_non_malignant
    )
    c(
        seq(2.0, 0.6, length.out = n_plateau),
        seq(0.45, -0.45, length.out = n_straddle),
        NA_real_,
        seq(-0.6, -2.0, length.out = n_tail)
    )
}

## Thresholds applied to a score. This is the STUB's behaviour, so assertions
## must not call it on both sides -- use the design constants instead.
##
## Boundaries belong to the OUTER state: 0.5 is malignant, -0.5 is
## non_malignant, and anything left (including NA) is undetermined. Stated here
## because the three conditions below have to be read together to see that every
## input lands in exactly one state.
sclet_ai_test_e2e_malignant_state <- function(scores) {
    state <- rep("undetermined", length(scores))
    state[!is.na(scores) & scores >= 0.5] <- "malignant"
    state[!is.na(scores) & scores < 0.5 & scores > -0.5] <- "likely_malignant"
    state[!is.na(scores) & scores <= -0.5] <- "non_malignant"
    factor(state,
        levels = c("malignant", "likely_malignant", "undetermined", "non_malignant"))
}

sclet_ai_test_e2e_malignant_object <- function() {
    design <- sclet_ai_test_e2e_malignant_design()
    n_cells <- design$n_malignant + design$n_likely_malignant +
        design$n_undetermined + design$n_non_malignant
    set.seed(1L)
    counts <- matrix(stats::rnbinom(50L * n_cells, mu = 5, size = 2),
        nrow = 50L, ncol = n_cells,
        dimnames = list(paste0("g", seq_len(50L)), paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    ## Reference cells are the normal population, so they are the LAST cells:
    ## in a real tumour the reference must be the low-burden tail, not the
    ## malignant plateau. Marking the high-score cells as reference strung the
    ## CNV threshold above the whole distribution.
    SummarizedExperiment::colData(sce)$reference <- seq_len(n_cells) >
        (n_cells - design$n_reference)
    SummarizedExperiment::colData(sce)$malignant_score <-
        sclet_ai_test_e2e_malignant_scores()
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
        scores <- cdata[["malignant_score"]]
        if (is.null(scores)) scores <- seq_len(ncol(object))

        ## Burden tracks the score so the two lines of evidence agree on which
        ## cells are abnormal, and it is high on the malignant tail, low on the
        ## reference tail. A cell with an NA score keeps a mid-range burden: it
        ## is still measured, it is the call that cannot classify it.
        burden <- rep(0.05, ncol(object))
        finite <- is.finite(scores)
        finite_scores <- scores[finite]
        if (length(finite_scores)) {
            score_range <- range(finite_scores)
            if (diff(score_range) > 0) {
                burden[finite] <- 0.02 + 0.10 *
                    (finite_scores - score_range[[1L]]) / diff(score_range)
            }
        }
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
        state <- sclet_ai_test_e2e_malignant_state(scores)
        cdata[[paste0(actual_name, "_state")]] <- state
        cdata[[paste0(actual_name, "_signature_score")]] <- scores
        SummarizedExperiment::colData(object) <- cdata

        burden_col <- paste0(cnv_id, "_burden")
        reference_col <- paste0(cnv_id, "_is_reference")
        burden <- cdata[[burden_col]]
        reference_mask <- cdata[[reference_col]]
        cnv_reference_threshold <- NA_real_
        if (!is.null(burden) && !is.null(reference_mask)) {
            ref_burden <- burden[as.logical(reference_mask)]
            if (any(is.finite(ref_burden))) {
                cnv_reference_threshold <- as.numeric(stats::quantile(
                    ref_burden,
                    probs = reference_quantile,
                    na.rm = TRUE
                ))
            }
        }
        if (is.null(burden)) {
            cnv_positive <- rep(FALSE, length(scores))
        } else if (is.na(cnv_reference_threshold)) {
            cnv_positive <- is.finite(burden)
        } else {
            cnv_positive <- is.finite(burden) & burden > cnv_reference_threshold
        }

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
                    n_undetermined = as.integer(unname(state_counts[["undetermined"]])),
                    n_non_malignant = as.integer(unname(state_counts[["non_malignant"]])),
                    cnv_available = !is.null(burden),
                    identity_conflict = FALSE,
                    signature_cut = 0,
                    signature_cut_is_median = 0L,
                    n_signature_positive = as.integer(sum(scores >= 0, na.rm = TRUE)),
                    n_cnv_positive = as.integer(sum(cnv_positive, na.rm = TRUE)),
                    cnv_reference_threshold = cnv_reference_threshold,
                    n_pre_demote_supported = as.integer(sum(state != "non_malignant", na.rm = TRUE))
                ),
                created_at = Sys.time()
            )
        )
        object
    }
}

test_that("run_cci survives RunAIAnalysis end to end with bounded evidence", {
    sclet_ai_test_cci_pin_engines()
    object <- sclet_ai_test_cci_object()
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
        RunCCI = sclet_ai_test_cci_engine_stub(id = "e2e_cci", method = "CellChat"),
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

    expect_false(sclet_ai_test_has_per_cell(malignant_node$values, ncol(result$object)))

    ## The four-state call really is four states. `undetermined` is only
    ## reachable when a cell's score is missing, so the fixture now carries one
    ## such cell: a mapping that silently collapses to three states fails here.
    design <- sclet_ai_test_e2e_malignant_design()
    observed_state <- SummarizedExperiment::colData(result$object)$e2e_malignant_state
    expect_s3_class(observed_state, "factor")
    expect_equal(levels(observed_state),
        c("malignant", "likely_malignant", "undetermined", "non_malignant"))
    expect_length(observed_state, ncol(result$object))
    expect_setequal(as.character(observed_state), levels(observed_state))

    ## Counts against the DESIGN, not against the function the stub used to
    ## write the column: same-source assertions cannot fail.
    summary <- sclet:::sclet_get_state_record(result$object, "malignancy", "e2e_malignant")$summary
    expect_equal(summary$n_malignant, design$n_malignant)
    expect_equal(summary$n_likely_malignant, design$n_likely_malignant)
    expect_equal(summary$n_undetermined, design$n_undetermined)
    expect_equal(summary$n_non_malignant, design$n_non_malignant)
    expect_equal(
        summary$n_malignant + summary$n_likely_malignant +
            summary$n_undetermined + summary$n_non_malignant,
        ncol(result$object)
    )

    ## The CNV line is a real second observation, not an empty set: normal
    ## cells are the reference, so cells must sit above its threshold, and the
    ## abnormal ones must be among them. (The exact count is not pinned: the
    ## threshold is the reference quantile, which by construction leaves some
    ## normal cells just above it.)
    expect_true(isTRUE(summary$cnv_available))
    expect_false(is.na(summary$cnv_reference_threshold))
    expect_gt(summary$n_cnv_positive, 0L)
    expect_gte(summary$n_cnv_positive, summary$n_malignant)

    ## Product-side self-check (computed by the package from the ledger, not
    ## written by a stub): this round's report verified.
    expect_equal(result$report$verification_status, "passed")
})
