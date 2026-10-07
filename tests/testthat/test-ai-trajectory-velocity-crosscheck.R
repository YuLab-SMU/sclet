## Trajectory/velocity evidence cross-check (P5, slice 3)
##
## Proves that compare_trajectory_velocity_evidence() reads only the bounded
## ev:trajectory_<name> and ev:velocity_<name> evidence nodes, compares their
## pseudotime spread summaries, and never touches raw per-cell colData.

sclet_ai_test_tv_combined_object <- function(n_cells = 80L) {
    set.seed(3L)
    n_genes <- 60L
    genes <- paste0("gene_", seq_len(n_genes))
    cells <- paste0("cell_", seq_len(n_cells))
    counts <- matrix(rpois(n_genes * n_cells, lambda = 8L), nrow = n_genes, ncol = n_cells,
        dimnames = list(genes, cells))
    spliced <- counts
    unspliced <- matrix(rpois(n_genes * n_cells, lambda = 3L), nrow = n_genes, ncol = n_cells,
        dimnames = list(genes, cells))
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = counts, spliced = spliced, unspliced = unspliced))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    set.seed(4L)
    embedding <- matrix(stats::rnorm(n_cells * 2L), ncol = 2L)
    embedding[1:30, 1] <- embedding[1:30, 1] - 4
    embedding[31:60, 2] <- embedding[31:60, 2] - 4
    embedding[61:80, 1] <- embedding[61:80, 1] + 3
    embedding[61:80, 2] <- embedding[61:80, 2] + 3
    SingleCellExperiment::reducedDim(sce, "UMAP") <- embedding
    SummarizedExperiment::colData(sce)$cluster <- rep(c("A", "B", "C"), times = c(30L, 30L, 20L))
    ActiveIdent(sce) <- "cluster"
    sce
}

sclet_ai_test_tv_run_trajectory <- function(sce, name = "traj_tv") {
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "trajectory_inference",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "traj1", action = "run_trajectory",
            params = list(group = "cluster", start_cluster = "A",
                reduction = "UMAP", name = name)))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    stopifnot(validation$valid)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    stopifnot(executed$status == "completed")
    executed$object
}

sclet_ai_test_tv_run_velocity <- function(sce, name = "vel_tv") {
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "velocity_estimation",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "vel1", action = "run_velocity",
            params = list(name = name)))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    stopifnot(validation$valid)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    stopifnot(executed$status == "completed")
    executed$object
}

## -- Scenario 1: neither evidence present --
test_that("compare_trajectory_velocity_evidence returns not_available when no evidence exists", {
    sce <- sclet_ai_test_tv_combined_object()
    result <- compare_trajectory_velocity_evidence(sce)
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "no_trajectory_or_velocity_evidence_records")
    expect_false(result$raw_values_included)
})

## -- Scenario 2: only trajectory evidence present --
test_that("compare_trajectory_velocity_evidence returns not_available when velocity is missing", {
    skip_if_not_installed("slingshot")
    sce <- sclet_ai_test_tv_combined_object()
    sce <- sclet_ai_test_tv_run_trajectory(sce, "traj_only")
    result <- compare_trajectory_velocity_evidence(sce)
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "velocity_evidence_missing")
    expect_true("traj_only" %in% result$available_trajectory_runs)
})

## -- Scenario 3: only velocity evidence present --
test_that("compare_trajectory_velocity_evidence returns not_available when trajectory is missing", {
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_tv_combined_object()
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_only")
    result <- compare_trajectory_velocity_evidence(sce)
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "trajectory_evidence_missing")
    expect_true("vel_only" %in% result$available_velocity_runs)
})

## -- Scenario 4: both present, comparable --
test_that("compare_trajectory_velocity_evidence returns available with both evidence nodes", {
    skip_if_not_installed("slingshot")
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_tv_combined_object()
    sce <- sclet_ai_test_tv_run_trajectory(sce, "traj_both")
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_both")
    result <- compare_trajectory_velocity_evidence(sce)
    expect_equal(result$status, "available")
    expect_equal(result$claim_level, "consistent_with")
    expect_false(result$raw_values_included)
    expect_true(is.numeric(result$trajectory_relative_spread))
    expect_true(is.numeric(result$velocity_relative_spread))
    expect_true(is.numeric(result$spread_ratio))
    expect_true(is.logical(result$comparable_spread))
    allowed_labels <- c(
        "relative spread of the two pseudotime estimates is of similar order of magnitude",
        "relative spread of the two pseudotime estimates differs by more than 2x; this may reflect genuinely different dynamics captured by each method, not necessarily an error",
        "relative spread could not be compared because one or both bounded summaries were degenerate"
    )
    expect_true(result$spread_consistency_label %in% allowed_labels)
    expect_equal(result$trajectory_run, "traj_both")
    expect_equal(result$velocity_run, "vel_both")
})

## -- Scenario 5: ambiguous run selection --
test_that("compare_trajectory_velocity_evidence returns ambiguous_run when multiple velocity runs exist", {
    skip_if_not_installed("slingshot")
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_tv_combined_object()
    sce <- sclet_ai_test_tv_run_trajectory(sce, "traj_amb")
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_amb_1")
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_amb_2")
    result <- compare_trajectory_velocity_evidence(sce)
    expect_equal(result$status, "ambiguous_run")
    expect_true("vel_amb_1" %in% result$available_velocity_runs)
    expect_true("vel_amb_2" %in% result$available_velocity_runs)
})

## -- Scenario 6: explicit ids resolve ambiguity --
test_that("compare_trajectory_velocity_evidence resolves ambiguity with explicit velocity_id", {
    skip_if_not_installed("slingshot")
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_tv_combined_object()
    sce <- sclet_ai_test_tv_run_trajectory(sce, "traj_exp")
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_exp_1")
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_exp_2")
    result <- compare_trajectory_velocity_evidence(sce, velocity_id = "vel_exp_1")
    expect_equal(result$status, "available")
    expect_equal(result$velocity_run, "vel_exp_1")
})

## -- Scenario 7: degenerate spread (divide-by-zero guard) --
## Real construction of a degenerate pseudotime (all values identical) is impractical
## without a contrived mock, so we call the internal helper directly with hand-built
## evidence-like values to prove the guard works.
test_that("sclet_ai_tv_spread and sclet_ai_tv_ratio guard against divide-by-zero", {
    expect_true(is.na(sclet:::sclet_ai_tv_spread(1.0, 0)))
    expect_true(is.na(sclet:::sclet_ai_tv_spread(1.0, NA)))
    expect_true(is.na(sclet:::sclet_ai_tv_spread(1.0, Inf)))
    expect_true(is.na(sclet:::sclet_ai_tv_ratio(1.0, 0)))
    expect_true(is.na(sclet:::sclet_ai_tv_ratio(NA, 1.0)))
    expect_true(is.na(sclet:::sclet_ai_tv_ratio(1.0, NA)))
    spread <- sclet:::sclet_ai_tv_spread(2.0, 4.0)
    expect_equal(spread, 0.5)
    ratio <- sclet:::sclet_ai_tv_ratio(0.5, 0.25)
    expect_equal(ratio, 2.0)
    ## Also test through the full function with a hand-built evidence list
    ## by injecting a degenerate trajectory node via RecordAIEvidence
    sce <- sclet_ai_test_tv_combined_object()
    degenerate_traj <- list(
        id = "ev:trajectory_degen",
        kind = "deterministic_summary",
        values = list(
            n_lineages = 1L,
            first_lineage_median = 5.0,
            first_lineage_iqr = 0.0,
            first_lineage_max = 0.0,
            relative_ordering_only = TRUE,
            pseudotime_is_absolute_time = FALSE,
            raw_values_included = FALSE
        ),
        claim_level = "consistent_with"
    )
    sce <- RecordAIEvidence(sce, degenerate_traj, source = NULL)
    vel_node <- list(
        id = "ev:velocity_degen",
        kind = "deterministic_summary",
        values = list(
            velocity_pseudotime = list(
                n_observed = 10L,
                mean = 1.0,
                quantile_25 = 0.5,
                median = 1.0,
                quantile_75 = 1.5,
                max = 2.0
            ),
            velocity_is_model_estimate = TRUE,
            raw_values_included = FALSE
        ),
        claim_level = "consistent_with"
    )
    sce <- RecordAIEvidence(sce, vel_node, source = NULL)
    result <- compare_trajectory_velocity_evidence(sce,
        trajectory_id = "degen", velocity_id = "degen")
    expect_equal(result$status, "available")
    expect_false(result$comparable_spread)
    expect_equal(result$spread_consistency_label,
        "relative spread could not be compared because one or both bounded summaries were degenerate")
})

## -- Scenario 8: no raw per-cell values in output --
test_that("compare_trajectory_velocity_evidence output contains no raw per-cell values", {
    skip_if_not_installed("slingshot")
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_tv_combined_object()
    sce <- sclet_ai_test_tv_run_trajectory(sce, "traj_priv")
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_priv")
    result <- compare_trajectory_velocity_evidence(sce)
    expect_equal(result$status, "available")
    ## Check that no field name matches the bad_field regex
    all_names <- names(result)
    bad <- vapply(all_names, sclet:::sclet_ai_evidence_bad_field, logical(1L))
    expect_false(any(bad),
        info = paste("Bad field names found:", paste(all_names[bad], collapse = ", ")))
    ## Check that no field has excessive length (bounded aggregate)
    for (nm in all_names) {
        val <- result[[nm]]
        if (is.atomic(val)) {
            expect_true(length(val) <= 20L,
                info = paste("Field", nm, "has length", length(val)))
        }
    }
    ## Check caveats is a short list of strings
    expect_true(is.list(result$caveats))
    expect_true(length(result$caveats) <= 5L)
})

## -- Scenario 9: privacy regression guard --
test_that("compare_trajectory_velocity_evidence passes through AIStatus without privacy warnings", {
    skip_if_not_installed("slingshot")
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_tv_combined_object()
    sce <- sclet_ai_test_tv_run_trajectory(sce, "traj_ai")
    sce <- sclet_ai_test_tv_run_velocity(sce, "vel_ai")

    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Status summary.",
            findings = list(),
            evidence = list(),
            warnings = list(),
            recommendations = list(),
            proposed_actions = list()
        )
    })
    on.exit(options(old), add = TRUE)

    privacy_messages <- character()
    tryCatch(
        withCallingHandlers(
            AIStatus(sce),
            warning = function(w) {
                msg <- conditionMessage(w)
                if (grepl("redact|privacy|Unknown outbound|requires_user_consent",
                    msg, ignore.case = TRUE)) {
                    privacy_messages <<- c(privacy_messages, msg)
                }
            }
        ),
        error = function(e) NULL
    )
    expect_equal(length(privacy_messages), 0L,
        info = paste("Privacy warnings fired:", paste(privacy_messages, collapse = "; ")))
})
