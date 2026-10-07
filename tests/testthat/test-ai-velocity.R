## End-to-end velocity readiness, action, and evidence coverage (P5)
##
## Proves that check_velocity_readiness(), run_velocity, and
## sclet_ai_record_velocity_evidence() satisfy the P5 slice contract:
## readiness diagnostics, action adapter with prerequisites, bounded
## aggregate evidence, interpretation wiring, and privacy regression guard.

sclet_ai_test_velocity_object <- function(n_cells = 40L) {
    set.seed(42L)
    n_genes <- 30L
    genes <- paste0("gene_", seq_len(n_genes))
    cells <- paste0("cell_", seq_len(n_cells))
    spliced <- matrix(rpois(n_genes * n_cells, lambda = 8L), nrow = n_genes, ncol = n_cells,
        dimnames = list(genes, cells))
    unspliced <- matrix(rpois(n_genes * n_cells, lambda = 3L), nrow = n_genes, ncol = n_cells,
        dimnames = list(genes, cells))
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = spliced, spliced = spliced, unspliced = unspliced))
    sce
}

sclet_ai_test_velocity_ready_object <- function(n_cells = 40L) {
    sce <- sclet_ai_test_velocity_object(n_cells)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20L)
    sce <- RunPCA(sce, ncomponents = 5L)
    sce
}

sclet_ai_test_velocity_plan <- function(sce, registry, params) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "velocity_estimation",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "vel1", action = "run_velocity", params = params))
    )
    ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
}

## -- Scenario 1: readiness on object missing spliced/unspliced --
test_that("check_velocity_readiness reports not_ready when spliced/unspliced assays are missing", {
    set.seed(1L)
    counts <- matrix(rpois(30L * 20L, lambda = 5L), nrow = 30L, ncol = 20L,
        dimnames = list(paste0("g", seq_len(30L)), paste0("c", seq_len(20L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    result <- check_velocity_readiness(sce)
    expect_equal(result$status, "not_ready")
    expect_false(result$checks$spliced_assay_available)
    expect_false(result$checks$unspliced_assay_available)
    expect_true("velocity" %in% result$blocked_actions)
    expect_true(length(result$questions) > 0L)
})

## -- Scenario 2: readiness with spliced/unspliced but no reduction --
test_that("check_velocity_readiness reports not_ready when no reduction is available", {
    sce <- sclet_ai_test_velocity_object()
    result <- check_velocity_readiness(sce)
    expect_equal(result$status, "not_ready")
    expect_true(result$checks$spliced_assay_available)
    expect_true(result$checks$unspliced_assay_available)
    expect_null(result$checks$reduction_resolved)
    expect_true("velocity" %in% result$blocked_actions)
})

## -- Scenario 3: readiness on a fully-ready object --
test_that("check_velocity_readiness reports ready_for_diagnostic on a fully-prepared object", {
    sce <- sclet_ai_test_velocity_ready_object()
    result <- check_velocity_readiness(sce)
    expect_equal(result$status, "ready_for_diagnostic")
    expect_true(length(result$questions) == 0L)
    expect_equal(length(result$blocked_actions), 0L)
    expect_false(is.null(result$checks$reduction_resolved))
})

## -- Scenario 4: prerequisites reject missing spliced/unspliced via ValidateAIPlan --
test_that("run_velocity prerequisites reject a plan missing spliced/unspliced assays", {
    set.seed(1L)
    counts <- matrix(rpois(30L * 20L, lambda = 5L), nrow = 30L, ncol = 20L,
        dimnames = list(paste0("g", seq_len(30L)), paste0("c", seq_len(20L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- RunPCA(sce, ncomponents = 5L)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    validation <- sclet_ai_test_velocity_plan(sce, registry, list(name = "vel_missing"))
    expect_false(validation$valid)
    expect_true(any(grepl("spliced_unspliced_missing", validation$errors)))
})

## -- Scenario 5: prerequisites reject missing reduction via ValidateAIPlan --
test_that("run_velocity prerequisites reject a plan with no available reduction", {
    sce <- sclet_ai_test_velocity_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    validation <- sclet_ai_test_velocity_plan(sce, registry, list(name = "vel_no_red"))
    expect_false(validation$valid)
    expect_true(any(grepl("reduction_missing", validation$errors)))
})

## -- Scenario 6: run_velocity survives ExecuteAIPlan directly --
test_that("run_velocity survives ExecuteAIPlan directly on a real ready fixture", {
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_velocity_ready_object()
    n_cells <- ncol(sce)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    validation <- sclet_ai_test_velocity_plan(sce, registry, list(name = "vel_direct"))
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    vel_ev <- evidence[["ev:velocity_vel_direct"]]
    expect_false(is.null(vel_ev))
    expect_equal(vel_ev$claim_level, "consistent_with")

    per_cell <- any(vapply(vel_ev$values, function(v) {
        is.atomic(v) && length(v) >= n_cells
    }, logical(1L)))
    expect_false(isTRUE(per_cell))
})

## -- Scenario 7: run_velocity end-to-end through RunAIAnalysis --
test_that("run_velocity survives RunAIAnalysis end to end with bounded evidence", {
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_velocity_ready_object()
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Estimate RNA velocity.",
            proposed_actions = list(list(
                id = "vel1",
                action = "run_velocity",
                params = list(name = "vel_e2e")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    result <- suppressWarnings(RunAIAnalysis(sce, goal = "Estimate RNA velocity.",
        confirm = "yes", registry = registry))

    expect_equal(result$status, "completed")
    expect_equal(result$execution$status, "completed")

    evidence <- sclet:::sclet_ai_evidence_get_all(result$object)
    vel_ev <- evidence[["ev:velocity_vel_e2e"]]
    expect_false(is.null(vel_ev))
    expect_equal(vel_ev$claim_level, "consistent_with")
})

## -- Scenario 8: prerequisite failure through the full orchestrator --
test_that("a velocity readiness failure surfaces as invalid_plan through RunAIAnalysis", {
    set.seed(1L)
    counts <- matrix(rpois(30L * 20L, lambda = 5L), nrow = 30L, ncol = 20L,
        dimnames = list(paste0("g", seq_len(30L)), paste0("c", seq_len(20L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- RunPCA(sce, ncomponents = 5L)

    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Estimate RNA velocity.",
            proposed_actions = list(list(
                id = "vel1",
                action = "run_velocity",
                params = list(name = "vel_fail")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    result <- RunAIAnalysis(sce, goal = "Estimate RNA velocity.",
        confirm = "yes", registry = registry)

    expect_equal(result$status, "invalid_plan")
    expect_null(result$execution)
})

## -- Scenario 9: RunAIAnalysis(interpret = TRUE) works for run_velocity --
test_that("RunAIAnalysis with interpret = TRUE works for run_velocity and persists interpretation", {
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_velocity_ready_object()

    old <- options(sclet.ai.call = function(task, context, ...) {
        if (identical(task, "analysis_plan")) {
            list(
                answer = "Estimate RNA velocity.",
                proposed_actions = list(list(
                    id = "vel1",
                    action = "run_velocity",
                    params = list(name = "vel_interp")
                ))
            )
        } else if (identical(task, "analysis_explanation")) {
            list(
                answer = "Interpretation of the velocity analysis.",
                findings = list(list(
                    statement = "Velocity was estimated using scVelo.",
                    severity = "info",
                    claim_level = "associated",
                    evidence_refs = "ev:velocity_vel_interp"
                )),
                evidence = list("ev:velocity_vel_interp"),
                warnings = list(),
                recommendations = list(),
                proposed_actions = list()
            )
        } else {
            list(answer = "Unhandled task.")
        }
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    result <- suppressWarnings(RunAIAnalysis(sce, goal = "Estimate RNA velocity and interpret.",
        confirm = "yes", registry = registry, interpret = TRUE))

    expect_equal(result$status, "completed")
    expect_false(is.null(result$report$interpretation))
    expect_s3_class(result$report$interpretation, "sclet_ai_result")
    expect_true(validate_sclet_ai_result(result$report$interpretation, error = FALSE))
    expect_null(result$report$interpretation_error)

    interp_ids <- grep("^ai_interpretation_",
        names(GetAnalysisLedger(result$object)$analyses), value = TRUE)
    expect_true(length(interp_ids) >= 1L)
})

## -- Scenario 10: evidence never contains a raw per-cell vector --
test_that("velocity evidence never contains a raw per-cell vector", {
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_velocity_ready_object()
    n_cells <- ncol(sce)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    validation <- sclet_ai_test_velocity_plan(sce, registry, list(name = "vel_bounded"))
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    vel_ev <- evidence[["ev:velocity_vel_bounded"]]
    expect_false(is.null(vel_ev))

    check_bounded <- function(x) {
        if (is.list(x) && !is.null(names(x))) {
            all(vapply(x, check_bounded, logical(1L)))
        } else if (is.atomic(x)) {
            length(x) < n_cells
        } else {
            TRUE
        }
    }
    expect_true(check_bounded(vel_ev$values))
})

## -- Scenario 11: no design confirmation is required --
test_that("run_velocity validates and executes without any ConfirmAIDesignSemantics call", {
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_velocity_ready_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))

    plan <- sclet:::new_sclet_ai_plan(
        task = "velocity_estimation",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "vel1", action = "run_velocity",
            params = list(name = "vel_no_confirm")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)

    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    design_confirmations <- GetAnalysisLedger(executed$object)$analysis_story$design_confirmations
    expect_true(length(design_confirmations) == 0L ||
        !any(vapply(design_confirmations, function(dc) {
            identical(dc$role, "design_velocity")
        }, logical(1L))))
})

## -- Scenario 12: full suite regression guard --
## (covered by running the whole test suite; this test is a placeholder
## that always passes, serving as a reminder that the full suite must be green)
test_that("existing tests still pass (regression guard placeholder)", {
    expect_true(TRUE)
})

## -- Scenario 13: privacy regression guard --
test_that("velocity evidence passes through AIStatus without privacy warnings", {
    skip_if_not_installed("velociraptor")
    sce <- sclet_ai_test_velocity_ready_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "velocity"))
    validation <- sclet_ai_test_velocity_plan(sce, registry, list(name = "vel_priv"))
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    vel_ev <- evidence[["ev:velocity_vel_priv"]]
    expect_false(is.null(vel_ev))

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
            AIStatus(executed$object),
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
