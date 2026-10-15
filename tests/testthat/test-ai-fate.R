## Fate (CellRank) contract: readiness, action prerequisites, execution with a
## mocked Python backend, and bounded evidence (P2 dynamic-process slice).
##
## No real Python / basilisk environment is ever built here: the CellRank
## Python bridge `sclet_run_cellrank_backend()` is mocked, exactly the way
## test-scenic.R mocks `sclet_basilisk_run()`. Everything asserted below is
## produced by sclet's own R code (colData columns, trajectory state record,
## evidence node), not by the mock.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

## Minimal object: counts + spliced + unspliced, used for the "no velocity
## record" readiness / prerequisite cases.
sclet_ai_test_fate_object <- function(n_cells = 24L, n_genes = 20L) {
    set.seed(20251014L)
    genes <- paste0("gene_", seq_len(n_genes))
    cells <- paste0("cell_", seq_len(n_cells))
    spliced <- matrix(rpois(n_genes * n_cells, lambda = 8L),
        nrow = n_genes, ncol = n_cells, dimnames = list(genes, cells))
    unspliced <- matrix(rpois(n_genes * n_cells, lambda = 3L),
        nrow = n_genes, ncol = n_cells, dimnames = list(genes, cells))
    velocity <- matrix(rnorm(n_genes * n_cells),
        nrow = n_genes, ncol = n_cells, dimnames = list(genes, cells))
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = spliced, spliced = spliced, unspliced = unspliced, velocity = velocity)
    )
}

## Same object plus the pieces fate inference needs: a velocity state record
## (built without running velociraptor), a "PCA" reducedDim, and a cluster
## colData column wired up as the active identity.
sclet_ai_test_fate_ready_object <- function(n_cells = 24L, velocity_id = "vel_fixture",
                                            with_velocity_assay = TRUE, n_clusters = 3L) {
    sce <- sclet_ai_test_fate_object(n_cells)
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(
        stats::rnorm(n_cells * 5L), nrow = n_cells, ncol = 5L,
        dimnames = list(colnames(sce), paste0("PC", seq_len(5L)))
    )
    sce$cluster <- factor(paste0("cluster_", seq_len(n_clusters))[seq_len(n_cells) %% n_clusters + 1L])
    ActiveIdent(sce) <- "cluster"

    artifacts <- list(
        analysis_key = velocity_id,
        colData = character(),
        reductions = character()
    )
    if (isTRUE(with_velocity_assay)) {
        artifacts$velocity_assay <- "velocity"
    }
    ## Mirrors the record RunVelocity() writes (R/velocity.R), minus the
    ## velociraptor call: a velocity state record is the only thing fate
    ## readiness needs on top of the assays.
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "velocity",
        id = velocity_id,
        method = "velociraptor::scvelo",
        inputs = list(assays = c("spliced", "unspliced"), reduction = "PCA"),
        artifacts = artifacts,
        params = list(mode = "deterministic", subset.row = NULL),
        summary = list(n_features = nrow(sce), n_cells = ncol(sce)),
        active = TRUE
    )
    sce
}

## What the Python side would return: per-cell terminal states (character, one
## per cell), fate probabilities (n_cells x 2), and no driver table.
sclet_ai_test_fate_backend_mock <- function(n_cells, n_lineages = 2L, seed = 7L) {
    set.seed(seed)
    terminal_states <- rep(c("cluster_1", "cluster_2"), length.out = n_cells)
    probs <- matrix(stats::runif(n_cells * n_lineages), nrow = n_cells, ncol = n_lineages)
    probs <- probs / rowSums(probs)
    function(backend, python, mtx_dir) {
        list(
            terminal_states = terminal_states,
            absorption_probs = probs,
            absorption_prob_names = c("A", "B")[seq_len(n_lineages)],
            lineage_drivers = NULL
        )
    }
}

## Validate a one-step fate plan through ValidateAIPlan, the way
## test-ai-velocity.R validates run_velocity prerequisites. `cluster_key` is
## always passed explicitly rather than relying on the ActiveIdent fallback.
sclet_ai_test_fate_plan <- function(sce, registry, params) {
    params <- utils::modifyList(list(cluster_key = "cluster"), as.list(params))
    plan <- sclet:::new_sclet_ai_plan(
        task = "fate_analysis",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "fate1", action = "run_fate_analysis", params = params))
    )
    ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
}

## Guard: the Python CellRank bridge must never actually run. If the mock is
## ever bypassed, this turns a silent env build into a loud test failure.
sclet_ai_test_fate_no_python <- function(mtx_dir) {
    stop("Real Python CellRank bridge reached: sclet_run_cellrank_python must stay mocked.")
}

## -- Scenario 1: readiness is not_ready without a velocity record ------------
test_that("check_fate_readiness reports not_ready when no velocity record exists", {
    sce <- sclet_ai_test_fate_object()
    result <- check_fate_readiness(sce)

    expect_equal(result$status, "not_ready")
    expect_true("fate" %in% result$blocked_actions)
    expect_true(length(result$questions) > 0L)
    expect_false(isTRUE(result$checks$velocity_record_available))
})

## -- Scenario 2: readiness turns ready once velocity + PCA + clusters exist ---
test_that("check_fate_readiness reports ready_for_diagnostic with a velocity record, PCA and clusters", {
    sce <- sclet_ai_test_fate_ready_object()
    result <- check_fate_readiness(sce)

    expect_equal(result$status, "ready_for_diagnostic")
    expect_equal(length(result$blocked_actions), 0L)
    expect_equal(length(result$questions), 0L)
    expect_true(isTRUE(result$checks$velocity_record_available))
    expect_equal(result$checks$velocity_id_resolved, "vel_fixture")
    expect_true(result$checks$cluster_key_resolvable)
    expect_false(is.null(result$checks$reduction_resolved))
})

## -- Scenario 3: prerequisites reject a missing velocity record --------------
test_that("run_fate_analysis prerequisites reject a missing velocity record", {
    sce <- sclet_ai_test_fate_object()
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(
        stats::rnorm(ncol(sce) * 5L), nrow = ncol(sce), ncol = 5L,
        dimnames = list(colnames(sce), paste0("PC", seq_len(5L)))
    )
    sce$cluster <- factor(rep(c("cluster_1", "cluster_2"), length.out = ncol(sce)))
    ActiveIdent(sce) <- "cluster"

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
    validation <- sclet_ai_test_fate_plan(sce, registry, list(name = "fate_no_velocity"))

    expect_false(validation$valid)
    ## ValidateAIPlan wraps the prerequisite verdict as
    ## "prerequisites not met for action <id>: <verdict>", so the token is
    ## matched anywhere in the error, not anchored at the start (same
    ## convention as test-ai-velocity.R's "spliced_unspliced_missing").
    expect_true(any(grepl("velocity_record_missing:", validation$errors)))
})

## -- Scenario 4: the mocked backend drives a full fate run -------------------
test_that("run_fate_analysis executes end to end on a mocked CellRank backend", {
    sce <- sclet_ai_test_fate_ready_object()
    n_cells <- ncol(sce)

    backend_calls <- list()
    mock <- function(backend, python, mtx_dir) {
        backend_calls <<- c(backend_calls, list(list(backend = backend, mtx_dir = mtx_dir)))
        expected <- c("velocity.mtx", "spliced.mtx", "unspliced.mtx",
            "clusters.txt", "cell_names.txt", "gene_names.txt", "reduction.txt", "emb.csv")
        expect_true(all(expected %in% list.files(mtx_dir)))
        sclet_ai_test_fate_backend_mock(n_cells)(backend, python, mtx_dir)
    }

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
    validation <- sclet_ai_test_fate_plan(sce, registry, list(name = "fate_e2e"))
    expect_true(validation$valid)

    executed <- with_mocked_bindings(
        ExecuteAIPlan(
            sce, validation$plan, registry,
            validation = validation,
            dry_run = FALSE,
            confirmation = validation$confirmation_token
        ),
        sclet_run_cellrank_backend = mock,
        sclet_run_cellrank_python = sclet_ai_test_fate_no_python,
        .package = "sclet"
    )
    expect_equal(executed$status, "completed")
    ## The Python bridge really was reached -- and only through the mock.
    expect_equal(length(backend_calls), 1L)

    out <- executed$object

    ## A trajectory state record exists and is attributed to CellRank.
    trajectory <- sclet:::sclet_get_state_record(out, "trajectory", id = "fate_e2e")
    expect_false(is.null(trajectory))
    expect_equal(trajectory$method, "CellRank")
    expect_equal(trajectory$type, "trajectory")

    ## No cell was dropped or reordered by the fate action.
    expect_equal(ncol(out), n_cells)
    expect_equal(colnames(out), colnames(sce))

    ## RunCellRank's own colData contract survived the mocked backend.
    cd <- SummarizedExperiment::colData(out)
    expect_true("cellrank_terminal_states" %in% colnames(cd))
    expect_equal(length(cd[["cellrank_terminal_states"]]), n_cells)
    expect_true(all(c("cellrank_fate_A", "cellrank_fate_B") %in% colnames(cd)))
    expect_equal(sum(as.numeric(cd[["cellrank_fate_A"]]) + as.numeric(cd[["cellrank_fate_B"]])),
        n_cells, tolerance = 1e-6)

    ## And the bounded evidence node was registered.
    evidence <- sclet:::sclet_ai_evidence_get_all(out)
    fate_ev <- evidence[["ev:fate_fate_e2e"]]
    expect_false(is.null(fate_ev))
})

## -- Scenario 5: fate evidence stays bounded and privacy-clean ---------------
test_that("fate evidence is bounded, anonymous and carries no per-cell vector", {
    sce <- sclet_ai_test_fate_ready_object()
    n_cells <- ncol(sce)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
    validation <- sclet_ai_test_fate_plan(sce, registry, list(name = "fate_bounded"))
    expect_true(validation$valid)

    executed <- with_mocked_bindings(
        ExecuteAIPlan(
            sce, validation$plan, registry,
            validation = validation,
            dry_run = FALSE,
            confirmation = validation$confirmation_token
        ),
        sclet_run_cellrank_backend = sclet_ai_test_fate_backend_mock(n_cells),
        sclet_run_cellrank_python = sclet_ai_test_fate_no_python,
        .package = "sclet"
    )

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    fate_ev <- evidence[["ev:fate_fate_bounded"]]
    expect_false(is.null(fate_ev))

    values <- fate_ev$values
    expect_true(length(values) > 0L)

    ## No field name leaks identifiers, paths or per-cell keys.
    bad <- vapply(names(values), sclet:::sclet_ai_evidence_bad_field, logical(1L))
    expect_false(any(bad))

    ## Terminal states enter evidence only as anonymized cluster codes, never
    ## as the raw label the backend returned.
    expect_true("terminal_state_codes" %in% names(values))
    expect_true(length(values$terminal_state_codes) > 0L)
    expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$",
        values$terminal_state_codes)))

    ## Every character value is an anonymized group/cluster code.
    char_values <- unlist(lapply(values, function(v) {
        if (is.character(v)) v else NULL
    }))
    if (length(char_values)) {
        expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$", char_values)))
    }

    expect_false(isTRUE(values$raw_values_included))
    ## Fate probabilities are flagged as model estimates, not measured fate.
    expect_true(isTRUE(values$absorption_probability_is_model_estimate))
    expect_false(isTRUE(values$pseudotime_is_absolute_time))

    ## No per-cell vector -- nothing with one entry per cell.
    has_per_cell <- function(x) {
        if (is.list(x)) {
            any(vapply(x, has_per_cell, logical(1L)))
        } else {
            is.atomic(x) && length(x) >= n_cells
        }
    }
    expect_false(has_per_cell(values))
    ## The structure passes the shared evidence value guard as a whole.
    expect_true(sclet:::sclet_ai_evidence_value_ok(fate_ev$values))
})

## -- Scenario 6: readiness never nominates a terminal state ------------------
test_that("check_fate_readiness never suggests or ranks a terminal state", {
    ready <- check_fate_readiness(sclet_ai_test_fate_ready_object())
    blocked <- check_fate_readiness(sclet_ai_test_fate_object())

    for (result in list(ready, blocked)) {
        ## Nothing that reads like a suggestion or a ranking, at any depth.
        flat_names <- names(result)
        stack <- list(result)
        while (length(stack)) {
            node <- stack[[1L]]
            stack <- stack[-1L]
            if (is.list(node) && !is.null(names(node))) {
                flat_names <- c(flat_names, names(node))
                stack <- c(stack, node)
            }
        }
        flat_names <- unique(flat_names[!is.na(flat_names)])
        expect_equal(length(grep("suggest|recommend|rank|terminal_state_.*(best|top)",
            flat_names, ignore.case = TRUE)), 0L)
    }
})
