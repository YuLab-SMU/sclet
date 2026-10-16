## Cell-cell interaction contract (batch B of the single-cell re-scope):
## readiness, the run_cci action, bounded group-flow evidence, and the shared
## dependency group that stops two engines from counting as independent lines.
##
## Why the engine is mocked and everything else is real. A real CellChat run on
## a test-scale fixture is not usable: with random counts it reports "The number
## of highly variable ligand-receptor pairs used for signaling inference is 0"
## and then fails inside CellChat with `subscript out of bounds`, so a test that
## drove the real engine would be asserting against CellChat's private database
## and thresholds rather than against this contract. Only RunCCI is therefore
## stubbed. The interaction table it writes is built in exactly the shape
## RunCellChat()/RunCellPhoneDB() persist -- analysis-layer `interactions` plus
## a type="communication" state record -- so get_cci(), the recorder, the
## dependency grouping, check_claim_ceiling() and the read-only summary are all
## real code paths over that fixture.
##
## The leak test is the heart of this file because in this layer the payload IS
## gene identity: ligand and receptor columns hold symbols, and CellPhoneDB's
## pathway column holds "LIGAND_RECEPTOR" pair strings. Population labels are
## deliberately non-code-looking here ("beta cells" / "alpha cells") so the
## anonymization is actually exercised; test-ai-fate.R keeps the same rule of
## never building fixtures whose labels already satisfy the code pattern.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

sclet_ai_test_cci_genes <- function() {
    c("CXCL12", "ACKR1", "HLA-A", "CD8A", "APP", "CD74", "ACTB", "GAPDH")
}

sclet_ai_test_cci_object <- function(n_cells = 40L) {
    genes <- sclet_ai_test_cci_genes()
    counts <- matrix(rpois(length(genes) * n_cells, lambda = 8L),
        nrow = length(genes), ncol = n_cells,
        dimnames = list(genes, paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    ## Labels that are NOT already anonymization-shaped, so a leak of a raw
    ## label and a legitimate code are distinguishable in the assertions below.
    sce$cluster <- factor(rep(c("beta cells", "alpha cells"), length.out = n_cells),
        levels = c("beta cells", "alpha cells"))
    ActiveIdent(sce) <- "cluster"
    sce
}

## The standardized interaction table the CCI runners persist. Columns and
## their order mirror R/cci.R (source_group, target_group, ligand, receptor,
## pathway, probability, p_value, method).
sclet_ai_test_cci_table <- function(method = "CellChat") {
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

## Persist a communication result the way the real runners do: the interaction
## table into the analysis layer (where get_cci() reads it from) and a
## type="communication" state record that carries only counts.
sclet_ai_test_cci_store <- function(object, id = "cci1", method = "CellChat",
                                   analysis_key = "cellchat", table = sclet_ai_test_cci_table(method),
                                   group = "cluster") {
    object <- sclet_set_analysis(object, analysis_key, list(
        id = id,
        method = method,
        object = NULL,
        interactions = table,
        group = group
    ))
    sclet_set_analysis_state(
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

## Stub for RunCCI(). It writes through the real storage path so everything
## downstream is exercised; `writes = FALSE` reproduces the silent no-record
## case the action has to detect.
sclet_ai_test_cci_engine_stub <- function(id = "cci1", method = "CellChat",
                                          analysis_key = "cellchat", writes = TRUE,
                                          table = sclet_ai_test_cci_table(method)) {
    function(object, method = NULL, ...) {
        if (!writes) return(object)
        sclet_ai_test_cci_store(object, id = id, method = method,
            analysis_key = analysis_key, table = table)
    }
}

## Pin the CCI engine gate for the duration of the calling test. The gate asks
## the real library, so an unpinned test asserts whatever the machine happens to
## have installed: CI installs nichenetr through its own setup step, while many
## developer machines do not have it, and CellChat is absent from some. The gate
## is what is under test here, not this environment, so both branches are forced
## deterministically: CellChat and basilisk present, nichenetr absent. Every
## other package keeps its real answer, and the pin is undone when the test ends.
sclet_ai_test_cci_pin_engines <- function(env = parent.frame()) {
    real_req_ns <- base::requireNamespace
    testthat::local_mocked_bindings(
        requireNamespace = function(package, ...) {
            if (identical(package, "CellChat")) return(TRUE)
            if (identical(package, "basilisk")) return(TRUE)
            if (identical(package, "nichenetr")) return(FALSE)
            real_req_ns(package, ...)
        },
        .package = "base",
        .env = env
    )
}

sclet_ai_test_cci_registry <- function(object) {
    AIDefaultExecutionRegistry(object, include = c("read", "communication"))
}

sclet_ai_test_cci_plan <- function(object, registry, params) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "communication",
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        actions = list(list(id = "cc1", action = "run_cci", params = as.list(params)))
    )
    ValidateAIPlan(plan, object = object, registry = registry, strict = FALSE)
}

## Every character value anywhere in a nested object.
sclet_ai_test_cci_characters <- function(x) {
    if (is.character(x)) return(x)
    if (is.list(x)) return(unlist(lapply(x, sclet_ai_test_cci_characters), use.names = FALSE))
    character()
}

## Every piece of text an object can carry, as one string. capture.output(str())
## is NOT used for the leak scan: str() elides long character vectors and
## truncates wide data frames, so a value that really leaked could sit in an
## omitted section and the scan would pass vacuously. Names are included
## because a field called `ligand` is itself a leak signal.
sclet_ai_test_cci_all_text <- function(x) {
    pieces <- character()
    collect <- function(value, name) {
        if (!is.null(name) && length(name)) pieces[[length(pieces) + 1L]] <<- as.character(name)
        if (is.character(value)) {
            pieces[[length(pieces) + 1L]] <<- paste(value, collapse = "\r")
        } else if (is.list(value)) {
            if (!is.null(names(value))) {
                pieces[[length(pieces) + 1L]] <<- paste(names(value), collapse = "\r")
            }
            for (i in seq_along(value)) collect(value[[i]], names(value)[i])
        } else if (is.data.frame(value)) {
            for (column in colnames(value)) collect(value[[column]], column)
        }
    }
    collect(x, NULL)
    paste(unlist(pieces, use.names = FALSE), collapse = "\n")
}

## The identity strings this layer must never let escape.
sclet_ai_test_cci_identity_strings <- function(table) {
    unique(c(
        "CXCL12", "ACKR1", "HLA-A", "CD8A", "APP", "CD74",
        table$ligand, table$receptor, table$pathway,
        "beta cells", "alpha cells"
    ))
}

## -- Block 1: readiness is blocked on an empty object ---------------------
test_that("check_communication_readiness is blocked on a bare object", {
    bare <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 3L, 4L)))
    result <- check_communication_readiness(bare)
    expect_equal(result$status, "not_ready")
    expect_true("communication" %in% result$blocked_actions)
    expect_false(result$checks$group_column_resolvable)
    expect_true(length(result$questions) > 0L)
    expect_false(isTRUE(result$checks$raw_values_included))
})

## -- Block 2: readiness passes and refuses to pick an engine --------------
test_that("readiness reports the engine count without recommending one", {
    sclet_ai_test_cci_pin_engines()
    object <- sclet_ai_test_cci_object()
    result <- check_communication_readiness(object)
    expect_equal(result$status, "ready_for_diagnostic")
    expect_equal(length(result$blocked_actions), 0L)
    expect_true(result$checks$group_column_resolvable)
    expect_equal(as.character(result$checks$group_column_resolved), "cluster")
    expect_true(result$checks$expression_source_available)
    expect_gt(length(result$checks$available_methods), 1L)

    ## Several engines available, none named: the user must choose, and the
    ## question carries only the COUNT, not an order or a preference.
    choice <- grep("engines are available", result$questions, value = TRUE)
    expect_length(choice, 1L)
    expect_true(grepl("must choose", choice))
    expect_match(choice, as.character(length(result$checks$available_methods)), fixed = TRUE)
    expect_false(any(grepl("recommend|prefer|best|use CellChat|suggest", result$questions)))
    ## Group labels never appear in the readiness output.
    dumped <- sclet_ai_test_cci_all_text(result)
    expect_false(grepl("beta cells", dumped, fixed = TRUE))
    expect_false(grepl("alpha cells", dumped, fixed = TRUE))
})

## -- Block 3: a named engine's package gate is per-engine -----------------
## The engine gate is pinned so both branches are exercised for real rather than
## skipped, and so the assertion holds on a machine where CellChat is absent.
## The point is that a missing engine blocks WITHOUT steering the user to a
## different engine.
test_that("readiness blocks a named engine whose package is missing", {
    sclet_ai_test_cci_pin_engines()
    object <- sclet_ai_test_cci_object()

    ready <- check_communication_readiness(object, method = "CellChat")
    expect_equal(ready$status, "ready_for_diagnostic")
    expect_true(isTRUE(ready$checks$method_packages_available[["CellChat"]]))

    blocked <- check_communication_readiness(object, method = "NicheNet")
    expect_equal(blocked$status, "not_ready")
    expect_true("communication" %in% blocked$blocked_actions)
    expect_false(isFALSE(blocked$checks$group_column_resolvable))
    expect_true(any(grepl("nichenetr", blocked$questions)))
    ## Reporting the gap is not recommending a substitute.
    expect_false(any(grepl("instead|prefer|use CellChat|switch", blocked$questions)))
})

## -- Block 4: the recorder leaks nothing ----------------------------------
test_that("recorded communication evidence carries no ligand, receptor, pathway or label", {
    object <- sclet_ai_test_cci_store(sclet_ai_test_cci_object())
    recorded <- sclet:::sclet_ai_record_communication_evidence(object, name = "cci1",
        method = "CellChat")
    node <- sclet:::sclet_ai_evidence_get_all(recorded)[["ev:communication_cci1"]]
    expect_false(is.null(node))
    expect_equal(node$claim_level, "consistent_with")

    node_dump <- sclet_ai_test_cci_all_text(node)
    table <- sclet_ai_test_cci_table()
    for (token in sclet_ai_test_cci_identity_strings(table)) {
        expect_false(grepl(token, node_dump, fixed = TRUE))
    }
    ## POSITIVE CONTROL: the very same scan DOES find those tokens in the
    ## analysis-layer table they were withheld from. Without this, "absent from
    ## the node" could pass because the patterns never match anything anywhere.
    table_dump <- sclet_ai_test_cci_all_text(table)
    expect_true(all(vapply(sclet_ai_test_cci_identity_strings(table), function(token) {
        grepl(token, table_dump, fixed = TRUE)
    }, logical(1L))))
    ## Character values are codes only...
    codes <- sclet_ai_test_cci_characters(node$values)
    expect_gt(length(codes), 0L)
    expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$", codes)))
    ## ...and the populations really were anonymized, not merely copied.
    expect_true(all(node$values$source_group_codes == "cluster_1" |
        node$values$source_group_codes == "cluster_2"))
    expect_false(isTRUE(node$values$raw_values_included))
    expect_true(isTRUE(node$values$communication_is_model_estimate))
    expect_equal(node$values$n_interactions, 5L)
    expect_equal(node$values$n_significant_pairs, 2L)
    ## Row-major flow grid: one count per (source, target) pair.
    expect_equal(length(node$values$pair_flow_counts_rowmajor),
        node$values$n_source_groups * node$values$n_target_groups)
    expect_equal(sum(node$values$pair_flow_counts_rowmajor), node$values$n_interactions)

    ## The detector has teeth: the same scan fires on a tampered node.
    tampered <- node
    tampered$values$decoy <- "CXCL12"
    expect_true(grepl("CXCL12", sclet_ai_test_cci_all_text(tampered), fixed = TRUE))
})

## -- Block 5: two engines over the same inputs are one line ---------------
test_that("two communication engines on shared inputs are not independent evidence", {
    object <- sclet_ai_test_cci_store(sclet_ai_test_cci_object(), id = "cc_cellchat",
        method = "CellChat", analysis_key = "cellchat")
    object <- sclet:::sclet_ai_record_communication_evidence(object, name = "cc_cellchat",
        method = "CellChat")
    object <- sclet_ai_test_cci_store(object, id = "cc_cpdb", method = "CellPhoneDB",
        analysis_key = "cellphonedb", table = sclet_ai_test_cci_table("CellPhoneDB"))
    object <- sclet:::sclet_ai_record_communication_evidence(object, name = "cc_cpdb",
        method = "CellPhoneDB")

    both <- check_claim_ceiling(object,
        evidence_ids = c("ev:communication_cc_cellchat", "ev:communication_cc_cpdb"))
    expect_equal(both$n_independent_groups, 1L)

    ## CONTROL: an unrelated second source really is a second line. Without
    ## this, the assertion above could pass because grouping never splits.
    unrelated <- RecordAIEvidence(object, list(
        id = "ev:unrelated_control", kind = "deterministic_summary",
        values = list(n_observed = 4L, raw_values_included = FALSE),
        claim_level = "associated"), source = NULL, parents = character(), scope = NULL)
    control <- check_claim_ceiling(unrelated,
        evidence_ids = c("ev:communication_cc_cellchat", "ev:unrelated_control"))
    expect_equal(control$n_independent_groups, 2L)
})

## -- Block 6: the action runs end to end ----------------------------------
test_that("run_cci executes with a stubbed engine and records bounded evidence", {
    sclet_ai_test_cci_pin_engines()
    object <- sclet_ai_test_cci_object()
    registry <- sclet_ai_test_cci_registry(object)
    validation <- sclet_ai_test_cci_plan(object, registry, list(method = "CellChat"))
    expect_true(validation$valid)

    executed <- with_mocked_bindings(
        ExecuteAIPlan(object, validation$plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        RunCCI = sclet_ai_test_cci_engine_stub(id = "cci_action", method = "CellChat"),
        .package = "sclet"
    )
    expect_equal(executed$status, "completed")
    expect_true(has_cci(executed$object, "cci_action"))
    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true("ev:communication_cci_action" %in% names(nodes))
    expect_false(isTRUE(nodes[["ev:communication_cci_action"]]$values$raw_values_included))
    expect_true(isTRUE(nodes[["ev:communication_cci_action"]]$values$engine_was_named_by_caller))

    ## The real accessor reads the fixture table the stub stored.
    table <- get_cci(executed$object, id = "cci_action")
    expect_true(all(c("source_group", "target_group", "ligand", "receptor") %in% colnames(table)))
})

## -- Block 7: a run that records nothing must not report success ----------
test_that("a silent engine run surfaces communication_no_result", {
    sclet_ai_test_cci_pin_engines()
    object <- sclet_ai_test_cci_object()
    registry <- sclet_ai_test_cci_registry(object)
    validation <- sclet_ai_test_cci_plan(object, registry, list(method = "CellChat"))
    expect_true(validation$valid)

    executed <- with_mocked_bindings(
        ExecuteAIPlan(object, validation$plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        RunCCI = sclet_ai_test_cci_engine_stub(writes = FALSE),
        .package = "sclet"
    )
    expect_equal(executed$status, "failed")
    messages <- unlist(lapply(executed$results, function(step) step$error))
    expect_true(any(grepl("communication_no_result", messages)))
    expect_false("ev:communication_cci1" %in%
        names(sclet:::sclet_ai_evidence_get_all(executed$object)))
})

## -- Block 8: typed prerequisites -----------------------------------------
test_that("run_cci refuses to pick an engine and reports typed reasons", {
    sclet_ai_test_cci_pin_engines()
    object <- sclet_ai_test_cci_object()
    registry <- sclet_ai_test_cci_registry(object)
    ## ValidateAIPlan prefixes errors with "prerequisites not met for action
    ## cc1: ", so tokens are matched anywhere and never with a ^ anchor.
    unnamed <- sclet_ai_test_cci_plan(object, registry, list())
    expect_false(unnamed$valid)
    expect_true(any(grepl("method_missing", unnamed$errors)))
    ## The refusal explains itself and offers no default.
    expect_true(any(grepl("No engine is chosen", unnamed$errors)))

    invalid <- sclet_ai_test_cci_plan(object, registry, list(method = "StringApper"))
    expect_false(invalid$valid)
    expect_true(any(grepl("method_invalid", invalid$errors)))

    bare <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 3L, 4L)))
    no_group <- sclet_ai_test_cci_plan(bare, sclet_ai_test_cci_registry(bare),
        list(method = "CellChat"))
    expect_false(no_group$valid)
    expect_true(any(grepl("group_column_missing", no_group$errors)))

    no_package <- sclet_ai_test_cci_plan(object, registry, list(method = "NicheNet"))
    expect_false(no_package$valid)
    expect_true(any(grepl("communication_package_missing", no_package$errors)))
})

## -- Block 9: the read-only summary ---------------------------------------
test_that("summarize_communication_evidence is bounded and recommends nothing", {
    bare <- sclet_ai_test_cci_object()
    absent <- summarize_communication_evidence(bare)
    expect_equal(absent$status, "not_available")
    expect_equal(absent$reason, "no_communication_record")
    expect_null(absent$recommendation)

    object <- sclet_ai_test_cci_store(bare)
    object <- sclet:::sclet_ai_record_communication_evidence(object, name = "cci1",
        method = "CellChat")
    summary <- summarize_communication_evidence(object, id = "cci1")
    expect_equal(summary$status, "available")
    expect_null(summary$recommendation)
    expect_false(isTRUE(summary$raw_values_included))
    expect_true(summary$evidence_node_recorded)
    expect_equal(summary$n_interactions, 5L)
    expect_equal(summary$method_code, 1L)

    dumped <- sclet_ai_test_cci_all_text(summary)
    for (token in sclet_ai_test_cci_identity_strings(sclet_ai_test_cci_table())) {
        expect_false(grepl(token, dumped, fixed = TRUE))
    }
})

## -- Block 10: everything here is read-only -------------------------------
test_that("the communication diagnostics never mutate the object", {
    object <- sclet_ai_test_cci_store(sclet_ai_test_cci_object())
    object <- sclet:::sclet_ai_record_communication_evidence(object, name = "cci1",
        method = "CellChat")
    fingerprint <- GetAnalysisLedger(object)$fingerprint
    evidence <- paste(sort(names(sclet:::sclet_ai_evidence_get_all(object))), collapse = ",")

    check_communication_readiness(object)
    check_communication_readiness(object, method = "CellChat")
    summarize_communication_evidence(object, id = "cci1")
    check_claim_ceiling(object)

    expect_equal(GetAnalysisLedger(object)$fingerprint, fingerprint)
    expect_equal(paste(sort(names(sclet:::sclet_ai_evidence_get_all(object))), collapse = ","),
        evidence)
})

## -- Block 11: the registry surface is opt-in -----------------------------
test_that("communication is an opt-in group and run_cci is a mutating action", {
    object <- sclet_ai_test_cci_object()
    expect_false("run_cci" %in% names(AIDefaultExecutionRegistry(object)))
    registry <- AIDefaultExecutionRegistry(object, include = "all")
    expect_true("run_cci" %in% names(registry))

    descriptor <- registry$run_cci
    expect_true(descriptor$mutates_object)
    expect_true(descriptor$requires_confirmation)
    expect_equal(descriptor$returns, "sce")
    expect_equal(descriptor$estimated_cost, "high")
    expect_true(all(c("communication", "ai_evidence") %in% descriptor$allowed_state_types))
})
