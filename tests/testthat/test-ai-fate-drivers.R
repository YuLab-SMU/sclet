## Lineage-driver contract: readiness, the bounded per-lineage summary, the
## across-lineage comparison, and the two new fate actions.
##
## The whole point of this contract is *withholding*: the stored CellRank driver
## table carries 40 real gene symbols and two real lineage labels ("L1", "L2"),
## and not one of them may leave these functions. The privacy regression below
## is therefore the most valuable test in this file, and it walks the returned
## structures recursively for both character values and field names.
##
## No real Python / basilisk environment is ever built here: the fixture is an
## in-memory SingleCellExperiment carrying a hand-written trajectory state
## record, exactly the way test-ai-fate.R builds its velocity record without
## running velociraptor.

## ---------------------------------------------------------------------------
## Fixture
## ---------------------------------------------------------------------------

## One fixture builder for the whole file. `n_lineages` selects whether the
## synthetic CellRank record describes two lineages (the default, which is what
## every interesting branch needs) or a single one (for the "a comparison needs
## at least two" case); `with_drivers = FALSE` records a completed run whose
## driver computation did not run.
##
## The record mirrors what RunCellRank() writes (R/cellrank.R): a "trajectory"
## state record with method "CellRank", the driver table under
## artifacts$lineage_drivers, and the fate probability column names under
## artifacts$fate_probability_names. `method` is what marks the record as
## CellRank -- the artifact's own `analysis_key` is only a secondary signal.
sclet_ai_test_drivers_object <- function(n_cells = 30L, n_drivers_per_lineage = 20L,
                                        n_lineages = 2L, with_drivers = TRUE) {
    set.seed(20251015L)
    genes <- paste0("gene_", seq_len(60L))
    cells <- paste0("cell_", seq_len(n_cells))
    counts <- matrix(rpois(length(genes) * n_cells, lambda = 6L),
        nrow = length(genes), ncol = n_cells,
        dimnames = list(genes, cells))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(
        stats::rnorm(n_cells * 5L), nrow = n_cells, ncol = 5L,
        dimnames = list(cells, paste0("PC", seq_len(5L)))
    )
    sce$cluster <- factor(paste0("cluster_", seq_len(3L))[seq_len(n_cells) %% 3L + 1L])
    ActiveIdent(sce) <- "cluster"

    lineage_labels <- c("L1", "L2")[seq_len(n_lineages)]
    fate_cols <- paste0("cellrank_fate_", c("A", "B")[seq_len(n_lineages)])

    ## Terminal states are labelled with real cluster names, so the driver
    ## summary has something legitimate to anonymize.
    sce[["cellrank_terminal_states"]] <- rep(c("cluster_1", "cluster_2"), length.out = n_cells)
    for (i in seq_along(fate_cols)) {
        values <- stats::runif(n_cells)
        sce[[fate_cols[[i]]]] <- values
    }

    ## 40 drivers with mixed signs: the sign is exactly what must NOT be
    ## readable from the bounded summary, which reports |metric| aggregates.
    drivers <- NULL
    if (isTRUE(with_drivers)) {
        n_drivers <- n_drivers_per_lineage * length(lineage_labels)
        corr <- rep(c(-0.42, 0.87, -0.13, 0.65, -0.91), length.out = n_drivers)
        drivers <- data.frame(
            gene = paste0("gene_", seq_len(n_drivers)),
            lineage = rep(lineage_labels, each = n_drivers_per_lineage),
            corr = corr,
            stringsAsFactors = FALSE
        )
    }

    artifacts <- list(
        analysis_key = "cr",
        terminal_state_col = "cellrank_terminal_states",
        fate_probability_cols = fate_cols,
        fate_probability_names = lineage_labels,
        lineage_drivers = drivers
    )
    if (isTRUE(with_drivers)) {
        artifacts$has_lineage_drivers <- TRUE
    }

    sclet_set_analysis_state(
        object = sce,
        type = "trajectory",
        id = "cr",
        method = "CellRank",
        inputs = list(reduction = "PCA", cluster_key = "cluster", velocity_id = "vel_fixture"),
        artifacts = artifacts,
        params = list(),
        summary = list(
            n_terminal_states = 2L,
            n_fate_dimensions = as.integer(length(lineage_labels)),
            n_lineage_drivers = if (is.null(drivers)) 0L else as.integer(nrow(drivers))
        ),
        active = TRUE
    )
}

## A bare object with no trajectory record at all, for the not_ready branch.
sclet_ai_test_drivers_bare_object <- function(n_cells = 12L, n_genes = 10L) {
    set.seed(20251015L)
    SingleCellExperiment::SingleCellExperiment(assays = list(
        counts = matrix(rpois(n_genes * n_cells, lambda = 6L),
            nrow = n_genes, ncol = n_cells,
            dimnames = list(paste0("gene_", seq_len(n_genes)), paste0("cell_", seq_len(n_cells))))
    ))
}

## Every real gene symbol and lineage label the fixture stores. The privacy
## test asserts none of these appear in any returned structure.
sclet_ai_test_drivers_secrets <- function() {
    c(paste0("gene_", seq_len(60L)), "L1", "L2")
}

## Walk a nested list and hand every leaf to `visit`. Recursion covers unnamed
## list elements too, so a structure that buries a gene symbol one level deeper
## than expected is still caught.
sclet_ai_test_drivers_walk <- function(node, visit) {
    if (is.list(node)) {
        for (element in node) {
            sclet_ai_test_drivers_walk(element, visit)
        }
    } else {
        visit(node)
    }
    invisible(NULL)
}

## All field names anywhere in a nested list, at every depth.
sclet_ai_test_drivers_field_names <- function(node) {
    names <- character()
    stack <- list(node)
    while (length(stack)) {
        current <- stack[[1L]]
        stack <- stack[-1L]
        if (is.list(current) && !is.null(names(current))) {
            names <- c(names, names(current))
            stack <- c(stack, unname(current))
        }
    }
    unique(names[!is.na(names) & nzchar(names)])
}

## Every character value anywhere in a nested list.
sclet_ai_test_drivers_char_values <- function(node) {
    values <- character()
    sclet_ai_test_drivers_walk(node, function(leaf) {
        if (is.character(leaf)) values <<- c(values, leaf)
    })
    values
}

## Collect the values stored under one field name, at any depth and in any
## number, so a code-bearing field is checked wherever it appears.
sclet_ai_test_drivers_values_for <- function(node, field) {
    found <- list()
    stack <- list(node)
    while (length(stack)) {
        current <- stack[[1L]]
        stack <- stack[-1L]
        if (!is.list(current)) next
        current_names <- names(current)
        if (!is.null(current_names)) {
            for (position in which(current_names == field)) {
                found[[length(found) + 1L]] <- current[[position]]
            }
            stack <- c(stack, unname(current))
        }
    }
    found
}

## TRUE when any leaf of the structure has one entry per cell -- the shape a
## per-cell fate probability vector would take if it ever leaked.
sclet_ai_test_drivers_has_per_cell <- function(node, n_cells) {
    found <- FALSE
    sclet_ai_test_drivers_walk(node, function(leaf) {
        if (is.atomic(leaf) && length(leaf) == n_cells) found <<- TRUE
    })
    found
}

## The anonymized lineage codes the result exposes, whichever field carries
## them (the comparison uses `lineage_codes`, the summary a single
## `lineage_code`).
sclet_ai_test_drivers_codes <- function(result) {
    codes <- character()
    for (field in c("lineage_code", "lineage_codes")) {
        codes <- c(codes, unlist(sclet_ai_test_drivers_values_for(result, field),
            use.names = FALSE))
    }
    codes <- as.character(codes)
    codes[!is.na(codes) & nzchar(codes)]
}

## Validate one one-step plan the way test-ai-fate.R does, with the plan-level
## fingerprint taken from the object's own ledger.
sclet_ai_test_drivers_plan <- function(sce, registry, action, params) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "fate_drivers",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "drv1", action = action, params = as.list(params)))
    )
    ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
}

## -- Scenario 1: readiness ---------------------------------------------------
test_that("check_fate_drivers_readiness is blocked on a bare object and ready on the fixture", {
    bare <- check_fate_drivers_readiness(sclet_ai_test_drivers_bare_object())
    expect_equal(bare$status, "not_ready")
    expect_true("fate_drivers" %in% bare$blocked_actions)
    expect_true(length(bare$questions) > 0L)
    expect_false(isTRUE(bare$checks$cellrank_record_available))
    expect_null(bare$checks$driver_metric_column)

    ready <- check_fate_drivers_readiness(sclet_ai_test_drivers_object())
    expect_equal(ready$status, "ready_for_diagnostic")
    expect_equal(length(ready$blocked_actions), 0L)
    expect_true(isTRUE(ready$checks$cellrank_record_available))
    expect_true(isTRUE(ready$checks$lineage_driver_table_available))
    expect_true(isTRUE(ready$checks$driver_table_has_lineage_column))
    ## The metric column NAME is resolved (never the metric values).
    expect_false(is.null(ready$checks$driver_metric_column))
    expect_equal(ready$checks$driver_metric_column, "corr")
    expect_equal(ready$checks$n_lineages, 2L)
    expect_equal(ready$checks$n_drivers, 40L)
    expect_false(isTRUE(ready$checks$raw_values_included))

    ## A completed CellRank run whose driver computation did not run is still
    ## record-available, but is not silently reported as ready to summarize.
    no_drivers <- check_fate_drivers_readiness(
        sclet_ai_test_drivers_object(with_drivers = FALSE))
    expect_true(isTRUE(no_drivers$checks$cellrank_record_available))
    expect_false(isTRUE(no_drivers$checks$lineage_driver_table_available))
    expect_false(isTRUE(no_drivers$checks$driver_metric_column))
    expect_true(length(no_drivers$questions) > 0L)
})

## -- Scenario 2: an unspecified lineage must be asked about ------------------
test_that("summarize_fate_lineage_drivers asks which lineage instead of choosing one", {
    sce <- sclet_ai_test_drivers_object()
    result <- summarize_fate_lineage_drivers(sce, lineage = NULL)

    expect_equal(result$status, "clarification_required")
    expect_equal(result$n_available_lineages, 2L)
    expect_true(is.character(result$question) && nzchar(result$question))
    expect_false(result$raw_values_included)

    ## The question must not name the lineages it is asking about: asking
    ## "which of L1/L2?" would leak the very labels this contract withholds.
    expect_false(grepl("L1", result$question, fixed = TRUE))
    expect_false(grepl("L2", result$question, fixed = TRUE))
})

## -- Scenario 3: an explicitly named lineage is summarized -------------------
test_that("summarize_fate_lineage_drivers summarizes one named lineage as route codes", {
    sce <- sclet_ai_test_drivers_object()
    result <- summarize_fate_lineage_drivers(sce, lineage = "L1", top_n = 7L)

    expect_equal(result$status, "available")
    expect_true(grepl("^(group|cluster|condition|sample|route)_[0-9]+$", result$lineage_code))
    expect_equal(result$n_drivers_total, 20L)
    expect_equal(result$n_drivers_summarized, 7L)
    expect_equal(result$driver_metric, "corr")
    expect_equal(length(result$driver_rank_codes), 7L)
    expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$",
        result$driver_rank_codes)))
    expect_equal(result$claim_level, "associated")
    expect_false(result$raw_values_included)

    ## The metric summary is bounded and reports |metric| magnitudes only.
    summary <- result$driver_abs_metric_summary
    expect_equal(summary$n_observed, 7L)
    expect_true(is.numeric(summary$mean) && is.finite(summary$mean))
    expect_true(summary$mean >= 0)

    ## top_n is capped by the number of drivers that lineage actually has.
    all_rows <- summarize_fate_lineage_drivers(sce, lineage = "L1", top_n = 1000L)
    expect_equal(all_rows$n_drivers_summarized, 20L)

    ## An unknown lineage is refused rather than silently defaulted.
    unknown <- summarize_fate_lineage_drivers(sce, lineage = "not_a_lineage")
    expect_equal(unknown$status, "not_available")
    expect_equal(unknown$reason, "lineage_unknown")
})

## -- Scenario 4: the privacy regression --------------------------------------
test_that("driver summaries and comparisons never leak a gene symbol or lineage label", {
    sce <- sclet_ai_test_drivers_object()
    secrets <- sclet_ai_test_drivers_secrets()
    results <- list(
        readiness = check_fate_drivers_readiness(sce),
        summary = summarize_fate_lineage_drivers(sce, lineage = "L1", top_n = 10L),
        summary_clarify = summarize_fate_lineage_drivers(sce, lineage = NULL),
        comparison = compare_fate_lineages(sce)
    )

    for (name in names(results)) {
        result <- results[[name]]
        field_names <- sclet_ai_test_drivers_field_names(result)

        ## (a) No gene symbol or lineage label anywhere -- not in values, not
        ## in field names, not in a printed rendering of the whole structure.
        text <- paste(capture.output(str(result)), collapse = "\n")
        for (secret in secrets) {
            if (grepl("^gene_", secret)) {
                ## gene_1 is a prefix of gene_10..gene_19, so match on a
                ## delimiter rather than a bare substring.
                expect_false(grepl(paste0(secret, "([^0-9]|$)"), text),
                    info = paste(name, "leaked", secret))
            } else {
                expect_false(grepl(secret, text, fixed = TRUE),
                    info = paste(name, "leaked", secret))
            }
            expect_false(secret %in% sclet_ai_test_drivers_field_names(result),
                info = paste(name, "leaked", secret, "as a field name"))
        }

        ## (b) Every IDENTITY field is an anonymized code (EX4). The rule cannot
        ## apply to the whole structure: the contract requires free-text
        ## `questions`, prose `caveats` and `driver_metric` (a column name such
        ## as "corr"). The already-green check_fate_readiness() emits prose
        ## notes/questions and would fail a whole-structure version of this
        ## rule, so it is scoped to the fields that carry lineage identity.
        identity_fields <- c("lineage_code", "lineage_codes", "driver_rank_codes",
            "terminal_state_codes")
        for (field in intersect(identity_fields, field_names)) {
            values <- unlist(sclet_ai_test_drivers_values_for(result, field), use.names = FALSE)
            values <- as.character(values)
            values <- values[!is.na(values) & nzchar(values)]
            if (!length(values)) next
            expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$", values)),
                info = paste(name, "returned a non-anonymized value in", field))
        }
        ## And the lineage identity the result exposes really is code-shaped.
        codes <- sclet_ai_test_drivers_codes(result)
        if (length(codes)) {
            expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$", codes)),
                info = paste(name, "returned a non-anonymized lineage code"))
        }

        ## (c) Raw values are always withheld.
        if (!is.null(result$raw_values_included)) {
            expect_false(isTRUE(result$raw_values_included))
        }

        ## (d) No per-cell vector: nothing in the structure has one entry per
        ## cell, which is how a per-cell fate probability would leak.
        expect_false(sclet_ai_test_drivers_has_per_cell(result, ncol(sce)))
    }

    ## The comparison hides gene identity inside per-lineage summaries too.
    comparison <- results$comparison
    expect_equal(comparison$status, "available")
    for (lineage in comparison$lineages) {
        expect_false(isTRUE(lineage$raw_values_included))
    }
})

## -- Scenario 5: across-lineage comparison -----------------------------------
test_that("compare_fate_lineages compares lineages without selecting one", {
    sce <- sclet_ai_test_drivers_object()
    result <- compare_fate_lineages(sce)

    expect_equal(result$status, "available")
    expect_equal(result$n_lineages, 2L)
    expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$",
        result$lineage_codes)))
    expect_equal(length(result$lineage_codes), 2L)
    expect_equal(length(unique(result$lineage_codes)), 2L)
    expect_equal(sort(names(result$lineages)), sort(result$lineage_codes))
    expect_equal(result$claim_level, "associated")
    expect_false(result$raw_values_included)

    ## Distinctness is a bounded numeric fraction in [0, 1], and the fixture's
    ## two lineages share no driver, so the overlap is exactly 0.
    expect_true(is.numeric(result$distinctness$jaccard_mean))
    expect_true(result$distinctness$jaccard_mean >= 0 &&
        result$distinctness$jaccard_mean <= 1)
    expect_equal(result$distinctness$jaccard_min, 0)
    expect_equal(result$distinctness$jaccard_max, 0)
    expect_equal(result$distinctness$n_pairs, 1L)
    expect_false(isTRUE(result$distinctness$raw_values_included))

    ## A single-lineage record cannot be compared against anything.
    single <- compare_fate_lineages(sclet_ai_test_drivers_object(n_lineages = 1L))
    expect_equal(single$status, "not_available")
    expect_equal(single$reason, "at_least_two_lineages_required")
    expect_equal(single$n_lineages, 1L)

    ## Nor can an object with no CellRank record at all.
    none <- compare_fate_lineages(sclet_ai_test_drivers_bare_object())
    expect_equal(none$status, "not_available")
    expect_equal(none$reason, "at_least_two_lineages_required")
})

## -- Scenario 6: nothing ranks, suggests or nominates ------------------------
test_that("no driver function ranks, suggests or nominates a lineage", {
    sce <- sclet_ai_test_drivers_object()
    results <- list(
        check_fate_drivers_readiness(sce),
        check_fate_drivers_readiness(sclet_ai_test_drivers_bare_object()),
        summarize_fate_lineage_drivers(sce, lineage = "L1"),
        summarize_fate_lineage_drivers(sce, lineage = NULL),
        compare_fate_lineages(sce)
    )

    for (result in results) {
        field_names <- sclet_ai_test_drivers_field_names(result)

        ## Nothing that reads like a nomination, at any depth (EX1). The pattern
        ## is token-delimited rather than a bare substring: a plain "rank" also
        ## matches "cellrank", and `cellrank_record_available` /
        ## `cellrank_id_resolved` are contract-mandated check names.
        ##
        ## The delimiter class is [^a-z0-9] rather than \b because R's \b
        ## treats "_" as a word character, so a \b pattern silently misses
        ## "suggest_lineage" -- exactly the field the guard exists to catch.
        ## `[a-z]*` after the verb group covers inflections ("suggested",
        ## "recommending"), which a bare stem-plus-delimiter alternation drops.
        ## `driver_rank_codes` and the explicit-null `recommendation`
        ## placeholder are required by the contract, so all four names are
        ## legitimate rather than violations.
        offending <- grep("(^|[^a-z0-9])(suggest|recommend|best|winner)[a-z]*([^a-z0-9]|$)|(^|[^a-z0-9])rank(ing)?([^a-z0-9]|$)|top_lineage",
            field_names, ignore.case = TRUE, value = TRUE)
        offending <- setdiff(offending, c("cellrank_record_available", "cellrank_id_resolved",
            "driver_rank_codes", "recommendation"))
        expect_equal(unique(offending), character(0))

        ## EX2: `recommendation` must exist and be NULL on the comparison,
        ## mirroring compare_trajectory_roots() -- the refusal is explicit
        ## rather than merely silent.
        if ("recommendation" %in% field_names) {
            expect_null(result$recommendation)
        }
    }

    ## The comparison carries the explicit null placeholder.
    comparison <- compare_fate_lineages(sce)
    expect_true("recommendation" %in% names(comparison))
    expect_null(comparison$recommendation)

    ## The guard still catches a real nomination, in bare, underscore-joined
    ## and inflected form, so it has not been hollowed out.
    nomination_pattern <- "(^|[^a-z0-9])(suggest|recommend|best|winner)[a-z]*([^a-z0-9]|$)|(^|[^a-z0-9])rank(ing)?([^a-z0-9]|$)|top_lineage"
    expect_equal(length(grep(nomination_pattern,
        c("suggested_lineage", "suggest_lineage", "suggest", "suggestion", "ranking",
            "best_lineage", "winner", "top_lineage"),
        ignore.case = TRUE)), 8L)
    ## ...while the contract-mandated cellrank fields are genuinely untouched.
    expect_equal(length(grep(nomination_pattern,
        c("cellrank_record_available", "cellrank_id_resolved"),
        ignore.case = TRUE)), 0L)
})

## -- Scenario 7: confirm_fate_terminal_states prerequisites ------------------
test_that("confirm_fate_terminal_states rejects missing and unknown terminal states", {
    sce <- sclet_ai_test_drivers_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
    expect_true("confirm_fate_terminal_states" %in% names(registry))

    ## A missing terminal_state_codes is rejected with a typed reason.
    missing_codes <- sclet_ai_test_drivers_plan(
        sce, registry, "confirm_fate_terminal_states", list())
    expect_false(missing_codes$valid)
    ## ValidateAIPlan prefixes prerequisite errors with
    ## "prerequisites not met for action <step_id>: ", so the token is matched
    ## anywhere in the error and never with a leading ^ anchor.
    expect_true(any(grepl("terminal_state_codes_missing", missing_codes$errors)))

    ## An unknown code is rejected as well.
    unknown_code <- sclet_ai_test_drivers_plan(
        sce, registry, "confirm_fate_terminal_states",
        list(terminal_state_codes = "cluster_99"))
    expect_false(unknown_code$valid)
    expect_true(any(grepl("terminal_state_unknown", unknown_code$errors)))

    ## A code that really is recorded as a terminal state is accepted. The
    ## fixture's terminal-state column holds exactly two labels, which the
    ## action exposes positionally as cluster_1 / cluster_2.
    known_code <- sclet_ai_test_drivers_plan(
        sce, registry, "confirm_fate_terminal_states",
        list(terminal_state_codes = "cluster_1"))
    expect_true(known_code$valid)

    ## Without a CellRank record the action is blocked before any code check,
    ## so the missing-record reason wins over the missing-codes one.
    bare <- sclet_ai_test_drivers_bare_object()
    bare_registry <- AIDefaultExecutionRegistry(bare, include = c("read", "fate"))
    no_record <- sclet_ai_test_drivers_plan(
        bare, bare_registry, "confirm_fate_terminal_states",
        list(terminal_state_codes = "cluster_1"))
    expect_false(no_record$valid)
    expect_true(any(grepl("cellrank_record_missing", no_record$errors)))
})

## -- Scenario 8: summarize_fate_drivers is a read-only action ----------------
test_that("summarize_fate_drivers is registered as a read-only action needing no confirmation", {
    sce <- sclet_ai_test_drivers_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
    expect_true("summarize_fate_drivers" %in% names(registry))

    action <- registry$summarize_fate_drivers
    expect_s3_class(action, "sclet_ai_action")
    expect_equal(action$name, "summarize_fate_drivers")
    expect_false(action$mutates_object)
    expect_false(action$requires_confirmation)
    ## A read-only action returns its bounded summary directly rather than a
    ## modified object. "value" is the established kind for "the handler
    ## returns a non-SCE value": the three existing read-only actions in the
    ## `read` group (inspect_status, inspect_ledger, check_qc) all use it, and
    ## every consumer only tests identical(returns, "sce").
    expect_equal(action$returns, "value")
    expect_true(action$estimated_cost %in% c("low", "medium", "high"))

    ## `allowed_state_types` declares which state types the action *reads*; a
    ## read-only action still reads the trajectory and evidence records. The
    ## read-only guarantee is carried by mutates_object, which is FALSE.
    expect_true(all(c("trajectory", "ai_evidence") %in% action$allowed_state_types))

    ## The handler really is the bounded driver summary, and it is reachable
    ## without a confirmation token.
    expect_true(is.function(action$handler))
    result <- action$handler(sce, list(lineage = "L1", top_n = 5L))
    expect_equal(result$status, "available")
    expect_equal(result$n_drivers_summarized, 5L)
    expect_false(result$raw_values_included)
})
