# The execution trace as a graph: the record must survive the round.
#
# Every case here reads what ExecuteAIPlan() WROTE, never what it returned.
# That is the whole point of the change: a trace that only exists in the return
# value dies with the R session, and the next round cannot see it.

tg_sce <- function() {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(60, 5), nrow = 10L, ncol = 6L))
    )
    SummarizedExperiment::colData(sce)$condition <- rep(c("ctrl", "stim"), each = 3L)
    sce
}

# One action that writes one named record of `type`.
tg_writer <- function(name, record_id, type = "preprocess", summary = list(n = 1)) {
    AIAction(
        name = name,
        handler = function(object, params) {
            sclet:::sclet_set_analysis_state(
                object,
                type = type,
                id = record_id,
                method = "tg_mock",
                summary = summary
            )
        },
        returns = "sce",
        allowed_state_types = type
    )
}

tg_registry <- function(...) {
    AIExecutionRegistry(list(...))
}

tg_plan <- function(sce, actions) {
    sclet:::new_sclet_ai_plan(
        task = "trace_case",
        actions = actions,
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
}

tg_run <- function(sce, actions, registry) {
    plan <- tg_plan(sce, actions)
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    if (!isTRUE(validation$valid)) {
        stop("fixture plan is invalid: ", paste(unlist(validation$errors), collapse = "; "))
    }
    list(
        plan = plan,
        validation = validation,
        executed = ExecuteAIPlan(
            sce, plan, registry,
            validation = validation,
            dry_run = FALSE,
            confirmation = validation$confirmation_token
        )
    )
}

tg_edges <- function(graph, kind) {
    graph$edges[graph$edges$kind == kind, , drop = FALSE]
}

test_that("the plan's topology survives execution", {
    registry <- tg_registry(
        a = tg_writer("a", "rec_a"),
        b = tg_writer("b", "rec_b"),
        c = tg_writer("c", "rec_c"),
        d = tg_writer("d", "rec_d"),
        e = tg_writer("e", "rec_e")
    )
    run <- tg_run(tg_sce(), list(
        list(id = "s1", action = "a"),
        list(id = "s2", action = "b", depends_on = "s1"),
        list(id = "s3", action = "c", depends_on = "s2"),
        list(id = "s4", action = "d", depends_on = "s2"),
        list(id = "s5", action = "e", depends_on = c("s3", "s4"))
    ), registry)

    graph <- AnalysisGraph(run$executed$object)
    expect_equal(graph$meta$n_rounds, 1L)
    expect_equal(graph$meta$legacy_rounds, 0L)

    deps <- tg_edges(graph, "depends_on")
    pairs <- paste(deps$from, deps$to, sep = " -> ")
    round_id <- paste0("round:", run$executed$execution_id)
    step <- function(id) paste0("step:", run$executed$execution_id, ":", id)
    expect_setequal(
        pairs,
        c(paste(step("s2"), step("s1"), sep = " -> "),
          paste(step("s3"), step("s2"), sep = " -> "),
          paste(step("s4"), step("s2"), sep = " -> "),
          paste(step("s5"), step("s3"), sep = " -> "),
          paste(step("s5"), step("s4"), sep = " -> "))
    )
    expect_true(all(tg_edges(graph, "step_of")$to == round_id))

    expect_setequal(graph_consumers(graph, step("s2"))$ids, c(step("s3"), step("s4")))
    expect_true(graph_consumers(graph, step("s2"))$complete)
    expect_setequal(graph_producers(graph, "record:preprocess::rec_b")$ids, step("s2"))
})

test_that("a step that writes two records owns both, and only it", {
    two <- AIAction(
        name = "two",
        handler = function(object, params) {
            object <- sclet:::sclet_set_analysis_state(object, type = "preprocess",
                id = "rec_two_1", method = "tg_mock", summary = list(n = 1))
            sclet:::sclet_set_analysis_state(object, type = "preprocess",
                id = "rec_two_2", method = "tg_mock", summary = list(n = 2))
        },
        returns = "sce",
        allowed_state_types = "preprocess"
    )
    registry <- tg_registry(one = tg_writer("one", "rec_one"), two = two)
    run <- tg_run(tg_sce(), list(
        list(id = "s1", action = "one"),
        list(id = "s2", action = "two", depends_on = "s1")
    ), registry)

    graph <- AnalysisGraph(run$executed$object)
    produced <- tg_edges(graph, "produced")
    two_step <- paste0("step:", run$executed$execution_id, ":s2")
    expect_setequal(
        produced$to[produced$from == two_step],
        c("record:preprocess::rec_two_1", "record:preprocess::rec_two_2")
    )
    expect_equal(sum(produced$to == "record:preprocess::rec_two_1"), 1L)
    expect_equal(sum(produced$to == "record:preprocess::rec_two_2"), 1L)
    expect_setequal(graph_producers(graph, "record:preprocess::rec_one")$ids,
        paste0("step:", run$executed$execution_id, ":s1"))
})

test_that("a re-run of the same record id is refreshed, not produced", {
    rewrite <- AIAction(
        name = "rewrite",
        handler = function(object, params) {
            Sys.sleep(0.02) # created_at must advance for "refreshed" to be visible
            sclet:::sclet_set_analysis_state(object, type = "preprocess",
                id = "rec_same", method = "tg_mock",
                summary = list(round = if (is.null(params$round)) 0L else as.integer(params$round)))
        },
        returns = "sce",
        allowed_state_types = "preprocess"
    )
    registry <- tg_registry(rewrite = rewrite)
    run <- tg_run(tg_sce(), list(
        list(id = "s1", action = "rewrite", params = list(round = 1L)),
        list(id = "s2", action = "rewrite", params = list(round = 2L), depends_on = "s1")
    ), registry)

    graph <- AnalysisGraph(run$executed$object)
    produced <- tg_edges(graph, "produced")
    refreshed <- tg_edges(graph, "refreshed")
    expect_equal(nrow(produced[produced$to == "record:preprocess::rec_same", , drop = FALSE]), 1L)
    expect_equal(nrow(refreshed[refreshed$to == "record:preprocess::rec_same", , drop = FALSE]), 1L)
    expect_identical(
        produced$from[produced$to == "record:preprocess::rec_same"],
        paste0("step:", run$executed$execution_id, ":s1")
    )
    expect_identical(
        refreshed$from[refreshed$to == "record:preprocess::rec_same"],
        paste0("step:", run$executed$execution_id, ":s2")
    )
})

test_that("a legacy round keeps its nodes and evidence edges, and says it is incomplete", {
    sce <- tg_sce()
    sce <- sclet:::sclet_set_analysis_state(sce, type = "preprocess", id = "rec_legacy",
        method = "tg_mock", summary = list(n = 1))
    legacy <- list(
        id = "ai_execution_legacy",
        type = "ai_execution",
        method = "ExecuteAIPlan",
        inputs = list(plan_id = "plan_legacy", context_fingerprint = "sclet-ai-00000000",
            action_ids = "s1"),
        params = list(dry_run = FALSE, requires_confirmation = TRUE, steps = list(
            list(id = "s1", action = "write_legacy", status = "completed",
                duration_sec = 0.1, params = list())
        )),
        artifacts = list(results = list()),
        summary = list(status = "completed", n_actions = 1L, n_completed = 1L, n_failed = 0L),
        created_at = Sys.time()
    )
    sce <- sclet:::sclet_set_analysis(sce, "ai_execution_legacy", legacy)
    sce <- sclet:::sclet_set_analysis_state(sce, type = "ai_evidence", id = "ev_legacy",
        method = "sclet_ai_evidence", active = FALSE,
        inputs = list(source = "rec_legacy", parents = character()),
        summary = list(id = "ev_legacy", kind = "test", source = "rec_legacy",
            parents = character(), claim_level = "measured", status = "completed"))

    graph <- AnalysisGraph(sce)
    expect_equal(graph$meta$legacy_rounds, 1L)
    expect_equal(graph$meta$n_rounds, 1L)
    round <- graph$nodes[graph$nodes$id == "round:ai_execution_legacy", , drop = FALSE]
    expect_equal(nrow(round), 1L)
    expect_false(isTRUE(round$complete))

    # Evidence edges are the ones that already worked: a legacy round must not
    # lose them (rev 1 would have returned empty edge tables here).
    source_of <- tg_edges(graph, "source_of")
    expect_equal(nrow(source_of), 1L)
    expect_identical(source_of$from, "evidence:ev_legacy")
    expect_identical(source_of$to, "record:preprocess::rec_legacy")

    expect_false(graph_producers(graph, "record:preprocess::rec_legacy")$complete)
})

test_that("the graph survives saveRDS", {
    registry <- tg_registry(a = tg_writer("a", "rec_a"))
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "a")), registry)
    graph <- AnalysisGraph(run$executed$object)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(graph, path)
    restored <- readRDS(path)
    expect_identical(restored$nodes, graph$nodes)
    expect_identical(restored$edges, graph$edges)
    expect_identical(restored$meta, graph$meta)
})

test_that("node and edge tables are a whitelist, not a copy of the records", {
    leaky <- AIAction(
        name = "leaky",
        handler = function(object, params) {
            sclet:::sclet_set_analysis_state(
                object,
                type = "preprocess",
                id = "rec_leaky",
                method = "tg_mock",
                params = list(barcode = "AAACCTGAGA"),
                summary = list(top_cell = "AAACCTGAGA", n = 1)
            )
        },
        returns = "sce",
        allowed_state_types = "preprocess"
    )
    registry <- tg_registry(leaky = leaky)
    sce <- tg_sce()
    # A value that lives in colData must not reach the graph either: the
    # projection reads the record store, and this pins that it stays that way.
    SummarizedExperiment::colData(sce)$cell_barcode <- rep("AAACCTGAGA", ncol(sce))
    run <- tg_run(sce, list(list(id = "s1", action = "leaky")), registry)
    graph <- AnalysisGraph(run$executed$object)

    expect_identical(
        names(graph$nodes),
        c("id", "kind", "label", "type", "status", "created_at", "complete",
          "action", "duration_sec", "attempts")
    )
    expect_identical(
        names(graph$edges),
        c("from", "to", "kind", "parent_kind", "declared_input")
    )
    flat <- c(unlist(graph$nodes, use.names = FALSE), unlist(graph$edges, use.names = FALSE))
    expect_false(any(grepl("AAACCTGAGA", as.character(flat), fixed = TRUE)))
    # The record itself still holds the value; the graph is the redacted view.
    ledger <- GetAnalysisLedger(run$executed$object, detail = "full", include_artifacts = TRUE)
    expect_true("rec_leaky" %in% names(ledger$state_records$preprocess))
})

test_that("dry runs and unrecorded executions leave no trace", {
    registry <- tg_registry(a = tg_writer("a", "rec_a"))
    sce <- tg_sce()
    plan <- tg_plan(sce, list(list(id = "s1", action = "a")))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)

    dry <- ExecuteAIPlan(sce, plan, registry, validation = validation, dry_run = TRUE)
    expect_true(dry$dry_run)
    unrecorded <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token, record = FALSE)
    expect_false(unrecorded$recorded)

    graph <- AnalysisGraph(unrecorded$object)
    expect_equal(graph$meta$n_rounds, 0L)
    # The action still wrote `rec_a`; what must not exist is a round or a step.
    expect_equal(sum(graph$nodes$kind %in% c("round", "step")), 0L)
    expect_equal(nrow(graph$edges), 0L)
    # No round is "no trace", and a query must not call that complete.
    expect_false(graph_producers(graph, "record:preprocess::rec_a")$complete)
})

test_that("a value-returning action is a step without record edges", {
    registry <- tg_registry(
        value_action = AIAction(
            name = "value_action",
            handler = function(object, params) list(n = ncol(object)),
            returns = "value"
        )
    )
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "value_action")), registry)

    graph <- AnalysisGraph(run$executed$object)
    step <- graph$nodes[graph$nodes$kind == "step", , drop = FALSE]
    expect_equal(nrow(step), 1L)
    expect_equal(step$status, "completed")
    expect_equal(nrow(tg_edges(graph, "produced")), 0L)
    expect_equal(nrow(tg_edges(graph, "refreshed")), 0L)
})

test_that("a failed step has no record edges and stops the round", {
    failing <- AIAction(
        name = "failing",
        handler = function(object, params) stop("tg fixture failure", call. = FALSE),
        returns = "sce"
    )
    registry <- tg_registry(
        a = tg_writer("a", "rec_a"),
        failing = failing,
        c = tg_writer("c", "rec_c")
    )
    sce <- tg_sce()
    plan <- tg_plan(sce, list(
        list(id = "s1", action = "a"),
        list(id = "s2", action = "failing", depends_on = "s1"),
        list(id = "s3", action = "c", depends_on = "s2")
    ))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)

    graph <- AnalysisGraph(executed$object)
    steps <- graph$nodes[graph$nodes$kind == "step", , drop = FALSE]
    expect_equal(nrow(steps), 2L)
    failed <- steps[steps$label == "failing", , drop = FALSE]
    expect_equal(failed$status, "failed")
    produced <- tg_edges(graph, "produced")
    expect_equal(sum(produced$from == failed$id), 0L)
    expect_false(paste0("step:", executed$execution_id, ":s3") %in% graph$nodes$id)
})

test_that("a retried step is one node that records its attempts", {
    flaky <- AIAction(
        name = "flaky",
        handler = local({
            seen <- 0L
            function(object, params) {
                seen <<- seen + 1L
                if (seen < 2L) stop("first attempt fails", call. = FALSE)
                sclet:::sclet_set_analysis_state(object, type = "preprocess",
                    id = "rec_flaky", method = "tg_mock", summary = list(n = 1))
            }
        }),
        returns = "sce",
        allowed_state_types = "preprocess",
        idempotent = TRUE
    )
    registry <- tg_registry(flaky = flaky)
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "flaky", max_retries = 1L)), registry)

    graph <- AnalysisGraph(run$executed$object)
    steps <- graph$nodes[graph$nodes$kind == "step", , drop = FALSE]
    expect_equal(nrow(steps), 1L)
    expect_equal(steps$attempts, 2L)
    expect_setequal(tg_edges(graph, "produced")$to, "record:preprocess::rec_flaky")
})

test_that("record-to-record edges come only from declared inputs", {
    chained <- AIAction(
        name = "chained",
        handler = function(object, params) {
            object <- sclet:::sclet_set_analysis_state(object, type = "graph",
                id = "knn_graph", method = "tg_mock", summary = list(n = 1))
            sclet:::sclet_set_analysis_state(object, type = "clustering",
                id = "louvain", method = "tg_mock",
                inputs = list(graph = "knn_graph"),
                summary = list(n = 1))
        },
        returns = "sce",
        allowed_state_types = c("graph", "clustering")
    )
    registry <- tg_registry(chained = chained)
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "chained")), registry)

    graph <- AnalysisGraph(run$executed$object)
    consumed <- tg_edges(graph, "consumed")
    expect_equal(nrow(consumed), 1L)
    expect_identical(consumed$from, "record:clustering::louvain")
    expect_identical(consumed$to, "record:graph::knn_graph")
    expect_identical(consumed$declared_input, "graph")
})

test_that("an evidence parent that is a record id is typed as a record", {
    evidential <- AIAction(
        name = "evidential",
        handler = function(object, params) {
            object <- sclet:::sclet_set_analysis_state(object, type = "preprocess",
                id = "rec_source", method = "tg_mock", summary = list(n = 1))
            sclet:::sclet_set_analysis_state(object, type = "ai_evidence",
                id = "ev_child", method = "sclet_ai_evidence", active = FALSE,
                inputs = list(source = "rec_source", parents = "rec_source"),
                summary = list(id = "ev_child", kind = "test", source = "rec_source",
                    parents = "rec_source", claim_level = "measured", status = "completed"))
        },
        returns = "sce",
        allowed_state_types = c("preprocess", "ai_evidence")
    )
    registry <- tg_registry(evidential = evidential)
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "evidential")), registry)

    graph <- AnalysisGraph(run$executed$object)
    derived <- tg_edges(graph, "derived_from")
    expect_equal(nrow(derived), 1L)
    expect_identical(derived$from, "evidence:ev_child")
    expect_identical(derived$to, "record:preprocess::rec_source")
    expect_identical(derived$parent_kind, "record")
    # The evidence node's `inputs` mirror source/parents; `consumed` must not
    # restate those edges under a vaguer kind.
    expect_equal(nrow(tg_edges(graph, "consumed")), 0L)
    expect_equal(graph$meta$unresolved_edges, 0L)
    expect_setequal(graph_producers(graph, "record:preprocess::rec_source")$ids,
        paste0("step:", run$executed$execution_id, ":s1"))
})

test_that("continue_on_error keeps the planner edge and records the failure", {
    failing <- AIAction(
        name = "failing",
        handler = function(object, params) stop("tg continue_on_error", call. = FALSE),
        returns = "sce"
    )
    registry <- tg_registry(
        a = tg_writer("a", "rec_a"),
        failing = failing,
        b = tg_writer("b", "rec_b")
    )
    run <- tg_run(tg_sce(), list(
        list(id = "s1", action = "a"),
        list(id = "s2", action = "failing", depends_on = "s1", continue_on_error = TRUE),
        list(id = "s3", action = "b", depends_on = "s2")
    ), registry)

    graph <- AnalysisGraph(run$executed$object)
    steps <- graph$nodes[graph$nodes$kind == "step", , drop = FALSE]
    expect_equal(nrow(steps), 3L)
    failed <- steps[steps$label == "failing", , drop = FALSE]
    expect_equal(failed$status, "failed")
    expect_equal(nrow(tg_edges(graph, "produced")[tg_edges(graph, "produced")$from == failed$id, , drop = FALSE]), 0L)
    dependents <- tg_edges(graph, "depends_on")
    expect_true(paste0("step:", run$executed$execution_id, ":s3") %in% dependents$from)
})

test_that("a re-run that writes a NEW id is produced, and a repeated parent is not a duplicate edge", {
    rewriter <- AIAction(
        name = "rewriter",
        handler = function(object, params) {
            id <- if (identical(params$round, 1L)) "rec_new_1" else "rec_new_2"
            sclet:::sclet_set_analysis_state(object, type = "preprocess",
                id = id, method = "tg_mock", summary = list(n = 1))
        },
        returns = "sce",
        allowed_state_types = "preprocess"
    )
    registry <- tg_registry(rewriter = rewriter)
    run <- tg_run(tg_sce(), list(
        list(id = "s1", action = "rewriter", params = list(round = 1L)),
        list(id = "s2", action = "rewriter", params = list(round = 2L), depends_on = "s1")
    ), registry)

    graph <- AnalysisGraph(run$executed$object)
    produced <- tg_edges(graph, "produced")
    expect_setequal(produced$to,
        c("record:preprocess::rec_new_1", "record:preprocess::rec_new_2"))
    expect_equal(nrow(tg_edges(graph, "refreshed")), 0L)

    # A repeated parent must not become two identical edges.
    evidential <- AIAction(
        name = "evidential",
        handler = function(object, params) {
            object <- sclet:::sclet_set_analysis_state(object, type = "preprocess",
                id = "rec_parent", method = "tg_mock", summary = list(n = 1))
            sclet:::sclet_set_analysis_state(object, type = "ai_evidence",
                id = "ev_dup", method = "sclet_ai_evidence", active = FALSE,
                inputs = list(source = "rec_parent", parents = c("rec_parent", "rec_parent")),
                summary = list(id = "ev_dup", kind = "test", source = "rec_parent",
                    parents = c("rec_parent", "rec_parent"), claim_level = "measured",
                    status = "completed"))
        },
        returns = "sce",
        allowed_state_types = c("preprocess", "ai_evidence")
    )
    run2 <- tg_run(tg_sce(), list(list(id = "s1", action = "evidential")),
        tg_registry(evidential = evidential))
    graph2 <- AnalysisGraph(run2$executed$object)
    derived <- tg_edges(graph2, "derived_from")
    expect_equal(nrow(derived), 1L)
})

test_that("an evidence node with no source does not break the projection", {
    # RecordAIEvidence(source = NULL) is a supported call; the projection used
    # to abort on it because `as.character(NULL)` has length 0.
    sce <- tg_sce()
    sce <- sclet:::sclet_set_analysis_state(sce, type = "ai_evidence", id = "ev_sourceless",
        method = "sclet_ai_evidence", active = FALSE,
        inputs = list(parents = character()),
        summary = list(id = "ev_sourceless", kind = "test", source = NULL,
            parents = character(), claim_level = "measured", status = "completed"))

    graph <- expect_silent(AnalysisGraph(sce))
    expect_equal(nrow(graph$edges), 0L)
    expect_equal(graph$meta$unresolved_edges, 0L)
    expect_true("evidence:ev_sourceless" %in% graph$nodes$id)
})

test_that("a round with no recorded steps is not called complete", {
    sce <- tg_sce()
    empty <- list(
        id = "ai_execution_empty",
        type = "ai_execution",
        method = "ExecuteAIPlan",
        inputs = list(plan_id = "plan_empty", context_fingerprint = "sclet-ai-00000000",
            action_ids = character()),
        params = list(dry_run = FALSE, steps = list()),
        artifacts = list(results = list()),
        summary = list(status = "failed", n_actions = 0L, n_completed = 0L, n_failed = 0L),
        created_at = Sys.time()
    )
    sce <- sclet:::sclet_set_analysis(sce, "ai_execution_empty", empty)
    graph <- AnalysisGraph(sce)
    round <- graph$nodes[graph$nodes$id == "round:ai_execution_empty", , drop = FALSE]
    expect_equal(nrow(round), 1L)
    expect_false(isTRUE(round$complete))
    expect_equal(graph$meta$legacy_rounds, 1L)
    # No rounds at all is also "no trace", not "complete".
    expect_false(graph_neighbors(graph, "round:ai_execution_empty")$complete)
})

test_that("record stamps degrade to unknown instead of aborting the round", {
    expect_true(is.na(sclet:::sclet_ai_record_stamp(NULL)))
    expect_true(is.na(sclet:::sclet_ai_record_stamp("not a time at all")))
    expect_true(is.na(sclet:::sclet_ai_record_stamp(NA_character_)))
    expect_equal(sclet:::sclet_ai_record_stamp(1717171717), 1717171717)
    iso <- sclet:::sclet_ai_record_stamp("2026-10-22T22:46:00")
    expect_equal(iso, as.numeric(as.POSIXct("2026-10-22 22:46:00")))
    expect_true(is.finite(iso))
})

test_that("a consumed edge comes from the declared field, including detest_id", {
    chained <- AIAction(
        name = "chained",
        handler = function(object, params) {
            object <- sclet:::sclet_set_analysis_state(object, type = "detest",
                id = "markers", method = "tg_mock", summary = list(n = 1))
            sclet:::sclet_set_analysis_state(object, type = "enrichment",
                id = "enr", method = "tg_mock",
                inputs = list(detest_id = "markers"),
                summary = list(n = 1))
        },
        returns = "sce",
        allowed_state_types = c("detest", "enrichment")
    )
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "chained")),
        tg_registry(chained = chained))
    graph <- AnalysisGraph(run$executed$object)
    consumed <- tg_edges(graph, "consumed")
    expect_equal(nrow(consumed), 1L)
    expect_identical(consumed$from, "record:enrichment::enr")
    expect_identical(consumed$to, "record:detest::markers")
    expect_identical(consumed$declared_input, "detest_id")
})

test_that("an unrecorded execution leaves no ai_execution record at all", {
    registry <- tg_registry(a = tg_writer("a", "rec_a"))
    sce <- tg_sce()
    plan <- tg_plan(sce, list(list(id = "s1", action = "a")))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token, record = FALSE)
    analyses <- sclet:::sclet_get_state(executed$object)$analyses
    expect_false(any(vapply(analyses, function(rec) {
        identical(as.character(rec$type %||% ""), "ai_execution")
    }, logical(1L))))
})

test_that("a typed _id field resolves to its own record type", {
    chained <- AIAction(
        name = "chained",
        handler = function(object, params) {
            object <- sclet:::sclet_set_analysis_state(object, type = "graph",
                id = "markers", method = "tg_mock", summary = list(n = 1))
            object <- sclet:::sclet_set_analysis_state(object, type = "detest",
                id = "markers", method = "tg_mock", summary = list(n = 2))
            sclet:::sclet_set_analysis_state(object, type = "enrichment",
                id = "enr_typed", method = "tg_mock",
                inputs = list(detest_id = "markers"),
                summary = list(n = 1))
        },
        returns = "sce",
        allowed_state_types = c("graph", "detest", "enrichment")
    )
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "chained")),
        tg_registry(chained = chained))
    graph <- AnalysisGraph(run$executed$object)
    consumed <- tg_edges(graph, "consumed")
    expect_equal(nrow(consumed), 1L)
    expect_identical(consumed$to, "record:detest::markers")
    expect_identical(consumed$declared_input, "detest_id")
})

test_that("provenance bookkeeping is not mistaken for a record reference", {
    noisy <- AIAction(
        name = "noisy",
        handler = function(object, params) {
            object <- sclet:::sclet_set_analysis_state(object, type = "preprocess",
                id = "report", method = "tg_mock", summary = list(n = 1))
            sclet:::sclet_set_analysis_state(object, type = "enrichment",
                id = "enr_meta", method = "tg_mock",
                inputs = list(task = "report", plan_id = "report",
                    context_fingerprint = "report"),
                summary = list(n = 1))
        },
        returns = "sce",
        allowed_state_types = c("preprocess", "enrichment")
    )
    run <- tg_run(tg_sce(), list(list(id = "s1", action = "noisy")),
        tg_registry(noisy = noisy))
    graph <- AnalysisGraph(run$executed$object)
    consumed <- tg_edges(graph, "consumed")
    expect_equal(nrow(consumed), 0L)
})

test_that("unresolvable and multiple evidence references are counted, not dropped", {
    sce <- tg_sce()
    sce <- sclet:::sclet_set_analysis_state(sce, type = "preprocess", id = "rec_a",
        method = "tg_mock", summary = list(n = 1))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "preprocess", id = "rec_b",
        method = "tg_mock", summary = list(n = 1))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "ai_evidence", id = "ev_dangling",
        method = "sclet_ai_evidence", active = FALSE,
        inputs = list(source = c("rec_a", "rec_b"), parents = "does_not_exist"),
        summary = list(id = "ev_dangling", kind = "test",
            source = c("rec_a", "rec_b"), parents = "does_not_exist",
            claim_level = "measured", status = "completed"))
    graph <- AnalysisGraph(sce)
    source_of <- tg_edges(graph, "source_of")
    # Both sources resolve, so both are drawn; the missing parent is counted.
    expect_setequal(source_of$to, c("record:preprocess::rec_a", "record:preprocess::rec_b"))
    expect_equal(graph$meta$unresolved_edges, 1L)
})

test_that("a malformed step store degrades to incomplete, not to an error", {
    sce <- tg_sce()
    malformed <- list(
        id = "ai_execution_malformed",
        type = "ai_execution",
        method = "ExecuteAIPlan",
        inputs = list(plan_id = "plan_malformed"),
        params = list(dry_run = FALSE, steps = "not a list"),
        artifacts = list(results = list()),
        summary = list(status = "completed"),
        created_at = Sys.time()
    )
    sce <- sclet:::sclet_set_analysis(sce, "ai_execution_malformed", malformed)
    graph <- expect_silent(AnalysisGraph(sce))
    round <- graph$nodes[graph$nodes$id == "round:ai_execution_malformed", , drop = FALSE]
    expect_equal(nrow(round), 1L)
    expect_false(isTRUE(round$complete))
    expect_equal(graph$meta$legacy_rounds, 1L)
})

test_that("malformed step elements are skipped, not fatal, and the round is incomplete", {
    sce <- tg_sce()
    mixed <- list(
        id = "ai_execution_mixed", type = "ai_execution", method = "ExecuteAIPlan",
        inputs = list(plan_id = "plan_mixed"),
        params = list(steps = list(list(id = "s1", action = "a", status = "completed",
            attempts = 1L, duration_sec = 0.1, params = list(),
            depends_on = character(), record_ids = character(),
            refreshed_ids = character()), "garbage")),
        artifacts = list(results = list()), summary = list(status = "completed"),
        created_at = Sys.time()
    )
    sce <- sclet:::sclet_set_analysis(sce, "ai_execution_mixed", mixed)
    graph <- expect_silent(AnalysisGraph(sce))
    steps <- graph$nodes[graph$nodes$kind == "step", , drop = FALSE]
    expect_equal(nrow(steps), 1L)
    round <- graph$nodes[graph$nodes$id == "round:ai_execution_mixed", , drop = FALSE]
    expect_false(isTRUE(round$complete))
})

test_that("evidence sources and parents are accounted for order-independently", {
    ev <- function(parents, source = "rec_b") {
        list(id = "ev_multi", kind = "test", source = source, parents = parents,
            claim_level = "measured", status = "completed")
    }
    sce <- tg_sce()
    sce <- sclet:::sclet_set_analysis_state(sce, type = "preprocess", id = "rec_a",
        method = "tg_mock", summary = list(n = 1))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "preprocess", id = "rec_b",
        method = "tg_mock", summary = list(n = 1))

    build <- function(node) {
        sclet:::sclet_set_analysis_state(sce, type = "ai_evidence", id = "ev_multi",
            method = "sclet_ai_evidence", active = FALSE,
            inputs = list(source = node$source, parents = node$parents),
            summary = node)
    }

    # "ghost" first must not hide the resolvable source.
    g1 <- AnalysisGraph(build(ev(parents = character(), source = c("ghost", "rec_b"))))
    expect_identical(tg_edges(g1, "source_of")$to, "record:preprocess::rec_b")
    expect_equal(g1$meta$unresolved_edges, 1L)

    # A repeated, unresolvable parent counts once.
    g2 <- AnalysisGraph(build(ev(parents = c("ghost", "ghost"))))
    expect_equal(g2$meta$unresolved_edges, 1L)

    # A list of parents is still read (not silently dropped).
    g3 <- AnalysisGraph(build(ev(parents = list("rec_a"))))
    expect_identical(tg_edges(g3, "derived_from")$to, "record:preprocess::rec_a")
    expect_equal(g3$meta$unresolved_edges, 0L)
})

test_that("a step with an empty id is not silently counted as recorded", {
    sce <- tg_sce()
    bad <- list(
        id = "ai_execution_noid", type = "ai_execution", method = "ExecuteAIPlan",
        inputs = list(plan_id = "plan_noid"),
        params = list(steps = list(list(status = "completed", action = "a",
            depends_on = character(), record_ids = character(),
            refreshed_ids = character()))),
        artifacts = list(results = list()), summary = list(status = "completed"),
        created_at = Sys.time()
    )
    sce <- sclet:::sclet_set_analysis(sce, "ai_execution_noid", bad)
    graph <- AnalysisGraph(sce)
    expect_equal(sum(graph$nodes$kind == "step"), 0L)
    round <- graph$nodes[graph$nodes$id == "round:ai_execution_noid", , drop = FALSE]
    expect_false(isTRUE(round$complete))
    expect_equal(graph$meta$legacy_rounds, 1L)
})
