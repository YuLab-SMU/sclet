sclet_ai_claim_test_node <- function(object, id, source = NULL, claim = "associated",
    dependency_group = "g1", parents = character(), kind = "deterministic_summary") {
    evidence <- list(
        id = id,
        kind = kind,
        values = list(flag = 1L),
        claim_level = claim,
        source = source,
        dependency_group = dependency_group,
        parents = parents
    )
    RecordAIEvidence(object, evidence, source = NULL, parents = parents, scope = NULL)
}

sclet_ai_claim_test_object <- function() {
    SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 3L, 4L)))
}

test_that("claim ceiling reports unsupported when no evidence exists", {
    sce <- sclet_ai_claim_test_object()
    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_equal(result$status, "unsupported")
    expect_null(result$ceiling_claim_level)
    expect_equal(result$n_independent_groups, 0L)
    expect_equal(result$reasons, "no evidence has been recorded on this object")
})

test_that("a single line of support is capped at associated and cannot be upgraded", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "consistent_with")
    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_equal(result$n_independent_groups, 1L)
    expect_equal(result$ceiling_claim_level, "associated")

    downgraded <- sclet:::sclet_ai_claim_ceiling(sce, proposed_claim_level = "consistent_with")
    expect_equal(downgraded$status, "downgraded")
    expect_true(any(grepl("exceeds the ceiling", downgraded$reasons)))
    expect_equal(sclet:::sclet_ai_claim_ceiling(sce, proposed_claim_level = "associated")$status, "allowed")
})

test_that("many nodes from one run count as a single line of support", {
    sce <- sclet_ai_claim_test_object()
    # same source, deliberately different dependency groups: dependency_group
    # alone would over-count these as independent
    for (i in 1:8) {
        sce <- sclet_ai_claim_test_node(sce, paste0("ev:n", i),
            source = "srcA", claim = "associated", dependency_group = paste0("g", i))
    }
    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_equal(result$n_nodes, 8L)
    expect_equal(result$n_independent_groups, 1L)
    expect_true(all(vapply(result$dependent_pairs, function(p) {
        "shared_source" %in% p$reasons
    }, logical(1L))))
    expect_equal(result$ceiling_claim_level, "associated")
})

test_that("two independent strong sources reach the consistent_with ceiling", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA",
        claim = "consistent_with", dependency_group = "g1")
    sce <- sclet_ai_claim_test_node(sce, "ev:b", source = "srcB",
        claim = "consistent_with", dependency_group = "g2")
    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_equal(result$n_independent_groups, 2L)
    expect_equal(result$ceiling_claim_level, "consistent_with")
    expect_equal(
        sclet:::sclet_ai_claim_ceiling(sce, proposed_claim_level = "consistent_with")$status,
        "allowed"
    )
    expect_setequal(vapply(result$groups, function(g) g$group, integer(1L)), c(1L, 2L))
})

test_that("the ceiling is bounded by the weakest supporting node", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:strong", source = "srcA",
        claim = "consistent_with", dependency_group = "g1")
    sce <- sclet_ai_claim_test_node(sce, "ev:weak", source = "srcB",
        claim = "observed", dependency_group = "g2")
    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_equal(result$n_independent_groups, 2L)
    expect_equal(result$ceiling_claim_level, "associated")
})

test_that("a shared dependency group merges nodes even across different sources", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA",
        claim = "consistent_with", dependency_group = "same")
    sce <- sclet_ai_claim_test_node(sce, "ev:b", source = "srcB",
        claim = "consistent_with", dependency_group = "same")
    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_equal(result$n_independent_groups, 1L)
    expect_true(any(vapply(result$dependent_pairs, function(p) {
        "shared_dependency_group" %in% p$reasons
    }, logical(1L))))
})

test_that("user_decision evidence never supports an AI-generated claim", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:human", source = NULL,
        claim = "observed", dependency_group = "h1", kind = "user_decision")
    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_equal(result$status, "unsupported")
    expect_equal(result$n_independent_groups, 0L)
    expect_true(any(grepl("human decision cannot support", result$reasons)))

    # and it is ignored, not counted, when real evidence is also present
    sce2 <- sclet_ai_claim_test_node(sce, "ev:ai", source = "srcA",
        claim = "associated", dependency_group = "g1")
    mixed <- sclet:::sclet_ai_claim_ceiling(sce2)
    expect_equal(mixed$n_independent_groups, 1L)
    expect_equal(mixed$user_decision_nodes, "ev:human")
    expect_true(any(grepl("user_decision nodes were excluded", mixed$reasons)))
})

test_that("claim ceiling is read-only and validates its inputs", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "associated")
    before_nodes <- length(sclet:::sclet_ai_evidence_get_all(sce))
    before_fp <- GetAnalysisLedger(sce)$fingerprint

    result <- sclet:::sclet_ai_claim_ceiling(sce)
    expect_false(result$raw_values_included)

    # auditing adds no evidence and changes nothing
    expect_equal(length(sclet:::sclet_ai_evidence_get_all(sce)), before_nodes)
    expect_equal(GetAnalysisLedger(sce)$fingerprint, before_fp)

    expect_error(sclet:::sclet_ai_claim_ceiling(sce, evidence_ids = "ev:nope"), "unknown or incomplete")
    expect_error(
        sclet:::sclet_ai_claim_ceiling(sce, proposed_claim_level = "certain"),
        "proposed_claim_level must be one of"
    )
    expect_error(sclet:::sclet_ai_claim_ceiling(list(a = 1)), "SingleCellExperiment")
})

test_that("evidence from a real two-source pipeline yields two independent lines", {
    skip_if_not_installed("BiocNeighbors")
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(rpois(50 * 40, 5),
        nrow = 50L, ncol = 40L, dimnames = list(paste0("g", seq_len(50L)), paste0("c", seq_len(40L))))))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 1.2)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation", "rare_cell"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "evidence_chain",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "d1", action = "run_de_test", params = list(name = "det1")),
            list(id = "r1", action = "run_rare_cell_detection",
                params = list(name = "rare1", k = 5L, dims = c(1L, 2L, 3L), rare_threshold = 1000L))
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)

    ids <- names(sclet:::sclet_ai_evidence_get_all(executed$object))
    expect_true(length(ids) > 2L)
    result <- sclet:::sclet_ai_claim_ceiling(executed$object)
    # many nodes, but only two underlying analyses produced them
    expect_true(result$n_nodes > 2L)
    expect_equal(result$n_independent_groups, 2L)
    expect_setequal(
        unlist(lapply(result$groups, function(g) g$sources)),
        c("det1", "rare1")
    )
})