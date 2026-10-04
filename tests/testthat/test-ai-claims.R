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

test_that("a causal finding is rejected regardless of how much evidence it cites", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "consistent_with")
    sce <- sclet_ai_claim_test_node(sce, "ev:b", source = "srcB", claim = "consistent_with")
    result <- list(findings = list(list(
        claim_level = "causal", evidence_refs = c("ev:a", "ev:b")
    )))
    audit <- sclet:::sclet_ai_audit_result_claims(sce, result)
    expect_equal(audit$status, "overclaimed")
    expect_true(any(grepl("not a claim level the AI may assert", audit$problems)))
    expect_equal(audit$findings[[1L]]$allowed_claim_level, "hypothesis")
})

test_that("a consistent_with finding needs evidence that reaches the ceiling", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "consistent_with")
    # only one line of support, so the ceiling is capped at associated
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:a"))
    ))
    expect_equal(audit$status, "overclaimed")
    expect_true(any(grepl("exceeds the evidence ceiling 'associated'", audit$problems)))

    # add a second independent strong source and the same finding is allowed
    sce2 <- sclet_ai_claim_test_node(sce, "ev:b", source = "srcB",
        claim = "consistent_with", dependency_group = "g2")
    audit2 <- sclet:::sclet_ai_audit_result_claims(sce2, list(
        findings = list(list(claim_level = "consistent_with", evidence_refs = c("ev:a", "ev:b")))
    ))
    expect_equal(audit2$status, "ok")
    expect_equal(audit2$n_problems, 0L)
})

test_that("suggestive is held to the consistent_with bar and hypothesis needs nothing", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "associated")
    suggestive <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "suggestive", evidence_refs = "ev:a"))
    ))
    expect_equal(suggestive$status, "overclaimed")
    expect_equal(suggestive$findings[[1L]]$allowed_claim_level, "associated")

    # a hypothesis is explicitly speculative and asserts nothing
    hypothesis <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "hypothesis"))
    ))
    expect_equal(hypothesis$status, "ok")
})

test_that("a consistent_with finding with no evidence_refs is rejected", {
    sce <- sclet_ai_claim_test_object()
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with"))
    ))
    expect_equal(audit$status, "overclaimed")
    expect_true(any(grepl("cites no evidence_refs", audit$problems)))
})

test_that("record_claims auditing refuses to write an over-claimed result to the ledger", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "associated")
    overclaimed <- new_sclet_ai_result(
        "review",
        answer = "something",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "causal", evidence_refs = "ev:a"))
    )
    expect_error(
        RecordAIResult(sce, overclaimed, id = "ai_over"),
        "refusing to record an over-claimed AI result"
    )
    # nothing was written
    expect_false("ai_over" %in% names(GetAnalysisLedger(sce)$analyses))

    # the audit can be turned off explicitly, and then it records
    updated <- RecordAIResult(sce, overclaimed, id = "ai_raw", audit_claims = FALSE)
    expect_true("ai_raw" %in% names(GetAnalysisLedger(updated)$analyses))
})

test_that("a supported result still records normally with auditing enabled", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "consistent_with")
    sce <- sclet_ai_claim_test_node(sce, "ev:b", source = "srcB",
        claim = "consistent_with", dependency_group = "g2")
    supported <- new_sclet_ai_result(
        "review",
        answer = "something",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "consistent_with", evidence_refs = c("ev:a", "ev:b")))
    )
    updated <- RecordAIResult(sce, supported, id = "ai_ok")
    expect_true("ai_ok" %in% names(GetAnalysisLedger(updated)$analyses))
})

test_that("result claim auditing is read-only", {
    sce <- sclet_ai_claim_test_object()
    sce <- sclet_ai_claim_test_node(sce, "ev:a", source = "srcA", claim = "associated")
    before_fp <- GetAnalysisLedger(sce)$fingerprint
    before_n <- length(sclet:::sclet_ai_evidence_get_all(sce))
    result <- list(findings = list(list(claim_level = "causal", evidence_refs = "ev:a")))
    audit <- sclet:::sclet_ai_audit_result_claims(sce, result)
    expect_equal(audit$status, "overclaimed")
    expect_equal(GetAnalysisLedger(sce)$fingerprint, before_fp)
    expect_equal(length(sclet:::sclet_ai_evidence_get_all(sce)), before_n)
    # the audit reports a ceiling but never rewrites the finding it inspected
    expect_equal(result$findings[[1L]]$claim_level, "causal")
})

test_that("results without findings are unaffected by claim auditing", {
    sce <- sclet_ai_claim_test_object()
    expect_equal(sclet:::sclet_ai_audit_result_claims(sce, list())$status, "ok")
    plain <- new_sclet_ai_result("status_review",
        answer = "all inputs are visible",
        context = list(schema_version = "1.0", fingerprint = "test"),
        metadata = list(provider = "mock"))
    updated <- RecordAIResult(sce, plain, id = "ai_plain")
    expect_true("ai_plain" %in% names(GetAnalysisLedger(updated)$analyses))
})

test_that("a finding citing an unknown evidence ref is reported, not crashed on", {
    sce <- sclet_ai_claim_test_object()
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:does_not_exist"))
    ))
    expect_equal(audit$status, "overclaimed")
    expect_true(any(grepl("is not supported by the cited evidence", audit$problems)))
    expect_equal(audit$findings[[1L]]$allowed_claim_level, "observed")
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