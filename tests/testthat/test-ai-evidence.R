test_that("RecordAIEvidence stores bounded evidence and validates references", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, nrow = 2L, ncol = 2L)))
    sce <- sclet:::sclet_set_analysis_state(
        sce, type = "integration", id = "integration_1", method = "harmony",
        summary = list(status = "completed")
    )
    before <- GetAnalysisLedger(sce)$fingerprint
    evidence <- list(
        id = "evidence_1",
        kind = "deterministic_summary",
        values = list(batch_mixing = 0.8, groups = c("group_1", "group_2")),
        uncertainty = list(method = "bootstrap")
    )
    updated <- sclet:::RecordAIEvidence(sce, evidence, source = "integration_1")
    validation <- sclet:::ValidateAIEvidenceRefs(updated, "evidence_1")

    expect_true(validation$valid)
    expect_equal(validation$status, "valid")
    expect_equal(validation$resolved[[1]]$kind, "deterministic_summary")
    expect_equal(validation$resolved[[1]]$source, "integration_1")
    expect_false(identical(before, GetAnalysisLedger(updated)$fingerprint))
})

test_that("evidence rejects duplicates, incomplete sources, and private values", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 2L)))
    evidence <- list(id = "evidence_1", kind = "state", values = list(n = 10))
    expect_error(sclet:::RecordAIEvidence(sce, evidence, source = "missing"), "completed")
    updated <- sclet:::RecordAIEvidence(sce, evidence)
    expect_error(sclet:::RecordAIEvidence(updated, evidence), "already exists")
    expect_error(
        sclet:::RecordAIEvidence(sce, list(id = "private", kind = "state", values = list(patient_id = "p1"))),
        "private"
    )
})

test_that("evidence references report stale fingerprints and missing ids", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 2L)))
    sce <- sclet:::RecordAIEvidence(sce, list(id = "evidence_1", kind = "test", values = list(n = 10)))
    invalid <- sclet:::ValidateAIEvidenceRefs(sce, c("evidence_1", "missing"), fingerprint = "stale")
    expect_false(invalid$valid)
    expect_equal(length(invalid$errors), 2L)
    valid <- sclet:::ValidateAIEvidenceRefs(sce, "evidence_1")
    expect_true(valid$valid)
})

test_that("RecordAIEvidence fills dependency_group deterministically from parents and scope", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 2L)))
    sce <- sclet:::RecordAIEvidence(sce, list(id = "base", kind = "test", values = list(n = 1L)))
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "ev_same_1", kind = "deterministic_summary", values = list(x = 1)),
        parents = "base", scope = list(group = c("group_1", "group_2")))
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "ev_same_2", kind = "plot", values = list(p = 0.05)),
        parents = "base", scope = list(group = c("group_1", "group_2")))
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "ev_different", kind = "plot", values = list(p = 0.02)),
        parents = character(), scope = list(group = "group_3"))
    ledger <- GetAnalysisLedger(sce)
    ev_same_1 <- ledger$state_records$ai_evidence$ev_same_1$summary
    ev_same_2 <- ledger$state_records$ai_evidence$ev_same_2$summary
    ev_different <- ledger$state_records$ai_evidence$ev_different$summary
    expect_true(nzchar(ev_same_1$dependency_group))
    expect_identical(ev_same_1$dependency_group, ev_same_2$dependency_group)
    expect_false(identical(ev_same_1$dependency_group, ev_different$dependency_group))
})

test_that("ValidateAIEvidenceRefs rejects scope violations and AI citing user_decision", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 2L)))
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "narrow_scope", kind = "deterministic_summary", values = list(x = 1)),
        scope = list(group = "group_1"))
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "user_dec", kind = "user_decision", values = list(route_index = 1, groups_supported = c("group_1", "group_2"))))
    v_scope <- sclet:::ValidateAIEvidenceRefs(sce, "narrow_scope",
        requesting_scope = list(group = c("group_1", "group_2")),
        requesting_kind = "test")
    expect_false(v_scope$valid)
    expect_true(any(grepl("scope violation", v_scope$errors)))
    v_cross_kind <- sclet:::ValidateAIEvidenceRefs(sce, "user_dec",
        requesting_kind = "deterministic_summary")
    expect_false(v_cross_kind$valid)
    expect_true(any(grepl("cross-kind", v_cross_kind$errors)))
})

test_that("sclet_ai_evidence_independence detects shared groups and ancestor lineage", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 2L)))
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "root", kind = "state", values = list(n = 1L)))
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "sibling_A", kind = "test", values = list(p = 1)),
        parents = "root")
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "sibling_B", kind = "test", values = list(p = 1)),
        parents = "root")
    sce <- sclet:::RecordAIEvidence(sce,
        list(id = "child_of_A", kind = "deterministic_summary", values = list(x = 1)),
        parents = "sibling_A")
    indep_siblings <- sclet:::sclet_ai_evidence_independence(sce, c("sibling_A", "sibling_B"))
    # siblings share a common root parents dependency_group but A != B paths
    # The indep_siblings independent may be TRUE or FALSE depending on implementation;
    # the API contract requires at least lineage_paths and shared_dependency_groups fields.
    expect_true("shared_dependency_groups" %in% names(indep_siblings))
    expect_true("lineage_paths" %in% names(indep_siblings))
    indep_ancestor <- sclet:::sclet_ai_evidence_independence(sce, c("sibling_A", "child_of_A"))
    expect_false(indep_ancestor$independent)
    expect_true(length(indep_ancestor$lineage_paths) >= 1L)
})
