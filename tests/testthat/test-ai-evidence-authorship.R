# Tests for engine evidence authorship (v2), per
# .dev/spec-ai-evidence-authorship-v2.md.
#
# v2 replaces v1's "attach the sole node to any cite-nothing finding" rule,
# which was measured defective: it attached `ev:findallmarkers` to findings
# asserting that evidence was ABSENT and raised them off the floor, with a clean
# audit (`.dev/analysis-authorship-attachment-risk.md`). v2 never attributes --
# it transcribes: attach only when the finding's own prose names exactly one
# evidence-node id as a whole token. That makes the multi-node object the
# POSITIVE case (selection comes from the model's words) instead of a forced
# decline.
#
# The v1 contract this file replaced asserted the attach-and-raise behaviour for
# cite-nothing, name-nothing findings. Those are now the floor cases and are
# titled so, because a future reader must not mistake the new floor for a
# regression of the old raise.

authorship_fixture <- function(seed = 1L) {
    set.seed(seed)
    g <- 40L
    n <- 20L
    m <- matrix(rpois(g * n, 5), nrow = g)
    rownames(m) <- paste0("g", 1:g)
    colnames(m) <- paste0("c", 1:n)
    SingleCellExperiment::SingleCellExperiment(list(counts = m))
}

add_evidence <- function(object, id, claim_level = "associated",
                         kind = "deterministic_summary", scope = NULL) {
    RecordAIEvidence(object, list(id = id, kind = kind, claim_level = claim_level,
        values = list(n = 1L)), scope = scope)
}

authorship_finding <- function(claim_level = "hypothesis", evidence_refs = NULL,
                               statement = "the clusters are separated") {
    finding <- list(id = "f1", statement = statement, severity = "info",
        claim_level = claim_level)
    if (!is.null(evidence_refs)) finding$evidence_refs <- evidence_refs
    finding
}

# A finding that names the node in prose is the only thing v2 transcribes.
named_finding <- function(id, claim_level = "hypothesis") {
    authorship_finding(claim_level = claim_level,
        statement = paste("as recorded by", id, "the result is stable"))
}

# A finding asserting evidence is MISSING -- the v1 false-attachment case (C5).
absence_finding <- function(claim_level = "hypothesis") {
    authorship_finding(claim_level = claim_level,
        statement = "no doublet detection or mitochondrial check was recorded")
}

test_that("a cite-nothing finding that NAMES the sole node is transcribed and raised", {
    sce <- add_evidence(authorship_fixture(), "ev:only")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, named_finding("ev:only"))
    expect_true(attachment$attached)
    expect_identical(attachment$rule, "prose_named_unique")
    expect_identical(attachment$finding$evidence_refs, "ev:only")
    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, attachment$finding)
    expect_identical(derivation$claim_level, "associated")
    expect_true(derivation$derived_from_ceiling)
})

test_that("a cite-nothing finding that names NOTHING keeps the floor (v1 C5 defect gone)", {
    # Criterion 1 + criterion 8: the exact case v1 wrongly attached.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, absence_finding())
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "no_prose_reference")
    expect_null(attachment$finding$evidence_refs)
    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, attachment$finding)
    expect_identical(derivation$claim_level, "hypothesis")
    expect_identical(derivation$ceiling_status, "no_refs")
})

test_that("the whole-token scan does not match a prefix inside a longer id", {
    expect_true(sclet:::sclet_ai_prose_names_token("see ev:markers here", "ev:markers"))
    expect_false(sclet:::sclet_ai_prose_names_token("see ev:markers_extra", "ev:markers"))
    expect_false(sclet:::sclet_ai_prose_names_token("xev:markers not this", "ev:markers"))
    expect_true(sclet:::sclet_ai_prose_names_token("node (ev:markers), end", "ev:markers"))
    expect_true(sclet:::sclet_ai_prose_names_token("ev:markers", "ev:markers"))
    expect_false(sclet:::sclet_ai_prose_names_token("", "ev:markers"))
})

test_that("multi-node: a finding naming exactly ONE node selects it", {
    # v1 refused whenever >= 2 nodes existed; v2 selects from the model's own
    # words, so a two-node object is now the positive case.
    sce <- authorship_fixture()
    sce <- add_evidence(sce, "ev:a", claim_level = "associated", scope = list(sample = "group_1"))
    sce <- add_evidence(sce, "ev:b", claim_level = "measured", scope = list(sample = "group_2"))
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, named_finding("ev:b"))
    expect_true(attachment$attached)
    expect_identical(attachment$finding$evidence_refs, "ev:b")
    # The ceiling is computed from the SELECTED node (measured), not from the
    # pair -- citing one node cannot be attenuated by an unrelated second one.
    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, attachment$finding)
    expect_identical(derivation$claim_level, "measured")
})

test_that("multi-node: a finding naming BOTH nodes is ambiguous and declined", {
    sce <- authorship_fixture()
    sce <- add_evidence(sce, "ev:a", scope = list(sample = "group_1"))
    sce <- add_evidence(sce, "ev:b", scope = list(sample = "group_2"))
    f <- authorship_finding(
        statement = "cross-referencing ev:a and ev:b together implies a trajectory")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, f)
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "ambiguous_prose_reference")
})

test_that("the model's own citation is never overwritten", {
    sce <- add_evidence(authorship_fixture(), "ev:only")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(
        sce, authorship_finding(evidence_refs = "ev:other",
            statement = "ev:only corroborates this"))
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "model_cited")
    expect_identical(attachment$finding$evidence_refs, "ev:other")
})

test_that("causal and estimated never gain a reference, even naming it in prose", {
    sce <- add_evidence(authorship_fixture(), "ev:only")
    for (level in c("causal", "estimated")) {
        attachment <- sclet:::sclet_ai_attach_finding_evidence(
            sce, named_finding("ev:only", claim_level = level))
        expect_false(attachment$attached, info = level)
        expect_identical(attachment$rule, "requested_level_excluded", info = level)
        expect_null(attachment$finding$evidence_refs, info = level)
    }
    causal <- sclet:::sclet_ai_derive_finding_claim_level(
        sce, authorship_finding(claim_level = "causal"))
    expect_true(causal$refused)
    expect_identical(causal$ceiling_status, "banned")
    # estimated stays provenance and is not clamped or raised.
    est <- sclet:::sclet_ai_derive_finding_claim_level(
        sce, authorship_finding(claim_level = "estimated", evidence_refs = "ev:only"))
    expect_identical(est$claim_level, "estimated")
    expect_identical(est$ceiling_status, "provenance")
})

test_that("a user decision is never treated as support", {
    sce <- add_evidence(authorship_fixture(), "ev:decision", kind = "user_decision")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, named_finding("ev:decision"))
    expect_false(attachment$attached)
    expect_true(attachment$rule %in% c("none", "no_prose_reference"))
})

test_that("an object with no evidence attaches nothing", {
    attachment <- sclet:::sclet_ai_attach_finding_evidence(
        authorship_fixture(), named_finding("ev:none-present"))
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "none")
})

test_that("the derivation CLAMPS: an observed request is not raised to associated", {
    # C1 / criterion 5. v1 overwrote this into a relational claim.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    audit <- sclet_ai_audit_result_claims(sce, list(findings = list(
        named_finding("ev:only", claim_level = "observed"))))
    entry <- audit$findings[[1L]]
    expect_identical(entry$claim_level, "observed")
    expect_identical(entry$status, "ok")
    expect_equal(audit$n_raised, 0L)
    expect_equal(audit$n_lowered, 0L)
})

test_that("the derivation still RAISES from the floor when a ref is present", {
    # Criterion 6: min() with no floor exemption would make attachment inert
    # (the 0-of-23 result). hypothesis ranks 0, so it takes the ceiling.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    audit <- sclet_ai_audit_result_claims(sce, list(findings = list(
        authorship_finding(claim_level = "hypothesis", evidence_refs = "ev:only"))))
    entry <- audit$findings[[1L]]
    expect_identical(entry$claim_level, "associated")
    expect_identical(entry$status, "raised")
    expect_equal(audit$n_raised, 1L)
    expect_equal(audit$status, "raised")
})

test_that("a downward move reports 'corrected' and the two counts are split", {
    # Criterion 3 / vocabulary: consistent_with above an associated ceiling is
    # a lowering, distinct from an upgrade.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    audit <- sclet_ai_audit_result_claims(sce, list(findings = list(
        authorship_finding(claim_level = "consistent_with", evidence_refs = "ev:only"))))
    expect_identical(audit$findings[[1L]]$status, "corrected")
    expect_equal(audit$n_lowered, 1L)
    expect_equal(audit$n_raised, 0L)
    # n_downgraded kept as an alias so existing readers see a number, not NULL.
    expect_equal(audit$n_downgraded, 1L)
})

test_that("the audit reports attachment without laundering it", {
    sce <- add_evidence(authorship_fixture(), "ev:only")
    result <- new_sclet_ai_result(task = "qc_review", answer = "a",
        findings = list(named_finding("ev:only"), absence_finding()))
    audit <- sclet_ai_audit_result_claims(sce, result)
    expect_true(audit$findings[[1L]]$evidence_refs_attached)
    expect_identical(audit$findings[[1L]]$attachment_rule, "prose_named_unique")
    expect_false(audit$findings[[2L]]$evidence_refs_attached)
    expect_identical(audit$findings[[2L]]$attachment_rule, "no_prose_reference")
})

test_that("the recorded finding carries the attached reference, its level and its rule", {
    sce <- add_evidence(authorship_fixture(), "ev:only")
    result <- new_sclet_ai_result(task = "qc_review", answer = "a",
        findings = list(named_finding("ev:only")))
    recorded <- RecordAIResult(sce, result, id = "rec1")
    finding <- sclet:::sclet_get_state(recorded)$analyses[["rec1"]]$artifacts$findings[[1L]]
    expect_identical(finding$claim_level, "associated")
    expect_identical(finding$requested_claim_level, "hypothesis")
    expect_identical(as.character(finding$evidence_refs), "ev:only")
    expect_true(finding$evidence_refs_attached)
    expect_identical(finding$evidence_attachment_rule, "prose_named_unique")
})

test_that("an over-claimed causal result is still refused", {
    sce <- add_evidence(authorship_fixture(), "ev:only")
    result <- new_sclet_ai_result(task = "qc_review", answer = "a",
        findings = list(authorship_finding(claim_level = "causal", evidence_refs = "ev:only")))
    expect_error(RecordAIResult(sce, result, id = "r"), "over-claimed")
})

test_that("attachment never lifts a finding above the node's own level", {
    # The ceiling is the min over cited nodes' own levels: an `observed` node
    # cannot carry a finding to `associated` even when transcribed.
    sce <- add_evidence(authorship_fixture(), "ev:weak", claim_level = "observed")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, named_finding("ev:weak"))
    expect_true(attachment$attached)
    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, attachment$finding)
    expect_identical(derivation$claim_level, "observed")
})

test_that("the v1 sole-node fallback is gone (regression guard)", {
    # The single most load-bearing assertion here: a cite-nothing,
    # name-nothing finding on a one-node object must keep the floor. If the
    # sole-node rule is ever re-added, this fails.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, authorship_finding())
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "no_prose_reference")
})
