# Tests for engine-attached evidence authorship, per
# .dev/spec-ai-evidence-authorship.md.
#
# The problem this guards: `sclet_ai_derive_finding_claim_level()` is correct,
# and with no `evidence_refs` on a finding there is nothing for the ceiling to
# reason about, so every finding lands at `hypothesis` through the "no refs"
# path. Measured on pbmc3k (sclet-benchmark-ai/results/derived_claim_levels.jsonl,
# 2026-10-13): the two `with_evidence` conditions had derivation_available=TRUE,
# the object carried ev:findallmarkers and the ceiling was `allowed`, and still
# 0 of 23 and 0 of 13 findings moved, because none cited anything.
#
# Authorship closes that by attaching the evidence the engine knows the step
# created. The whole risk is that attaching RAISES a claim, so these tests are
# mostly about when it must REFUSE: two candidates is a guess, the model's own
# citation is authoritative, and `causal`/`estimated` must never gain a ref.

authorship_fixture <- function() {
    set.seed(1)
    g <- 40L
    n <- 20L
    m <- matrix(rpois(g * n, 5), nrow = g)
    rownames(m) <- paste0("g", 1:g)
    colnames(m) <- paste0("c", 1:n)
    SingleCellExperiment::SingleCellExperiment(list(counts = m))
}

add_evidence <- function(object, id, claim_level = "associated", kind = "deterministic_summary") {
    RecordAIEvidence(
        object,
        list(id = id, kind = kind, claim_level = claim_level, values = list(n = 1L))
    )
}

authorship_finding <- function(claim_level = "hypothesis", evidence_refs = NULL) {
    finding <- list(
        id = "f1",
        statement = "the clusters are separated",
        severity = "info",
        claim_level = claim_level
    )
    if (!is.null(evidence_refs)) finding$evidence_refs <- evidence_refs
    finding
}

test_that("a sole evidence node is attached and raises the finding", {
    sce <- add_evidence(authorship_fixture(), "ev:only")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, authorship_finding())

    expect_true(attachment$attached)
    expect_identical(attachment$rule, "sole_evidence_node")
    expect_identical(attachment$finding$evidence_refs, "ev:only")

    # The derivation must now consult the ceiling instead of taking the
    # "no refs" floor -- this is the whole point of the change.
    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, attachment$finding)
    expect_identical(derivation$claim_level, "associated")
    expect_identical(derivation$ceiling_status, "allowed")
    expect_true(derivation$derived_from_ceiling)
})

test_that("two candidates is a guess, so nothing is attached", {
    # The anti-guessing criterion, and the most important assertion in this
    # file. An attachment raises a claim, so choosing between two nodes could
    # report an unsupported claim as clean -- strictly worse than the honest
    # floor. Under-attaching is the safe direction.
    sce <- authorship_fixture()
    sce <- add_evidence(sce, "ev:a", claim_level = "associated")
    sce <- add_evidence(sce, "ev:b", claim_level = "observed")

    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, authorship_finding())
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "ambiguous_multiple_candidates")
    expect_null(attachment$finding$evidence_refs)

    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, attachment$finding)
    expect_identical(derivation$claim_level, "hypothesis")
    expect_identical(derivation$ceiling_status, "no_refs")
})

test_that("the model's own citation is never overwritten", {
    sce <- authorship_fixture()
    sce <- add_evidence(sce, "ev:only")
    sce <- add_evidence(sce, "ev:other")

    # Even with a single node the engine could attach, an existing citation
    # stands: the model is authoritative for what it says it cited.
    single <- add_evidence(authorship_fixture(), "ev:only")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(
        single, authorship_finding(evidence_refs = "ev:other")
    )
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "model_cited")
    expect_identical(attachment$finding$evidence_refs, "ev:other")

    multi <- sclet:::sclet_ai_attach_finding_evidence(
        sce, authorship_finding(evidence_refs = "ev:b")
    )
    expect_false(multi$attached)
    expect_identical(multi$finding$evidence_refs, "ev:b")
})

test_that("causal and estimated never gain a reference", {
    # `causal` is banned, and `validate_sclet_ai_result()` demands a reference
    # for a causal claim -- so attaching one would turn a result the engine
    # rejects into one it accepts, defeating the ban by helpfulness.
    # `estimated` is a provenance statement, orthogonal to the evidence ladder:
    # the derivation returns that verdict without consulting evidence at all.
    sce <- add_evidence(authorship_fixture(), "ev:only")

    for (level in c("causal", "estimated")) {
        attachment <- sclet:::sclet_ai_attach_finding_evidence(
            sce, authorship_finding(claim_level = level)
        )
        expect_false(attachment$attached, info = level)
        expect_identical(attachment$rule, "requested_level_excluded", info = level)
        expect_null(attachment$finding$evidence_refs, info = level)
    }

    causal <- sclet:::sclet_ai_derive_finding_claim_level(
        sce, authorship_finding(claim_level = "causal")
    )
    expect_true(causal$refused)
    expect_identical(causal$ceiling_status, "banned")

    estimated <- sclet:::sclet_ai_derive_finding_claim_level(
        sce, authorship_finding(claim_level = "estimated")
    )
    expect_identical(estimated$claim_level, "estimated")
    expect_identical(estimated$ceiling_status, "provenance")
})

test_that("a user decision is never treated as support", {
    sce <- add_evidence(
        authorship_fixture(), "ev:decision",
        kind = "user_decision"
    )
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, authorship_finding())
    expect_false(attachment$attached)
    expect_identical(attachment$rule, "no_candidates")

    # A real node alongside a decision is still a single candidate.
    sce2 <- add_evidence(sce, "ev:real")
    attachment2 <- sclet:::sclet_ai_attach_finding_evidence(sce2, authorship_finding())
    expect_true(attachment2$attached)
    expect_identical(attachment2$finding$evidence_refs, "ev:real")
})

test_that("an object with no evidence attaches nothing", {
    attachment <- sclet:::sclet_ai_attach_finding_evidence(
        authorship_fixture(), authorship_finding()
    )
    expect_false(attachment$attached)
    # "none" (no evidence at all on the object) is distinct from
    # "no_candidates" (evidence exists, but none of it can support a finding,
    # e.g. the only node is a user decision).
    expect_identical(attachment$rule, "none")
})

test_that("the audit reports the attachment rather than laundering it", {
    # A reader must be able to tell a reference the model cited from one the
    # engine attached, because only the former is evidence the model reasoned
    # about. The rule name travels with the entry.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    result <- new_sclet_ai_result(
        task = "qc_review", answer = "a",
        findings = list(authorship_finding())
    )
    audit <- sclet_ai_audit_result_claims(sce, result)
    entry <- audit$findings[[1L]]
    expect_true(entry$evidence_refs_attached)
    expect_identical(entry$attachment_rule, "sole_evidence_node")
    expect_identical(entry$claim_level, "associated")

    # A finding that cited for itself reports no attachment.
    cited <- new_sclet_ai_result(
        task = "qc_review", answer = "a",
        findings = list(authorship_finding(evidence_refs = "ev:only"))
    )
    cited_entry <- sclet_ai_audit_result_claims(sce, cited)$findings[[1L]]
    expect_false(cited_entry$evidence_refs_attached)
    expect_identical(cited_entry$attachment_rule, "model_cited")
})

test_that("the recorded finding carries the attached reference and its level", {
    # The reference has to travel with the recorded finding: a raised
    # `claim_level` whose stated basis is absent from the ledger would be
    # exactly the unsupported claim the audit exists to prevent.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    result <- new_sclet_ai_result(
        task = "qc_review", answer = "a",
        findings = list(authorship_finding())
    )
    recorded <- RecordAIResult(sce, result, id = "rec1")
    finding <- sclet:::sclet_get_state(recorded)$analyses[["rec1"]]$artifacts$findings[[1L]]

    expect_identical(finding$claim_level, "associated")
    expect_identical(finding$requested_claim_level, "hypothesis")
    expect_identical(as.character(finding$evidence_refs), "ev:only")
    expect_true(finding$evidence_refs_attached)
    expect_true(finding$derived_from_ceiling)
})

test_that("an over-claimed causal result is still refused", {
    # The ban must survive authorship, including when the model supplied the
    # reference itself.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    result <- new_sclet_ai_result(
        task = "qc_review", answer = "a",
        findings = list(authorship_finding(claim_level = "causal", evidence_refs = "ev:only"))
    )
    expect_error(
        RecordAIResult(sce, result, id = "rec_causal"),
        "over-claimed"
    )
})

test_that("a result that already cited its evidence is unchanged by authorship", {
    # The change must not alter any verdict on a result that cited for itself.
    sce <- add_evidence(authorship_fixture(), "ev:only")
    finding <- authorship_finding(claim_level = "consistent_with", evidence_refs = "ev:only")
    result <- new_sclet_ai_result(task = "qc_review", answer = "a", findings = list(finding))

    audit <- sclet_ai_audit_result_claims(sce, result)
    direct <- sclet:::sclet_ai_derive_finding_claim_level(sce, finding)

    expect_identical(audit$findings[[1L]]$claim_level, direct$claim_level)
    expect_false(audit$findings[[1L]]$evidence_refs_attached)
    # One node at "associated" means "consistent_with" is above the ceiling.
    expect_identical(audit$findings[[1L]]$claim_level, "associated")
    expect_identical(audit$findings[[1L]]$status, "downgraded")
})

test_that("the attached flag is never mistaken for an evidence id", {
    # Regression guard for R's `$` partial matching. A recorded finding carries
    # BOTH `evidence_refs` (character) and `evidence_refs_attached` (logical),
    # and `finding$evidence_refs` partial-matches the latter when the former is
    # absent -- yielding `FALSE`, which `as.character()` turns into the literal
    # id "FALSE". That id then failed `ValidateAIEvidenceRefs()`, so every
    # result whose findings cited nothing was REJECTED at record time, and the
    # derivation consulted the ceiling with a bogus id for the same reason.
    # `[[` does not partial-match; these assertions pin that down.
    recorded <- list(list(
        statement = "s", severity = "info", claim_level = "hypothesis",
        evidence_refs_attached = FALSE
    ))
    expect_identical(sclet:::sclet_ai_finding_refs(recorded), character())

    sce <- add_evidence(authorship_fixture(), "ev:only")
    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, recorded[[1L]])
    # "no_refs" is the honest verdict; "unsupported" would mean the derivation
    # tried to resolve a phantom id.
    expect_identical(derivation$ceiling_status, "no_refs")

    # And a finding that cites nothing at all still records successfully.
    bare <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "no refs", severity = "info", claim_level = "hypothesis"
        ))
    )
    recorded_object <- RecordAIResult(bare, result, id = "ai_no_refs")
    expect_true(!is.null(
        sclet:::sclet_get_state(recorded_object)$analyses[["ai_no_refs"]]
    ))
})

test_that("authorship does not raise a claim beyond what the node declares", {
    # The ceiling is the minimum over cited nodes' own claim levels, so
    # authorship cannot manufacture strength -- it can only surface strength a
    # recorded node already carries. A node at `observed` must not lift a
    # finding to `associated`.
    sce <- add_evidence(authorship_fixture(), "ev:weak", claim_level = "observed")
    attachment <- sclet:::sclet_ai_attach_finding_evidence(sce, authorship_finding())
    expect_true(attachment$attached)
    derivation <- sclet:::sclet_ai_derive_finding_claim_level(sce, attachment$finding)
    expect_identical(derivation$claim_level, "observed")
})
