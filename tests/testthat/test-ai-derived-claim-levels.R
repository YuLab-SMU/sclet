# Tests for `spec-ai-derived-claim-levels.md` -- derive finding claim levels
# from evidence instead of asking the model. Each test corresponds to a row of
# the spec's acceptance table, plus the causal invariant and the "same model
# result, different evidence" property.

# A minimal SCE for tests that do not need real evidence.
sclet_ai_derived_test_object <- function() {
    SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 3L, 4L)))
}

sclet_ai_derived_test_node <- function(object, id, source = NULL, claim = "associated",
    dependency_group = "g1", parents = character(), kind = "deterministic_summary") {
    evidence <- list(
        id = id, kind = kind, values = list(flag = 1L),
        claim_level = claim, source = source,
        dependency_group = dependency_group, parents = parents
    )
    RecordAIEvidence(object, evidence, source = NULL, parents = parents, scope = NULL)
}

# An SCE carrying the kind of evidence `ExecuteAIPlan()` would normally
# produce. Uses a small local fixture (the `pbmc3k` benchmark fixture is not
# in this repo, but the test needs an object with at least one completed
# evidence node so the ceiling can be computed).
sclet_ai_derived_test_object_with_evidence <- function() {
    sce <- sclet_ai_derived_test_object()
    sce <- sclet_ai_derived_test_node(sce, "ev:findallmarkers", source = "findallmarkers",
        claim = "associated", dependency_group = "g_fm")
    sce
}

# Build an SCE that carries evidence from a real two-source pipeline, the
# same kind of object the spec's "associated" / "consistent_with" cases are
# written against. This is the recipe the spec authors used.
sclet_ai_derived_test_object_with_pipeline <- function() {
    # Raw counts with two real populations, then the FULL pipeline through the
    # registry, because evidence is created only by AIAction handlers.
    #
    # Do NOT pre-process before the plan and do NOT use a 40-cell toy: running
    # find_clusters a second time on an already-clustered small object collapses
    # it to a single cluster, so run_de_test fails with "Must have at least 2
    # groups defined" and no evidence is ever created. That was the original
    # form of this helper and it silently produced zero evidence nodes.
    skip_if_not_installed("BiocNeighbors")
    set.seed(1)
    n_per_group <- 20L
    n_genes <- 40L
    n <- n_per_group * 2L
    mu <- matrix(stats::runif(n_genes, 1, 5), nrow = n_genes, ncol = n)
    mu[1:10, seq_len(n_per_group)] <- mu[1:10, seq_len(n_per_group)] * 6
    mu[11:20, (n_per_group + 1L):n] <- mu[11:20, (n_per_group + 1L):n] * 6
    counts <- matrix(stats::rpois(n_genes * n, lambda = as.vector(mu)), nrow = n_genes)
    dimnames(counts) <- list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n)))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))

    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- sclet:::new_sclet_ai_plan(
        task = "evidence_chain",
        actions = list(
            list(id = "s1", action = "normalize_data", params = list()),
            list(id = "s2", action = "find_variable_features", params = list(nfeatures = 20L)),
            list(id = "s3", action = "scale_data", params = list()),
            list(id = "s4", action = "run_pca", params = list(ncomponents = 10L)),
            list(id = "s5", action = "find_neighbors",
                params = list(dims = as.integer(1:10)),
                depends_on = "s4"),
            list(id = "s6", action = "find_clusters", params = list()),
            list(id = "s7", action = "run_de_test", params = list(name = "markers"))
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    ex <- ExecuteAIPlan(sce, plan, registry, validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token, record = FALSE)
    if (!identical(ex$status, "completed")) {
        failed <- Filter(function(r) identical(r$status, "failed"), ex$results)
        stop("pipeline did not complete: ",
            paste(vapply(failed, function(r) paste0(r$id, "=", r$error %||% "?"), character(1)),
                collapse = "; "),
            call. = FALSE)
    }
    if (!length(sclet:::sclet_ai_evidence_get_all(ex$object))) {
        stop("pipeline completed but produced no evidence; the fixture no longer ",
            "exercises the evidence path", call. = FALSE)
    }
    ex$object
}

# --- Spec row 1: no claim_level -> hypothesis, status ok. ---
#
# v2 contract (`.dev/spec-ai-evidence-authorship-v2.md`). Under v1 this object
# (exactly ONE node) got the node attached to a cite-nothing finding and the
# ceiling raised it to `associated`. v1's rule is measured defective -- it
# attached the node even to findings asserting evidence was ABSENT -- so v2
# attaches only when the finding's own prose names a node id. These rows now
# assert that: naming -> transcribed and raised; not naming -> floor, regardless
# of how many nodes exist.
test_that("spec row 1: a finding that names the sole node is transcribed and raised", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(statement = "as ev:findallmarkers records, the clusters separate"))
    ))
    expect_true(audit$findings[[1L]]$evidence_refs_attached)
    expect_equal(audit$findings[[1L]]$attachment_rule, "prose_named_unique")
    expect_equal(audit$findings[[1L]]$claim_level, "associated")
    expect_equal(audit$findings[[1L]]$requested_claim_level, "hypothesis")
    expect_true(audit$findings[[1L]]$derived_from_ceiling)
    expect_equal(audit$findings[[1L]]$status, "raised")
})

test_that("spec row 1: a cite-nothing finding that names nothing keeps the floor", {
    # The v1 regression guard: a sole node is NO LONGER attached automatically.
    sce <- sclet_ai_derived_test_object_with_evidence()
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(statement = "x"))
    ))
    expect_false(audit$findings[[1L]]$evidence_refs_attached)
    expect_equal(audit$findings[[1L]]$attachment_rule, "no_prose_reference")
    expect_equal(audit$status, "ok")
    expect_equal(audit$findings[[1L]]$claim_level, "hypothesis")
    expect_equal(audit$findings[[1L]]$requested_claim_level, "hypothesis")
    expect_false(audit$findings[[1L]]$derived_from_ceiling)
})

# --- Spec row 2: estimated -> passes untouched. ---
test_that("spec row 2: an estimated finding is provenance and survives the audit", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "estimated", evidence_refs = "ev:findallmarkers"))
    ))
    expect_equal(audit$status, "ok")
    expect_equal(audit$findings[[1L]]$claim_level, "estimated")
    expect_equal(audit$findings[[1L]]$requested_claim_level, "estimated")
    expect_false(audit$findings[[1L]]$derived_from_ceiling)
})

# --- Spec row 3: consistent_with with refs, ceiling associated -> downgrade
#     to associated. ---
test_that("spec row 3: consistent_with with refs and ceiling associated is downgraded to associated", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:findallmarkers"))
    ))
    expect_equal(audit$status, "corrected")
    expect_equal(audit$findings[[1L]]$claim_level, "associated")
    expect_equal(audit$findings[[1L]]$requested_claim_level, "consistent_with")
    expect_equal(audit$findings[[1L]]$allowed_claim_level, "associated")
    expect_true(audit$findings[[1L]]$derived_from_ceiling)
    expect_true(any(grepl("exceeds the ceiling", audit$findings[[1L]]$ceiling_reasons)))
})

# --- Spec row 4: causal is a hard ban; never silently downgraded. ---
test_that("spec row 4: causal stays overclaimed and RecordAIResult still refuses", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "causal", evidence_refs = "ev:findallmarkers"))
    ))
    # the audit must not silently downgrade causal; the status is overclaimed
    # and the problem text says so.
    expect_equal(audit$status, "overclaimed")
    expect_equal(audit$findings[[1L]]$status, "overclaimed")
    expect_true(any(grepl("not a claim level the AI may assert", audit$problems)))
    # the audit's per-finding verdict preserves the request as
    # `requested_claim_level` so the reason for refusal is visible.
    expect_equal(audit$findings[[1L]]$requested_claim_level, "causal")
    # and the entry reports the ban as its basis.
    expect_equal(audit$findings[[1L]]$ceiling_status, "banned")

    # RecordAIResult still refuses to write a causal finding.
    causal_result <- new_sclet_ai_result(
        "review",
        answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "causal", evidence_refs = "ev:findallmarkers"))
    )
    expect_error(
        RecordAIResult(sce, causal_result, id = "ai_causal"),
        "refusing to record an over-claimed AI result"
    )
    expect_false("ai_causal" %in% names(GetAnalysisLedger(sce)$analyses))
})

# --- Spec row 5: consistent_with with no refs -> downgraded to hypothesis. ---
test_that("spec row 5: consistent_with naming no ref lands at hypothesis as a correction", {
    # v2: attachment is name-based, not count-based. A cite-nothing, name-nothing
    # finding on a two-node object declines to attach and exercises the genuine
    # "no refs" path: the request exceeds a finding with no basis, floor holds.
    sce <- sclet_ai_derived_test_object_with_evidence()
    sce <- sclet_ai_derived_test_node(sce, "ev:second", source = "second",
        claim = "associated", dependency_group = "g_second")
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with"))
    ))
    expect_equal(audit$status, "corrected")
    expect_equal(audit$findings[[1L]]$claim_level, "hypothesis")
    expect_equal(audit$findings[[1L]]$requested_claim_level, "consistent_with")
    expect_equal(audit$findings[[1L]]$allowed_claim_level, "hypothesis")
    expect_equal(audit$findings[[1L]]$ceiling_status, "no_refs")
    # the basis travels with the entry
    expect_true(any(grepl("no evidence_refs", audit$findings[[1L]]$ceiling_reasons)))
})

test_that("spec row 5: on a sole-node object a named ref lifts a floor request", {
    # v2: the sole node is only transcribed when the finding names it. A
    # `consistent_with` request that names the node is clamped DOWN to the
    # ceiling (associated), not raised -- the raise-from-floor only applies to a
    # floor request. Naming nothing keeps the hypothesis floor.
    sce <- sclet_ai_derived_test_object_with_evidence()
    named <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with",
            statement = "as ev:findallmarkers records, clusters separate"))
    ))
    expect_true(named$findings[[1L]]$evidence_refs_attached)
    # consistent_with (rank 4) > associated (rank 3) -> clamped to associated, corrected.
    expect_equal(named$findings[[1L]]$claim_level, "associated")
    expect_equal(named$findings[[1L]]$status, "corrected")
    expect_equal(named$findings[[1L]]$ceiling_status, "downgraded")
    # Same request without naming -> the sole node is NOT attached (v1 regression guard).
    unnamed <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with"))
    ))
    expect_false(unnamed$findings[[1L]]$evidence_refs_attached)
    expect_equal(unnamed$findings[[1L]]$claim_level, "hypothesis")
})

# --- Spec row 6: unsupported ceiling -> downgraded to hypothesis. ---
test_that("spec row 6: an unsupported ceiling lands the finding at hypothesis", {
    sce <- sclet_ai_derived_test_object()
    sce <- sclet_ai_derived_test_node(sce, "ev:a", source = "srcA", claim = "associated")
    audit <- sclet:::sclet_ai_audit_result_claims(sce, list(
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:does_not_exist"))
    ))
    expect_equal(audit$status, "corrected")
    expect_equal(audit$findings[[1L]]$claim_level, "hypothesis")
    expect_equal(audit$findings[[1L]]$allowed_claim_level, "hypothesis")
    expect_equal(audit$findings[[1L]]$ceiling_status, "unsupported")
    expect_true(any(grepl("unknown or incomplete", audit$findings[[1L]]$ceiling_reasons)))
})

# --- Spec: same model result, different evidence -> different claim levels. ---
test_that("the same model result against two objects with different evidence yields different levels", {
    # No evidence -> the finding lands at hypothesis. It cites NO refs here on
    # purpose: citing a ref the object does not carry trips the separate,
    # pre-existing "unresolvable evidence" gate in RecordAIResult(), which is
    # correct behaviour (asserting evidence that does not exist is an error)
    # and is covered by its own test above. Using an unresolvable ref here would
    # test that gate rather than the derivation.
    sce_empty <- sclet_ai_derived_test_object()
    res_empty <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "consistent_with"))
    )
    rec_empty <- RecordAIResult(sce_empty, res_empty, id = "ai_empty")
    expect_equal(
        sclet:::sclet_get_state(rec_empty)$analyses$ai_empty$artifacts$findings[[1L]]$claim_level,
        "hypothesis"
    )
    expect_equal(
        sclet:::sclet_get_state(rec_empty)$analyses$ai_empty$artifacts$findings[[1L]]$requested_claim_level,
        "consistent_with"
    )

    # with one strong evidence node -> the finding is downgraded to associated
    sce_with <- sclet_ai_derived_test_object_with_evidence()
    res_with <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:findallmarkers"))
    )
    rec_with <- RecordAIResult(sce_with, res_with, id = "ai_with")
    expect_equal(
        sclet:::sclet_get_state(rec_with)$analyses$ai_with$artifacts$findings[[1L]]$claim_level,
        "associated"
    )
    expect_equal(
        sclet:::sclet_get_state(rec_with)$analyses$ai_with$artifacts$findings[[1L]]$requested_claim_level,
        "consistent_with"
    )
})

# --- Spec: recording against an object with no evidence yields hypothesis. ---
test_that("recording against an object with no evidence yields hypothesis", {
    sce <- sclet_ai_derived_test_object()
    res <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "consistent_with"))
    )
    rec <- RecordAIResult(sce, res, id = "ai_noev")
    rec_finding <- sclet:::sclet_get_state(rec)$analyses$ai_noev$artifacts$findings[[1L]]
    expect_equal(rec_finding$claim_level, "hypothesis")
    expect_equal(rec_finding$requested_claim_level, "consistent_with")
})

# --- estimated survives the recording pipeline. ---
test_that("an estimated finding is recorded as estimated", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    res <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "estimated", evidence_refs = "ev:findallmarkers"))
    )
    rec <- RecordAIResult(sce, res, id = "ai_est")
    rec_finding <- sclet:::sclet_get_state(rec)$analyses$ai_est$artifacts$findings[[1L]]
    expect_equal(rec_finding$claim_level, "estimated")
    expect_equal(rec_finding$requested_claim_level, "estimated")
    expect_false(isTRUE(rec_finding$derived_from_ceiling))
})

# --- no finding is dropped by derivation (preserves `result$findings`). ---
test_that("derivation does not drop any finding from the recorded result", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    res <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(
            list(statement = "a", claim_level = "hypothesis"),
            list(statement = "b", claim_level = "consistent_with", evidence_refs = "ev:findallmarkers"),
            list(statement = "c", claim_level = "estimated", evidence_refs = "ev:findallmarkers")
        )
    )
    rec <- RecordAIResult(sce, res, id = "ai_multi")
    findings <- sclet:::sclet_get_state(rec)$analyses$ai_multi$artifacts$findings
    expect_length(findings, 3L)
    expect_equal(findings[[1L]]$statement, "a")
    expect_equal(findings[[2L]]$statement, "b")
    expect_equal(findings[[3L]]$statement, "c")
})

# --- the audit does not mutate the caller's `result`. ---
test_that("the audit does not rewrite the caller's `result` list", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    finding <- list(claim_level = "consistent_with", evidence_refs = "ev:findallmarkers")
    result <- list(findings = list(finding))
    audit <- sclet:::sclet_ai_audit_result_claims(sce, result)
    # the caller's list is untouched
    expect_equal(result$findings[[1L]]$claim_level, "consistent_with")
    # the audit's entry is a different object
    expect_equal(audit$findings[[1L]]$claim_level, "associated")
})

# --- recording without audit_claims preserves the model's word. ---
test_that("audit_claims = FALSE records the model-asserted level unchanged", {
    sce <- sclet_ai_derived_test_object_with_evidence()
    res <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:findallmarkers"))
    )
    rec <- RecordAIResult(sce, res, id = "ai_off", audit_claims = FALSE)
    rec_finding <- sclet:::sclet_get_state(rec)$analyses$ai_off$artifacts$findings[[1L]]
    expect_equal(rec_finding$claim_level, "consistent_with")
    expect_null(rec_finding$requested_claim_level)
})

# --- the `hypothesis` requested_claim_level is also preserved. ---
test_that("a hypothesis finding records with no derivation", {
    # Two nodes: no attachment, so the finding records exactly as requested and
    # nothing is derived. (On the single-node fixture authorship would attach
    # and raise it -- see `.dev/spec-ai-evidence-authorship.md`.)
    sce <- sclet_ai_derived_test_object_with_evidence()
    sce <- sclet_ai_derived_test_node(sce, "ev:second", source = "second",
        claim = "associated", dependency_group = "g_second")
    res <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "hypothesis"))
    )
    rec <- RecordAIResult(sce, res, id = "ai_hyp")
    rec_finding <- sclet:::sclet_get_state(rec)$analyses$ai_hyp$artifacts$findings[[1L]]
    expect_equal(rec_finding$claim_level, "hypothesis")
    expect_equal(rec_finding$requested_claim_level, "hypothesis")
    expect_false(isTRUE(rec_finding$derived_from_ceiling))
    expect_false(isTRUE(rec_finding$evidence_refs_attached))
})

# --- The `RecordAIResult` exit path for unresolvable evidence refs: the
#     audit reports a downgrade, but `RecordAIResult` itself still refuses
#     via `ValidateAIEvidenceRefs`. ---
test_that("unresolvable evidence refs are still refused at the RecordAIResult level", {
    sce <- sclet_ai_derived_test_object()
    sce <- sclet_ai_derived_test_node(sce, "ev:a", source = "srcA", claim = "associated")
    res <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:does_not_exist"))
    )
    expect_error(
        RecordAIResult(sce, res, id = "ai_unres"),
        "refusing to record an AI result citing unresolvable evidence"
    )
    expect_false("ai_unres" %in% names(GetAnalysisLedger(sce)$analyses))
})

# --- A real two-source pipeline is the spec's "different evidence" object. ---
test_that("a real pipeline object with `ev:findallmarkers` records consistent_with as associated", {
    skip_if_not_installed("BiocNeighbors")
    sce <- sclet_ai_derived_test_object_with_pipeline()
    expect_true(length(sclet:::sclet_ai_evidence_get_all(sce)) > 0L)
    res <- new_sclet_ai_result(
        "review", answer = "x",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(list(claim_level = "consistent_with", evidence_refs = "ev:findallmarkers"))
    )
    rec <- RecordAIResult(sce, res, id = "ai_pipeline")
    rec_finding <- sclet:::sclet_get_state(rec)$analyses$ai_pipeline$artifacts$findings[[1L]]
    expect_equal(rec_finding$claim_level, "associated")
    expect_equal(rec_finding$requested_claim_level, "consistent_with")
    expect_true(isTRUE(rec_finding$derived_from_ceiling))
})
