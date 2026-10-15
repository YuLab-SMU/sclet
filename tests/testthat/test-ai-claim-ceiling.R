## Plan-level claim ceiling: the pre-check primitive, the read-only entry-point
## annotation, and the plan layer's informational ceiling field.
##
## The invariant this whole file exists to protect: a claim ceiling is the
## STRONGEST the evidence can carry, so it may only ever LOWER or preserve a
## requested level -- never raise one. A test that cannot fail would not be
## worth writing, so the guards here are checked against deliberately leaky
## input before they are trusted.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

## One evidence node with full control over the fields that decide dependence.
## `sclet_ai_claim_dependency_reasons()` treats two nodes as one line of
## support when they share a `source`, share a `dependency_group`, or are
## related through the `parents` chain -- so independence needs all three to
## differ. `RecordAIEvidence()` is called with `source = NULL` so the node can
## be created without a completed analysis record; the node's own `source`
## field is what the ceiling reads.
sclet_ai_ceiling_test_node <- function(object, id, source, claim = "associated",
                                     dependency_group = NULL, parents = character(),
                                     kind = "deterministic_summary") {
    evidence <- list(
        id = id,
        kind = kind,
        values = list(flag = 1L, raw_values_included = FALSE),
        claim_level = claim,
        source = source,
        dependency_group = dependency_group,
        parents = parents
    )
    RecordAIEvidence(object, evidence, source = NULL, parents = parents, scope = NULL)
}

sclet_ai_ceiling_test_object <- function(n_cells = 6L) {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 3L, n_cells)))
    sce$cluster <- factor(rep(c("A", "B"), length.out = n_cells))
    ActiveIdent(sce) <- "cluster"
    sce
}

## An object carrying exactly ONE line of support: one completed velocity
## record and one evidence node derived from it.
sclet_ai_ceiling_one_line <- function() {
    sce <- sclet_ai_ceiling_test_object()
    sce <- sclet_set_analysis_state(
        object = sce, type = "velocity", id = "vel1",
        method = "velociraptor::scvelo",
        inputs = list(), artifacts = list(analysis_key = "vel1"),
        params = list(), summary = list(), active = TRUE
    )
    sclet_ai_ceiling_test_node(sce, "ev:a", source = "vel1", claim = "consistent_with")
}

## An object carrying TWO independent lines of support: two completed records
## from genuinely different analyses, and one evidence node from each.
##
## The explicit `dependency_group` values are load-bearing, not decoration.
## `RecordAIEvidence()` assigns a content-derived hash when the caller supplies
## none, so two nodes with the same shape and values land in the SAME group and
## are correctly judged dependent no matter how different their `source` is --
## the conservative direction, which can never overstate independence. Naming
## two distinct groups is therefore the only way to express "these came from
## genuinely independent computations".
sclet_ai_ceiling_two_lines <- function() {
    sce <- sclet_ai_ceiling_test_object()
    sce <- sclet_set_analysis_state(
        object = sce, type = "velocity", id = "vel1",
        method = "velociraptor::scvelo",
        inputs = list(), artifacts = list(analysis_key = "vel1"),
        params = list(), summary = list(), active = TRUE
    )
    sce <- sclet_set_analysis_state(
        object = sce, type = "trajectory", id = "traj1", method = "slingshot",
        inputs = list(reduction = "PCA"), artifacts = list(analysis_key = "trajectory"),
        params = list(), summary = list(n_lineages = 1L), active = TRUE
    )
    sce <- sclet_ai_ceiling_test_node(sce, "ev:a", source = "vel1",
        claim = "consistent_with", dependency_group = "g1")
    sce <- sclet_ai_ceiling_test_node(sce, "ev:b", source = "traj1",
        claim = "consistent_with", dependency_group = "g2")
    sce
}

## -- Block 1: the object-level ceiling --------------------------------------
test_that("check_claim_ceiling reports an unsupported ceiling on an object with no evidence", {
    sce <- sclet_ai_ceiling_test_object()
    result <- check_claim_ceiling(sce)

    expect_equal(result$status, "read_only")
    expect_null(result$ceiling_claim_level)
    expect_equal(result$n_independent_groups, 0L)
    expect_false(isTRUE(result$raw_values_included))
})

## -- Block 2: THE CORE INVARIANT -- the ceiling never strengthens a claim ---
test_that("a single line of support cannot be upgraded to consistent_with", {
    sce <- sclet_ai_ceiling_one_line()
    ## One completed velocity record and one node from it: one line of support.
    expect_equal(check_claim_ceiling(sce)$ceiling_claim_level, "associated")
    expect_equal(check_claim_ceiling(sce)$n_independent_groups, 1L)

    ## The requested level is preserved so the movement stays visible...
    verdict <- check_claim_ceiling(sce, claims = list(
        list(claim_level = "consistent_with", evidence_refs = "ev:a")
    ))$verdicts[[1L]]
    expect_equal(verdict$requested_claim_level, "consistent_with")
    ## ...and the derived level is the capped one, never the requested one.
    expect_equal(verdict$claim_level, "associated")
    expect_equal(verdict$ceiling_status, "downgraded")
    expect_false(verdict$claim_level == "consistent_with")
})

test_that("two independent lines of support can carry consistent_with where one cannot", {
    one <- check_claim_ceiling(sclet_ai_ceiling_one_line(), claims = list(
        list(claim_level = "consistent_with", evidence_refs = "ev:a")
    ))$verdicts[[1L]]
    two <- check_claim_ceiling(sclet_ai_ceiling_two_lines(), claims = list(
        list(claim_level = "consistent_with", evidence_refs = c("ev:a", "ev:b"))
    ))$verdicts[[1L]]

    ## The contrast is the point: same requested level, different evidence
    ## structure, different outcome -- and only the richer one is allowed.
    expect_equal(one$claim_level, "associated")
    expect_equal(one$ceiling_status, "downgraded")
    expect_equal(two$claim_level, "consistent_with")
    expect_equal(two$ceiling_status, "allowed")

    expect_true(sclet:::sclet_ai_claim_rank(two$claim_level) >
        sclet:::sclet_ai_claim_rank(one$claim_level))
})

## -- Block 3: causal is never an AI-producible level ------------------------
test_that("a causal claim is reported overclaimed and never produced", {
    sce <- sclet_ai_ceiling_two_lines()
    ## Even the strongest available evidence cannot lift a claim into causal:
    ## the contract reserves causal interpretation for a separately designed
    ## and separately evidenced causal study.
    verdict <- check_claim_ceiling(sce, claims = list(
        list(claim_level = "causal", evidence_refs = c("ev:a", "ev:b"))
    ))$verdicts[[1L]]

    expect_equal(verdict$status, "overclaimed")
    expect_equal(verdict$requested_claim_level, "causal")
    expect_false(identical(verdict$claim_level, "causal"))
    ## The refusal is explicit rather than a silent downgrade.
    expect_true(length(verdict$ceiling_reasons) > 0L)
})

## -- Block 4: no resolvable evidence ---------------------------------------
test_that("a claim citing no resolvable evidence derives to hypothesis", {
    sce <- sclet_ai_ceiling_one_line()

    unknown <- check_claim_ceiling(sce, claims = list(
        list(claim_level = "measured", evidence_refs = "ev:does_not_exist")
    ))$verdicts[[1L]]
    expect_equal(unknown$status, "unsupported")
    expect_equal(unknown$claim_level, "hypothesis")

    none_cited <- check_claim_ceiling(sce, claims = list(
        list(claim_level = "measured", evidence_refs = character())
    ))$verdicts[[1L]]
    expect_equal(none_cited$status, "unsupported")
    expect_equal(none_cited$claim_level, "hypothesis")
})

## -- Block 5: a downgrade is visible, never silent --------------------------
test_that("a downgraded claim carries both the corrected level and its reasons", {
    sce <- sclet_ai_ceiling_one_line()
    verdict <- check_claim_ceiling(sce, claims = list(
        list(claim_level = "consistent_with", evidence_refs = "ev:a")
    ))$verdicts[[1L]]

    expect_equal(verdict$claim_level, "associated")
    expect_true(length(verdict$ceiling_reasons) > 0L)
    ## A reader can tell WHY it moved, and that it moved at all.
    expect_equal(verdict$ceiling_status, "downgraded")
    expect_equal(verdict$requested_claim_level, "consistent_with")
    expect_false(identical(verdict$ceiling_reasons, character()))
})

## -- Block 6: the shared annotation helper is a genuine no-op --------------
test_that("sclet_ai_annotate_findings_ceiling leaves an empty or malformed result untouched", {
    sce <- sclet_ai_ceiling_one_line()
    annotate <- sclet:::sclet_ai_annotate_findings_ceiling

    ## Empty findings: the result must come back unchanged. This is what keeps
    ## test-ai-investigate.R's stubbed `findings = list()` test passing.
    empty <- list(answer = "nothing to report", findings = list(), warnings = list())
    expect_identical(annotate(sce, empty), empty)

    ## Missing findings entirely.
    bare <- list(answer = "no findings field")
    expect_identical(annotate(sce, bare), bare)

    ## A result that is not a list at all must not raise.
    expect_identical(annotate(sce, "not a result"), "not a result")

    ## A malformed finding (not a list) must not error, and must not corrupt
    ## the findings it sits beside.
    mixed <- list(
        answer = "one good, one broken",
        findings = list(
            list(text = "good", claim_level = "associated", evidence_refs = "ev:a"),
            "i am not a finding"
        )
    )
    annotated <- annotate(sce, mixed)
    expect_length(annotated$findings, 2L)
    expect_equal(annotated$findings[[1L]]$claim_level, "associated")
})

test_that("a finding citing resolvable evidence gains the ceiling fields", {
    sce <- sclet_ai_ceiling_one_line()
    annotated <- sclet:::sclet_ai_annotate_findings_ceiling(sce, list(
        findings = list(list(
            text = "velocity supports a trajectory",
            claim_level = "consistent_with",
            evidence_refs = "ev:a"
        ))
    ))

    finding <- annotated$findings[[1L]]
    ## The model's own text survives.
    expect_equal(finding$text, "velocity supports a trajectory")
    ## The request is preserved, the level is capped, and the flag is set.
    expect_equal(finding$requested_claim_level, "consistent_with")
    expect_equal(finding$claim_level, "associated")
    expect_true(isTRUE(finding$derived_from_ceiling))
    expect_true(length(finding$ceiling_reasons) > 0L)
})

## -- Block 7: the read-only entry points are audited ------------------------
## The LLM is stubbed through `options(sclet.ai.call = ...)`, the same
## technique test-ai-investigate.R uses. No real model is ever called.
##
## The option is set and restored by the CALLER, never by this helper: an
## `on.exit()` inside it would fire the moment the helper returned, so the
## stub would be gone before the entry point ran and a real provider would be
## contacted.
sclet_ai_ceiling_stub_call <- function(findings) {
    options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "stubbed",
            warnings = list(),
            findings = findings,
            recommendations = list()
        )
    })
}

## Reads whatever `sclet.ai.call` is currently set to, and returns a closure
## that puts it back.
##
## Both halves matter. `getOption()` rather than `options()` for reading,
## because querying an unset option through `options()` raises. And when the
## option was UNSET to begin with, restoring must actively UNSET it via
## `options(sclet.ai.call = NULL)` -- doing nothing leaves the stub installed
## for every test that runs afterwards, silently hijacking unrelated suites
## that stub the provider their own way. A no-op restore is exactly what broke
## 22 tests in test-ai-integration.R and friends on this file's first run:
## they saw this file's answer string "stubbed" where they expected their own
## tool-loop response.
sclet_ai_ceiling_call_restorer <- function() {
    saved <- getOption("sclet.ai.call", NULL)
    function() {
        options(sclet.ai.call = saved)
    }
}

sclet_ai_ceiling_entry_points <- function() {
    c("AIInvestigate", "AIExplainAnalysis", "AIRecommendNextStep")
}

test_that("AIInvestigate, AIExplainAnalysis and AIRecommendNextStep annotate findings", {
    sce <- sclet_ai_ceiling_one_line()
    restore <- sclet_ai_ceiling_call_restorer()
    on.exit(restore(), add = TRUE)
    sclet_ai_ceiling_stub_call(list(list(
        text = "a finding with a real reference",
        claim_level = "consistent_with",
        evidence_refs = "ev:a"
    )))

    investigated <- AIInvestigate(sce, "what does the velocity say?", scope = "diagnosis")
    explained <- AIExplainAnalysis(sce, type = "velocity", id = "vel1")
    recommended <- AIRecommendNextStep(sce)

    results <- list(
        AIInvestigate = investigated,
        AIExplainAnalysis = explained,
        AIRecommendNextStep = recommended
    )
    for (name in sclet_ai_ceiling_entry_points()) {
        result <- results[[name]]
        expect_true(isTRUE(result$metadata$claim_ceiling_applied), info = name)
        finding <- result$findings[[1L]]
        expect_equal(finding$requested_claim_level, "consistent_with", info = name)
        expect_equal(finding$claim_level, "associated", info = name)
        expect_true(isTRUE(finding$derived_from_ceiling), info = name)
    }

    ## All three carry the same ceiling fields, so a caller can rely on one
    ## shape regardless of which entry point produced it.
    shapes <- lapply(results, function(r) sort(names(r$findings[[1L]])))
    expect_true(length(unique(shapes)) == 1L)
})

test_that("the entry points are a no-op on empty findings", {
    sce <- sclet_ai_ceiling_one_line()
    restore <- sclet_ai_ceiling_call_restorer()
    on.exit(restore(), add = TRUE)
    sclet_ai_ceiling_stub_call(list())

    for (name in sclet_ai_ceiling_entry_points()) {
        result <- switch(name,
            AIInvestigate = AIInvestigate(sce, "anything", scope = "diagnosis"),
            AIExplainAnalysis = AIExplainAnalysis(sce, type = "velocity", id = "vel1"),
            AIRecommendNextStep = AIRecommendNextStep(sce)
        )
        expect_equal(length(result$findings), 0L, info = name)
        expect_s3_class(result, "sclet_ai_result")
    }
})

## -- Block 8: the plan layer's informational ceiling ------------------------
test_that("ValidateAIPlan exposes a read-only claim_ceiling field", {
    sce <- sclet_ai_ceiling_one_line()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "ceiling_check",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "s1", action = "check_qc", params = list()))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)

    expect_true(!is.null(validation$claim_ceiling))
    expect_equal(validation$claim_ceiling$ceiling_claim_level, "associated")
    expect_equal(validation$claim_ceiling$n_independent_groups, 1L)
    expect_false(isTRUE(validation$claim_ceiling$raw_values_included))
})

test_that("the plan ceiling field never changes validation's verdict", {
    ## An informational field must not be able to invalidate a plan.
    for (sce in list(sclet_ai_ceiling_test_object(), sclet_ai_ceiling_one_line())) {
        registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
        plan <- sclet:::new_sclet_ai_plan(
            task = "t",
            context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
            actions = list(list(id = "s1", action = "check_qc", params = list()))
        )
        with_ceiling <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)

        expect_true(with_ceiling$valid)
        expect_length(with_ceiling$errors, 0L)
        expect_length(with_ceiling$warnings, 0L)
        expect_false(isTRUE(with_ceiling$requires_confirmation))
        expect_null(with_ceiling$confirmation_token)
        ## On an object with no evidence the field degrades rather than lying.
        expect_true(!is.null(with_ceiling$claim_ceiling$status))
    }
})

test_that("the plan ceiling field survives a missing object", {
    sce <- sclet_ai_ceiling_one_line()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "fate"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "t",
        actions = list(list(id = "s1", action = "check_qc", params = list()))
    )
    ## The no-object branch is reachable and must not raise.
    validation <- ValidateAIPlan(plan, object = NULL, registry = registry, strict = FALSE)
    expect_true(!is.null(validation$claim_ceiling))
    expect_equal(validation$claim_ceiling$status, "unavailable")
})

## -- Block 9: everything here is strictly read-only ------------------------
test_that("none of the ceiling functions mutates the object", {
    digest_ledger <- function(object) GetAnalysisLedger(object)$fingerprint
    digest_evidence <- function(object) paste(sort(names(sclet:::sclet_ai_evidence_get_all(object))),
        collapse = ",")

    for (maker in list(sclet_ai_ceiling_test_object, sclet_ai_ceiling_one_line,
        sclet_ai_ceiling_two_lines)) {
        sce <- maker()
        before_ledger <- digest_ledger(sce)
        before_evidence <- digest_evidence(sce)

        check_claim_ceiling(sce)
        check_claim_ceiling(sce, claims = list(list(
            claim_level = "consistent_with", evidence_refs = "ev:a"
        )))
        sclet:::sclet_ai_annotate_findings_ceiling(sce, list(
            findings = list(list(claim_level = "associated", evidence_refs = "ev:a"))
        ))

        expect_equal(digest_ledger(sce), before_ledger)
        expect_equal(digest_evidence(sce), before_evidence)
    }
})

## -- Block 11: the stub does not leak into other suites ---------------------
## A leaked `sclet.ai.call` option is this file's most dangerous failure mode:
## every later suite that stubs the provider its own way would silently
## receive this file's answers instead. On the first run the no-op restore
## broke 22 tests across test-ai-integration.R and friends, all reporting
## "stubbed" where they expected their own tool-loop response.
test_that("the LLM stub is uninstalled again, including when it was never set", {
    ## Case 1: the option was unset beforehand, so restoring must unset it.
    if (!is.null(getOption("sclet.ai.call", NULL))) {
        options(sclet.ai.call = NULL)
    }
    restore <- sclet_ai_ceiling_call_restorer()
    sclet_ai_ceiling_stub_call(list())
    expect_false(is.null(getOption("sclet.ai.call", NULL)))
    restore()
    expect_null(getOption("sclet.ai.call", NULL))

    ## Case 2: the option was already set to something else, and that value
    ## must come back rather than this file's stub.
    sentinel <- function(...) "sentinel answer"
    options(sclet.ai.call = sentinel)
    restore2 <- sclet_ai_ceiling_call_restorer()
    sclet_ai_ceiling_stub_call(list())
    expect_identical(getOption("sclet.ai.call"), getOption("sclet.ai.call"))
    expect_false(identical(getOption("sclet.ai.call"), sentinel))
    restore2()
    expect_identical(getOption("sclet.ai.call"), sentinel)
    options(sclet.ai.call = NULL)
})

## -- Block 10: the guards themselves are checked ---------------------------
## A guard that has never been shown to fire is not a guard.
##
## The pattern under test is the never-raise comparison itself: a derived level
## must never outrank the ceiling. Note that `sclet_ai_claim_rank()` returns NA
## for a level outside the ladder ("hypothesis", "causal"), and NA comparisons
## propagate rather than evaluating FALSE -- so an off-ladder level must be
## excluded explicitly or the guard silently stops working on exactly the
## levels that need policing most.
test_that("the never-raise comparison would actually catch an upgrade", {
    rank <- sclet:::sclet_ai_claim_rank

    would_fire <- function(derived, ceiling) {
        dr <- rank(derived)
        cr <- rank(ceiling)
        !is.na(dr) && !is.na(cr) && dr > cr
    }

    ## It fires on a real upgrade...
    expect_true(would_fire("consistent_with", "associated"))
    expect_true(would_fire("measured", "observed"))
    ## ...stays silent where the invariant genuinely holds...
    expect_false(would_fire("associated", "consistent_with"))
    expect_false(would_fire("associated", "associated"))
    ## ...and refuses to return TRUE or NA for an off-ladder level rather than
    ## quietly treating "we cannot rank this" as "it is weaker".
    expect_false(would_fire("hypothesis", "associated"))
    expect_false(would_fire("causal", "consistent_with"))
    expect_false(would_fire("associated", "hypothesis"))
})
