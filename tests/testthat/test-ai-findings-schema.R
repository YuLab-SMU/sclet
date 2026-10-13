# Tests for the structured-output finding schema and the evidence_refs
# citation vocabulary, per .dev/spec-interpretive-evidence-context.md.
#
# The measured problem (2026-10-12, pbmc3k): the schema declared findings as
# `z_array(z_any())`, so `evidence_refs` had no output channel -- the model
# described `ev:findallmarkers` correctly in prose but cited it in 0 of 12
# (AIReviewQC) and 0 of 19 (AIStatus) findings. These tests pin the fix: the
# field is named in the schema, and the instruction that tells the model to
# fill it reaches both the tool-loop prompt and the structured-summary prompt.

test_that("finding schema names evidence_refs and the claim-level ladder", {
    skip_if_not_installed("aisdk")
    schema <- sclet:::sclet_ai_native_result_schema()
    expect_s3_class(schema, "z_schema")

    finding_schema <- schema$properties$findings$items
    expect_s3_class(finding_schema, "z_object")
    expect_true("evidence_refs" %in% names(finding_schema$properties))
    expect_true("claim_level" %in% names(finding_schema$properties))
    expect_true("severity" %in% names(finding_schema$properties))
    # Tasks name the finding text either way; both stay declared.
    expect_true("statement" %in% names(finding_schema$properties))
    expect_true("description" %in% names(finding_schema$properties))
    # Task-specific extras must not be rejected by the provider.
    expect_true(isTRUE(finding_schema$additionalProperties))

    # The enum must not silently drop a level the local validator accepts.
    declared <- finding_schema$properties$claim_level$enum
    expect_setequal(declared, sclet:::sclet_ai_finding_claim_levels())
    expect_setequal(
        finding_schema$properties$severity$enum,
        sclet:::sclet_ai_finding_severities()
    )
    expect_true("causal" %in% declared)
    expect_true("hypothesis" %in% declared)
})

test_that("the local validator accepts every level the schema declares", {
    # The two vocabularies are read from one helper (ai-adapter.R), so a level
    # the provider is allowed to emit can never be rejected by the validator.
    for (level in sclet:::sclet_ai_finding_claim_levels()) {
        result <- new_sclet_ai_result(
            "test_task",
            findings = list(list(
                statement = "s", severity = "info", claim_level = level,
                evidence_refs = if (identical(level, "causal")) "ev:x" else character()
            ))
        )
        expect_s3_class(result, "sclet_ai_result")
    }
    # ...and a level outside the ladder is still rejected.
    expect_error(
        new_sclet_ai_result(
            "test_task",
            findings = list(list(statement = "s", claim_level = "proven"))
        ),
        class = "sclet_ai_invalid_output"
    )
    for (sev in sclet:::sclet_ai_finding_severities()) {
        expect_s3_class(
            new_sclet_ai_result(
                "test_task",
                findings = list(list(statement = "s", severity = sev))
            ),
            "sclet_ai_result"
        )
    }
})

test_that("the citation instruction is carried by the structured-summary prompt", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    sce <- RecordAIEvidence(sce, list(
        id = "ev:findallmarkers",
        kind = "deterministic_summary",
        values = list(n_rows = 100L)
    ))

    agent <- new.env(parent = emptyenv())
    agent$run <- function(task, max_steps = 10) {
        list(text = "the only evidence record ev:findallmarkers is associated level")
    }
    captured <- new.env(parent = emptyenv())
    old <- options(
        sclet.ai.create_agent = function(...) agent,
        sclet.ai.generate_object = function(model, prompt, schema, schema_name,
            system, mode, max_retries, ...) {
            captured$system <- system
            captured$schema <- schema
            list(object = list(
                answer = "structured",
                findings = list(list(
                    statement = "only evidence is associated level",
                    severity = "info",
                    claim_level = "hypothesis",
                    evidence_refs = "ev:findallmarkers"
                )),
                evidence = list(), warnings = list(),
                recommendations = list(), proposed_actions = list()
            ))
        }
    )
    on.exit(options(old), add = TRUE)

    result <- AIStatus(sce)
    expect_s3_class(result, "sclet_ai_result")
    expect_true(isTRUE(result$metadata$structured_output))

    # The instruction itself reached the provider's system prompt...
    expect_true(grepl("evidence_refs", captured$system, fixed = TRUE))
    expect_true(grepl("ai_evidence", captured$system, fixed = TRUE))
    # ...and the schema it was handed names the field.
    expect_true(
        "evidence_refs" %in% names(captured$schema$properties$findings$items$properties)
    )
    # The round trip preserves the citation.
    expect_equal(result$findings[[1L]]$evidence_refs, "ev:findallmarkers")
})

test_that("the citation instruction is also carried by the tool-loop prompt", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    agent <- new.env(parent = emptyenv())
    captured <- new.env(parent = emptyenv())
    agent$run <- function(task, max_steps = 10) {
        list(text = "tool loop response")
    }
    old <- options(
        sclet.ai.create_agent = function(name, description, system_prompt, ...) {
            captured$system_prompt <- system_prompt
            agent
        },
        sclet.ai.generate_object = function(...) {
            list(object = list(
                answer = "structured",
                findings = list(), evidence = list(), warnings = list(),
                recommendations = list(), proposed_actions = list()
            ))
        }
    )
    on.exit(options(old), add = TRUE)

    AIStatus(sce)
    # The mock factory only receives the named args, so read it back through
    # the option the adapter used.
    expect_true(grepl("evidence_refs", captured$system_prompt, fixed = TRUE))
    expect_true(grepl("ai_evidence", captured$system_prompt, fixed = TRUE))
})

test_that("the instruction forbids citing ids that are not in the ledger", {
    instruction <- sclet:::sclet_ai_evidence_citation_instruction()
    expect_true(grepl("never invent", instruction, fixed = TRUE))
    expect_true(grepl("verbatim", instruction, fixed = TRUE))
})
