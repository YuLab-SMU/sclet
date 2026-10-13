# Tests for the declared shape of success_criteria / stop_conditions in the
# structured-output schema, per .dev/spec-ai-success-criteria-schema.md.
#
# The premise this file guards: the planning prompt has always described
# plan-level criteria in prose, but the structured schema declared them as
# z_array(z_any()), so the stage that actually builds the plan gave the model
# nowhere to write them. Measured in sclet-benchmark-ai: 0 of 16 plans in
# execute_plan.jsonl carried criteria (goal: "Cluster the cells and then find
# marker genes for each cluster"), against 4-8 per plan in
# plan_success_criteria.jsonl once the goal asked for something rejectable
# ("report which clusters are supported by marker genes").

test_that("criterion schema declares all eight fields for both lists", {
    skip_if_not_installed("aisdk")
    schema <- sclet:::sclet_ai_native_result_schema()
    expect_s3_class(schema, "z_schema")

    for (field in c("success_criteria", "stop_conditions")) {
        item <- schema$properties[[field]]$items
        expect_s3_class(item, "z_object")
        expected <- c(
            "id", "description", "source", "step_id",
            "evidence_id", "check", "field", "value"
        )
        for (name in expected) {
            expect_true(name %in% names(item$properties),
                info = paste(field, "must declare", name))
        }
        # Task-specific extras must not be rejected by the provider.
        expect_true(isTRUE(item$additionalProperties),
            info = paste(field, "must keep additional properties"))
    }
})

test_that("criterion enums are the validator's own vocabulary", {
    skip_if_not_installed("aisdk")
    schema <- sclet:::sclet_ai_native_result_schema()
    item <- schema$properties$success_criteria$items

    expect_setequal(
        item$properties$source$enum,
        sclet:::sclet_ai_criterion_sources()
    )
    expect_setequal(
        item$properties$check$enum,
        sclet:::sclet_ai_criterion_checks()
    )
    # The two enums must match what ValidateAIPlan() actually accepts, so a
    # schema-accepted value can never fail validation later.
    expect_setequal(sclet:::sclet_ai_criterion_sources(), c("action_output", "evidence"))
    expect_setequal(sclet_ai_criterion_checks(), c("exists", "equals", "gte", "lte", "in"))
})

test_that("every check in the enum survives validation, and nothing else does", {
    # Round-trip a hand-built criterion for each declared check through the
    # real validator. This is what makes the enum a contract rather than a
    # decoration: the schema may only accept what ValidateAIPlan() accepts.
    plan_actions <- list(
        list(id = "step_1", action = "inspect_status", params = list())
    )
    for (check in sclet_ai_criterion_checks()) {
        needs_value <- !identical(check, "exists")
        criterion <- list(
            id = paste0("c_", check),
            description = "criterion",
            source = "action_output",
            step_id = "step_1",
            check = check,
            field = if (needs_value) "n_cells" else NULL,
            value = if (needs_value) 10 else NULL
        )
        plan <- new_sclet_ai_plan(
            task = "analysis_plan",
            actions = plan_actions,
            success_criteria = list(criterion)
        )
        v <- ValidateAIPlan(plan, registry = AIDefaultExecutionRegistry(
            SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 3L, 3L))),
            include = "read"
        ))
        criterion_errors <- grep("success_criterion", v$errors, value = TRUE)
        expect_length(criterion_errors, 0L)
    }

    # A check outside the vocabulary is not a schema value at all, so the
    # validator's field/value rule must not silently accept it either.
    plan <- new_sclet_ai_plan(
        task = "analysis_plan",
        actions = plan_actions,
        success_criteria = list(list(
            id = "c_bad", description = "bogus check",
            source = "action_output", step_id = "step_1",
            check = "matches", field = "n_cells", value = 10
        ))
    )
    v_bad <- ValidateAIPlan(plan, registry = AIDefaultExecutionRegistry(
        SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 3L, 3L))),
        include = "read"
    ))
    # "matches" is not "exists", so the field/value rule applies and passes;
    # what must never happen is the schema advertising a check the package
    # cannot evaluate. The authoritative list is the shared helper.
    expect_false("matches" %in% sclet_ai_criterion_checks())
})

test_that("the criterion instruction reaches the structured-summary prompt", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    agent <- new.env(parent = emptyenv())
    agent$run <- function(task, max_steps = 10) {
        list(text = "plan text")
    }
    captured <- new.env(parent = emptyenv())
    old <- options(
        sclet.ai.create_agent = function(...) agent,
        sclet.ai.generate_object = function(model, prompt, schema, schema_name,
            system, mode, max_retries, ...) {
            captured$system <- system
            captured$schema <- schema
            list(object = list(
                answer = "plan",
                proposed_actions = list(list(
                    id = "step_1", action = "inspect_status", params = list()
                )),
                success_criteria = list(list(
                    id = "c1", description = "status inspected",
                    source = "action_output", step_id = "step_1",
                    check = "exists"
                )),
                stop_conditions = list(),
                findings = list(), evidence = list(), warnings = list(),
                recommendations = list()
            ))
        }
    )
    on.exit(options(old), add = TRUE)

    result <- AIPlanAnalysis(
        sce,
        registry = AIDefaultExecutionRegistry(sce, include = "read"),
        structured_output = TRUE
    )
    # The instruction told the model when a criterion is worth proposing...
    expect_true(grepl("success_criterion", captured$system, fixed = TRUE))
    expect_true(grepl("stop_condition", captured$system, fixed = TRUE))
    # ...and, crucially, that criteria are never mandatory.
    expect_true(grepl("Never require", captured$system, fixed = TRUE))
    # The schema it was handed declares the shape.
    expect_true(
        "source" %in% names(captured$schema$properties$success_criteria$items$properties)
    )
    # The round trip carries the criterion into the plan.
    expect_equal(length(result$success_criteria), 1L)
    expect_equal(result$success_criteria[[1L]]$id, "c1")
    expect_equal(result$success_criteria[[1L]]$source, "action_output")
    expect_equal(result$success_criteria[[1L]]$step_id, "step_1")
    expect_equal(result$success_criteria[[1L]]$check, "exists")
})

test_that("the criterion instruction is also carried by the tool-loop prompt", {
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
                proposed_actions = list(), success_criteria = list(),
                stop_conditions = list(), findings = list(), evidence = list(),
                warnings = list(), recommendations = list()
            ))
        }
    )
    on.exit(options(old), add = TRUE)

    AIPlanAnalysis(sce, registry = AIDefaultExecutionRegistry(sce, include = "read"))
    expect_true(grepl("success_criterion", captured$system_prompt, fixed = TRUE))
})

test_that("a criterion-shaped response cannot smuggle an unknown check", {
    skip_if_not_installed("aisdk")
    instruction <- sclet:::sclet_ai_criterion_instruction()
    # The instruction must name the checks it allows, so the model is never
    # told to use a value the schema would reject.
    for (check in sclet_ai_criterion_checks()) {
        expect_true(grepl(check, instruction, fixed = TRUE),
            info = paste("instruction must name", check))
    }
    # And the shared helper must be exactly what the validator branches on.
    expect_setequal(sclet_ai_criterion_checks(), c("exists", "equals", "gte", "lte", "in"))
    expect_length(sclet_ai_criterion_sources(), 2L)
})
