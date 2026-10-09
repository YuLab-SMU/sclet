# A prerequisite verdict must not be derived from parameters the schema
# already rejected.
#
# Regression: ValidateAIPlan() used to run the per-action prerequisite probe
# even for a step whose params had just failed input_schema validation. The
# probe therefore received a value the action never expects, which produced a
# second, misleading error alongside the one actionable schema error:
#
#   - invalid params for action step_1: parameter n has invalid type
#   - prerequisite check failed for step_1: missing value where TRUE/FALSE needed
#
# ...or, worse, a confident but meaningless verdict (a numeric comparison
# silently becoming a lexicographic one on a character value). Since these
# errors are fed back to the model by the repair loop, the model was being
# pointed at the wrong thing to fix.
#
# The same class of defect had a second source: a malformed value reference.
# The plan syntax is ${step_id.output.field}, but a model writing
# ${step_4.reduction} (no `output` segment) was handled inconsistently --
# sclet_ai_collect_plan_references() accepted the token, so the step was
# charged with an undeclared dependency, while sclet_ai_project_plan_value()
# silently passed the raw template through, so the action's own prerequisite
# reported "reduction is not available: ${step_4.reduction}". The syntax
# mistake itself was never named.

prereq_test_sce <- function() {
    SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:12, nrow = 4))
    )
}

test_that("a step with invalid params does not get a derived prerequisite error", {
    action <- AIAction(
        name = "needs_int",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) {
            n <- as.integer(params$n)
            if (n > 2) "n_too_big" else TRUE
        },
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(n = "integer")
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "needs_int", params = list(n = "5x")))
    )

    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )

    # The actionable error is kept...
    expect_true(any(grepl("invalid params for action step_1", validation$errors)))
    # ...and the derived one is gone.
    expect_false(any(grepl("prerequisite", validation$errors)))
    expect_false(any(grepl("n_too_big", validation$errors)))
    expect_false(any(grepl("missing value where TRUE/FALSE needed", validation$errors)))
})

test_that("the skip is reported rather than silent", {
    action <- AIAction(
        name = "needs_int",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) TRUE,
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(n = "integer")
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "needs_int", params = list(n = "5x")))
    )

    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )

    expect_true(any(grepl("prerequisite check skipped for action step_1", validation$warnings)))
})

test_that("probing an invalid-param step emits no conditions of its own", {
    action <- AIAction(
        name = "needs_int",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) {
            n <- as.integer(params$n)
            if (n > 2) "n_too_big" else TRUE
        },
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(n = "integer")
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "needs_int", params = list(n = "5x")))
    )

    # The old probe leaked "NAs introduced by coercion" out of the validator.
    expect_no_warning(
        ValidateAIPlan(plan, object = prereq_test_sce(), registry = registry, strict = FALSE)
    )
})

test_that("a valid-param step still runs its prerequisite (skip is not blanket)", {
    action <- AIAction(
        name = "needs_int",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) "distinctive_prereq_verdict",
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(n = "integer")
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "needs_int", params = list(n = 5L)))
    )

    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )

    expect_true(any(grepl("distinctive_prereq_verdict", validation$errors)))
    expect_false(any(grepl("prerequisite check skipped", validation$warnings)))
})

test_that("skipping the probe cannot flip an invalid plan to valid", {
    action <- AIAction(
        name = "needs_int",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) TRUE,
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(n = "integer")
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "needs_int", params = list(n = "5x")))
    )

    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )

    # The schema error alone already marks the plan invalid, so dropping the
    # derived error can only remove noise, never change the verdict.
    expect_false(validation$valid)
    expect_true(length(validation$errors) >= 1L)
})

test_that("a skipped step still propagates its declared output to later steps", {
    producer <- AIAction(
        name = "producer",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) TRUE,
        returns = "sce",
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(n = "integer"),
        output_schema = list(active_reduction = "PCA", required_reductions = "PCA")
    )
    consumer <- AIAction(
        name = "consumer",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) {
            if (!identical(planned$active_reduction, "PCA")) "no_pca_planned" else TRUE
        },
        mutates_object = TRUE,
        requires_confirmation = TRUE
    )
    registry <- AIExecutionRegistry(list(producer, consumer))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(
            list(id = "step_1", action = "producer", params = list(n = "5x")),
            list(id = "step_2", action = "consumer", params = list(), depends_on = "step_1")
        )
    )

    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )

    # Skipping the *probe* for step_1 must not drop the output it declares:
    # a later step's prerequisite still needs to see the planned state.
    expect_false(any(grepl("no_pca_planned", validation$errors)))
    expect_true(any(grepl("prerequisite check skipped for action step_1", validation$warnings)))
})

test_that("malformed value references are detected, correct ones are not", {
    expect_identical(
        sclet:::sclet_ai_malformed_plan_references("${step_4.reduction}"),
        "${step_4.reduction}"
    )
    expect_identical(
        sclet:::sclet_ai_malformed_plan_references("${step_4.output.reduction}"),
        character()
    )
    # A bare ${a.b} has no field segment at all.
    expect_identical(
        sclet:::sclet_ai_malformed_plan_references("${step_4.output}"),
        "${step_4.output}"
    )
    # ...and it recurses into nested parameter lists.
    expect_identical(
        sclet:::sclet_ai_malformed_plan_references(
            list(a = 1L, b = list(c = "${step_1.reduction}"))
        ),
        "${step_1.reduction}"
    )
    expect_identical(sclet:::sclet_ai_malformed_plan_references(30L), character())
    expect_identical(sclet:::sclet_ai_malformed_plan_references("plain"), character())
})

reference_test_actions <- function() {
    producer <- AIAction(
        name = "producer",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) TRUE,
        returns = "sce",
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        output_schema = list(reduction = "PCA", required_reductions = "PCA")
    )
    consumer <- AIAction(
        name = "consumer",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) {
            if (!params$reduction %in% (planned$reductions %||% character())) {
                paste0("reduction is not available: ", params$reduction)
            } else {
                TRUE
            }
        },
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(reduction = "character")
    )
    AIExecutionRegistry(list(producer, consumer))
}

reference_test_plan <- function(ref) {
    sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(
            list(id = "step_1", action = "producer", params = list()),
            list(id = "step_2", action = "consumer", params = list(reduction = ref), depends_on = "step_1")
        )
    )
}

test_that("a malformed reference is named, not passed through to the action", {
    validation <- ValidateAIPlan(
        reference_test_plan("${step_1.reduction}"),
        object = prereq_test_sce(),
        registry = reference_test_actions(),
        strict = FALSE
    )

    expect_true(any(grepl("value_reference_malformed", validation$errors)))
    expect_true(any(grepl("step_1.output.reduction", validation$errors)))
    # The old symptom -- the raw template surfacing from the prerequisite --
    # must be gone.
    expect_false(any(grepl("reduction is not available", validation$errors)))
    expect_false(validation$valid)
})

test_that("a well-formed reference still resolves (control)", {
    validation <- ValidateAIPlan(
        reference_test_plan("${step_1.output.reduction}"),
        object = prereq_test_sce(),
        registry = reference_test_actions(),
        strict = FALSE
    )

    expect_true(validation$valid)
    expect_false(any(grepl("value_reference_malformed", validation$errors)))
})

test_that("the malformed-reference error suggests the corrected spelling", {
    expect_identical(
        sclet:::sclet_ai_suggest_reference_fix("${step_4.reduction}"),
        "${step_4.output.reduction}"
    )
    expect_identical(
        sclet:::sclet_ai_suggest_reference_fix("${step_4}"),
        "${step_4.output.<field>}"
    )
    # Already namespaced but missing a field: keep `output`, flag the field.
    expect_identical(
        sclet:::sclet_ai_suggest_reference_fix("${step_4.output}"),
        "${step_4.output.<field>}"
    )
    # A deeper path is preserved.
    expect_identical(
        sclet:::sclet_ai_suggest_reference_fix("${step_4.summary.n"),
        "${step_4.output.summary.n}"
    )
})

test_that("the reported suggestion is copy-pasteable back into the plan", {
    bad <- "${step_1.reduction}"
    fixed <- sclet:::sclet_ai_suggest_reference_fix(bad)
    validation <- ValidateAIPlan(
        reference_test_plan(fixed),
        object = prereq_test_sce(),
        registry = reference_test_actions(),
        strict = FALSE
    )
    # The suggestion must actually be a valid reference, not just prose.
    expect_true(validation$valid)
})

test_that("an absent parameter still gets its prerequisite diagnostic", {
    # Regression guard on the fix itself. Narrowing "unusable params" to
    # present-but-wrong-type must not swallow the prerequisite's own message
    # about a *missing* required parameter: run_trajectory's
    # start_cluster_missing, and the clarification question derived from it,
    # depend on the probe still running.
    action <- AIAction(
        name = "needs_start",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) {
            if (is.null(params$start_cluster) || !nzchar(params$start_cluster)) {
                return("start_cluster_missing: provide an explicit start_cluster")
            }
            TRUE
        },
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(
            group = "character",
            start_cluster = list(type = "character", required = TRUE)
        )
    )
    registry <- AIExecutionRegistry(list(action))

    for (params in list(
        list(group = "cluster"),
        list(group = "cluster", start_cluster = NULL)
    )) {
        plan <- sclet:::new_sclet_ai_plan(
            task = "analysis_plan",
            actions = list(list(id = "step_1", action = "needs_start", params = params))
        )
        validation <- ValidateAIPlan(
            plan,
            object = prereq_test_sce(),
            registry = registry,
            strict = FALSE
        )
        expect_true(any(grepl("start_cluster_missing", validation$errors)))
        expect_false(any(grepl("prerequisite check skipped", validation$warnings)))
    }
})

test_that("an enumeration violation still probes (only type errors skip)", {
    action <- AIAction(
        name = "enumerated",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) "enum_probe_ran",
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(method = c("scran", "seurat"))
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "enumerated", params = list(method = "banana")))
    )
    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )
    expect_true(any(grepl("must be one of", validation$errors)))
    expect_true(any(grepl("enum_probe_ran", validation$errors)))
})

test_that("param type validity has one definition shared by both callers", {
    expect_true(sclet:::sclet_ai_param_type_valid(30, "integer"))
    expect_true(sclet:::sclet_ai_param_type_valid(30L, "integer"))
    expect_false(sclet:::sclet_ai_param_type_valid("30", "integer"))
    expect_false(sclet:::sclet_ai_param_type_valid(30.5, "integer"))
    expect_true(sclet:::sclet_ai_param_type_valid("x", "character"))
    expect_true(sclet:::sclet_ai_param_type_valid(TRUE, "logical"))

    # Only present, wrongly-typed parameters are named.
    schema <- list(a = "integer", b = list(type = "character", required = TRUE), c = c("x", "y"))
    expect_identical(
        sclet:::sclet_ai_action_param_type_problems(list(a = "no", b = "yes"), schema),
        "a"
    )
    expect_identical(
        sclet:::sclet_ai_action_param_type_problems(list(b = "yes"), schema),
        character()
    )
    expect_identical(
        sclet:::sclet_ai_action_param_type_problems(list(c = "z"), schema),
        character()
    )
})

test_that("a type error names the value that was actually received", {
    action <- AIAction(
        name = "needs_int",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) TRUE,
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(n = "integer", k = "character")
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "needs_int", params = list(n = "30", k = 5L)))
    )
    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )
    joined <- paste(validation$errors, collapse = "\n")
    # Naming the received shape is what makes the error actionable: "invalid
    # type" alone left both the model and a benchmark reader guessing between
    # the number as text, a placeholder, and something else.
    expect_true(grepl("expected integer", joined, fixed = TRUE))
    expect_true(grepl("got the text \"30\"", joined, fixed = TRUE))
    expect_true(grepl("expected character", joined, fixed = TRUE))
    expect_true(grepl("got the number 5", joined, fixed = TRUE))
})

test_that("an enumeration error names the value that was actually received", {
    action <- AIAction(
        name = "enumerated",
        handler = function(object, params) object,
        prerequisites = function(object, params, planned = NULL) TRUE,
        mutates_object = TRUE,
        requires_confirmation = TRUE,
        input_schema = list(method = c("scran", "seurat"))
    )
    registry <- AIExecutionRegistry(list(action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "step_1", action = "enumerated", params = list(method = "logreg")))
    )
    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = registry,
        strict = FALSE
    )
    joined <- paste(validation$errors, collapse = "\n")
    expect_true(grepl("must be one of: scran, seurat", joined, fixed = TRUE))
    expect_true(grepl("got the text \"logreg\"", joined, fixed = TRUE))
})

test_that("parameter value descriptions are bounded and type-aware", {
    describe <- sclet:::sclet_ai_describe_param_value
    expect_identical(describe("30"), "the text \"30\"")
    expect_identical(describe(30L), "the number 30")
    expect_identical(describe(30.5), "the number 30.5")
    expect_identical(describe(TRUE), "the logical value TRUE")
    expect_identical(describe(c("a", "b")), "a character vector of length 2")
    expect_identical(describe(list(a = 1L, b = 2L)), "a list with 2 element(s) named a, b")
    expect_identical(describe(NA_character_), "the text NA (a missing value)")
    # Long payloads are truncated: the point is the shape, not the content.
    long <- describe(strrep("z", 500L))
    expect_true(nchar(long) < 100L)
    expect_true(grepl("...", long, fixed = TRUE))
})

reference_registry <- function() {
    AIExecutionRegistry(list(
        AIAction(
            name = "producer",
            handler = function(object, params) object,
            prerequisites = function(object, params, planned = NULL) TRUE,
            returns = "sce",
            mutates_object = TRUE,
            requires_confirmation = TRUE,
            output_schema = list(x = "X", y = "Y", z = "Z")
        ),
        AIAction(
            name = "consumer",
            handler = function(object, params) object,
            prerequisites = function(object, params, planned = NULL) TRUE,
            mutates_object = TRUE,
            requires_confirmation = TRUE,
            input_schema = list(a = "character", b = "character", c = "character")
        )
    ))
}

test_that("references are collected per reference, not as one flat vector", {
    # A flat vector only lines up with references when every reference has the
    # same number of segments. Walking it with a fixed stride read a field
    # name as a step id ("references unknown step(s): y") or read an empty
    # first segment and reported "references unknown step(s): " with no name.
    refs <- sclet:::sclet_ai_collect_plan_references(list(
        a = "${step_2.output.x}",
        b = "${step_3}",
        c = "${step_4.output.y}"
    ))
    expect_length(refs, 3L)
    expect_identical(refs[[1L]], c("step_2", "output", "x"))
    expect_identical(refs[[2L]], "step_3")
    expect_identical(refs[[3L]], c("step_4", "output", "y"))
})

test_that("mixed reference shapes produce no bogus unknown-step error", {
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(
            list(id = "step_1", action = "producer", params = list()),
            list(
                id = "step_2", action = "consumer",
                params = list(
                    a = "${step_1.output.x}",
                    b = "${step_1.output.y}",
                    c = "${step_1.output.z}"
                ),
                depends_on = "step_1"
            )
        )
    )
    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = reference_registry(),
        strict = FALSE
    )
    expect_true(validation$valid)
    expect_false(any(grepl("references unknown step", validation$errors)))
})

test_that("a malformed reference yields exactly one error, not two", {
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(
            list(id = "step_1", action = "producer", params = list()),
            list(id = "step_2", action = "consumer", params = list(a = "${step_1.x}"))
        )
    )
    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = reference_registry(),
        strict = FALSE
    )
    expect_true(any(grepl("value_reference_malformed", validation$errors)))
    # Deriving dependencies from a token we cannot parse only adds noise.
    expect_false(any(grepl("references undeclared dependency", validation$errors)))
    expect_false(any(grepl("references unknown step", validation$errors)))
})

test_that("the well-formedness grammar has one definition", {
    well_formed <- sclet:::sclet_ai_reference_is_well_formed
    expect_true(well_formed(c("step_1", "output", "x")))
    expect_true(well_formed(c("step_1", "output", "a", "b")))
    expect_false(well_formed(c("step_1", "x")))
    expect_false(well_formed(c("step_1", "output")))
    expect_false(well_formed(c("", "output", "x")))
    expect_false(well_formed(character()))
    # The detector and the grammar must agree on every one of these.
    templates <- c("${step_1.output.x}", "${step_1.output}", "${step_1.x}", "${.output.x}")
    expected_malformed <- c(FALSE, TRUE, TRUE, TRUE)
    for (i in seq_along(templates)) {
        got <- length(sclet:::sclet_ai_malformed_plan_references(templates[[i]])) > 0L
        expect_identical(got, expected_malformed[[i]], info = templates[[i]])
    }
})

test_that("an undeclared dependency error says how to fix it", {
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(
            list(id = "step_1", action = "producer", params = list()),
            list(id = "step_2", action = "consumer", params = list(a = "${step_1.output.x}"))
        )
    )
    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = reference_registry(),
        strict = FALSE
    )
    expect_true(any(grepl("add them to this step's depends_on", validation$errors)))
})

test_that("a placeholder action name is told to drop the step instead", {
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(
            list(id = "step_1", action = "producer", params = list()),
            list(id = "step_2", action = "none", params = list())
        )
    )
    validation <- ValidateAIPlan(
        plan,
        object = prereq_test_sce(),
        registry = reference_registry(),
        strict = FALSE
    )
    expect_true(any(grepl("placeholder such as \"none\"", validation$errors)))
    expect_true(any(grepl("remove the step from the plan", validation$errors)))
})
