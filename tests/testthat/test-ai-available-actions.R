# Tests for the action vocabulary exposed to planning models
# (`context$available_actions`), including its privacy behaviour.
#
# Why this file exists: `available_actions` was added to the planning context
# without any test, and the privacy classifier's handling of it is not
# predictable from reading its regex -- two candidate field shapes were
# silently destroyed before a working one was found (see
# .dev/spec-ai-ident-groups.md section 2.2). These tests pin the working shape
# and the classifier guard so a future refactor cannot silently re-break it.

sclet_ai_test_aa_object <- function() {
    SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
}

sclet_ai_test_aa_context <- function(sce = sclet_ai_test_aa_object()) {
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    sclet_ai_context(sce, target = "planning", detail = "full", registry = registry)
}

test_that("available_actions lists every registered action with its schemas", {
    ctx <- sclet_ai_test_aa_context()
    expect_true(is.list(ctx$available_actions))
    expect_setequal(
        names(ctx$available_actions),
        names(AIDefaultExecutionRegistry(sclet_ai_test_aa_object(), include = "all"))
    )
    for (entry in ctx$available_actions) {
        expect_true(all(c("name", "description", "input_schema", "output_schema") %in% names(entry)))
    }
})

test_that("output_schema exposes the concrete names a step produces", {
    ctx <- sclet_ai_test_aa_context()
    # The value a downstream action needs when it declares reduction = "PCA".
    expect_equal(ctx$available_actions$run_pca$output_schema$reduction, "PCA")
    # And the assay name normalize_data registers.
    expect_equal(ctx$available_actions$normalize_data$output_schema$assay, "logcounts")
})

test_that("output_schema never carries a `note` field", {
    # Every output_schema `note` is prose, and the privacy classifier hard-denies
    # the field name `note`: leaving it in produced 1 redaction + 1 spurious
    # privacy warning for each of the five actions that carry one.
    ctx <- sclet_ai_test_aa_context()
    leaked <- vapply(ctx$available_actions, function(a) !is.null(a$output_schema$note), logical(1L))
    expect_false(any(leaked))
})

test_that("available_actions survives sanitizing with no redaction or warning", {
    ctx <- sclet_ai_test_aa_context()
    san <- sclet_ai_sanitize_context(
        object = NULL, context = ctx, privacy = "standard", user_consent = FALSE
    )
    expect_length(grep("available_actions", san$policy$redacted_fields, value = TRUE), 0L)
    expect_length(san$policy$warnings, 0L)
    # The whole vocabulary must come through, not a truncated subset.
    expect_setequal(names(san$payload$available_actions), names(ctx$available_actions))
    expect_equal(san$payload$available_actions$run_pca$output_schema$reduction, "PCA")
})

test_that("schema metadata named after clusters is not treated as group label data", {
    # Regression guard for a `"group"`-classification false positive: the field
    # path available_actions.find_clusters.output_schema.required_states.clustering
    # used to be routed into the group k-anonymization path, where its value
    # ("any") failed the aggregation threshold and was redacted with a warning.
    # Schema descriptors are metadata, never label data, so paths under
    # input_schema/output_schema must skip the "group" classification.
    san <- sclet_ai_sanitize_context(
        object = NULL,
        context = list(
            available_actions = list(
                find_clusters = list(
                    name = "find_clusters",
                    input_schema = list(clustering = "character"),
                    output_schema = list(required_states = list(clustering = "any"))
                )
            )
        ),
        privacy = "standard", user_consent = FALSE
    )
    expect_length(san$policy$warnings, 0L)
    expect_length(grep("clustering", san$policy$redacted_fields, value = TRUE), 0L)
    expect_equal(
        san$payload$available_actions$find_clusters$output_schema$required_states$clustering,
        "any"
    )
})

test_that("the group guard still protects real label data outside schemas", {
    # The guard added above must be path-scoped: actual cluster label values
    # elsewhere are still anonymized or redacted as before.
    san <- sclet_ai_sanitize_context(
        object = NULL,
        context = list(
            available_actions = list(
                run_pca = list(output_schema = list(reduction = "PCA"))
            ),
            active_view = list(cluster_labels = c("a", "a", "b"))
        ),
        privacy = "standard", user_consent = FALSE
    )
    # Three labels with counts under the aggregation threshold -> withheld.
    expect_true(length(grep("cluster_labels", san$policy$redacted_fields, value = TRUE)) >= 1L)
    # ...while the schema value of the same name is untouched.
    expect_equal(san$payload$available_actions$run_pca$output_schema$reduction, "PCA")
})
