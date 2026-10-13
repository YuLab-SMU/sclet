# Tests for the variable-feature set as a referenceable plan output,
# per .dev/spec-ai-hvg-reference.md.
#
# The measured failure this guards: on pbmc3k all four planning runs tried to
# reference step 2's feature set and all four were rejected, in four different
# spellings (${step_2.output.required_hvg}, ${step_2.output.hvg},
# ${step_2.output}, ${step_2}) -- sclet-benchmark-ai, execute_plan. The reason
# was that `find_variable_features`' output_schema declared `required_hvg` (a
# flag meaning "the object has an HVG set"), never the set itself.
#
# The design constraint these tests pin down: the gene names must NOT reach the
# planning context. `output_schema` is exposed statically through
# `available_actions` before any step runs, so the static value is the marker
# "hvg" and only the post-execution summary carries the real names.

hvg_reference_fixture <- function(seed = 1L) {
    set.seed(seed)
    n <- 40L
    g <- 60L
    mu <- matrix(runif(g, 1, 5), nrow = g, ncol = n)
    mu[1:15, 1:20] <- mu[1:15, 1:20] * 6
    counts <- matrix(rpois(g * n, as.vector(mu)), nrow = g)
    rownames(counts) <- paste0("g", 1:g)
    colnames(counts) <- paste0("c", 1:n)
    SingleCellExperiment::SingleCellExperiment(list(counts = counts))
}

hvg_pipeline_plan <- function(features_value) {
    new_sclet_ai_plan(
        task = "hvg_reference",
        actions = list(
            list(id = "s1", action = "normalize_data", params = list()),
            list(
                id = "s2", action = "find_variable_features",
                params = list(nfeatures = 20L), depends_on = "s1"
            ),
            list(
                id = "s3", action = "scale_data",
                params = list(features = features_value), depends_on = "s2"
            )
        )
    )
}

test_that("the static declaration is a marker, not the gene set", {
    sce <- hvg_reference_fixture()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    schema <- reg$find_variable_features$output_schema

    expect_true("variable_features" %in% names(schema))
    # A scalar marker. A character vector of gene names here would push 2000
    # names into every planning context through available_actions.
    expect_true(is.character(schema$variable_features))
    expect_length(schema$variable_features, 1L)
    expect_identical(schema$variable_features, "hvg")
    # The flag it used to offer is still there: it means something different
    # and other code depends on it.
    expect_true(isTRUE(schema$required_hvg))
})

test_that("the planning context carries the field name but no gene name", {
    sce <- hvg_reference_fixture()
    # Give the object a real HVG set so the names exist and could leak.
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20L)
    real_hvgs <- VariableFeatures(sce)
    expect_true(length(real_hvgs) > 0L)

    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    descriptors <- sclet:::sclet_ai_describe_registry(reg)
    exposed <- descriptors$find_variable_features$output_schema

    # The name is advertised...
    expect_true("variable_features" %in% names(exposed))
    # ...but no actual gene name is anywhere in the serialised registry view.
    rendered <- paste(utils::capture.output(dput(descriptors)), collapse = "\n")
    for (gene in real_hvgs) {
        expect_false(grepl(gene, rendered, fixed = TRUE),
            info = paste("gene", gene, "must not leak into the planning context"))
    }
})

test_that("the reference resolves to the real set and scales exactly it", {
    sce <- hvg_reference_fixture()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- hvg_pipeline_plan("${s2.output.variable_features}")

    validation <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_true(validation$valid, info = paste(validation$errors, collapse = "; "))

    execution <- ExecuteAIPlan(
        sce, plan, reg,
        validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token, record = FALSE
    )
    expect_identical(execution$status, "completed")

    after <- execution$object
    hvg <- VariableFeatures(after)
    scaled <- as.matrix(SummarizedExperiment::assay(after, "scaled"))
    logcounts <- as.matrix(SummarizedExperiment::assay(after, "logcounts"))

    # The promise resolved to the real set, and the handler received it: this
    # is what the four rejected spellings could never do.
    expect_setequal(
        sclet:::sclet_get_layer(after, "scaled")$params$features,
        hvg
    )
    # Exactly the HVG rows were scaled -- check the values, not the dimensions.
    # ScaleData() with a subset keeps the full assay and writes only the
    # selected rows (R/preprocessing.R), so asserting nrow() would be wrong.
    untouched <- setdiff(rownames(scaled), hvg)
    expect_true(all(scaled[untouched, ] == logcounts[untouched, ]))
    expect_false(all(scaled[hvg, ] == logcounts[hvg, ]))
})

test_that("the execution summary carries the real names, downstream or not", {
    sce <- hvg_reference_fixture()
    sce <- NormalizeData(sce)
    descriptor <- AIDefaultExecutionRegistry(sce, include = "all")$find_variable_features
    after <- FindVariableFeatures(sce, nfeatures = 20L)

    summary <- sclet:::sclet_ai_execution_summary(after, descriptor)
    hvg <- VariableFeatures(after)
    expect_true("variable_features" %in% names(summary))
    # Real names, not the marker: the contract-fill loop must not have
    # overwritten the observation with the static declaration.
    expect_setequal(summary$variable_features, hvg)
    expect_false(identical(summary$variable_features, "hvg"))

    # An action that does not promise the field must not grow one.
    other <- AIDefaultExecutionRegistry(sce, include = "all")$normalize_data
    other_summary <- sclet:::sclet_ai_execution_summary(after, other)
    expect_false("variable_features" %in% names(other_summary))
})

test_that("the four spellings the model actually used are still rejected, and taught", {
    sce <- hvg_reference_fixture()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")

    for (spelling in c(
        "${s2.output.required_hvg}",
        "${s2.output.hvg}",
        "${s2.output}",
        "${s2}"
    )) {
        plan <- hvg_pipeline_plan(spelling)
        validation <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
        expect_false(validation$valid,
            info = paste(spelling, "must not validate"))
    }

    # The unresolvable-field rejection must name the field that does exist,
    # which is how a model repairs toward the right spelling.
    plan <- hvg_pipeline_plan("${s2.output.hvg}")
    validation <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    message <- paste(validation$errors, collapse = " ")
    expect_true(grepl("value_reference_unresolvable", message, fixed = TRUE))
    expect_true(grepl("variable_features", message, fixed = TRUE))
})

test_that("prose in features is still rejected at validation", {
    sce <- hvg_reference_fixture()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- hvg_pipeline_plan("the top HVGs from step 2")

    validation <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_false(validation$valid)
    expect_true(grepl("not usable", paste(validation$errors, collapse = " "), fixed = TRUE))
})

test_that("only the reference form is deferred; a literal marker is not", {
    # The deferred path is driven purely by the reference shape, so it cannot
    # be spoofed by a model typing the marker string. This matters because
    # validation and execution must agree: if the literal were accepted at
    # validation, the handler would still hand "hvg" to the gene selector and
    # fail after the plan had passed.
    is_deferred <- sclet:::sclet_ai_is_deferred_gene_selector
    expect_true(is_deferred("${s2.output.variable_features}"))
    expect_false(is_deferred("hvg"))
    # Other references are not gene-set promises; they keep their own errors.
    expect_false(is_deferred("${s2.output.n_cells}"))
    expect_false(is_deferred("${s2.output}"))
    expect_false(is_deferred("${s2}"))
    expect_false(is_deferred(c("g1", "g2")))
    expect_false(is_deferred("top HVGs from step 2"))
    expect_false(is_deferred(NULL))

    # And it is not just the string: the literal really is rejected by the
    # validator, so validate/execute cannot disagree.
    sce <- hvg_reference_fixture()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    validation <- ValidateAIPlan(
        hvg_pipeline_plan("hvg"), object = sce, registry = reg, strict = FALSE
    )
    expect_false(validation$valid)
    expect_true(grepl("not usable", paste(validation$errors, collapse = " "), fixed = TRUE))
})

test_that("the resolver turns the promise into the real set, or explains itself", {
    sce <- hvg_reference_fixture()
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 20L)
    hvg <- VariableFeatures(sce)

    resolved <- sclet:::sclet_ai_resolve_deferred_selector(
        sce, "${s2.output.variable_features}", "features"
    )
    expect_setequal(resolved, hvg)

    # A value that is not a promise passes through untouched, so ordinary gene
    # vectors, row names and range notation keep working.
    expect_identical(
        sclet:::sclet_ai_resolve_deferred_selector(sce, c("g1", "g2"), "features"),
        c("g1", "g2")
    )
    expect_identical(
        sclet:::sclet_ai_resolve_deferred_selector(sce, "1:20", "features"),
        "1:20"
    )

    # An object with no HVG set cannot satisfy the promise: the error must name
    # the cause rather than reaching the selector as a literal.
    bare <- hvg_reference_fixture()
    bare <- NormalizeData(bare)
    error <- tryCatch(
        sclet:::sclet_ai_resolve_deferred_selector(
            bare, "${s2.output.variable_features}", "features"
        ),
        error = function(e) e
    )
    expect_s3_class(error, "error")
    expect_true(grepl("no variable-feature set", conditionMessage(error), fixed = TRUE))
})

test_that("a plan-level criterion can name the field and be assessed", {
    # The second consumer of the declared name: `sclet_ai_referenceable_output_fields()`
    # offers it to criteria, and `sclet_ai_assess_success_criteria()` reads it out
    # of the step summary. Worth pinning, because a criterion naming a field the
    # summary does not carry reports `not_met` at execution -- the exact
    # validation-accepts-what-execution-cannot-do defect the field list exists to
    # prevent.
    sce <- hvg_reference_fixture()
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(
        task = "hvg_criterion",
        actions = list(
            list(id = "s1", action = "normalize_data", params = list()),
            list(
                id = "s2", action = "find_variable_features",
                params = list(nfeatures = 20L), depends_on = "s1"
            )
        ),
        success_criteria = list(list(
            id = "c1", description = "the HVG set was produced",
            source = "action_output", step_id = "s2",
            check = "exists", field = "variable_features"
        ))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_true(validation$valid, info = paste(validation$errors, collapse = "; "))

    execution <- ExecuteAIPlan(
        sce, plan, reg,
        validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token, record = FALSE
    )
    assessment <- execution$success_assessment[[1L]]
    expect_identical(assessment$status, "met")
    expect_setequal(
        as.character(unlist(assessment$observed_value)),
        VariableFeatures(execution$object)
    )
})

test_that("run_pca reads the same promise through subset_row", {
    sce <- hvg_reference_fixture()
    sce <- NormalizeData(sce)
    reg <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(
        task = "hvg_pca",
        actions = list(
            list(id = "s1", action = "normalize_data", params = list()),
            list(
                id = "s2", action = "find_variable_features",
                params = list(nfeatures = 20L), depends_on = "s1"
            ),
            list(
                id = "s3", action = "run_pca",
                # ncomponents is pinned to stay well under the fixture's 60
                # genes: the registry defaults to 50, and irlba warns "did not
                # converge" on an object that small. That warning is unrelated
                # to the reference being tested.
                params = list(
                    subset_row = "${s2.output.variable_features}",
                    ncomponents = 10L
                ),
                depends_on = "s2"
            )
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = reg, strict = FALSE)
    expect_true(validation$valid, info = paste(validation$errors, collapse = "; "))

    execution <- ExecuteAIPlan(
        sce, plan, reg,
        validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token, record = FALSE
    )
    expect_identical(execution$status, "completed")
    expect_true("PCA" %in% SingleCellExperiment::reducedDimNames(execution$object))
})
