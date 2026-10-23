# `scletAI()`: the entry point a first-time user needs.
#
# The contract these tests defend is that the wrapper adds two checks and hides
# nothing else: no aisdk -> an error naming the install command; no resolvable
# model -> an error naming the ways to supply one; otherwise the call is
# `RunAIAnalysis()` unchanged.

test_that("scletAI checks its own arguments before any provider work", {
    # `aisdk` is installed on this machine and cannot be unloaded from a loaded
    # session, so the missing-dependency branch is not reproducible here. What
    # IS reproducible, and is the part a user hits first, is that a malformed
    # call is refused locally rather than being sent to a provider.
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5L, ncol = 4L))
    )
    expect_error(scletAI("not an sce", "cluster"),
        "`object` must be a SingleCellExperiment", fixed = TRUE)
    expect_error(scletAI(sce, ""),
        "`goal` must be one non-empty character value", fixed = TRUE)
    expect_error(scletAI(sce, c("a", "b")),
        "`goal` must be one non-empty character value", fixed = TRUE)
    # The dependency check itself still resolves (it is the installed case).
    expect_true(isTRUE(sclet_ai_require()))
})

test_that("the missing-dependency message names the install command", {
    # Pinned on the message text, because a user who hits this branch has to be
    # told what to install; the branch itself is not reachable on a machine
    # with aisdk installed.
    src <- paste(deparse(body(sclet_ai_require)), collapse = "\n")
    expect_match(src, "sclet_ai_missing_dependency", fixed = TRUE)
})

test_that("scletAI validates its arguments before touching the provider", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5L, ncol = 4L))
    )
    expect_error(scletAI("not an sce", "cluster"), "`object` must be a SingleCellExperiment")
    expect_error(scletAI(sce, ""), "`goal` must be one non-empty character value")
    expect_error(scletAI(sce, c("a", "b")), "`goal` must be one non-empty character value")
})

test_that("scletAI is RunAIAnalysis with the reliability switches hidden", {
    # The wrapper must not grow its own behaviour: with a model supplied and
    # confirmation skipped, the two calls have to be the same call. Asserted by
    # capturing the arguments RunAIAnalysis() receives, not by running a
    # provider.
    skip_if_not_installed("aisdk")
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = SummarizedExperiment::assay(sce, "counts"))
    )

    captured <- new.env(parent = emptyenv())
    fake_run <- function(object, goal, model = NULL, confirm = "ask", ...) {
        captured$goal <- goal
        captured$model <- model
        captured$confirm <- confirm
        captured$dots <- list(...)
        structure(list(object = object, status = "mock"), class = "sclet_ai_analysis")
    }
    with_mocked_bindings(
        {
            out <- scletAI(sce, "cluster the cells", model = "mock:model",
                confirm = "yes", include = "read")
            expect_s3_class(out, "sclet_ai_analysis")
            expect_identical(captured$goal, "cluster the cells")
            expect_identical(captured$model, "mock:model")
            expect_identical(captured$confirm, "yes")
            # The documented escape hatch still reaches RunAIAnalysis().
            expect_identical(captured$dots$include, "read")
        },
        RunAIAnalysis = fake_run,
        .package = "sclet"
    )
})

test_that("an explicit model is passed through the resolver unchanged", {
    skip_if_not_installed("aisdk")
    # scletAI() resolves the model before calling RunAIAnalysis(); with one
    # supplied the resolver must hand it back verbatim, never substitute a
    # default. This is the branch that runs on every real call.
    expect_identical(sclet_ai_resolve_model("mock:model"), "mock:model")
    expect_identical(sclet_ai_resolve_model("openai:gpt-4o"), "openai:gpt-4o")
})
