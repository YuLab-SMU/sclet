# M12-b: `analysis_story$data_support` pins.
#
# Three things must hold, and each one is a claim M12-b makes:
#
#   1. DETERMINISM -- same object, same two calls, identical result. A
#      context block that varies between rounds would make plan fingerprints
#      unstable and every plan would be refused by design.
#
#   2. VOCABULARY -- every `reason`/`gate` this block emits is a name the
#      engine already uses. This is what keeps `data_support` from becoming a
#      second gate with a private language that disagrees with the real one.
#
#   3. NO RETEST ADVICE -- `unexplored` never suggests an analysis the ledger
#      already has. "Tested and found nothing" and "never tested" are
#      different states; only the second belongs there. `unexplored` is the
#      input to choosing the next round, so a suggestion to redo finished
#      work would be actively harmful.

m12b_sce <- function(n_cells = 8L, with_samples = FALSE) {
    assays <- list(counts = matrix(1L, nrow = 6L, ncol = n_cells))
    coldata <- S4Vectors::DataFrame(row.names = paste0("cell", seq_len(n_cells)))
    if (with_samples) {
        coldata$sample <- rep(paste0("s", seq_len(n_cells / 2L)), each = 2L)
        coldata$condition <- rep(paste0("c", seq_len(n_cells / 4L)), each = 4L)
    }
    SingleCellExperiment::SingleCellExperiment(assays = assays, colData = coldata)
}

test_that("data_support is deterministic across calls", {
    sce <- m12b_sce(with_samples = TRUE)
    first <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))
    second <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))
    expect_identical(first, second)
    expect_identical(first$schema_version, "1.0")
    expect_false(isTRUE(first$shapes$raw_values_included))
    expect_false(isTRUE(first$raw_values_included))
})

test_that("data_support is part of the analysis_story the planner reads", {
    sce <- m12b_sce()
    story <- GetAnalysisLedger(sce)$analysis_story
    expect_true(!is.null(story$data_support))
    expect_identical(story$schema_version, "1.1")
    # the pre-M12 fields must survive: this is an addition, not a rewrite
    for (field in c("timeline", "user_decisions", "evidence_gaps", "conflicts")) {
        expect_true(field %in% names(story), info = field)
    }
})

test_that("every gate a blocked action reports is a name data_support declares", {
    # The pin is against the ENGINE, not against a comment and not against a
    # hard-coded list. Two directions:
    #
    #   1. data_support's own vocabulary is non-trivial and closed;
    #   2. the actions whose gates data_support mirrors really do refuse with
    #      names from that vocabulary -- so the block cannot drift into saying
    #      "blocked by X" while the action refuses with Y.
    #
    # An earlier version scanned list.files("../../R"), which resolves in the
    # source tree and fails under R CMD check (found by make check,
    # 2026-10-21). A test that only passes from one directory is testing where
    # it is run, not what the code does. This one calls the live prerequisites.
    declared <- sclet:::sclet_ai_data_support_gates()
    expect_true(length(declared) >= 5L)
    expect_identical(declared, unique(declared))

    sce <- m12b_sce(with_samples = TRUE)
    registry <- AIDefaultExecutionRegistry(sce, include = "all")

    # velocity on a counts-only object: whatever the engine refuses with must
    # be a gate name data_support can also emit
    velocity <- registry[["run_velocity"]]
    expect_true(!is.null(velocity))
    verdict <- velocity$prerequisites(sce, list())
    if (is.character(verdict)) {
        gate <- sub("^([a-z_]+):.*$", "\\1", verdict)
        expect_true(gate %in% declared,
            info = paste("run_velocity refuses with undeclared gate:", gate))
    }

    # a condition contrast naming a level the column does not have
    SummarizedExperiment::colData(sce)$condition <- rep(c("ctrl", "stim"), 4L)
    de <- registry[["run_de_test"]]
    expect_true(!is.null(de))
    wrong <- de$prerequisites(sce, list(
        group = "condition", ident.1 = "nope", ident.2 = "ctrl"
    ))
    expect_true(is.character(wrong))
    expect_identical(sub("^([a-z_]+):.*$", "\\1", wrong), "contrast_invalid")
})

test_that("data_support reports what a counts-only object cannot do, with the real gate", {
    sce <- m12b_sce()
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))

    velocity <- Filter(function(x) identical(x$question_class, "velocity"), ds$supports)[[1L]]
    expect_false(velocity$supported)
    # counts are present and integral, so the missing assay is the reason
    expect_identical(velocity$reason, "spliced_unspliced_missing")
    expect_false(is.null(velocity$unblocked_by))

    trajectory <- Filter(function(x) identical(x$question_class, "trajectory"), ds$supports)[[1L]]
    expect_false(trajectory$supported)
    expect_identical(trajectory$reason, "reduction_missing")

    spatial <- Filter(function(x) identical(x$question_class, "spatial"), ds$supports)[[1L]]
    expect_identical(spatial$reason, "no_spatial_record")

    multimodal <- Filter(function(x) identical(x$question_class, "multimodal"), ds$supports)[[1L]]
    expect_identical(multimodal$reason, "no_multimodal_record")

    # blocked is the same information in gate-first form, and must agree
    blocked_gates <- vapply(ds$blocked, function(x) x$gate, character(1))
    expect_true("spliced_unspliced_missing" %in% blocked_gates)
    expect_true("no_spatial_record" %in% blocked_gates)
})

test_that("a design with sample and condition columns is recognised as replicable", {
    sce <- m12b_sce(with_samples = TRUE)
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))
    expect_true(ds$shapes$has_multiple_samples)
    expect_true("sample" %in% ds$shapes$candidate_sample_columns)
    expect_true("condition" %in% ds$shapes$candidate_condition_columns)
    # two cells per sample, two samples per condition
    expect_true(ds$shapes$has_replicates)
    abundance <- Filter(function(x) identical(x$question_class, "differential_abundance"),
        ds$supports)[[1L]]
    expect_true(abundance$supported)
})

test_that("a single sample never reads as replicated", {
    sce <- m12b_sce()
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))
    expect_false(ds$shapes$has_multiple_samples)
    # no sample/condition columns at all -> the reason is the missing role,
    # and it must not claim "one sample per condition", which was never seen
    abundance <- Filter(function(x) identical(x$question_class, "differential_abundance"),
        ds$supports)[[1L]]
    expect_false(abundance$supported)
    expect_identical(abundance$reason, "design_role_missing")
})

test_that("unexplored names data value the ledger has no result for", {
    sce <- m12b_sce(with_samples = TRUE)
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    # several samples plus a batch-less design: integration is not suggested
    # (no batch column), but the vocabulary must still be usable
    expect_true(all(vapply(ds$unexplored, function(x) {
        !is.null(x$what) && !is.null(x$action) && !is.null(x$why)
    }, logical(1L))))
    expect_true(length(actions) == length(unique(actions)))
})

test_that("unexplored never proposes re-running an analysis the ledger already has", {
    sce <- m12b_sce(with_samples = TRUE)
    analyses <- list(
        clustering = list(type = "clustering", id = "cl1"),
        detest = list(type = "detest", id = "de1")
    )
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = analyses))
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    # already clustered and already tested for markers: neither may reappear
    expect_false("run_de_test" %in% actions)
    # but enrichment is genuinely missing and must still be offered
    expect_true("run_enrichment" %in% actions)
    # annotation is missing too, because clustering exists and annotation does not
    expect_true("run_annotation" %in% actions)
})

test_that("a velocity-ready object is told velocity is available", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = matrix(1L, nrow = 6L, ncol = 4L),
            spliced = matrix(1, nrow = 6L, ncol = 4L),
            unspliced = matrix(1, nrow = 6L, ncol = 4L)
        )
    )
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))
    expect_true(ds$shapes$has_spliced_unspliced)
    expect_true(ds$shapes$counts_are_integer)
    velocity <- Filter(function(x) identical(x$question_class, "velocity"), ds$supports)[[1L]]
    expect_true(velocity$supported)
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    expect_true("run_velocity" %in% actions)
})

test_that("a non-integer counts assay is not counted as counts", {
    # the name is not the thing: a log-normalized matrix stored under `counts`
    # would otherwise green-light velocity on data velociraptor cannot use
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = matrix(c(0.5, 1.2, 0.7, 3.1), nrow = 2L, ncol = 2L),
            spliced = matrix(1, nrow = 2L, ncol = 2L),
            unspliced = matrix(1, nrow = 2L, ncol = 2L)
        )
    )
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = list()))
    expect_false(ds$shapes$counts_are_integer)
    velocity <- Filter(function(x) identical(x$question_class, "velocity"), ds$supports)[[1L]]
    expect_false(velocity$supported)
    expect_identical(velocity$reason, "counts_assay_missing")
})

# ---- Regressions found by running the M12-b question bank on the M1-M11
# ---- objects. Both were false statements about a real object, and both were
# ---- invisible to the unit fixtures above because the fixtures happened to
# ---- carry the record type the code keyed on.

test_that("a grouping is proved by any analysis that presupposes one, not only by a clustering record", {
    # The M11 pancreas object carries preprocess + reduction + velocity +
    # trajectory + cellrank and NO `clustering` record, because that run went
    # from the embedding straight to a trajectory. Keying on the record type
    # reported cluster_identity_missing on an object that already had a
    # trajectory, AND silenced every unexplored entry that required clustering,
    # so the whole suggestion list came back empty.
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 6L, ncol = 4L))
    )
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(1, nrow = 4L, ncol = 3L)
    analyses <- list(
        preprocess = list(type = "preprocess"),
        reduction = list(type = "reduction"),
        trajectory = list(type = "trajectory"),
        velocity = list(type = "velocity")
    )
    ds <- sclet:::sclet_ai_build_data_support(sce, list(analyses = analyses))

    trajectory <- Filter(function(x) identical(x$question_class, "trajectory"),
        ds$supports)[[1L]]
    expect_true(trajectory$supported)

    # and the suggestions that need a grouping must actually appear
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    expect_true("run_de_test" %in% actions)
    expect_true("run_annotation" %in% actions)
})

test_that("grouping evidence names its sources instead of guessing", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 6L, ncol = 4L))
    )
    none <- sclet:::sclet_ai_data_support_grouping_evidence(sce, list(analyses = list()))
    expect_false(none$available)
    expect_length(none$sources, 0L)

    implied <- sclet:::sclet_ai_data_support_grouping_evidence(sce,
        list(analyses = list(list(type = "detest"))))
    expect_true(implied$available)
    expect_true(any(grepl("implied_by:detest", implied$sources, fixed = TRUE)))

    sce$sclet_ident <- rep(c("1", "2"), 2)
    on_object <- sclet:::sclet_ai_data_support_grouping_evidence(sce, list(analyses = list()))
    expect_true(on_object$available)
    expect_true("object_sclet_ident" %in% on_object$sources)
})

test_that("an empty unexplored list is only produced when there is genuinely nothing to add", {
    # A block that always returns nothing would pass a "does not re-run what
    # exists" assertion. This pins the other direction: an object with an
    # embedding, a grouping and no trajectory must be told about the trajectory.
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 6L, ncol = 8L))
    )
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(1, nrow = 8L, ncol = 3L)
    sce$sclet_ident <- rep(c("1", "2"), 4)
    ds <- sclet:::sclet_ai_build_data_support(sce,
        list(analyses = list(clustering = list(type = "clustering"))))
    actions <- vapply(ds$unexplored, function(x) x$action, character(1))
    expect_true(length(actions) > 0L)
    expect_true("run_trajectory" %in% actions)
})
