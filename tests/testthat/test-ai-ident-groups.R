# Tests for identity group code exposure and resolution
# Spec: .dev/spec-ai-ident-groups.md

test_that("sclet_ai_ident_group_codes returns sorted deterministic mapping", {
    sce <- scuttle::mockSCE()
    sce$cluster <- sample(c("treat1", "treat2"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"

    mapping <- sclet_ai_ident_group_codes(sce)
    expect_equal(mapping$codes, c("group_1", "group_2"))
    expect_equal(sort(mapping$labels), c("treat1", "treat2"))
    expect_equal(sum(mapping$sizes), ncol(sce))

    tbl <- table(Idents(sce))
    expect_equal(mapping$sizes[mapping$labels == "treat1"], as.integer(tbl["treat1"]))
    expect_equal(mapping$sizes[mapping$labels == "treat2"], as.integer(tbl["treat2"]))

    sce_rev <- sce[, rev(seq_len(ncol(sce)))]
    mapping_rev <- sclet_ai_ident_group_codes(sce_rev)
    expect_identical(mapping, mapping_rev)
})

test_that("no active identity returns empty result without error", {
    sce <- scuttle::mockSCE()
    mapping <- sclet_ai_ident_group_codes(sce)
    expect_equal(mapping$n_groups, 0L)
    expect_equal(length(mapping$codes), 0L)
    expect_null(mapping$col)

    ctx <- sclet_ai_context(sce, target = "planning", detail = "full")
    expect_equal(ctx$capabilities$identity_levels$n_levels, 0L)
    expect_equal(length(ctx$capabilities$identity_levels$codes), 0L)
})

test_that("context exposes codes and counts but no real labels", {
    sce <- scuttle::mockSCE()
    set.seed(42)
    # The cluster labels use names that appear NOWHERE else on the object.
    # They used to be "treat1"/"treat2", which also happen to be the levels of
    # mockSCE()'s Treatment column -- so once `capabilities$contrast_levels`
    # legitimately exposed design levels (M12-a, 2026-10-21), the assertion
    # "the string treat1 appears nowhere" started failing on the DESIGN value,
    # not on a leaked cluster label. This test is about the clustering, so its
    # fixtures must not collide with a design column.
    sce$cluster <- sample(c("cladeA", "cladeB"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"

    ctx <- sclet_ai_context(sce, target = "planning", detail = "full")
    codes_list <- ctx$capabilities$identity_levels$codes
    expect_true(length(codes_list) >= 2L)
    code_values <- vapply(codes_list, function(x) x$code, character(1L))
    expect_true("group_1" %in% code_values)
    expect_true("group_2" %in% code_values)

    n_cells_values <- vapply(codes_list, function(x) x$n_cells, integer(1L))
    expect_true(all(n_cells_values > 0L))

    serialized <- paste(utils::capture.output(str(ctx)), collapse = " ")
    expect_false(grepl("cladeA", serialized, fixed = TRUE))
    expect_false(grepl("cladeB", serialized, fixed = TRUE))
})

test_that("contrast_levels never re-exposes the clustering column under real labels", {
    # The identity column holds the labels the analysis is meant to discover;
    # identity_levels anonymises them, so contrast_levels must not undo that by
    # listing the same column with its real values.
    sce <- scuttle::mockSCE()
    set.seed(42)
    sce$cluster <- sample(c("cladeA", "cladeB"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"

    ctx <- sclet_ai_context(sce, target = "planning", detail = "full")
    columns <- vapply(ctx$capabilities$contrast_levels$columns,
        function(x) x$column, character(1))
    expect_false("cluster" %in% columns)

    serialized <- paste(utils::capture.output(str(ctx$capabilities$contrast_levels)),
        collapse = " ")
    expect_false(grepl("cladeA", serialized, fixed = TRUE))

    # a genuine design column still comes through: that is the whole point
    expect_true("Treatment" %in% columns)
})

test_that("sanitized context preserves identity_levels and remains label-free", {
    sce <- scuttle::mockSCE()
    set.seed(42)
    # Same fixture rule as above: cluster labels must not collide with a
    # design column's levels, or the "label-free" assertion fails on a design
    # value that is supposed to be present.
    sce$cluster <- sample(c("cladeA", "cladeB"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"

    ctx <- sclet_ai_context(sce, target = "planning", detail = "full")
    sanitized <- sclet_ai_sanitize_context(sce, ctx, privacy = "standard")

    sanitized_levels <- sanitized$payload$capabilities$identity_levels
    expect_false(is.null(sanitized_levels))
    expect_equal(sanitized_levels$n_levels, ctx$capabilities$identity_levels$n_levels)
    expect_equal(length(sanitized_levels$codes), length(ctx$capabilities$identity_levels$codes))

    serialized <- paste(utils::capture.output(str(sanitized$payload)), collapse = " ")
    expect_false(grepl("cladeA", serialized, fixed = TRUE))
    expect_false(grepl("cladeB", serialized, fixed = TRUE))

    redacted <- sanitized$policy$redacted_fields
    identity_redacted <- redacted[grepl("identity_levels", redacted)]
    expect_equal(length(identity_redacted), 0L)

    message("=== Sanitized identity_levels subtree ===")
    print(sanitized_levels)
    message("=== Redacted fields matching identity_levels ===")
    print(identity_redacted)
})

test_that("small levels are omitted not redacted", {
    sce <- scuttle::mockSCE()
    set.seed(123)
    labels <- c(rep("large", 191L), rep("small", 9L))
    sce$cluster <- labels
    ActiveIdent(sce) <- "cluster"

    ctx <- sclet_ai_context(sce, target = "planning", detail = "full")
    levels <- ctx$capabilities$identity_levels
    expect_equal(levels$omitted_below_threshold, 1L)
    expect_false(levels$all_levels_available)
    expect_equal(levels$n_levels, 1L)

    code_values <- vapply(levels$codes, function(x) x$code, character(1L))
    expect_equal(length(code_values), 1L)
    expect_true("group_1" %in% code_values)
})

test_that("resolver maps codes and passes through other values", {
    sce <- scuttle::mockSCE()
    set.seed(42)
    sce$cluster <- sample(c("treat1", "treat2"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"

    mapping <- sclet_ai_ident_group_codes(sce)
    label1 <- mapping$labels[1L]

    resolved <- sclet_ai_resolve_ident_value(sce, "group_1")
    expect_equal(resolved, label1)

    resolved_vec <- sclet_ai_resolve_ident_value(sce, c("group_1", "group_2"))
    expect_equal(length(resolved_vec), 2L)
    expect_equal(resolved_vec[1L], mapping$labels[1L])
    expect_equal(resolved_vec[2L], mapping$labels[2L])

    expect_equal(sclet_ai_resolve_ident_value(sce, "treat1"), "treat1")
    expect_equal(sclet_ai_resolve_ident_value(sce, "group_99"), "group_99")
    expect_null(sclet_ai_resolve_ident_value(sce, NULL))
    expect_equal(sclet_ai_resolve_ident_value(sce, character(0)), character(0))
})

test_that("plan using code passes validation", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE()
    set.seed(42)
    sce$cluster <- sample(c("treat1", "treat2"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"

    registry <- AIDefaultExecutionRegistry(sce, include = "all")

    plan_code <- new_sclet_ai_plan(
        task = "test",
        actions = list(
            list(action = "run_de_test", params = list(ident.1 = "group_1"))
        ),
        requires_confirmation = TRUE
    )
    val_code <- ValidateAIPlan(plan_code, object = sce, registry = registry)
    expect_true(val_code$valid)

    mapping <- sclet_ai_ident_group_codes(sce)
    label1 <- mapping$labels[1L]
    plan_real <- new_sclet_ai_plan(
        task = "test",
        actions = list(
            list(action = "run_de_test", params = list(ident.1 = label1))
        ),
        requires_confirmation = TRUE
    )
    val_real <- ValidateAIPlan(plan_real, object = sce, registry = registry)
    expect_true(val_real$valid)
})

test_that("execution resolves code to real label", {
    # Fixture note: scuttle::mockSCE() with hand-sampled two-level labels makes
    # RunDEtest() fail with "duplicate 'row.names' are not allowed" on this
    # sclet version. That reproduces with a plain REAL label and no code
    # mapping involved at all, so it is a pre-existing DE-path issue unrelated
    # to this feature (recorded as a finding in sclet-benchmark-ai's
    # docs/action-coverage.md). Use the same clustered fixture that the
    # existing run_de_test coverage in test-ai-domain-adapter-e2e.R uses,
    # which is known to execute cleanly.
    set.seed(1L)
    counts <- matrix(rpois(50L * 60L, lambda = 5L), nrow = 50L, ncol = 60L,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(60L))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    # A fixed seed for the fixture, and one whose partition has no singleton
    # cluster: Louvain is seeded now (FindClusters default seed = 1), and the
    # default partition of THIS random fixture contains a one-cell cluster, which
    # the marker test correctly refuses. The product default is unchanged; the
    # fixture picks a partition the DE prerequisite accepts.
    sce <- FindClusters(sce, resolution = 1.2, seed = 42L)

    registry <- AIDefaultExecutionRegistry(sce, include = "all")

    plan <- new_sclet_ai_plan(
        task = "test",
        actions = list(
            list(action = "run_de_test", params = list(ident.1 = "group_1"))
        ),
        requires_confirmation = TRUE
    )

    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)

    # NOTE: no tryCatch(skip) around this. An earlier version wrapped the call
    # in tryCatch(error = function(e) NULL) and then skip_if(is.null(result)),
    # which silently skipped this scenario -- the single most load-bearing one
    # in the file -- when the call errored for an unrelated reason (arg
    # matching plus a missing confirmation token). A failure here must fail,
    # not skip.
    executed <- ExecuteAIPlan(
        sce, plan, registry,
        validation = validation,
        dry_run = FALSE,
        confirmation = validation$confirmation_token
    )
    expect_equal(executed$status, "completed")

    mapping <- sclet_ai_ident_group_codes(sce)
    label1 <- mapping$labels[mapping$codes == "group_1"]
    expect_length(label1, 1L)

    evidence <- sclet_ai_evidence_get_all(executed$object)
    expect_true(
        any(grepl(paste0("ev:detest_", label1), names(evidence), fixed = TRUE)),
        info = paste0(
            "handler must resolve group_1 to the real label (", label1,
            ") before calling RunDEtest(); evidence ids were: ",
            paste(names(evidence), collapse = ", ")
        )
    )
})

test_that("out-of-range code produces normal error", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE()
    set.seed(42)
    sce$cluster <- sample(c("treat1", "treat2"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"

    registry <- AIDefaultExecutionRegistry(sce, include = "all")

    plan <- new_sclet_ai_plan(
        task = "test",
        actions = list(
            list(action = "run_de_test", params = list(ident.1 = "group_99"))
        ),
        requires_confirmation = TRUE
    )
    val <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(val$valid)
    expect_true(any(grepl("unknown group labels", val$errors, ignore.case = TRUE)))
})

test_that("start_cluster resolves in run_trajectory", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("slingshot")
    sce <- scuttle::mockSCE()
    SummarizedExperiment::assay(sce, "logcounts") <- SummarizedExperiment::assay(sce, "counts")
    set.seed(42)
    sce$cluster <- sample(c("treat1", "treat2"), ncol(sce), replace = TRUE)
    ActiveIdent(sce) <- "cluster"
    sce <- runPCA(sce)
    sce <- runUMAP(sce)

    registry <- AIDefaultExecutionRegistry(sce, include = "all")

    plan <- new_sclet_ai_plan(
        task = "test",
        actions = list(
            list(action = "run_trajectory", params = list(
                group = "cluster",
                start_cluster = "group_1",
                reduction = "UMAP"
            ))
        ),
        requires_confirmation = TRUE
    )
    val <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(val$valid)
})

## Scenario 11 ("the full suite still passes") is intentionally not a test
## here. The spec requires it to be verified by running the suite itself
## (make check), not by nesting a whole-suite run inside this file -- see
## .dev/spec-ai-plan-repair.md section 4 scenario 9 for the measured reasons
## (nested runs were ~20x slower and downgraded file-level errors to skips).
##
## It was previously written as test_that(..., skip("documentation
## placeholder")), which added a spurious skip to every run without testing
## anything. A comment is the honest representation.
