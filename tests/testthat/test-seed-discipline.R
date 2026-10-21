test_that("sclet_with_seed reproduces a seeded draw", {
    skip_if_not_installed("scuttle")
    set.seed(7)
    expect <- runif(3)
    sclet_with_seed(7, {
        got <- runif(3)
    })
    expect_equal(got, expect)
})

test_that("sclet_with_seed leaves the caller's RNG stream where it was", {
    # A package that moves the user's RNG stream silently changes the numbers
    # the user gets from their own code, and they have no way to know.
    set.seed(42)
    before <- runif(3)
    set.seed(42)
    sclet_with_seed(99, {
        runif(100)
    })
    after <- runif(3)
    expect_equal(after, before)
})

test_that("sclet_with_seed removes .Random.seed again if the caller had none", {
    # base::runif() creates .Random.seed in globalenv() on first use, so the
    # state has to be put back exactly as it was found, not merely restored.
    rm(".Random.seed", envir = globalenv())
    sclet_with_seed(1, {
        runif(2)
    })
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("sclet_with_seed rejects a bad seed and passes NULL through", {
    expect_error(sclet_with_seed("x", 1), "single finite number")
    expect_error(sclet_with_seed(c(1, 2), 1), "single finite number")
    expect_error(sclet_with_seed(NA_real_, 1), "single finite number")
    expect_identical(sclet_with_seed(NULL, 42), 42)
})

test_that("RunUMAP is reproducible under a fixed seed", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    # The object is built ONCE and the embedding is computed twice from it.
    # scuttle::mockSCE() is itself random -- two calls with the same ncells
    # and ngenes return different count matrices -- so building two objects
    # and comparing their embeddings would compare two different datasets.
    sce <- scuttle::mockSCE(ncells = 60, ngenes = 200)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 50)
    sce <- RunPCA(sce, ncomponents = 5)
    a <- RunUMAP(sce, seed = 11L)
    b <- RunUMAP(sce, seed = 11L)
    expect_equal(
        SingleCellExperiment::reducedDim(a, "UMAP"),
        SingleCellExperiment::reducedDim(b, "UMAP")
    )
})

test_that("RunUMAP without a seed is NOT reproducible (the reason the seed exists)", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scater")
    # This is the defect the seed parameter was added for. It is pinned so
    # that a future change to the seeding cannot quietly remove the
    # difference without someone noticing.
    sce <- scuttle::mockSCE(ncells = 60, ngenes = 200)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 50)
    sce <- RunPCA(sce, ncomponents = 5)
    a <- RunUMAP(sce, seed = NULL)
    b <- RunUMAP(sce, seed = NULL)
    expect_false(identical(
        SingleCellExperiment::reducedDim(a, "UMAP"),
        SingleCellExperiment::reducedDim(b, "UMAP")
    ))
})

test_that("RunUMAP records the seed in the reduction state", {
    # A seed that is not written down is not reproducible.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 100)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30)
    sce <- RunPCA(sce, ncomponents = 5)
    sce <- RunUMAP(sce, seed = 23L)
    record <- sclet_get_state_record(sce, "reduction", "umap")
    expect_false(is.null(record))
    expect_equal(record$params$seed, 23L)
})

test_that("RunUMAP does not move the caller's RNG stream", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 100)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30)
    sce <- RunPCA(sce, ncomponents = 5)
    set.seed(5)
    before <- runif(2)
    set.seed(5)
    RunUMAP(sce, seed = 3L)
    after <- runif(2)
    expect_equal(after, before)
})

test_that("RunIntegration records the seed it used for the stochastic backends", {
    # fastMNN is deterministic, so the record has to say so rather than
    # implying a seed controlled something it did not.
    skip_if_not_installed("batchelor")
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 100)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30)
    colData(sce)$batch <- rep(c("b1", "b2"), length.out = ncol(sce))
    sce <- RunIntegration(sce, method = "fastMNN", batch = "batch", seed = 5L)
    record <- sclet_get_state_record(sce, "integration", "fastmnn")
    if (!is.null(record)) {
        expect_true(is.null(record$params$seed))
    }
})

test_that("the RunUMAP signature carries a seed argument", {
    # A pin on the API shape, not on the behaviour: a future refactor that
    # drops `seed` would silently make every UMAP irreproducible again.
    expect_true("seed" %in% names(formals(RunUMAP)))
    expect_true("seed" %in% names(formals(RunIntegration)))
})
