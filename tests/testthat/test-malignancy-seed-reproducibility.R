test_that("the module score is reproducible under a fixed seed", {
    # The control pool is drawn with `sample()`, so two calls on the same input
    # returned different scores (measured 2026-10-22: max abs difference
    # 0.0348). The scoring is what the malignant/not call is built on, so an
    # unseeded score makes the call itself unreproducible.
    set.seed(11)
    expr <- matrix(abs(stats::rnorm(300 * 40)) + 0.1, nrow = 300,
                   dimnames = list(paste0("G", 1:300), paste0("C", 1:40)))
    genes <- paste0("G", 1:5)
    a <- sclet_add_module_score(expr, genes, n_ctrl = 20, seed = 7L)
    b <- sclet_add_module_score(expr, genes, n_ctrl = 20, seed = 7L)
    expect_equal(a, b)
    # A different seed has to be able to give a different answer, or the seed
    # is not reaching the sampling at all.
    c <- sclet_add_module_score(expr, genes, n_ctrl = 20, seed = 8L)
    expect_false(identical(a, c))
})

test_that("seeding the score does not move the caller's RNG stream", {
    set.seed(11)
    expr <- matrix(abs(stats::rnorm(300 * 40)) + 0.1, nrow = 300,
                   dimnames = list(paste0("G", 1:300), paste0("C", 1:40)))
    set.seed(5)
    before <- stats::runif(2)
    set.seed(5)
    invisible(sclet_add_module_score(expr, paste0("G", 1:5), n_ctrl = 20, seed = 99L))
    expect_equal(stats::runif(2), before)
})

test_that("the k-means split is reproducible under a fixed seed", {
    set.seed(3)
    x <- stats::rnorm(60)
    y <- stats::rnorm(60)
    a <- sclet_kmeans_two_scores(x, y, nstart = 10, seed = 4L)
    b <- sclet_kmeans_two_scores(x, y, nstart = 10, seed = 4L)
    expect_equal(a$labels, b$labels)
})

test_that("ClassifyMalignantEpithelium separates the gene signature from the RNG seed", {
    # The trap: the argument named `seed` is a GENE SIGNATURE (the output of
    # FindMalignantSeedGenes), not an RNG seed. A caller who read `seed =` as
    # "make this reproducible" was wrong, and there was no other seed to use.
    f <- formals(ClassifyMalignantEpithelium)
    expect_true("seed" %in% names(f))
    expect_true("rng_seed" %in% names(f))
    expect_false(is.null(f$rng_seed))
})

test_that("ClassifyMalignantEpithelium is reproducible under a fixed rng_seed", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 60, ngenes = 200)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 50)
    sce <- RunPCA(sce, ncomponents = 5)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce, resolution = 0.5)
    cd <- SummarizedExperiment::colData(sce)
    cd$sample <- rep(paste0("s", 1:6), length.out = ncol(sce))
    SummarizedExperiment::colData(sce) <- cd
    signature <- list(malignant = rownames(sce)[1:5], normal = rownames(sce)[6:12])
    run <- function(seed) {
        suppressWarnings(ClassifyMalignantEpithelium(
            sce, seed = signature, sample_col = "sample",
            rng_seed = seed, verbose = FALSE
        ))
    }
    a <- run(1L)
    b <- run(1L)
    expect_equal(a$labels, b$labels)
})

test_that("CallPhenotypeCells takes an rng_seed too", {
    f <- formals(CallPhenotypeCells)
    expect_true("rng_seed" %in% names(f))
    expect_false(is.null(f$rng_seed))
})
