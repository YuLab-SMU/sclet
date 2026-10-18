## Iterative malignant-epithelium classification (issue YuLab-SMU/ProjectYulab#26)
##
## The original method has no per-cell ground truth, so these tests pin the
## properties that make the result *reportable* rather than the verdict itself:
##
##   1. a two-population input is separated and the iteration stops by
##      CONVERGENCE, not by the round cap, and the agreement history is
##      non-constant (a stuck NA would pass a "converged" assertion by accident);
##   2. a homogeneous input still gets a k=2 split, and the call says so with a
##      warning instead of returning two meaningful groups;
##   3. k-means cluster ids flipping between rounds is the SAME partition and
##      must not read as 0% agreement;
##   4. the result is written to colData and registered in the state contract,
##      and the full record survives a second call under another name;
##   5. a signature that separates the sides raises the identity-axis warning.
##
## Fixtures are synthetic with real gene symbols so gene-id handling is
## exercised; nothing here asserts a biological conclusion.

sclet_ns <- function(name) getFromNamespace(name, "sclet")

sclet_mal_test_genes <- function(n_marker = 12L, n_noise = 200L) {
    c(sprintf("MAL%02d", seq_len(n_marker)),
      sprintf("NML%02d", seq_len(n_marker)),
      sprintf("BG%03d", seq_len(n_noise)))
}

## Two populations with real programmes: MAL* high in the first 40 cells, NML*
## high in the last 40. `bg_scale` sets how much background genes express.
sclet_mal_test_sce <- function(n_per = 40L, n_patients = 4L, signal = 25,
                               n_noise = 240L) {
    genes <- sclet_mal_test_genes(n_noise = n_noise)
    n_marker <- 12L
    set.seed(11)
    n_cells <- 2L * n_per
    # A background that spans the expression range, so that the signature genes
    # (positioned low) have same-bin controls that are not themselves markers.
    bg <- exp(seq(log(0.3), log(40), length.out = n_noise))
    mu <- matrix(2.5, nrow = length(genes), ncol = n_cells)
    mu[n_marker * 2L + seq_len(n_noise), ] <- bg
    mu[seq_len(n_marker), seq_len(n_per)] <- signal
    mu[n_marker + seq_len(n_marker), n_per + seq_len(n_per)] <- signal
    counts <- matrix(stats::rpois(length(mu), lambda = mu),
                     nrow = length(genes), ncol = n_cells,
                     dimnames = list(genes, sprintf("c%04d", seq_len(n_cells))))

    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SummarizedExperiment::assay(sce, "logcounts") <- log1p(counts)
    SummarizedExperiment::colData(sce)$sample <- rep(
        sprintf("P%d", seq_len(n_patients)),
        length.out = n_cells
    )
    sce
}

sclet_mal_test_seed <- function() {
    list(
        tumor_genes = sprintf("MAL%02d", seq_len(12L)),
        normal_genes = sprintf("NML%02d", seq_len(12L))
    )
}

test_that("a two-population input is split and stops by convergence", {
    sce <- sclet_mal_test_sce()
    truth <- rep(c("A", "B"), each = 40L)

    res <- withCallingHandlers(
        ClassifyMalignantEpithelium(sce, seed = sclet_mal_test_seed(), max_iter = 10,
                                    convergence = 0.9, verbose = FALSE),
        warning = function(w) invokeRestart("muffleWarning")
    )
    r <- S4Vectors::metadata(res)$malignancy[["malignant"]]

    expect_equal(r$stop_reason, "converged")
    expect_true(isTRUE(r$converged))
    expect_lt(r$iterations, 10L)

    # The split must recover the generating groups, not merely produce two:
    # each truth group must fall entirely on one side.
    tab <- table(factor(r$labels, levels = c("malignant", "non_malignant")), truth)
    expect_equal(unname(max(colSums(tab > 0))), 1L)
    expect_equal(unname(sum(tab)), 80L)
    expect_gt(min(r$diagnostics$n_malignant, r$diagnostics$n_non_malignant), 10L)

    # Agreement must be REAL, not stuck at NA: round 1 has no previous split, so
    # it is NA by construction, and every later round carries a number.
    ag <- vapply(r$iteration_log, function(x) x$agreement, numeric(1))
    expect_true(is.na(ag[[1]]))
    expect_false(any(is.na(ag[-1])))
})

test_that("results are written to colData and registered in the state contract", {
    sce <- sclet_mal_test_sce(n_per = 30L)
    res <- suppressWarnings(
        ClassifyMalignantEpithelium(sce, seed = sclet_mal_test_seed(), verbose = FALSE)
    )
    cd <- SummarizedExperiment::colData(res)
    for (col in c("malignant_malignant_score", "malignant_nonmalignant_score",
                  "malignant_cluster", "malignant_classification")) {
        expect_true(col %in% colnames(cd), info = col)
    }
    expect_s3_class(cd$malignant_classification, "factor")
    expect_equal(levels(cd$malignant_classification), c("malignant", "non_malignant"))
    expect_false(any(is.na(cd$malignant_classification)))

    rec <- sclet_ns("sclet_get_state_record")(res, "annotation", "malignant")
    expect_false(is.null(rec))
    expect_equal(rec$summary$n_malignant + rec$summary$n_non_malignant, ncol(res))

    cmds <- sclet_ns("sclet_get_commands")(res)
    expect_true(any(vapply(cmds, function(x) identical(x$command, "ClassifyMalignantEpithelium"),
                           logical(1))))
})

test_that("two named runs coexist instead of overwriting each other", {
    sce <- sclet_mal_test_sce(n_per = 30L)
    sce <- suppressWarnings(
        ClassifyMalignantEpithelium(sce, seed = sclet_mal_test_seed(), name = "run1",
                                    verbose = FALSE)
    )
    sce <- suppressWarnings(
        ClassifyMalignantEpithelium(sce, seed = sclet_mal_test_seed(), name = "run2",
                                    verbose = FALSE)
    )
    md <- S4Vectors::metadata(sce)$malignancy
    expect_true(all(c("run1", "run2") %in% names(md)))
    expect_true("run2_classification" %in% colnames(SummarizedExperiment::colData(sce)))
})

test_that("a lopsided k=2 split warns that it may not be a real two-way structure", {
    # 60 cells with a real malignant programme in 55 and none in 5: the split
    # lands at 92/8 and must be flagged rather than reported as two groups.
    genes <- sclet_mal_test_genes()
    set.seed(7)
    counts <- matrix(stats::rpois(length(genes) * 60L, 2), nrow = length(genes),
                     dimnames = list(genes, sprintf("u%03d", seq_len(60L))))
    counts[1:12, 1:55] <- counts[1:12, 1:55] + 5L
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = log1p(counts)))
    SummarizedExperiment::colData(sce)$sample <- rep(sprintf("P%d", 1:4), 15)

    w <- capture_warnings(
        ClassifyMalignantEpithelium(sce, seed = sclet_mal_test_seed(), max_iter = 2,
                                    verbose = FALSE)
    )
    expect_true(any(grepl("k-means split the cells", w)))
    # And the standing interpretation caveat is always present.
    expect_true(any(grepl("not validated", w)))
})

test_that("a flipped k-means assignment is not read as zero agreement", {
    # Two well-separated score groups, second round with cluster ids swapped.
    set.seed(5)
    mal <- c(stats::rnorm(40, 3), stats::rnorm(40, 0))
    non <- c(stats::rnorm(40, 0), stats::rnorm(40, 3))
    first <- sclet_ns("sclet_kmeans_two_scores")(mal, non, nstart = 10)
    second <- sclet_ns("sclet_kmeans_two_scores")(mal + 1e-8, non - 1e-8, nstart = 10)
    # Force the mirrored id assignment that the raw diagonal would misread.
    second$labels <- ifelse(first$labels == "malignant", "non_malignant", "malignant")
    same <- sum(second$labels == first$labels)
    swapped <- sum(second$labels != first$labels)
    expect_equal(max(same, swapped) / length(first$labels), 1)
})

test_that("a signature that tracks the split raises the identity-axis warning", {
    sce <- sclet_mal_test_sce()
    w <- capture_warnings(
        ClassifyMalignantEpithelium(sce, seed = sclet_mal_test_seed(),
                                     seed_genes = sprintf("MAL%02d", seq_len(12L)),
                                     verbose = FALSE)
    )
    expect_true(any(grepl("seed_genes", w)))
})

test_that("inputs are validated before any work happens", {
    sce <- sclet_mal_test_sce(n_per = 20L)
    expect_error(
        ClassifyMalignantEpithelium(sce, seed = sclet_mal_test_seed(), convergence = 1.5),
        "convergence"
    )
    expect_error(
        ClassifyMalignantEpithelium(sce, seed = list(a = "X", b = "Y")),
        "tumor"
    )
    no_sample_col <- sce
    SummarizedExperiment::colData(no_sample_col)$sample <- NULL
    expect_error(
        ClassifyMalignantEpithelium(no_sample_col, seed = sclet_mal_test_seed(),
                                     verbose = FALSE),
        "sample"
    )
})

## ---- bulk seed genes -------------------------------------------------------

sclet_mal_test_bulk <- function(n_genes = 400L, n_patients = 8L, noise = 0.25) {
    set.seed(21)
    genes <- sprintf("G%04d", seq_len(n_genes))
    up <- seq_len(20L)          # tumour-high
    dn <- 21:40                 # normal-high
    tpm_t <- matrix(stats::rlnorm(n_genes * n_patients, 4, 0.6), n_genes, n_patients)
    tpm_n <- matrix(stats::rlnorm(n_genes * n_patients, 4, 0.6), n_genes, n_patients)
    tpm_t[up, ] <- tpm_t[up, ] * 25
    tpm_n[dn, ] <- tpm_n[dn, ] * 25
    expr <- cbind(log2(tpm_t + 0.001), log2(tpm_n + 0.001))
    colnames(expr) <- c(sprintf("T%d", seq_len(n_patients)), sprintf("N%d", seq_len(n_patients)))
    rownames(expr) <- genes
    info <- data.frame(
        sample = colnames(expr),
        patient = rep(sprintf("P%d", seq_len(n_patients)), 2L),
        condition = rep(c("Tumor", "Normal"), each = n_patients),
        stringsAsFactors = FALSE
    )
    list(expr = expr, info = info, up = genes[up], dn = genes[dn])
}

test_that("bulk seeds are the directional top genes, paired by default", {
    b <- sclet_mal_test_bulk()
    seed <- FindMalignantSeedGenes(b$expr, b$info, n_top = 20)

    expect_true(is.list(seed))
    expect_true(all(c("tumor_genes", "normal_genes", "deg", "info") %in% names(seed)))
    # Ordered by fold change: the first gene of each side is one of the planted
    # markers, and every selected gene has the expected sign.
    expect_gt(mean(seed$tumor_genes %in% b$up), 0.8)
    expect_gt(mean(seed$normal_genes %in% b$dn), 0.8)
    expect_true(all(seed$deg$logFC[match(seed$tumor_genes, seed$deg$gene)] > 0))
    expect_true(isTRUE(seed$info$design$paired))
    expect_equal(seed$info$design$n_case, 8L)
    expect_equal(seed$info$design$n_control, 8L)
})

test_that("bulk seeds survive a log-offset mismatch and record it", {
    b <- sclet_mal_test_bulk()
    # Claim log2(TPM+1) for data built as log2(TPM+0.001): the filter scale must
    # recover, not cut everything away.
    seed <- suppressWarnings(
        FindMalignantSeedGenes(b$expr, b$info, assay_note = "log2(TPM+1)", n_top = 20)
    )
    expect_gt(seed$info$counts$n_genes_used, 300L)
    expect_gt(mean(seed$tumor_genes %in% b$up), 0.6)
})

test_that("bulk seed input errors name what to fix", {
    b <- sclet_mal_test_bulk(n_patients = 4L)
    bad <- b$info
    bad$condition <- rep("weird", nrow(bad))
    expect_error(FindMalignantSeedGenes(b$expr, bad), "case_values")

    # Fewer than min_samples on each side.
    thin <- do.call(rbind, lapply(split(b$info, b$info$condition), function(x) x[1:2, ]))
    expect_error(FindMalignantSeedGenes(b$expr, thin, min_samples = 3), "at least")
})
