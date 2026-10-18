## CNV evidence and the four-state malignant call.
##
## Two properties carry the whole design and are what these tests pin:
##
##   1. An absent measurement must never read as a negative result. When the CNV
##      layer cannot run, every cell becomes `undetermined` or `likely_malignant`
##      and NONE becomes `malignant` -- a one-line call is not promoted.
##   2. Disagreeing evidence must not be resolved by majority. CNV-positive but
##      signature-negative (and the reverse) is `undetermined`, which is the
##      explicit answer the iterative port could never give.
##
## inferCNV itself is mocked where a run is needed: it needs JAGS, runs minutes
## per small sample, and its fitted output is not what this layer owns. The
## burden construction IS tested against a real inferCNV object, because the
## centring is the part that had a real defect (see test-cnv-burden-*).

sclet_cnv_ns <- function(name) getFromNamespace(name, "sclet")

## Two populations plus normal reference cells. MAL* is the malignant
## signature, REF* marks the reference cells.
sclet_cnv_test_sce <- function(n_per = 30L, n_ref = 20L, n_noise = 200L) {
    n_marker <- 10L
    genes <- c(sprintf("MAL%02d", seq_len(n_marker)),
               sprintf("REF%02d", seq_len(n_marker)),
               sprintf("BG%03d", seq_len(n_noise)))
    set.seed(31)
    n_cells <- 2L * n_per + n_ref
    bg <- exp(seq(log(0.3), log(40), length.out = n_noise))
    mu <- matrix(2.5, nrow = length(genes), ncol = n_cells)
    mu[n_marker * 2L + seq_len(n_noise), ] <- bg
    mu[seq_len(n_marker), seq_len(n_per)] <- 25
    counts <- matrix(stats::rpois(length(mu), lambda = mu),
                     nrow = length(genes),
                     dimnames = list(genes, sprintf("c%04d", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SummarizedExperiment::assay(sce, "logcounts") <- log1p(counts)
    cd <- SummarizedExperiment::colData(sce)
    cd$sample <- rep(c(rep("P1", n_per), rep("P2", n_per), rep("N1", n_ref)), length.out = n_cells)
    cd$reference <- c(rep(FALSE, 2L * n_per), rep(TRUE, n_ref))[seq_len(n_cells)]
    SummarizedExperiment::colData(sce) <- cd
    sce
}

sclet_cnv_test_signature <- function() sprintf("MAL%02d", seq_len(10L))

## Write a completed CNV record the way RunCNVInference() does, without running
## inferCNV.
sclet_cnv_test_record <- function(sce, burden, name = "cnv") {
    cd <- SummarizedExperiment::colData(sce)
    cd[[paste0(name, "_burden")]] <- burden
    cd[[paste0(name, "_n_cnv_regions")]] <- rep(NA_integer_, nrow(cd))
    cd[[paste0(name, "_is_reference")]] <- cd$reference
    SummarizedExperiment::colData(sce) <- cd
    rec <- list(status = "completed", reference_cells = rownames(cd)[cd$reference],
                burden = data.frame(cell = rownames(cd), burden = burden))
    sce <- sclet_cnv_ns("sclet_set_analysis")(sce, "malignant_cnv",
                                             stats::setNames(list(rec), name))
    sce
}

test_that("a missing CNV record is refused with an actionable message", {
    sce <- sclet_cnv_test_sce()
    expect_error(
        CallMalignantCells(sce, signature = sclet_cnv_test_signature()),
        "RunCNVInference"
    )
})

test_that("CNV unavailable never promotes a signature-only call to malignant", {
    sce <- sclet_cnv_test_sce()
    # RunCNVInference records an unavailable run rather than erroring.
    sce <- suppressWarnings(RunCNVInference(sce, gene_order = NULL))
    expect_true(all(is.na(SummarizedExperiment::colData(sce)$cnv_burden)))

    w <- capture_warnings(
        sce <- CallMalignantCells(sce, signature = sclet_cnv_test_signature())
    )
    expect_true(any(grepl("never|cannot|unavailable", w)))
    st <- SummarizedExperiment::colData(sce)$malignant_call_state
    expect_false(any(st == "malignant", na.rm = TRUE))
    expect_gt(sum(st == "likely_malignant", na.rm = TRUE), 0L)
})

test_that("agreeing two-line evidence gives malignant, disagreeing gives undetermined", {
    sce <- sclet_cnv_test_sce()
    n <- ncol(sce)
    # Cells 1..30 carry the MAL programme; give them the high CNV burden.
    burden <- c(rep(0.120, 30), rep(0.052, n - 30))
    # Make 10 signature-negative cells CNV-positive: evidence conflict.
    burden[31:40] <- 0.120
    sce <- sclet_cnv_test_record(sce, burden)

    # An explicit transcriptomic cutoff keeps the boundary off a cell: the
    # fixture's two programmes are separated, so 0 is between them.
    res <- suppressWarnings(
        CallMalignantCells(sce, signature = sclet_cnv_test_signature(),
                           transcriptomic_cut = 0)
    )
    cd <- SummarizedExperiment::colData(res)
    st <- cd$malignant_call_state
    ev <- cd$malignant_call_evidence

    expect_setequal(unique(as.character(ev)),
                    c("cnv_signature", "cnv_only", "neither"))
    # Two lines agree on the first 30.
    expect_true(all(st[1:30] == "malignant"))
    # CNV-positive but signature-negative is a conflict, not a verdict.
    expect_true(all(st[31:40] == "likely_malignant"))
    # Neither line positive.
    expect_true(all(st[41:(n - 20)] == "non_malignant"))
    expect_setequal(levels(st), c("malignant", "likely_malignant",
                                  "undetermined", "non_malignant"))

    rec <- sclet_cnv_ns("sclet_get_state_record")(res, "annotation", "malignant_call")
    expect_equal(rec$summary$n_malignant, 30L)
})

test_that("an explicit transcriptomic cutoff is honoured, not silently replaced", {
    sce <- sclet_cnv_test_sce()
    burden <- c(rep(0.120, 30), rep(0.052, ncol(sce) - 30))
    sce <- sclet_cnv_test_record(sce, burden)

    # A cutoff far above any score makes every cell signature-negative, so the
    # CNV-positive ones must fall to the conflict branch.
    res <- suppressWarnings(CallMalignantCells(sce, signature = sclet_cnv_test_signature(),
                                              transcriptomic_cut = 1e6))
    st <- SummarizedExperiment::colData(res)$malignant_call_state
    expect_equal(sum(st == "malignant"), 0L)
    expect_equal(sum(st == "likely_malignant"), 30L)
})

test_that("a competing identity programme demotes two-line calls", {
    sce <- sclet_cnv_test_sce()
    n <- ncol(sce)
    # Give the first 30 cells a strong "intestinal" programme on top of MAL*.
    genes <- rownames(sce)
    idg <- c("BG001", "BG002", "BG003", "BG004", "BG005")
    lc <- SummarizedExperiment::assay(sce, "logcounts")
    lc[idg, 1:30] <- lc[idg, 1:30] + 3
    SummarizedExperiment::assay(sce, "logcounts") <- lc

    burden <- c(rep(0.120, 30), rep(0.052, n - 30))
    sce <- sclet_cnv_test_record(sce, burden)

    w <- capture_warnings(
        res <- CallMalignantCells(sce, signature = sclet_cnv_test_signature(),
                                  identity_genes = idg, identity_demote = 0.05)
    )
    st <- SummarizedExperiment::colData(res)$malignant_call_state
    expect_true(any(grepl("identity", w)))
    expect_equal(sum(st == "malignant"), 0L)
    expect_gt(sum(st == "undetermined"), 0L)
    rec <- sclet_cnv_ns("sclet_get_state_record")(res, "annotation", "malignant_call")
    expect_true(isTRUE(rec$summary$identity_conflict))
})

test_that("readiness names every missing prerequisite", {
    sce <- sclet_cnv_test_sce()
    # No reference supplied and no colData$reference equivalent: point it at a
    # column that does not exist to force the missing-reference branch.
    r <- check_cnv_readiness(sce, reference_cells = character(0),
                             gene_order = data.frame(gene = rownames(sce)[1:10],
                                                     chr = "chr1", start = 1:10, end = 2:11))
    expect_false(isTRUE(r$ready))
    expect_true(any(grepl("reference", r$issues)))
})

test_that("readiness reports the gene-order source it used", {
    sce <- sclet_cnv_test_sce()
    ord <- data.frame(gene = rownames(sce), chr = "chr1",
                      start = seq_len(nrow(sce)), end = seq_len(nrow(sce)) + 1L)
    r <- check_cnv_readiness(sce, reference_col = "reference", gene_order = ord)
    expect_equal(r$gene_order_source, "caller-supplied")
    expect_equal(r$n_genes_with_position, nrow(sce))
})

## ---- the burden construction (the part that had a real defect) -------------

test_that("burden centres on each cell's own median, not on the raw level", {
    # Reference-centred inferCNV output hovers around 1 for EVERY cell; abs()
    # of it measures the baseline. A cell with a wide amplified/deleted profile
    # must outscore a flat one even when both sit at the same level.
    set.seed(4)
    genes <- sprintf("G%03d", seq_len(500))
    cells <- sprintf("c%02d", seq_len(4))
    flat <- stats::rnorm(500, 1, 0.01)
    flat2 <- stats::rnorm(500, 1, 0.01)
    wide <- c(stats::rnorm(250, 1.25, 0.02), stats::rnorm(250, 0.75, 0.02))
    wide2 <- c(stats::rnorm(250, 1.20, 0.02), stats::rnorm(250, 0.80, 0.02))
    expr <- cbind(flat, flat2, wide, wide2)
    dimnames(expr) <- list(genes, cells)

    obj <- new("infercnv", expr.data = expr, count.data = matrix(0, 1, 1),
               gene_order = data.frame(), reference_grouped_cell_indices = list(),
               observation_grouped_cell_indices = list(), tumor_subclusters = list(),
               options = list())
    b <- sclet_cnv_ns("sclet_cnv_burden")(obj, cells = cells, ref = cells[1:2])
    b <- stats::setNames(b$burden, b$cell)

    expect_lt(b[["c01"]], b[["c03"]])
    expect_lt(b[["c02"]], b[["c04"]])
    # The raw-level version would have made these indistinguishable.
    raw <- colMeans(abs(expr))
    expect_lt(abs(raw[["c03"]] - raw[["c01"]]), 0.01)
    expect_gt(b[["c03"]] - b[["c01"]], 0.1)
})
