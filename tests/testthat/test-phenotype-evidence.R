## The phenotype evidence registry and the generic call.
##
## What these tests pin is the rule that makes the framework worth having:
##
##   1. A line counts only if it is INDEPENDENT of the signature. A provider
##      computed from the same expression matrix cannot turn a one-line call into
##      a two-line one, no matter how many providers are registered.
##   2. An unavailable provider is reported and does not count. It never becomes
##      a silent "negative".
##   3. Conflict is not resolved by majority: yes + no is undetermined.
##   4. A user can register their own provider without subclassing anything.
##
## Point 1 is the property the iterative method in R/malignancy.R lacks: it
## re-derives its own signature and counts the result as evidence.

sclet_ph_ns <- function(name) getFromNamespace(name, "sclet")

sclet_ph_test_sce <- function(n_per = 30L, n_noise = 200L) {
    n_marker <- 10L
    genes <- c(sprintf("MAL%02d", seq_len(n_marker)),
               sprintf("NML%02d", seq_len(n_marker)),
               sprintf("BG%03d", seq_len(n_noise)))
    set.seed(41)
    n_cells <- 2L * n_per
    bg <- exp(seq(log(0.3), log(40), length.out = n_noise))
    mu <- matrix(2.5, nrow = length(genes), ncol = n_cells)
    mu[n_marker * 2L + seq_len(n_noise), ] <- bg
    mu[seq_len(n_marker), seq_len(n_per)] <- 25
    counts <- matrix(stats::rpois(length(mu), lambda = mu),
                     nrow = length(genes),
                     dimnames = list(genes, sprintf("c%04d", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SummarizedExperiment::assay(sce, "logcounts") <- log1p(counts)
    SummarizedExperiment::colData(sce)$sample <- rep(c("P1", "P2"), each = n_per)
    sce
}

sclet_ph_test_sig <- function() sprintf("MAL%02d", seq_len(10L))

## A provider that reports a fixed state for every cell. `dep` controls whether
## it claims independence.
sclet_ph_stub_provider <- function(nm, state, dep = "independent") {
    list(
        name = nm,
        label = paste("Stub:", nm),
        dependence = dep,
        available = function(object, ...) list(available = TRUE, issues = character()),
        decide = function(object, ...) list(state = rep(state, ncol(object)))
    )
}

with_provider <- function(provider, code) {
    PhenotypeEvidenceProviders("register", provider = provider)
    on.exit(PhenotypeEvidenceProviders("unregister", name = provider$name), add = TRUE)
    force(code)
}

test_that("the built-in registry lists its providers and validates lookups", {
    nms <- PhenotypeEvidenceProviders("list")
    expect_true(all(c("cnv", "module") %in% nms))
    expect_error(PhenotypeEvidenceProviders("get", name = "nope"), "No evidence provider")
    expect_error(CallPhenotypeCells(sclet_ph_test_sce(), signature = sclet_ph_test_sig(),
                                    evidence = "nope"),
                 "Unknown evidence provider")
})

test_that("an unavailable provider is reported and never counts as negative", {
    sce <- sclet_ph_test_sce()
    # No CNV record exists, so the cnv provider cannot run.
    w <- capture_warnings(
        res <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(), evidence = "cnv")
    )
    expect_true(any(grepl("unavailable", w)))
    cd <- SummarizedExperiment::colData(res)
    expect_true(all(cd$phenotype_call_cnv_state == "unavailable"))
    # Signature-positive cells are 'probable' (one line), never 'supported'.
    expect_false(any(cd$phenotype_call_call == "supported"))
    rec <- sclet_ph_ns("sclet_get_state_record")(res, "annotation", "phenotype_call")
    expect_equal(rec$summary$unavailable_providers, "cnv")
})

test_that("a shared-matrix provider cannot turn one line into two", {
    sce <- sclet_ph_test_sce()
    # Three providers, all positive, all sharing the signature's matrix. They
    # must be registered together: `with_provider` unregisters on exit, so
    # nesting would leave only the last one live.
    nms <- c("shared_a", "shared_b", "shared_c")
    for (nm in nms) {
        PhenotypeEvidenceProviders("register",
                                   provider = sclet_ph_stub_provider(nm, "positive",
                                                                     dep = "shared_with_signature"))
    }
    on.exit(for (nm in nms) PhenotypeEvidenceProviders("unregister", name = nm), add = TRUE)
    res <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(), evidence = nms)
    cd <- SummarizedExperiment::colData(res)
    # Every cell should be at most 'probable': the signature is the only line.
    expect_equal(sum(cd$phenotype_call_call == "supported"), 0L)
    rec <- sclet_ph_ns("sclet_get_state_record")(res, "annotation", "phenotype_call")
    expect_equal(rec$summary$independent_lines, 0L)
})

test_that("two independent agreeing lines reach supported, disagreeing does not", {
    sce <- sclet_ph_test_sce()

    with_provider(sclet_ph_stub_provider("indep_yes", "positive"), {
        res <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                                  evidence = "indep_yes", transcriptomic_cut = 0)
    })
    cd <- SummarizedExperiment::colData(res)
    score <- cd$phenotype_call_signature_score
    call <- as.character(cd$phenotype_call_call)
    # Signature-positive cells now have two agreeing lines.
    expect_true(all(call[score >= 0] == "supported"))
    expect_gt(sum(call == "supported"), 0L)
    # The stub says yes everywhere. A signature-negative cell is NOT a conflict:
    # a low score is the absence of that line, not a negative one, so the cell
    # has exactly the stub's single line and is 'probable', not 'supported'.
    expect_true(all(call[score < 0] == "probable"))
    # A genuine conflict needs a provider that says NO.
    with_provider(sclet_ph_stub_provider("indep_no2", "negative"), {
        res2 <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                                   evidence = "indep_no2", transcriptomic_cut = 0)
    })
    cd2 <- SummarizedExperiment::colData(res2)
    s2 <- cd2$phenotype_call_signature_score
    expect_true(all(as.character(cd2$phenotype_call_call)[s2 >= 0] == "undetermined"))

    with_provider(sclet_ph_stub_provider("indep_no", "negative"), {
        res2 <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                                   evidence = "indep_no")
    })
    cd2 <- SummarizedExperiment::colData(res2)
    expect_equal(sum(cd2$phenotype_call_call == "supported"), 0L)
    expect_gt(sum(cd2$phenotype_call_call == "undetermined"), 0L)
})

test_that("a competing hypothesis demotes a supported call", {
    sce <- sclet_ph_test_sce()
    n <- ncol(sce)

    # The competing programme must be measurable on BOTH the called cells and
    # the rest. Give the "alt" block to the first population, and make the
    # provider say yes only there, so that the two groups really differ.
    alt <- sprintf("BG%03d", seq_len(5L))
    lc <- SummarizedExperiment::assay(sce, "logcounts")
    lc[alt, 1:30] <- lc[alt, 1:30] + 3
    SummarizedExperiment::assay(sce, "logcounts") <- lc

    half_provider <- sclet_ph_stub_provider("half_yes", "positive")
    half_provider$decide <- function(object, ...) {
        list(state = ifelse(seq_len(ncol(object)) <= 30, "positive", "negative"))
    }
    PhenotypeEvidenceProviders("register", provider = half_provider)
    on.exit(PhenotypeEvidenceProviders("unregister", name = "half_yes"), add = TRUE)

    res <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                              evidence = "half_yes",
                              competing = list(alt = alt),
                              competing_demote = 0.05,
                              transcriptomic_cut = 0)
    rec <- sclet_ph_ns("sclet_get_analysis")(res, "phenotype_call")$phenotype_call
    expect_equal(rec$competing$alt$status, "measured")
    expect_true(isTRUE(rec$competing$alt$demote))
    cd <- SummarizedExperiment::colData(res)
    expect_equal(sum(cd$phenotype_call_call == "supported"), 0L)
})

test_that("a competing hypothesis with no contrast is reported, not scored", {
    sce <- sclet_ph_test_sce()
    PhenotypeEvidenceProviders("register", provider = sclet_ph_stub_provider("indep_yes", "positive"))
    on.exit(PhenotypeEvidenceProviders("unregister", name = "indep_yes"), add = TRUE)

    # The stub calls every cell, so there is no "rest" to separate from.
    res <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                              evidence = "indep_yes",
                              competing = list(alt = sprintf("MAL%02d", seq_len(10L))),
                              transcriptomic_cut = 0)
    rec <- sclet_ph_ns("sclet_get_analysis")(res, "phenotype_call")$phenotype_call
    expect_true(is.na(rec$competing$alt$separation))
    expect_false(isTRUE(rec$competing$alt$demote))
})

test_that("claims use the package's existing ladder, not a new one", {
    sce <- sclet_ph_test_sce()
    with_provider(sclet_ph_stub_provider("indep_yes", "positive"), {
        res <- CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                                  evidence = "indep_yes", transcriptomic_cut = 0)
    })
    cd <- SummarizedExperiment::colData(res)
    claim <- cd$phenotype_call_claim
    expect_true(all(claim %in% c("consistent_with", "associated", "hypothesis")))
    expect_true(all(claim[cd$phenotype_call_call == "supported"] == "consistent_with"))
    expect_true(all(claim[cd$phenotype_call_call == "probable"] == "associated"))
    # Nothing here may reach 'observed': a per-cell label is never a direct
    # observation of cell state.
    expect_false(any(claim == "observed"))
})

test_that("a badly formed provider is refused at registration", {
    bad <- list(name = "bad", label = "no decide")
    expect_error(PhenotypeEvidenceProviders("register", provider = bad), "decide")

    bad2 <- sclet_ph_stub_provider("bad2", "positive")
    bad2$dependence <- "maybe"
    expect_error(PhenotypeEvidenceProviders("register", provider = bad2), "dependence")

    bad3 <- sclet_ph_stub_provider("bad3", "positive")
    bad3$decide <- function(object, ...) list(state = rep("positive", 3))
    expect_error(
        with_provider(bad3, CallPhenotypeCells(sclet_ph_test_sce(),
                                               signature = sclet_ph_test_sig(),
                                               evidence = "bad3")),
        "states for"
    )
})

test_that("the signature cutoff source is recorded", {
    sce <- sclet_ph_test_sce()
    res <- suppressWarnings(CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                                               evidence = character(0)))
    rec <- sclet_ph_ns("sclet_get_analysis")(res, "phenotype_call")$phenotype_call
    expect_equal(rec$signature_cut$source, "median")

    res2 <- suppressWarnings(CallPhenotypeCells(sce, signature = sclet_ph_test_sig(),
                                                evidence = character(0),
                                                transcriptomic_cut = 0.5))
    rec2 <- sclet_ph_ns("sclet_get_analysis")(res2, "phenotype_call")$phenotype_call
    expect_equal(rec2$signature_cut$source, "caller")
    expect_equal(rec2$signature_cut$value, 0.5)
})
