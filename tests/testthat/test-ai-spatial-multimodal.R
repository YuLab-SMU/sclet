## Read-only spatial and multimodal contracts: readiness, bounded evidence
## summaries, and the across-analysis / across-modality comparisons.
##
## The heart of this file is the privacy regression (block 7) plus the
## no-recommendation guard (block 8). The stored spatial records carry real
## column names from the fixture's deconvolution result and the multimodal
## records are registered through the REAL register_multimodal_adt() /
## register_multimodal_atac(), so the records these functions read have the
## genuine shape rather than a hand-rolled approximation.
##
## No real Python / basilisk environment is ever built here: the fixtures are
## in-memory SingleCellExperiments carrying hand-written analysis-state
## records, exactly the way test-ai-fate-drivers.R builds its CellRank record
## without running CellRank.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

## Bare object: counts only, no spatial and no multimodal record at all. Used
## for every "not_ready" / "no_*_record" branch.
sclet_ai_test_sm_bare_object <- function(n_cells = 12L, n_genes = 10L) {
    set.seed(20251016L)
    counts <- matrix(rpois(n_genes * n_cells, lambda = 6L),
        nrow = n_genes, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(n_genes)), paste0("cell_", seq_len(n_cells))))
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
}

## Spatial fixture: counts plus a spatial reducedDim plus TWO different spatial
## state records, mirroring the two shapes R/spatial.R actually writes.
##
##   * `spatial_deconv` -- RunSpatialDeconvolution() writes a "spatial" record
##     whose method is "cell2location", whose artifacts$analysis_key is
##     "spatial_deconv", and whose artifacts$columns are the deconvolution
##     column names added to colData (R/spatial.R:124-140).
##   * `spatial_niche` -- RunSpatialNiche() writes a "spatial" record whose
##     method is "SVP_LISA", whose artifacts$analysis_key is "spatial_niche",
##     whose artifacts$features are the real feature names, and which carries
##     summary$n_features (R/spatial.R:495-515).
##
## `with_second = FALSE` records only the deconvolution analysis, for the
## "a comparison needs at least two analyses" case.
sclet_ai_test_sm_spatial_object <- function(n_cells = 24L, n_genes = 30L,
                                            n_features = 4L, with_second = TRUE) {
    set.seed(20251016L)
    genes <- paste0("gene_", seq_len(n_genes))
    cells <- paste0("cell_", seq_len(n_cells))
    counts <- matrix(rpois(n_genes * n_cells, lambda = 6L),
        nrow = n_genes, ncol = n_cells, dimnames = list(genes, cells))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    ## The coordinate representation spatial consumers read.
    coords <- cbind(
        x = stats::runif(n_cells, 0, 1000),
        y = stats::runif(n_cells, 0, 1000)
    )
    rownames(coords) <- cells
    SingleCellExperiment::reducedDim(sce, "spatial") <- coords

    ## RunSpatialDeconvolution() appends one colData column per cell type. These
    ## names are the identities that must never leave the summary.
    deconv_columns <- paste0("c2l_celltype_", seq_len(n_features))
    for (column in deconv_columns) {
        SummarizedExperiment::colData(sce)[[column]] <- stats::runif(n_cells)
    }

    niche_features <- paste0("gene_", seq_len(n_features))

    sce <- sclet_set_analysis_state(
        object = sce,
        type = "spatial",
        id = "deconv_fixture",
        method = "cell2location",
        inputs = list(
            ref_group_key = "celltype",
            ref_batch_key = "batch",
            n_common_genes = n_genes
        ),
        artifacts = list(
            analysis_key = "spatial_deconv",
            columns = deconv_columns
        ),
        params = list(max_epochs = 30000L)
    )

    if (isTRUE(with_second)) {
        sce <- sclet_set_analysis_state(
            object = sce,
            type = "spatial",
            id = "niche_fixture",
            method = "SVP_LISA",
            inputs = list(
                features = niche_features,
                n_features = length(niche_features)
            ),
            artifacts = list(
                analysis_key = "spatial_niche",
                features = niche_features,
                n_features = length(niche_features)
            ),
            summary = list(
                n_features = length(niche_features)
            )
        )
    }

    sce
}

## A minimal altExp with its own features and cells, so the ADT/ATAC record
## carries genuine n_features / n_cells from register_multimodal_adt().
sclet_ai_test_sm_alt_exp <- function(n_cells = 24L, n_features = 8L, prefix = "adt") {
    set.seed(20251016L)
    features <- paste0(prefix, "_feat_", seq_len(n_features))
    counts <- matrix(rpois(n_features * n_cells, lambda = 5L),
        nrow = n_features, ncol = n_cells,
        dimnames = list(features, paste0("cell_", seq_len(n_cells))))
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
}

## Bare counts plus an "ADT" altExp that has NOT been registered. This is the
## "attached but not registered" branch, and it is deliberately distinct from
## the bare object: the fix the readiness question should ask for is a
## registration call, not "load the data".
sclet_ai_test_sm_adt_unregistered <- function(n_cells = 24L, n_genes = 10L, n_features = 8L) {
    sce <- sclet_ai_test_sm_bare_object(n_cells = n_cells, n_genes = n_genes)
    SingleCellExperiment::altExp(sce, "ADT") <- sclet_ai_test_sm_alt_exp(
        n_cells = n_cells, n_features = n_features, prefix = "adt")
    sce
}

## Multimodal fixture: an ADT altExp registered through the REAL
## register_multimodal_adt(), and optionally an ATAC altExp registered through
## the REAL register_multimodal_atac(). `with_atac = FALSE` yields a
## single-modality object, for the "a comparison needs at least two" case.
sclet_ai_test_sm_multimodal_object <- function(n_cells = 24L, n_genes = 10L,
                                               n_adt = 8L, n_atac = 6L,
                                               with_atac = TRUE) {
    sce <- sclet_ai_test_sm_bare_object(n_cells = n_cells, n_genes = n_genes)
    SingleCellExperiment::altExp(sce, "ADT") <- sclet_ai_test_sm_alt_exp(
        n_cells = n_cells, n_features = n_adt, prefix = "adt")
    if (isTRUE(with_atac)) {
        SingleCellExperiment::altExp(sce, "ATAC") <- sclet_ai_test_sm_alt_exp(
            n_cells = n_cells, n_features = n_atac, prefix = "atac")
    }

    sce <- register_multimodal_adt(sce)
    if (isTRUE(with_atac)) {
        sce <- register_multimodal_atac(sce)
    }
    sce
}

## Every real identifier the fixtures store. The privacy regression asserts that
## none of these appears anywhere in any returned structure: not a cell
## barcode, not a gene symbol, not an altExp feature name, not a deconvolution
## column name.
sclet_ai_test_sm_secrets <- function(n_cells = 24L, n_genes = 30L,
                                     n_features = 4L, n_adt = 8L, n_atac = 6L) {
    c(
        paste0("cell_", seq_len(n_cells)),
        paste0("gene_", seq_len(n_genes)),
        paste0("c2l_celltype_", seq_len(n_features)),
        paste0("adt_feat_", seq_len(n_adt)),
        paste0("atac_feat_", seq_len(n_atac))
    )
}

## Walk a nested list and hand every leaf to `visit`. Recursion covers unnamed
## list elements too, so a structure that buries an identifier one level deeper
## than expected is still caught.
sclet_ai_test_sm_walk <- function(node, visit) {
    if (is.list(node)) {
        for (element in node) {
            sclet_ai_test_sm_walk(element, visit)
        }
    } else {
        visit(node)
    }
    invisible(NULL)
}

## All field names anywhere in a nested list, at every depth.
sclet_ai_test_sm_field_names <- function(node) {
    names <- character()
    stack <- list(node)
    while (length(stack)) {
        current <- stack[[1L]]
        stack <- stack[-1L]
        if (is.list(current) && !is.null(names(current))) {
            names <- c(names, names(current))
            stack <- c(stack, unname(current))
        }
    }
    unique(names[!is.na(names) & nzchar(names)])
}

## Collect the values stored under one field name, at any depth and in any
## number, so a code-bearing field is checked wherever it appears.
sclet_ai_test_sm_values_for <- function(node, field) {
    found <- list()
    stack <- list(node)
    while (length(stack)) {
        current <- stack[[1L]]
        stack <- stack[-1L]
        if (!is.list(current)) next
        current_names <- names(current)
        if (!is.null(current_names)) {
            for (position in which(current_names == field)) {
                found[[length(found) + 1L]] <- current[[position]]
            }
            stack <- c(stack, unname(current))
        }
    }
    found
}

## TRUE when any leaf of the structure has one entry per cell -- the shape a
## per-cell deconvolution proportion or per-cell modality vector would take if
## it ever leaked.
sclet_ai_test_sm_has_per_cell <- function(node, n_cells) {
    found <- FALSE
    sclet_ai_test_sm_walk(node, function(leaf) {
        if (is.atomic(leaf) && length(leaf) == n_cells) found <<- TRUE
    })
    found
}

## The anonymized codes the result exposes, whichever field carries them: the
## spatial comparison uses `analysis_codes`, the multimodal comparison
## `modality_codes`, and a single-modality summary may carry just one.
sclet_ai_test_sm_codes <- function(result) {
    codes <- character()
    for (field in c("analysis_code", "analysis_codes", "modality_code", "modality_codes")) {
        codes <- c(codes, unlist(sclet_ai_test_sm_values_for(result, field), use.names = FALSE))
    }
    codes <- as.character(codes)
    codes[!is.na(codes) & nzchar(codes)]
}

## The anonymized-code shape every identity field must satisfy (EX4).
sclet_ai_test_sm_code_pattern <- "^(group|cluster|condition|sample|route)_[0-9]+$"

## The no-recommendation pattern (EX1).
##
## The delimiter class is [^a-z0-9] rather than \b because R's \b treats "_" as
## a word character, so a \b-anchored pattern silently misses underscore-joined
## names like "best_method" -- exactly the field the guard exists to catch.
## `[a-z]*` after the verb group covers inflections ("recommended",
## "suggesting"), which a bare stem-plus-delimiter alternation drops. The
## contract-mandated `recommendation` placeholder is excluded by the caller,
## since it must exist and be NULL rather than be absent.
sclet_ai_test_sm_recommendation_pattern <- paste0(
    "(^|[^a-z0-9])(suggest|recommend|best|winner|top|prefer|optimal|favor)[a-z]*([^a-z0-9]|$)")

## Field names that read as a nomination. `top_n` is a legitimate bounded-count
## parameter, so it is excluded explicitly rather than by weakening the pattern.
sclet_ai_test_sm_offending_names <- function(field_names) {
    offending <- grep(sclet_ai_test_sm_recommendation_pattern, field_names,
        ignore.case = TRUE, value = TRUE)
    setdiff(offending, c("recommendation", "top_n"))
}

## In the privacy regression block 7 the contract forbids applying
## sclet_ai_evidence_bad_field(): that gate guards EVIDENCE NODE field names and
## these six functions never write evidence, so applying it here is a category
## error.

## Validate one one-step plan the way test-ai-fate-drivers.R does, with the
## plan-level fingerprint taken from the object's own ledger.
sclet_ai_test_sm_plan <- function(sce, registry, action, params = list()) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "spatial_multimodal",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "sm1", action = action, params = as.list(params)))
    )
    ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
}

## The six actions under test, grouped the way the registry groups them.
##
## NAMING (settled -- do not "re-fix" this back): the readiness functions and
## their registered actions both use check_*_readiness. The action names were
## briefly inspect_spatial_readiness / inspect_multimodal_readiness by mistake;
## there was never a pre-existing `inspect_*_readiness` action -- the existing
## readiness-ish names are inspect_status, inspect_ledger and check_qc, and
## check_qc the action wraps the check_qc() function. So the convention pairs
## the action name to the function name, and R/ai-execution.R was renamed to
## check_spatial_readiness / check_multimodal_readiness to match. These helpers
## are the single place the action names are written down, so a rename on either
## side surfaces here as one clearly-labelled failure rather than as six
## confusing ones.
sclet_ai_test_sm_spatial_actions <- function() {
    c("check_spatial_readiness", "summarize_spatial_evidence", "compare_spatial_analyses")
}

sclet_ai_test_sm_multimodal_actions <- function() {
    c("check_multimodal_readiness", "summarize_multimodal_evidence", "compare_multimodal_evidence")
}

sclet_ai_test_sm_all_actions <- function() {
    c(sclet_ai_test_sm_spatial_actions(), sclet_ai_test_sm_multimodal_actions())
}

## -- Block 1: spatial readiness ----------------------------------------------
test_that("check_spatial_readiness is blocked on a bare object and ready on the spatial fixture", {
    bare <- check_spatial_readiness(sclet_ai_test_sm_bare_object())
    expect_equal(bare$status, "not_ready")
    expect_true("spatial" %in% bare$blocked_actions)
    expect_true(length(bare$questions) > 0L)

    ready <- check_spatial_readiness(sclet_ai_test_sm_spatial_object())
    expect_equal(ready$status, "ready_for_diagnostic")
    expect_equal(length(ready$blocked_actions), 0L)
    expect_true(isTRUE(ready$checks$spatial_record_available))
    expect_false(isTRUE(ready$checks$raw_values_included))

    ## Both recorded analyses are named by their real analysis_key. The
    ## available-analysis list lives under `checks`, like every other
    ## readiness diagnostic in this package.
    available <- as.character(unlist(ready$checks$available_analyses, use.names = FALSE))
    expect_true("spatial_deconv" %in% available)
    expect_true("spatial_niche" %in% available)
})

## -- Block 2: spatial evidence summary ---------------------------------------
test_that("summarize_spatial_evidence reports no_spatial_record and analysis_not_available", {
    bare <- sclet_ai_test_sm_bare_object()
    none <- summarize_spatial_evidence(bare)
    expect_equal(none$status, "not_available")
    expect_equal(none$reason, "no_spatial_record")
    expect_false(isTRUE(none$raw_values_included))

    ## A named analysis with no matching record is refused rather than
    ## silently defaulted to whatever record happens to exist.
    sce <- sclet_ai_test_sm_spatial_object()
    unknown <- summarize_spatial_evidence(sce, analysis = "spatial_colocalization")
    expect_equal(unknown$status, "not_available")
    expect_equal(unknown$reason, "analysis_not_available")

    ## A recorded analysis really is summarized.
    known <- summarize_spatial_evidence(sce, analysis = "spatial_deconv")
    expect_equal(known$status, "available")
    expect_false(isTRUE(known$raw_values_included))
})

## -- Block 3: across-analysis spatial comparison -----------------------------
test_that("compare_spatial_analyses needs two analyses and anonymizes their codes", {
    single <- compare_spatial_analyses(sclet_ai_test_sm_spatial_object(with_second = FALSE))
    expect_equal(single$status, "not_available")
    expect_equal(single$reason, "at_least_two_analyses_required")

    ## Nor can an object with no spatial record at all be compared.
    none <- compare_spatial_analyses(sclet_ai_test_sm_bare_object())
    expect_equal(none$status, "not_available")
    expect_equal(none$reason, "at_least_two_analyses_required")

    ## Two analyses: the fixture records spatial_deconv and spatial_niche.
    result <- compare_spatial_analyses(sclet_ai_test_sm_spatial_object())
    expect_equal(result$status, "available")
    expect_equal(result$n_analyses, 2L)
    expect_equal(length(result$analysis_codes), 2L)
    expect_equal(length(unique(result$analysis_codes)), 2L)
    expect_true(all(grepl(sclet_ai_test_sm_code_pattern, result$analysis_codes)))
    expect_false(isTRUE(result$raw_values_included))
})

## -- Block 4: multimodal readiness -------------------------------------------
test_that("check_multimodal_readiness distinguishes not attached from not registered", {
    ## (a) Bare object: no altExp, no record.
    bare <- check_multimodal_readiness(sclet_ai_test_sm_bare_object())
    expect_equal(bare$status, "not_ready")
    expect_true("multimodal" %in% bare$blocked_actions)
    expect_true(length(bare$questions) > 0L)
    expect_false(isTRUE(bare$checks$adt_available))
    expect_false(isTRUE(bare$checks$adt_registered))
    bare_questions <- paste(unlist(bare$questions, use.names = FALSE), collapse = " ")

    ## (b) ADT altExp attached but never registered: available yet not
    ## registered. This is the branch whose remedy is a registration call.
    attached <- check_multimodal_readiness(sclet_ai_test_sm_adt_unregistered())
    expect_true(isTRUE(attached$checks$adt_available))
    expect_false(isTRUE(attached$checks$adt_registered))
    expect_true(length(attached$questions) > 0L)
    attached_questions <- paste(unlist(attached$questions, use.names = FALSE), collapse = " ")

    ## (c) After the REAL registration the record exists.
    registered <- check_multimodal_readiness(sclet_ai_test_sm_multimodal_object())
    expect_true(isTRUE(registered$checks$adt_available))
    expect_true(isTRUE(registered$checks$adt_registered))
    expect_true(isTRUE(registered$checks$atac_available))
    expect_true(isTRUE(registered$checks$atac_registered))
    expect_equal(registered$status, "ready_for_diagnostic")
    expect_equal(length(registered$blocked_actions), 0L)

    ## The two questions must not be interchangeable. The distinguishing
    ## markers are the ones that only make sense in one branch:
    ##   * the bare branch says the modality is NOT ATTACHED and asks the user
    ##     to attach it first;
    ##   * the attached-but-unregistered branch says it IS attached (altExp
    ##     exists) but is NOT REGISTERED, and explicitly warns against
    ##     re-attaching data that is already present.
    ## A bare "register" test is NOT sufficient: the not-attached wording also
    ## mentions registration when explaining why no record can exist.
    expect_true(grepl("is not registered", attached_questions, fixed = TRUE),
        info = "attached-but-unregistered branch must name the missing registration")
    expect_true(grepl("do not re-attach", attached_questions, fixed = TRUE),
        info = "attached-but-unregistered branch must not ask for a re-attach")
    expect_false(grepl("is not registered", bare_questions, fixed = TRUE),
        info = "bare branch has nothing registered to be missing")
    expect_false(grepl("do not re-attach", bare_questions, fixed = TRUE),
        info = "bare branch has nothing attached to avoid re-attaching")
    ## The messages begin a sentence, so the leading letter is capitalised.
    expect_true(grepl("No ADT modality is attached", bare_questions, fixed = TRUE),
        info = "bare branch must report the modality as not attached")
    ## The two branches are not the same message: the attached one says the
    ## altExp exists, which the bare one must never claim.
    expect_true(grepl("altExp \"ADT\" exists", attached_questions, fixed = TRUE),
        info = "attached branch must report the altExp as present")
    expect_false(grepl("altExp \"ADT\" exists", bare_questions, fixed = TRUE),
        info = "bare branch must not claim an altExp exists")
})

## -- Block 5: multimodal evidence summary ------------------------------------
test_that("summarize_multimodal_evidence reports per-modality counts without raw values", {
    bare <- sclet_ai_test_sm_bare_object()
    none <- summarize_multimodal_evidence(bare)
    expect_equal(none$status, "not_available")
    expect_equal(none$reason, "no_multimodal_record")

    sce <- sclet_ai_test_sm_multimodal_object()
    result <- summarize_multimodal_evidence(sce)
    expect_equal(result$status, "available")
    expect_false(isTRUE(result$raw_values_included))

    ## Every modality the fixture registered is summarized with the counts
    ## register_multimodal_adt() / _atac() actually recorded.
    modalities <- result$modalities
    expect_true(length(modalities) >= 2L)
    for (modality in modalities) {
        expect_true(is.numeric(modality$n_features) && length(modality$n_features) == 1L)
        expect_true(is.numeric(modality$n_cells) && length(modality$n_cells) == 1L)
        expect_false(isTRUE(modality$raw_values_included))
    }
    counts <- vapply(modalities, function(modality) as.integer(modality$n_features), integer(1L))
    expect_true(8L %in% counts)
    expect_true(6L %in% counts)
})

## -- Block 6: across-modality comparison -------------------------------------
test_that("compare_multimodal_evidence needs two modalities and anonymizes their codes", {
    single <- compare_multimodal_evidence(sclet_ai_test_sm_multimodal_object(with_atac = FALSE))
    expect_equal(single$status, "not_available")
    expect_equal(single$reason, "at_least_two_modalities_required")

    none <- compare_multimodal_evidence(sclet_ai_test_sm_bare_object())
    expect_equal(none$status, "not_available")
    expect_equal(none$reason, "at_least_two_modalities_required")

    ## Two modalities: ADT and ATAC, both registered for real.
    result <- compare_multimodal_evidence(sclet_ai_test_sm_multimodal_object())
    expect_equal(result$status, "available")
    expect_equal(length(result$modality_codes), 2L)
    expect_equal(length(unique(result$modality_codes)), 2L)
    expect_true(all(grepl(sclet_ai_test_sm_code_pattern, result$modality_codes)))
    expect_false(isTRUE(result$raw_values_included))

    ## The refusal to nominate a modality is explicit, not merely silent.
    expect_true("recommendation" %in% names(result))
    expect_null(result$recommendation)
})

## -- Block 7: the privacy regression -----------------------------------------
test_that("no spatial or multimodal function leaks a barcode, symbol or per-cell vector", {
    spatial <- sclet_ai_test_sm_spatial_object()
    multimodal <- sclet_ai_test_sm_multimodal_object()
    bare <- sclet_ai_test_sm_bare_object()
    secrets <- sclet_ai_test_sm_secrets()

    results <- list(
        spatial_readiness = check_spatial_readiness(spatial),
        spatial_readiness_bare = check_spatial_readiness(bare),
        spatial_summary = summarize_spatial_evidence(spatial, analysis = "spatial_deconv"),
        spatial_summary_unknown = summarize_spatial_evidence(spatial, analysis = "spatial_colocalization"),
        spatial_comparison = compare_spatial_analyses(spatial),
        multimodal_readiness = check_multimodal_readiness(multimodal),
        multimodal_readiness_bare = check_multimodal_readiness(bare),
        multimodal_summary = summarize_multimodal_evidence(multimodal),
        multimodal_comparison = compare_multimodal_evidence(multimodal)
    )

    for (name in names(results)) {
        result <- results[[name]]
        field_names <- sclet_ai_test_sm_field_names(result)

        ## (a) No barcode, gene symbol, altExp feature or deconvolution column
        ## name anywhere -- not in values, not in field names, not in a printed
        ## rendering of the whole structure.
        text <- paste(capture.output(str(result)), collapse = "\n")
        for (secret in secrets) {
            if (grepl("^(gene|cell|adt_feat|atac_feat|c2l_celltype)_", secret)) {
                ## gene_1 is a prefix of gene_10..gene_19, and cell_1 of
                ## cell_10.., so match on a delimiter rather than a substring.
                expect_false(grepl(paste0(secret, "([^0-9]|$)"), text),
                    info = paste(name, "leaked", secret))
            } else {
                expect_false(grepl(secret, text, fixed = TRUE),
                    info = paste(name, "leaked", secret))
            }
            expect_false(secret %in% field_names,
                info = paste(name, "leaked", secret, "as a field name"))
        }

        ## (b) Every ANONYMIZED CODE field is code-shaped. The rule cannot apply
        ## to the whole structure: the contract requires free-text `questions`
        ## and prose `caveats`. It is scoped to the fields carrying identity.
        identity_fields <- c("analysis_code", "analysis_codes", "modality_code",
            "modality_codes")
        for (field in intersect(identity_fields, field_names)) {
            values <- unlist(sclet_ai_test_sm_values_for(result, field), use.names = FALSE)
            values <- as.character(values)
            values <- values[!is.na(values) & nzchar(values)]
            if (!length(values)) next
            expect_true(all(grepl(sclet_ai_test_sm_code_pattern, values)),
                info = paste(name, "returned a non-anonymized value in", field))
        }
        codes <- sclet_ai_test_sm_codes(result)
        if (length(codes)) {
            expect_true(all(grepl(sclet_ai_test_sm_code_pattern, codes)),
                info = paste(name, "returned a non-anonymized code"))
        }

        ## (c) No per-cell vector: nothing has one entry per cell, which is how
        ## a per-cell deconvolution proportion would leak.
        expect_false(sclet_ai_test_sm_has_per_cell(result, ncol(spatial)),
            info = paste(name, "carried a per-cell vector"))

        ## (d) Raw values are always withheld where the field is present.
        if (!is.null(result$raw_values_included)) {
            expect_false(isTRUE(result$raw_values_included),
                info = paste(name, "reported raw_values_included"))
        }
    }
})

## -- Block 8: nothing recommends, ranks or nominates -------------------------
test_that("no spatial or multimodal function recommends an analysis or a modality", {
    spatial <- sclet_ai_test_sm_spatial_object()
    multimodal <- sclet_ai_test_sm_multimodal_object()

    results <- list(
        check_spatial_readiness(spatial),
        check_spatial_readiness(sclet_ai_test_sm_bare_object()),
        summarize_spatial_evidence(spatial, analysis = "spatial_deconv"),
        summarize_spatial_evidence(spatial, analysis = NULL),
        compare_spatial_analyses(spatial),
        check_multimodal_readiness(multimodal),
        check_multimodal_readiness(sclet_ai_test_sm_bare_object()),
        summarize_multimodal_evidence(multimodal),
        compare_multimodal_evidence(multimodal)
    )

    for (result in results) {
        offending <- sclet_ai_test_sm_offending_names(sclet_ai_test_sm_field_names(result))
        expect_equal(unique(offending), character(0))

        ## `recommendation` must exist and be NULL on the comparison returns,
        ## mirroring compare_fate_lineages(): the refusal is explicit rather
        ## than merely silent.
        if ("recommendation" %in% names(result)) {
            expect_null(result$recommendation)
        }
    }

    ## Both compare functions carry the explicit null placeholder.
    for (comparison in list(compare_spatial_analyses(spatial),
        compare_multimodal_evidence(multimodal))) {
        expect_true("recommendation" %in% names(comparison))
        expect_null(comparison$recommendation)
    }

    ## SELF-TEST: the pattern really catches the names it exists to catch, in
    ## bare, underscore-joined and inflected form. A future weakening of the
    ## pattern fails here loudly instead of passing on a blind spot.
    must_catch <- c("best_method", "suggested_method", "recommended_method",
        "top_analysis", "best_method_name", "suggest_method", "recommending",
        "top_modality", "prefer_method", "optimal_analysis", "favor_group",
        "winner")
    expect_equal(
        sort(sclet_ai_test_sm_offending_names(must_catch)),
        sort(must_catch)
    )

    ## ...while the contract-mandated placeholder and legitimate names are
    ## genuinely untouched.
    must_pass <- c("recommendation", "top_n", "n_analyses", "analysis_codes",
        "modality_codes", "raw_values_included", "status", "reason", "questions",
        "caveats", "n_features", "n_cells", "available_analyses", "spatial_record_available")
    expect_equal(sclet_ai_test_sm_offending_names(must_pass), character(0))
})

## -- Block 9: the six actions are registered and read-only -------------------
test_that("the six spatial and multimodal actions are registered as read-only", {
    sce <- sclet_ai_test_sm_spatial_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "spatial", "multimodal"))

    for (action_name in sclet_ai_test_sm_all_actions()) {
        expect_true(action_name %in% names(registry),
            info = paste(action_name, "is not registered"))
        action <- registry[[action_name]]
        expect_s3_class(action, "sclet_ai_action")
        expect_equal(action$name, action_name)
        expect_false(action$mutates_object)
        expect_false(action$requires_confirmation)
        ## "value" is the established kind for "the handler returns a non-SCE
        ## value": the three read-only actions in the `read` group all use it.
        expect_equal(action$returns, "value")
        expect_true(is.function(action$handler),
            info = paste(action_name, "has no handler"))
    }

    ## The default registry must not carry them: they are not part of `read`.
    default_registry <- AIDefaultExecutionRegistry(sce)
    expect_equal(intersect(sclet_ai_test_sm_all_actions(), names(default_registry)),
        character(0))

    ## include = "read" plus the spatial group yields exactly the three spatial
    ## actions and none of the multimodal ones.
    spatial_registry <- AIDefaultExecutionRegistry(sce, include = c("read", "spatial"))
    expect_true(all(sclet_ai_test_sm_spatial_actions() %in% names(spatial_registry)))
    expect_equal(intersect(sclet_ai_test_sm_multimodal_actions(), names(spatial_registry)),
        character(0))
})

## -- Block 10: typed prerequisites on a bare object --------------------------
test_that("all six actions fail typed prerequisites on a bare object instead of erroring", {
    sce <- sclet_ai_test_sm_bare_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "spatial", "multimodal"))

    reasons <- c(
        check_spatial_readiness = "no_spatial_record",
        summarize_spatial_evidence = "no_spatial_record",
        compare_spatial_analyses = "no_spatial_record",
        check_multimodal_readiness = "no_multimodal_record",
        summarize_multimodal_evidence = "no_multimodal_record",
        compare_multimodal_evidence = "no_multimodal_record"
    )
    ## The reason table must cover exactly the six registered actions.
    expect_equal(sort(names(reasons)), sort(sclet_ai_test_sm_all_actions()))

    for (action_name in names(reasons)) {
        validation <- sclet_ai_test_sm_plan(sce, registry, action_name, list())
        expect_false(validation$valid, info = paste(action_name, "validated on a bare object"))
        ## ValidateAIPlan prefixes prerequisite errors with
        ## "prerequisites not met for action <step_id>: ", so the token is
        ## matched anywhere in the error and NEVER with a leading ^ anchor.
        expect_true(any(grepl(reasons[[action_name]], validation$errors)),
            info = paste(action_name, "did not report", reasons[[action_name]]))
    }

    ## A region-free sanity check on the convention itself: the error is a
    ## prefixed string, which is exactly why an anchored pattern would be wrong.
    validation <- sclet_ai_test_sm_plan(sce, registry, "check_spatial_readiness", list())
    expect_true(any(grepl("prerequisites not met for action",
        validation$errors, fixed = TRUE)))
})
