# sclet_ai_registry_malignant(): built-in AI actions for CNV inference and malignant-cell calling.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# These run inferCNV through the basilisk-backed wrappers.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_malignant <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    # Malignant-cell identification as actions (merged 2026-10-18).
    #
    # The package grew RunCNVInference()/CallMalignantCells(), but the AI layer
    # could not reach either -- a planner had no way to name a copy-number run or
    # the two-line call. These two actions mirror the functions, with the
    # invariants the feature documents kept in the actions:
    #
    #   * the CNV reference must be EXOGENOUS and named by the caller; the
    #     action refuses to fall back to clusters the analysis itself produced;
    #   * the call is a TWO-line call; running it without a CNV record would
    #     silently downgrade it to signature-only under the stronger name, so the
    #     prerequisite requires the record.
    if ("malignant" %in% include) {
    actions <- c(actions, list(
        run_cnv_inference = AIAction(
            name = "run_cnv_inference",
            description = paste0(
                "Infer per-cell copy-number burden with inferCNV against an ",
                "EXOGENOUS reference: name normal cells with `reference_ids` ",
                "or a colData column with `reference_col`. Never use a cluster ",
                "the analysis itself produced as the reference. ",
                "EXPENSIVE: this is the slowest action in the registry. ",
                "Measured 2026-10-19 on 1,400 cells with HMM on (the ",
                "default): 66 minutes and 22.6 GB of memory; HMM post-",
                "processing alone (Gibbs sampling) was 27 minutes. ",
                "Run it only when the goal asks for copy-number or ",
                "malignant evidence, and tell the user it is a long step ",
                "before starting it. `hmm = FALSE` skips the HMM stage and ",
                "is much faster, but `run_cnv_inference`'s burden is then ",
                "built from the non-HMM profile."
            ),
            handler = function(object, params) {
                name <- params$name %||% "cnv"
                reference_cells <- if (!is.null(params$reference_ids)) {
                    sclet_ai_vector_param(params$reference_ids, "character")
                } else {
                    NULL
                }
                object <- RunCNVInference(
                    object,
                    reference_cells = reference_cells,
                    reference_col = params$reference_col %||% NULL,
                    species = params$species %||% "human",
                    assay = params$assay %||% "counts",
                    hmm = if (is.null(params$hmm)) TRUE else isTRUE(params$hmm),
                    cutoff = params$cutoff %||% 1,
                    num_threads = params$num_threads %||% 1,
                    name = name
                )
                sclet_ai_record_cnv_evidence(object, name = name)
            },
            input_schema = list(
                reference_ids = "character_vector",
                reference_col = "character",
                species = c("human", "mouse"),
                assay = "character",
                name = "character",
                hmm = "logical",
                cutoff = "number",
                num_threads = "integer"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!requireNamespace("infercnv", quietly = TRUE)) {
                    return("optional_package_missing: infercnv is required for copy-number inference; install Bioconductor 'infercnv'")
                }
                assays <- SummarizedExperiment::assayNames(object)
                assay <- params$assay %||% "counts"
                if (!assay %in% assays) {
                    return(paste0("counts_assay_missing: assay '", assay,
                        "' not found; inferCNV needs raw counts (available: ",
                        paste(assays, collapse = ", "), ")"))
                }
                count_problem <- sclet_ai_assay_count_problem(object, assay)
                if (!is.null(count_problem)) return(count_problem)
                if (is.null(params$reference_col) && is.null(params$reference_ids)) {
                    return(paste0(
                        "reference_missing: copy-number inference needs an ",
                        "EXOGENOUS reference; name normal cells with ",
                        "reference_ids, or a colData column with reference_col. ",
                        "Do not use clusters the analysis itself produced."))
                }
                if (!is.null(params$reference_col)) {
                    cd <- colnames(SummarizedExperiment::colData(object))
                    if (!is.character(params$reference_col) ||
                            !params$reference_col %in% cd) {
                        return(paste0("reference_col_missing: '",
                            params$reference_col, "' is not a colData column (available: ",
                            paste(cd, collapse = ", "), ")"))
                    }
                }
                if (!is.null(params$reference_ids)) {
                    cells <- sclet_ai_vector_param(params$reference_ids, "character")
                    known <- intersect(cells, colnames(object))
                    if (length(known) < 5L) {
                        return(paste0("reference_too_small: ", length(known),
                            " supplied reference cell(s) exist on the object; at ",
                            "least 5 are needed for the CNV cutoff, otherwise the ",
                            "CNV line is unavailable and the call degrades silently"))
                    }
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(
                required_states = list(cnv = "any", ai_evidence = "any"),
                cnv_id = "cnv_record_id"
            ),
            allowed_state_types = c("cnv", "ai_evidence"),
            mutates_object = TRUE,
            requires_confirmation = TRUE,
            estimated_cost = "high",
            idempotent = FALSE
        ),
        call_malignant_cells = AIAction(
            name = "call_malignant_cells",
            description = paste0(
                "Combine copy-number and transcriptional evidence into the ",
                "four-state malignant call (malignant / likely_malignant / ",
                "undetermined / non_malignant). Requires a completed CNV run ",
                "and either a signature score column or signature genes; ",
                "`identity_markers` names a competing identity programme to ",
                "demote."
            ),
            handler = function(object, params) {
                name <- params$name %||% "malignant_call"
                signature <- NULL
                if (!is.null(params$signature)) {
                    genes <- sclet_ai_vector_param(params$signature, "character")
                    signature <- stats::setNames(list(genes), "malignant")
                }
                identity_genes <- if (!is.null(params$identity_markers)) {
                    sclet_ai_vector_param(params$identity_markers, "character")
                } else {
                    NULL
                }
                object <- CallMalignantCells(
                    object,
                    cnv_id = params$cnv_id %||% "cnv",
                    signature = signature,
                    signature_score = params$signature_score %||% NULL,
                    transcriptomic_cut = params$transcriptomic_cut %||% NULL,
                    cnv_cut = params$cnv_cut %||% NULL,
                    reference_quantile = params$reference_quantile %||% 0.99,
                    identity_genes = identity_genes,
                    identity_demote = params$identity_demote %||% 0.1,
                    assay = params$assay %||% NULL,
                    name = name
                )
                state_col <- SummarizedExperiment::colData(object)[[
                    paste0(name, "_state")]]
                state_counts <- table(as.character(state_col))
                count_of <- function(level) {
                    value <- unname(state_counts[level])
                    if (length(value) && !is.na(value)) as.integer(value) else 0L
                }
                state_rec <- tryCatch(
                    sclet_get_state_record(object, "annotation", name),
                    error = function(e) NULL)
                state_summary <- state_rec$summary %||% list()
                values <- list(
                    n_cells = as.integer(ncol(object)),
                    n_malignant = count_of("malignant"),
                    n_likely_malignant = count_of("likely_malignant"),
                    n_undetermined = count_of("undetermined"),
                    n_non_malignant = count_of("non_malignant"),
                    # Read from the call's own state record, where
                    # CallMalignantCells() computed them, instead of
                    # re-deriving a second answer here.
                    cnv_available = as.integer(isTRUE(state_summary$cnv_available)),
                    identity_conflict = as.integer(isTRUE(state_summary$identity_conflict)),
                    # What actually decided these counts. The signature
                    # cut defaults to the POOL MEDIAN, so "signature
                    # positive" is a relative statement that lands near
                    # half the cells by construction (measured: 22,006 of
                    # 44,011 in a user report, 2026-10-20); and the CNV
                    # threshold comes from the reference burden. Without
                    # both numbers the four-state counts cannot be
                    # re-derived or challenged from the ledger.
                    signature_cut_value = as.numeric(state_summary$signature_cut %||% NA_real_),
                    signature_cut_is_median = as.integer(isTRUE(
                        state_summary$signature_cut_is_median)),
                    cnv_reference_threshold = as.numeric(
                        state_summary$cnv_reference_threshold %||% NA_real_),
                    n_cnv_positive = as.integer(state_summary$n_cnv_positive %||% NA_integer_),
                    n_signature_positive = as.integer(
                        state_summary$n_signature_positive %||% NA_integer_),
                    n_pre_demote_supported = as.integer(
                        state_summary$n_pre_demote_supported %||% NA_integer_)
                )
                sclet_ai_record_evidence_for_run(object, list(
                        id = paste0("ev:malignant_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        values = values,
                        claim_level = "associated"
                    ),
                    source = name,
                    parents = character(),
                    scope = NULL,
                    replace = isTRUE(params$replace_evidence))
            },
            input_schema = list(
                cnv_id = "character",
                signature = "character_vector",
                signature_score = "character",
                transcriptomic_cut = "number",
                cnv_cut = "number",
                reference_quantile = "number",
                identity_markers = "character_vector",
                identity_demote = "number",
                assay = "character",
                name = "character",
                # Re-running the same call on an object whose evidence node
                # already exists: refresh it instead of refusing. Off by
                # default, because two runs that disagree should be visible
                # rather than silently reconciled.
                replace_evidence = "logical"
            ),
            prerequisites = function(object, params, planned = NULL) {
                cnv_id <- params$cnv_id %||% "cnv"
                records <- sclet_get_analysis(object, "malignant_cnv")
                planned_cnv <- planned$cnv_ids %||% character()
                # A well-formed output reference counts as pointing at the
                # promise: `cnv_id = "${<step>.output.cnv_id}"` resolves at
                # execution from the producing step's summary, which carries
                # the real record id. Demanding the LITERAL id here rejected
                # the exact spelling the 2026-10-19 benchmark measured as the
                # one a model actually writes (the same defect the enrichment
                # gate had for detest_id). The reference shape ties the id to
                # a step of this plan; execution stays the authority on
                # whether that step really registered a CNV run.
                is_output_reference <- is.character(cnv_id) &&
                    length(cnv_id) == 1L &&
                    grepl("^\\$\\{[A-Za-z0-9_.-]+\\.output\\.[A-Za-z0-9_.-]+\\}$",
                        cnv_id)
                if (is_output_reference) {
                    if (!length(planned_cnv)) {
                        return(paste0(
                            "cnv_record_missing: the call references ",
                            cnv_id, ", but no run_cnv_inference step ",
                            "precedes it in this plan; run run_cnv_inference ",
                            "first. The malignant call is a TWO-line call and ",
                            "will not be silently downgraded to ",
                            "signature-only."))
                    }
                    cnv_id <- planned_cnv[[length(planned_cnv)]]
                }
                have_real <- !is.null(records) && !is.null(records[[cnv_id]])
                have_planned <- length(planned_cnv) &&
                    cnv_id %in% planned_cnv
                if (!have_real && !have_planned) {
                    return(paste0(
                        "cnv_record_missing: no CNV run '", cnv_id,
                        "' on the object or earlier in this plan; run ",
                        "run_cnv_inference first. The malignant call is a ",
                        "TWO-line call and will not be silently downgraded ",
                        "to signature-only."))
                }
                if (is.null(params$signature) && is.null(params$signature_score)) {
                    return(paste0(
                        "signature_missing: the malignant call needs a ",
                        "transcriptional line -- pass signature_score (a colData ",
                        "column) or signature"))
                }
                if (!is.null(params$signature_score)) {
                    cd <- colnames(SummarizedExperiment::colData(object))
                    if (!is.character(params$signature_score) ||
                            !params$signature_score %in% cd) {
                        return(paste0("signature_score_missing: '",
                            params$signature_score, "' is not a colData column (available: ",
                            paste(cd, collapse = ", "), ")"))
                    }
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(
                required_states = list(malignancy = "any", ai_evidence = "any")
            ),
            allowed_state_types = c("malignancy", "ai_evidence"),
            mutates_object = TRUE,
            requires_confirmation = TRUE,
            estimated_cost = "low",
            idempotent = FALSE
        )
    ))
    }

    # ---- iterative signature classifier ----------------------------------
    #
    # A different analysis from `call_malignant_cells`: that one combines CNV
    # with a fixed signature; this one clusters an AddModuleScore difference and
    # refines it with pseudobulk limma. Two analyses, two records, two evidence
    # nodes -- the ids must not collide.
    if ("malignant" %in% include) {
    actions <- c(actions, list(
        classify_malignant_epithelium = AIAction(
            name = "classify_malignant_epithelium",
            description = paste(
                "Iteratively classify malignant epithelium from two marker-gene signatures",
                "(tumor vs normal) with AddModuleScore, k-means and pseudobulk limma,",
                "registering one bounded evidence node (malignant/non-malignant counts,",
                "fraction, iterations, convergence). The two gene vectors are the user's",
                "signature: this action never invents them. A signature-based call is",
                "recorded at claim level 'associated', never as a measurement."
            ),
            handler = function(object, params) {
                name <- params$name %||% "malignant"
                args <- list(
                    object = object,
                    seed = list(tumor = params$seed_tumor, normal = params$seed_normal),
                    name = name
                )
                if (!is.null(params$method)) args$method <- params$method
                if (!is.null(params$sample_col)) args$sample_col <- params$sample_col
                if (!is.null(params$n_top)) args$n_top <- as.integer(params$n_top)
                if (!is.null(params$padj)) args$padj <- as.numeric(params$padj)
                result <- do.call(ClassifyMalignantEpithelium, args)
                record <- sclet_ai_mainline_record(result, "malignancy", name,
                    "classify_malignant_epithelium")
                summary <- record$summary %||% list()
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:malignant_iterative_",
                            sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "associated",
                        values = list(
                            n_malignant = sclet_ai_mainline_number(summary$n_malignant),
                            n_non_malignant = sclet_ai_mainline_number(summary$n_non_malignant),
                            malignant_fraction = sclet_ai_mainline_number(summary$malignant_fraction),
                            iterations = sclet_ai_mainline_number(summary$iterations),
                            converged = isTRUE(summary$converged)
                        )
                    ),
                    source = name
                )
            },
            input_schema = list(
                seed_tumor = list(type = "character_vector", required = TRUE),
                seed_normal = list(type = "character_vector", required = TRUE),
                method = c("iterative", "cnv"),
                sample_col = "character",
                n_top = "integer",
                padj = "number",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                tumor <- params$seed_tumor
                normal <- params$seed_normal
                if (!length(tumor) || !length(normal)) {
                    return("seed_missing: classify_malignant_epithelium needs `seed_tumor` and `seed_normal` as two non-empty gene vectors (the user's signature, not a guess)")
                }
                if (!is.character(tumor) || !is.character(normal)) {
                    return("seed_invalid: `seed_tumor` and `seed_normal` must be character vectors of gene identifiers")
                }
                genes <- c(tumor, normal)
                if (!any(genes %in% rownames(object))) {
                    return("seed_not_found: none of the seed genes is in rownames(object); check the identifier type (symbol vs Ensembl)")
                }
                sample_col <- params$sample_col %||% "sample"
                if (!sample_col %in% colnames(SummarizedExperiment::colData(object))) {
                    return(paste0(
                        "sample_column_missing: '", sample_col,
                        "' is not a colData column; the pseudobulk step needs the sample unit. Available: ",
                        paste(colnames(SummarizedExperiment::colData(object)), collapse = ", ")
                    ))
                }
                if (ncol(object) < 20L) {
                    return("too_few_cells_for_split: the iterative classifier needs at least 2 * min_cells_per_group (20) cells")
                }
                if (!base::requireNamespace("limma", quietly = TRUE)) {
                    return("optional_package_missing: limma is required for classify_malignant_epithelium")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(malignancy = "any", ai_evidence = "any"),
                note = paste(
                    "Per-cell scores and the classification land in colData. The evidence",
                    "node carries counts and the convergence flag only."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("malignancy", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE
        )
    ))
    }

    actions
}
