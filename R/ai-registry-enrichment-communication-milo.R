# sclet_ai_registry_enrichment_communication_milo(): built-in AI actions for enrichment, cell-cell communication and differential abundance.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# The three expensive tail groups, order preserved.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_enrichment_communication_milo <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("enrichment" %in% include) {
    actions <- c(actions, list(
        run_enrichment = AIAction(
            name = "run_enrichment",
            description = paste(
                "Summarize a completed differential-expression marker table into enrichment categories",
                "with clusterProfiler, and register bounded aggregate evidence for it. The marker table",
                "and its cluster column are required: enrichment is an interpretation layer over the DE",
                "result, not a second measurement of it. The database, ontology, key type and the",
                "significance cutoffs are the user's biological modelling choices and are never picked",
                "or recommended here. No gene symbol and no category name or accession is recorded as",
                "evidence; the full result stays in the state record for get_enrichment()."
            ),
            handler = function(object, params) {
                name <- params$name %||% "enrichment"
                # A re-run that produces nothing must not be reported as a
                # success: RunEnrichment() warns and returns the object
                # unchanged when no marker passes the cutoffs, when the
                # marker table has no p-value column (whose comparison
                # yields logical(0), the same zero-row outcome), or when it
                # lacks a logFC column. Comparing the record's created_at
                # before and after distinguishes "this call wrote a result"
                # from "a stale record from an earlier call is still here".
                before <- tryCatch(
                    sclet_get_state_record(object, "enrichment", name),
                    error = function(e) NULL
                )
                before_created <- before$created_at %||% NA_character_
                # The marker-table gate may have been skipped at
                # validation time: in a chained plan the detest record was
                # only a promise then. Now the upstream step has really run,
                # so apply the same helper here. Without this, a chain whose
                # DE table has no cluster column would validate, run, write
                # a degenerate one-group record and report success, because
                # RunEnrichment() warns rather than errors in that case and
                # the record it writes is non-empty.
                shape_problem <- sclet_ai_enrichment_marker_problem(
                    object, params$detest_id)
                if (!is.null(shape_problem)) {
                    stop(shape_problem, call. = FALSE)
                }
                result <- RunEnrichment(
                    object,
                    id = params$detest_id,
                    p_val_cutoff = params$p_val_cutoff %||% 0.05,
                    logfc_cutoff = params$logfc_cutoff %||% 0.25,
                    db = params$db %||% "GO",
                    orgDb = params$orgDb %||% "org.Hs.eg.db",
                    keyType = params$keyType %||% "SYMBOL",
                    ont = params$ont %||% "BP",
                    name = name
                )
                after <- tryCatch(
                    sclet_get_state_record(result, "enrichment", name),
                    error = function(e) NULL
                )
                wrote_record <- !is.null(after) &&
                    !identical(as.character(after$created_at %||% NA_character_),
                        as.character(before_created))
                if (!wrote_record) {
                    stop(paste0(
                        "enrichment_no_result: run_enrichment produced no enrichment record '", name,
                        "'. Either no marker passed p_val_cutoff / logfc_cutoff, or the marker table",
                        " lacks a p-value or logFC column for the cutoffs to be applied. Lower the",
                        " cutoffs, or re-run RunDEtest() so the marker table carries them. Nothing",
                        " was registered as evidence and no claim is supported by this step."
                    ), call. = FALSE)
                }
                # A record can exist and still hold nothing. On PBMC 10k
                # (2026-10-19) the KEGG round wrote an enrichment record with
                # n_enriched = 0 and artifacts$result NULL -- clusterProfiler
                # had warned "No enrichment found in any gene cluster" -- and
                # the step still reported "completed" with no evidence node.
                # The guard above only asks whether a record appeared, which
                # an empty run satisfies. Empty is not absent: fail here, as
                # the cutoff path above does, and say which way out exists.
                n_enriched <- after$summary$n_enriched %||% NA_integer_
                has_result <- !is.null(after$artifacts$result)
                if ((!is.na(n_enriched) && n_enriched == 0L) || !has_result) {
                    stop(paste0(
                        "enrichment_no_result: run_enrichment recorded '", name,
                        "' with zero enriched terms (n_enriched = ", n_enriched, ", result stored: ",
                        has_result, "). clusterProfiler found nothing for db = ", params$db %||% "GO",
                        " / keyType = ", params$keyType %||% "SYMBOL",
                        "; KEGG expects ENTREZID, so symbol input commonly comes back empty. ",
                        "Nothing was registered as evidence and no claim is supported by this step. ",
                        "Try db = \"GO\" with keyType = \"SYMBOL\", or convert the ids first."
                    ), call. = FALSE)
                }
                sclet_ai_record_enrichment_evidence(
                    result,
                    name = name,
                    detest_id = params$detest_id
                )
            },
            input_schema = list(
                detest_id = "character",
                db = "character",
                orgDb = "character",
                keyType = "character",
                ont = "character",
                p_val_cutoff = "number",
                logfc_cutoff = "number",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                # `db` is a biological modelling choice: GO and KEGG use
                # different annotation universes and different statistics.
                # The handler still falls back to GO for a direct
                # RunEnrichment() call, but an AI plan must not inherit that
                # fallback silently -- every sibling action with such a
                # choice (run_cci `method`, run_milo design roles,
                # run_annotation `ref`) refuses to run without it, and this
                # one gate was missing, so "GO" could arrive without anyone
                # choosing it (found 2026-10-19 on PBMC 10k, where a plan
                # omitting `db` validated).
                db <- params$db
                if (is.null(db) || !length(db) || !nzchar(trimws(as.character(db)[[1L]]))) {
                    return(paste0(
                        "db_missing: run_enrichment requires an explicit `db`: \"GO\" or \"KEGG\". ",
                        "The annotation database is the user's modelling choice and the engine ",
                        "default (GO) is not applied to an AI plan; pass the value, or ask the ",
                        "user which one to use."
                    ))
                }
                if (!as.character(db)[[1L]] %in% c("GO", "KEGG")) {
                    return(paste0(
                        "db_invalid: `db` must be \"GO\" or \"KEGG\"; got \"",
                        as.character(db)[[1L]], "\". The annotation database is the user's ",
                        "modelling choice and is never selected automatically."
                    ))
                }
                # KEGG annotates ENTREZID. A SYMBOL table is accepted by
                # the handler and then comes back empty, which the M3 guard
                # turns into `enrichment_no_result` -- a true but late
                # diagnosis (measured 2026-10-19: 46.8 s spent to learn it).
                # Refuse here, where the caller can act on it.
                if (identical(as.character(params$db %||% ""), "KEGG") &&
                    identical(as.character(params$keyType %||% "SYMBOL"), "SYMBOL")) {
                    return(paste0(
                        "key_type_incompatible: KEGG annotates ENTREZID, not SYMBOL, so a SYMBOL ",
                        "marker table comes back empty. Convert the ids first (register the ",
                        "converted table as its own detest record and set keyType = \"ENTREZID\"), ",
                        "or use db = \"GO\" with keyType = \"SYMBOL\"."
                    ))
                }
                detest <- tryCatch(get_detest(object, params$detest_id),
                    error = function(e) NULL)
                # TRUE only when the record is on the object right now, so
                # its contents can be inspected. A chained plan has only a
                # promise (see below), which is not the same thing and must
                # not be treated as if it were.
                record_available <- !is.null(detest)
                if (!record_available) {
                    # A plan that chains run_de_test -> run_enrichment
                    # cannot have its upstream record at validation time,
                    # because nothing has executed yet. The projection in
                    # sclet_ai_planned_state_records() is the promise that
                    # an earlier step registers it, and the handler
                    # re-checks the marker table for real at execution.
                    planned_detest <- as.character(
                        planned$state_records$detest %||% character())
                    requested <- params$detest_id
                    # A well-formed output reference counts as pointing at
                    # the promise: `detest_id =
                    # "${<step>.output.detest_id}"` resolves at execution
                    # from the producing step's summary, which carries the
                    # real record id (run_de_test now exposes it).
                    # Demanding the LITERAL id here rejected every chained
                    # plan that used the reference -- the exact spelling
                    # the 2026-10-16 benchmark measured as the only
                    # legitimate one once detest_id became referenceable.
                    # The reference shape itself ties the id to a step of
                    # this plan; execution stays the authority on whether
                    # that step really registered a detest record.
                    is_output_reference <- is.character(requested) &&
                        length(requested) == 1L &&
                        grepl("^\\$\\{[A-Za-z0-9_.-]+\\.output\\.[A-Za-z0-9_.-]+\\}$", requested)
                    if (is_output_reference) {
                        requested <- planned_detest
                    }
                    if (length(requested) == 1L && nzchar(requested)) {
                        if (!requested %in% planned_detest) {
                            return(paste0(
                                "detest_state_missing: run_enrichment is computed from a completed ",
                                "differential-expression marker table; run RunDEtest() first. ",
                                "Requested detest id '", as.character(requested),
                                "' is neither an existing record nor produced by an earlier step of ",
                                "this plan. Planned detest records: ",
                                if (length(planned_detest)) paste(planned_detest, collapse = ", ") else "none", "."
                            ))
                        }
                    } else if (!length(planned_detest)) {
                        return(paste0(
                            "detest_state_missing: run_enrichment is computed from a completed ",
                            "differential-expression marker table; run RunDEtest() first."
                        ))
                    }
                }
                # The marker table's shape can only be judged once the
                # record exists. When it is merely planned, this gate is
                # skipped here and the handler applies the same helper
                # after the upstream step has really run.
                if (record_available) {
                    problem <- sclet_ai_enrichment_marker_problem(
                        object, params$detest_id)
                    if (!is.null(problem)) return(problem)
                }
                db <- params$db %||% "GO"
                if (!is.character(db) || length(db) != 1L || !db %in% c("GO", "KEGG")) {
                    return(paste0(
                        "db_invalid: run_enrichment accepts db = \"GO\" or \"KEGG\"; got '",
                        paste(as.character(db), collapse = ", "), "'. The database is a modelling ",
                        "choice and is never selected automatically."
                    ))
                }
                org_db <- params$orgDb %||% "org.Hs.eg.db"
                missing_packages <- c(
                    if (!base::requireNamespace("clusterProfiler", quietly = TRUE)) "clusterProfiler" else NULL,
                    if (is.character(org_db) && length(org_db) == 1L && nzchar(org_db) &&
                        !base::requireNamespace(org_db, quietly = TRUE)) org_db else NULL
                )
                if (length(missing_packages)) {
                    return(paste0(
                        "enrichment_package_missing: run_enrichment needs ",
                        paste(missing_packages, collapse = " and "),
                        "; install it via BiocManager::install()."
                    ))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(enrichment = "any", ai_evidence = "any"),
                note = paste(
                    "Enrichment categories are derived from the DE marker table and are model output,",
                    "not independent measurements: they inherit that table's thresholds and confounders,",
                    "and the recorded evidence declares the DE node as its parent so a claim ceiling",
                    "keeps counting the two as a single line of support. No gene symbol, category name",
                    "or category accession enters the evidence layer, and no cutoff, database or",
                    "ontology is recommended by this action."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("detest", "enrichment", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = TRUE
        )
    ))
    }

    if ("communication" %in% include) {
    actions <- c(actions, list(
        run_cci = AIAction(
            name = "run_cci",
            description = paste(
                "Infer cell-cell interactions with the engine the user names, and register bounded",
                "aggregate evidence for it. The method is required and never defaulted: CellChat,",
                "CellPhoneDB and NicheNet use different ligand-receptor databases, expression tests and",
                "sender/receiver framing, so choosing between them is a modelling decision, not a",
                "parameter to fill in. No ligand, receptor or pathway identity is recorded as",
                "evidence; the interaction table stays in the analysis layer for get_cci()."
            ),
            handler = function(object, params) {
                method <- params$method
                # Record ids are engine-specific and the caller may omit
                # name, so the newly written record is located by comparing
                # created_at across communication records rather than
                # assumed. This also detects a run that wrote no record.
                before <- sclet_get_state_records(object, "communication")
                before_created <- vapply(before, function(record) {
                    as.character(record$created_at %||% "")
                }, character(1L))
                result <- RunCCI(
                    object,
                    method = method,
                    keep_object = isTRUE(params$keep_object %||% FALSE),
                    group = params$group,
                    assay_name = params$assay_name,
                    layer = params$layer,
                    name = params$name
                )
                after <- sclet_get_state_records(result, "communication")
                created_after <- vapply(after, function(record) {
                    as.character(record$created_at %||% "")
                }, character(1L))
                written <- names(created_after)[vapply(names(created_after), function(id) {
                    previous <- if (id %in% names(before_created)) {
                        before_created[[id]]
                    } else {
                        ""
                    }
                    !identical(created_after[[id]], previous)
                }, logical(1L))]
                if (!length(written)) {
                    stop(paste0(
                        "communication_no_result: run_cci with method '", method,
                        "' wrote no communication record, so nothing can be reported or claimed from",
                        " it. Check that the group column resolves to at least two populations with",
                        " cells in the expression source, then re-run."
                    ), call. = FALSE)
                }
                sclet_ai_record_communication_evidence(
                    result,
                    name = written[[1L]],
                    method = method
                )
            },
            input_schema = list(
                method = "character",
                group = "character",
                assay_name = "character",
                layer = "character",
                name = "character",
                keep_object = "logical"
            ),
            prerequisites = function(object, params, planned = NULL) {
                engines <- sclet_ai_clarification_engines()
                method <- params$method
                if (is.null(method) || !is.character(method) || length(method) != 1L ||
                        !nzchar(trimws(method))) {
                    return(paste0(
                        "method_missing: run_cci requires an explicit method because the engines are",
                        " not interchangeable (different ligand-receptor databases and statistical tests).",
                        " Pass one of: ", paste(engines, collapse = ", "),
                        ". No engine is chosen for the user."
                    ))
                }
                if (!method %in% engines) {
                    return(paste0(
                        "method_invalid: run_cci accepts method = \"CellChat\", \"CellPhoneDB\" or",
                        " \"NicheNet\"; got '", method, "'. No engine is substituted automatically."
                    ))
                }
                columns <- colnames(SummarizedExperiment::colData(object))
                group_key <- params$group %||% tryCatch(ActiveIdent(object), error = function(e) NULL)
                if (!is.null(group_key) && !as.character(group_key) %in% columns &&
                        !identical(as.character(group_key), "colLabels")) {
                    group_key <- NULL
                }
                if (is.null(group_key)) {
                    candidates <- grep("cluster|leiden|louvain|ident", columns,
                        ignore.case = TRUE, value = TRUE)
                    if (length(candidates)) group_key <- candidates[[1L]]
                }
                if (is.null(group_key)) {
                    return(paste0(
                        "group_column_missing: run_cci groups cells into sender and receiver",
                        " populations and none is resolvable; set an active identity or pass 'group'.",
                        " Available colData columns: ", paste(columns, collapse = ", "), "."
                    ))
                }
                if (!length(SummarizedExperiment::assayNames(object))) {
                    return(
                        "expression_source_missing: run_cci needs an expression assay and this object has none."
                    )
                }
                required_packages <- c(CellChat = "CellChat", CellPhoneDB = "basilisk",
                    NicheNet = "nichenetr")
                if (identical(method, "CellPhoneDB")) {
                    # CellPhoneDB is the one engine that crosses into Python.
                    missing_bridge <- c("reticulate", "basilisk")[
                        !vapply(c("reticulate", "basilisk"),
                            base::requireNamespace, logical(1L), quietly = TRUE)]
                    if (length(missing_bridge)) {
                        return(paste0(
                            "communication_backend_unavailable: CellPhoneDB runs through ",
                            paste(missing_bridge, collapse = " and "),
                            "; install them to use this engine."
                        ))
                    }
                }
                needed <- required_packages[[method]]
                if (!base::requireNamespace(needed, quietly = TRUE)) {
                    return(paste0(
                        "communication_package_missing: run_cci with method '", method, "' needs the '",
                        needed, "' package; install it via BiocManager::install('", needed, "')."
                    ))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(communication = "any", ai_evidence = "any"),
                note = paste(
                    "Inferred interactions are model estimates computed from expression plus a curated",
                    "ligand-receptor database; a predicted pair is not evidence that cells signalled, and",
                    "an absent pair is uninformative because each engine filters against its own database.",
                    "Only group-flow counts and strength quantiles reach the evidence layer, and every",
                    "communication run on an object shares a dependency group, so agreement between two",
                    "engines over the same expression matrix and labels is not two independent lines of",
                    "support. This action never chooses an engine, and removes or filters no cell."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("communication", "ai_evidence"),
            estimated_cost = "high",
            idempotent = TRUE
        )
    ))
    }

    if ("milo" %in% include) {
    actions <- c(actions, list(
        run_milo = AIAction(
            name = "run_milo",
            description = paste(
                "Test differential abundance of neighbourhoods with Milo/miloR and register one",
                "bounded evidence node. The sample and condition design roles are mandatory: this",
                "action never guesses which colData column plays which role, and they must be",
                "declared through ConfirmAIDesignSemantics() naming the columns yourself. Only",
                "neighbourhood counts and FDR quantiles reach the evidence layer; the design frame,",
                "the model formula and every sample or condition label stay in the state record."
            ),
            handler = function(object, params) {
                name <- params$name %||% "milo"
                # Same silent-path discipline as enrichment and CCI: Milo
                # can reach testNhoods and write a record with nothing
                # usable in it, so created_at is compared before and after
                # and an absent record is a typed failure, not a success.
                before <- sclet_get_state_records(object, "milo")
                before_created <- vapply(before, function(record) {
                    as.character(record$created_at %||% "")
                }, character(1L))
                milo_out <- RunMilo(
                    object,
                    sample_col = params$sample_col,
                    condition_col = params$condition_col,
                    covariates = params$covariates,
                    contrasts = params$contrasts,
                    prop = params$prop %||% 0.1,
                    name = name
                )
                # RunMilo() returns list(sce=, milo=, da_results=), not an
                # SCE. The handler assumed an SCE, so every execution died
                # after the expensive part with `no metadata method for
                # class list` -- and no benchmark had ever run this action
                # (coverage: none), so nothing caught it (2026-10-20).
                result <- if (inherits(milo_out, "SingleCellExperiment")) {
                    milo_out
                } else {
                    milo_out$sce
                }
                after <- sclet_get_state_records(result, "milo")
                created_after <- vapply(after, function(record) {
                    as.character(record$created_at %||% "")
                }, character(1L))
                written <- names(created_after)[vapply(names(created_after), function(id) {
                    previous <- if (id %in% names(before_created)) before_created[[id]] else ""
                    !identical(created_after[[id]], previous)
                }, logical(1L))]
                if (!length(written)) {
                    stop(paste0(
                        "milo_no_result: run_milo wrote no milo record for name '", name,
                        "', so nothing can be reported or claimed from it. Check that the design ",
                        "columns resolve and that the neighbourhood graph was built, then re-run."
                    ), call. = FALSE)
                }
                sclet_ai_record_milo_evidence(result, name = written[[1L]])
            },
            input_schema = list(
                sample_col = list(type = "character", required = TRUE),
                condition_col = list(type = "character", required = TRUE),
                covariates = "character_vector",
                contrasts = "character_vector",
                prop = "number",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!base::requireNamespace("miloR", quietly = TRUE)) {
                    return(paste0(
                        "milo_package_missing: run_milo needs miloR; ",
                        "install it via BiocManager::install('miloR')."
                    ))
                }
                columns <- colnames(SummarizedExperiment::colData(object))
                # Deferred projection, same as every other consuming gate:
                # a plan that runs run_pca earlier in the SAME plan has not
                # executed yet, so `reducedDimNames(object)` is empty at
                # validation time and this gate refused every chained milo
                # plan (measured 2026-10-16, three records). The planned
                # view carries the reductions earlier steps declare, and
                # execution stays the authority through RunMilo()'s own
                # reducedDim check.
                reductions <- planned$reductions %||%
                    tryCatch(SingleCellExperiment::reducedDimNames(object),
                        error = function(e) character())
                if (!"PCA" %in% reductions) {
                    return(paste0(
                        "reduction_missing: run_milo builds its neighbourhood graph on a PCA ",
                        "embedding; run RunPCA() first. Available reductions: ",
                        paste(reductions, collapse = ", "), "."
                    ))
                }
                for (role in c("sample_col", "condition_col")) {
                    value <- params[[role]]
                    plain <- role
                    if (is.null(value) || !is.character(value) || length(value) != 1L ||
                            !nzchar(trimws(value))) {
                        return(paste0(
                            plain, "_missing: run_milo requires the '", plain,
                            "' column name to be supplied explicitly; the design role is never ",
                            "guessed from column names. Available colData columns: ",
                            paste(columns, collapse = ", "), "."
                        ))
                    }
                    if (!value %in% columns) {
                        return(paste0(
                            plain, "_not_in_colData: '", value, "' is not a colData column. ",
                            "Available colData columns: ", paste(columns, collapse = ", "), "."
                        ))
                    }
                }
                if (is.null(params$contrasts) || !length(params$contrasts)) {
                    return(paste0(
                        "contrasts_missing: run_milo needs explicit contrasts, e.g. ",
                        "\"trtA - ctrl\". The condition contrast is a statement about the design, ",
                        "not something to infer. Distinct condition values: ",
                        paste(utils::head(unique(as.character(
                            SummarizedExperiment::colData(object)[[params$condition_col]]
                        )), 8L), collapse = ", "), "."
                    ))
                }
                # Design semantics gate, shared with run_integration rather
                # than forked: the user's own declaration must cover both
                # roles, mapped onto the role names Milo uses.
                confirmation <- tryCatch(sclet_ai_find_design_confirmation(object, design = list(
                    sample = params$sample_col,
                    condition = params$condition_col
                )), error = function(e) NULL)
                if (is.null(confirmation)) {
                    return(paste0(
                        "design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, ",
                        "design = list(sample = '", params$sample_col,
                        "', condition = '", params$condition_col,
                        "')) before running Milo. Differential abundance attributes neighbourhoods ",
                        "to samples and contrasts between conditions, so both roles must be the ",
                        "user's own declaration; no column is inferred here."
                    ))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(milo = "any", ai_evidence = "any"),
                note = paste(
                    "Differential abundance results are associations between neighbourhoods and",
                    "conditions under a chosen graph and test; they are not causal and not",
                    "measured effects. Counts and FDR quantiles only reach the evidence layer;",
                    "the design frame, model formula and all sample or condition labels stay in",
                    "the state record. Every Milo node on an object shares a dependency group, so",
                    "a second contrast over the same graph is not a second line of support."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("reduction", "milo", "ai_evidence"),
            estimated_cost = "high",
            idempotent = TRUE
        )
    ))
    }

    actions
}
