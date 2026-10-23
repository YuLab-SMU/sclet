# sclet_ai_registry_annotation(): built-in AI actions for annotation and differential expression.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# The largest single group; its own file keeps the registry scannable.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_annotation <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("annotation" %in% include) {
    actions <- c(actions, list(
        run_de_test = AIAction(
            name = "run_de_test",
            description = "Run a differential-expression or all-cluster marker analysis and register evidence. For an experimental-group contrast (condition/treatment/disease) set group = the colData column and ident.1/ident.2 to its levels; add covariates = a colData column such as donor for a paired pseudobulk analysis.",
            handler = function(object, params) {
                ident.1 <- if (!is.null(params$ident.1)) sclet_ai_vector_param(params$ident.1, "character") else NULL
                ident.2 <- if (!is.null(params$ident.2)) sclet_ai_vector_param(params$ident.2, "character") else NULL
                group <- if (!is.null(params$group)) sclet_ai_vector_param(params$group, "character") else NULL
                covariates <- if (!is.null(params$covariates)) sclet_ai_vector_param(params$covariates, "character") else NULL
                # The thresholds this run actually uses. `params` carries the
                # caller's values and the defaults are applied once, here, so
                # the evidence node below can record what ran rather than a
                # constant: measured 2026-10-19 on PBMC 10k, a plan asking
                # logfc.threshold = 0.25 produced an evidence node reading
                # logfc_threshold = 0.1 -- a node describing a different
                # analysis than the one that executed.
                de_min_pct <- params$min.pct %||% 0.01
                de_logfc_threshold <- params$logfc.threshold %||% 0.1
                if (!is.null(group)) {
                    if (length(group) != 1L) stop("`group` must be one colData column name.")
                    ident.1 <- sclet_ai_resolve_group_value(object, group, ident.1)
                    ident.2 <- sclet_ai_resolve_group_value(object, group, ident.2)
                } else {
                    ident.1 <- sclet_ai_resolve_ident_value(object, ident.1)
                    ident.2 <- sclet_ai_resolve_ident_value(object, ident.2)
                }
                result <- RunDEtest(
                    object,
                    ident.1 = ident.1,
                    ident.2 = ident.2,
                    group = group,
                    covariates = covariates,
                    name = params$name %||% if (is.null(ident.1)) "findallmarkers" else "detest",
                    min.pct = de_min_pct,
                    logfc.threshold = de_logfc_threshold
                )
                rec_name <- params$name %||% if (is.null(ident.1)) "findallmarkers" else "detest"
                de_records <- sclet_get_state_records(result, "detest")
                rec <- de_records[[rec_name]]
                artifact <- if (!is.null(rec)) {
                    tryCatch(rec$artifacts$result, error = function(e) NULL)
                } else NULL
                if (!is.null(artifact) && is.data.frame(artifact)) {
                    logfc_col <- grep("logFC|avg_log2FC|avg_logFC|log2FC", colnames(artifact), value = TRUE, ignore.case = TRUE)[[1L]]
                    if (!is.na(logfc_col)) {
                        p_col <- grep("^p_val_adj$|^adj.P.Val$|padj|FDR|p_adj", colnames(artifact), value = TRUE, ignore.case = TRUE)[[1L]]
                        p_col <- if (!is.na(p_col)) p_col else grep("^p_val$|P.Value|^pvalue$", colnames(artifact), value = TRUE, ignore.case = TRUE)[[1L]]
                        # p is the engine's own p_val_adj filter inside
                        # FindMarkers, not a plan parameter; the logFC
                        # threshold is the one this run was given.
                        p_thresh <- 0.05
                        fc_thresh <- de_logfc_threshold
                        sig <- if (is.na(p_col)) {
                            rep(TRUE, nrow(artifact))
                        } else {
                            artifact[[p_col]] < p_thresh
                        }
                        sig[is.na(sig)] <- FALSE
                        up <- sig & artifact[[logfc_col]] >= fc_thresh
                        top_n <- 10L
                        ord <- order(ifelse(up, -artifact[[logfc_col]], NA_real_), na.last = NA)
                        top_genes <- if (length(ord)) as.character(rownames(artifact)[utils::head(ord[seq_len(min(top_n, length(ord)))], top_n)]) else character()
                        idents <- if (!is.null(group)) {
                            as.character(SummarizedExperiment::colData(result)[[group]])
                        } else {
                            Idents(result)
                        }
                        groups <- unique(as.character(idents))
                        # The marker table arrives ALREADY filtered by the
                        # caller's own thresholds (FindAllMarkers drops
                        # pval > return.thresh and |logFC| < logfc.threshold
                        # before returning), so recounting "padj < 0.05"
                        # restates the filter instead of measuring anything:
                        # measured 12,559 of 12,572 on PBMC 10k. `n_rows`
                        # already says how many rows came back.
                        aggregate_values <- list(
                            n_rows = as.integer(nrow(artifact)),
                            n_upregulated = as.integer(sum(up)),
                            n_groups = as.integer(length(groups)),
                            smallest_group_n = as.integer(min(as.integer(table(as.character(idents), useNA = "no")))),
                            pvalue_threshold = p_thresh,
                            logfc_threshold = fc_thresh,
                            min_pct = de_min_pct,
                            raw_values_included = FALSE
                        )
                        # `n_top_up` is a real count (how many of the top
                        # logFC genes were found), but the rank slots and the
                        # fraction slots were never measured -- `top_up_ranks`
                        # was 1:10 by construction and `top_up_fractions` was
                        # 10 x NA on every run. Promising a field and leaving
                        # it empty is worse than not promising it.
                        aggregate_values$n_top_up <- as.integer(length(top_genes))
                        evidence <- list(
                            id = if (is.null(ident.1)) {
                                "ev:findallmarkers"
                            } else if (!is.null(group)) {
                                paste0("ev:detest_", group, "_",
                                    paste(ident.1, collapse = "_"), "_vs_",
                                    paste(ident.2, collapse = "_"))
                            } else {
                                paste0("ev:detest_", paste(ident.1, collapse = "_"),
                                    if (length(ident.2)) paste0("_vs_", paste(ident.2, collapse = "_")))
                            },
                            kind = "deterministic_summary",
                            values = aggregate_values,
                            claim_level = "associated"
                        )
                        # Route through the re-run-aware recorder. The
                        # marker node id is FIXED ("ev:findallmarkers"), so
                        # re-running run_de_test with different thresholds
                        # hits the same id; swallowing that left the first
                        # run's n_rows/min_pct in the ledger while the
                        # record held the second run's table (measured
                        # 2026-10-21: record 0 rows, node still n_rows = 2,
                        # min_pct = 0.01 after a stricter re-run).
                        sclet_ai_record_evidence_for_run(result, evidence,
                            source = rec_name,
                            parents = if (!is.null(rec$id)) as.character(rec$id) else character(),
                            scope = NULL,
                            replace = isTRUE(params$replace_evidence))
                    } else result
                } else result
            },
            input_schema = list(
                ident.1 = "character_vector",
                ident.2 = "character_vector",
                group = "character",
                covariates = "character_vector",
                name = "character",
                min.pct = "number",
                logfc.threshold = "number",
                # Re-running marker detection with different thresholds on
                # an object that already has the node: refresh instead of
                # refusing. Off by default so two runs that disagree are
                # visible.
                replace_evidence = "logical"
            ),
            prerequisites = function(object, params, planned = NULL) {
                ident.1 <- if (is.null(params$ident.1)) NULL else sclet_ai_vector_param(params$ident.1, "character")
                ident.2 <- if (is.null(params$ident.2)) NULL else sclet_ai_vector_param(params$ident.2, "character")
                group <- if (is.null(params$group)) NULL else sclet_ai_vector_param(params$group, "character")
                covariates <- if (is.null(params$covariates)) NULL else sclet_ai_vector_param(params$covariates, "character")
                if (!is.null(group)) {
                    cd <- SummarizedExperiment::colData(object)
                    if (length(group) != 1L || !group %in% colnames(cd)) {
                        return(paste0("group_missing: 'group' must name one colData column; available: ",
                            paste(colnames(cd), collapse = ", ")))
                    }
                    values <- unique(as.character(cd[[group]]))
                    if (is.null(ident.1) || is.null(ident.2) ||
                            length(ident.1) != 1L || length(ident.2) != 1L) {
                        return(paste0("contrast_missing: a '", group,
                            "' contrast needs one ident.1 and one ident.2 (levels: ",
                            paste(values, collapse = ", "),
                            "). fix: pass both levels explicitly; they are",
                            " listed above and also in",
                            " context$capabilities$contrast_levels"))
                    }
                    if (!ident.1 %in% values || !ident.2 %in% values) {
                        # The valid levels are already in the message; what
                        # was missing is the instruction to COPY one of
                        # them rather than to invent a neighbouring label.
                        # Measured 2026-10-21: a round-2 plan wrote a level
                        # the object does not have, and the error named the
                        # right set without saying to read it.
                        bad <- unique(c(ident.1, ident.2)[
                            !c(ident.1, ident.2) %in% values])
                        return(paste0(
                            "contrast_invalid: ident.1/ident.2 must be levels of '",
                            group, "': ", paste(values, collapse = ", "),
                            " (you wrote: ", paste(bad, collapse = ", "),
                            "). fix: copy two of the listed levels verbatim,",
                            " or read them from",
                            " context$capabilities$contrast_levels -- do not",
                            " invent or translate a level name"))
                    }
                    if (identical(ident.1, ident.2)) {
                        return("contrast_invalid: ident.1 and ident.2 must differ")
                    }
                    if (length(covariates)) {
                        missing <- setdiff(covariates, colnames(cd))
                        if (length(missing)) {
                            return(paste0("covariates_missing: ",
                                paste(missing, collapse = ", ")))
                        }
                        if (!requireNamespace("limma", quietly = TRUE)) {
                            return("covariates_unavailable: covariate-adjusted DE needs the 'limma' package; install it or omit covariates (unadjusted cell-level test)")
                        }
                    }
                    return(TRUE)
                }
                active_ident <- ActiveIdent(object)
                planned_ident <- (planned %||% list())$active_ident
                if (!is.null(planned_ident) &&
                    !identical(as.character(planned_ident), as.character(active_ident))) {
                    # An earlier step in this plan will (re)set the identity
                    # -- typically find_clusters writing colLabels -- so the
                    # object's *current* identity is not the one this step
                    # will see, and neither its group sizes nor its label
                    # set can be checked before that step runs. Refusing the
                    # plan here would make the ordinary "cluster, then find
                    # markers" request unplannable from a raw object, which
                    # is the common case. Defer: the producing step's output
                    # contract (required_states/active_ident) and the
                    # execution-time check are the authority.
                    return(TRUE)
                }
                ident.1 <- sclet_ai_resolve_ident_value(object, ident.1)
                ident.2 <- sclet_ai_resolve_ident_value(object, ident.2)
                if (is.null(active_ident)) {
                    return("cluster identity is not set: run FindClusters() or set ActiveIdent(object) <- 'colname' before running marker/DE tests")
                }
                idents <- Idents(object)
                if (is.null(idents)) {
                    return("cluster identity is available at slot but has no values")
                }
                min_cells_per_group <- 2L
                ident_vec <- as.character(idents)
                grp_n <- table(ident_vec, useNA = "no")
                if (length(grp_n) < 1L || min(as.integer(grp_n)) < min_cells_per_group) {
                    return(paste0("at least one group/identity is too small; minimum required per group: ", min_cells_per_group))
                }
                if (is.null(ident.1) && length(grp_n) < 2L) {
                    return("FindAllMarkers requires at least two distinct identity groups; object has 1")
                }
                if (!is.null(ident.1) && !all(ident.1 %in% names(grp_n))) {
                    return(paste0("ident.1 contains unknown group labels: ", paste(setdiff(ident.1, names(grp_n)), collapse = ", ")))
                }
                if (!is.null(ident.2) && !all(ident.2 %in% names(grp_n))) {
                    return(paste0("ident.2 contains unknown group labels: ", paste(setdiff(ident.2, names(grp_n)), collapse = ", ")))
                }
                if (!is.null(ident.1) && !is.null(ident.2) && length(intersect(ident.1, ident.2))) {
                    return("ident.1 and ident.2 overlap; DE comparison groups must be disjoint")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(detest = "any"),
                # A referenceable name for the detest record this step
                # registers. Static value is a marker, not the id: the
                # schema is exposed through `available_actions` before any
                # step has run, and the real id is written into the
                # execution summary as an observation (see
                # `sclet_ai_execution_summary()`), which is also where the
                # `${step.output.detest_id}` projection resolves it from.
                # Needed because `run_enrichment` requires `detest_id` and
                # no registered action used to expose one -- measured as
                # the reason every model-proposed DE -> enrichment chain
                # died at validation (sclet-benchmark-ai, 2026-10-16,
                # three spellings across three legs). Same defect class as
                # `variable_features`, fixed there in 1bddedc.
                detest_id = "detest_record_id"
            ),
            mutates_object = TRUE,
            allowed_state_types = c("detest", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE
        ),
        run_annotation = AIAction(
            name = "run_annotation",
            description = "Run SingleR reference-based annotation with an explicitly supplied reference and record aggregate evidence.",
            handler = function(object, params) {
                ref_arg <- params$ref
                if (is.character(ref_arg) && length(ref_arg) == 1L && nzchar(ref_arg)) {
                    # Bare, `celldex::Name` and `<celldex::Name>` are one
                    # reference; the literal spelling a model writes must
                    # not decide whether an otherwise-planned run executes.
                    ref_name <- sclet_celldex_reference_name(ref_arg)
                    if (!is.null(ref_name) && !requireNamespace("celldex", quietly = TRUE)) {
                        # Name the missing package rather than calling the
                        # reference unsupported: the latter reads as the
                        # user's mistake and suggests no repair.
                        stop(paste0(
                            "celldex is required to resolve ref='", ref_arg, "'; install it via ",
                            "BiocManager::install('celldex') or pass a SummarizedExperiment ref directly."
                        ), call. = FALSE)
                    }
                    fun <- if (is.null(ref_name)) NULL else sclet_celldex_dataset(ref_name)
                    if (is.null(fun)) {
                        stop(paste0(
                            "reference '", ref_arg, "' is not a supported celldex dataset name. ",
                            "Provide a celldex dataset such as HumanPrimaryCellAtlasData ",
                            "or celldex::MonacoImmuneData, or a SummarizedExperiment ref directly."
                        ), call. = FALSE)
                    }
                    ref <- do.call(fun, list())
                } else {
                    ref <- ref_arg
                }
                labels_arg <- params$labels
                result <- RunSingleR(
                    object,
                    ref = ref,
                    labels = labels_arg,
                    layer = params$layer %||% NULL,
                    name = params$name %||% "singler"
                )
                name <- params$name %||% "singler"
                ann_records <- sclet_get_state_records(result, "annotation")
                ann_rec <- ann_records[[name]]
                if (!is.null(ann_rec)) {
                    artifacts <- if (!is.null(ann_rec)) tryCatch(ann_rec$artifacts, error = function(e) NULL) else NULL
                    labels_col <- if (!is.null(artifacts)) as.character(artifacts$labels_col) else NULL
                    score_col <- if (!is.null(artifacts)) as.character(artifacts$score_col) else NULL
                    cd <- SummarizedExperiment::colData(result)
                    lab_vec <- if (!is.null(labels_col) && labels_col %in% colnames(cd)) cd[[labels_col]] else NULL
                    score_vec <- if (!is.null(score_col) && score_col %in% colnames(cd)) cd[[score_col]] else NULL
                    idents <- if (!is.null(ActiveIdent(result))) as.character(Idents(result)) else NULL
                    n_cells <- if (!is.null(idents)) length(idents) else ncol(result)
                    idents_vec <- if (is.null(idents)) rep("group_0", n_cells) else {
                        u <- unique(as.character(idents))
                        map <- paste0("group_", seq_along(u))
                        names(map) <- u
                        as.character(map[as.character(idents)])
                    }
                    aggregate_values <- list(n_cells = as.integer(n_cells))
                    if (!is.null(lab_vec)) {
                        lab_vec_ch <- as.character(lab_vec)
                        lab_vec_ch[!nzchar(lab_vec_ch) | is.na(lab_vec_ch)] <- "NA"
                        tb <- table(lab_vec_ch, useNA = "no")
                        top_n <- 15L
                        ord <- order(as.integer(tb), decreasing = TRUE)
                        top_labels <- names(tb)[utils::head(ord, top_n)]
                        top_fractions <- round(as.integer(tb[top_labels]) / n_cells, 4L)
                        aggregate_values$n_labels <- as.integer(length(tb))
                        tl <- rep(NA_integer_, 15L)
                        if (length(top_labels) > 0L) tl[seq_len(min(15L, length(top_labels)))] <- seq_len(min(15L, length(top_labels)))
                        tf <- rep_len(NA_real_, 15L)
                        if (length(top_fractions) > 0L) tf[seq_len(min(15L, length(top_fractions)))] <- top_fractions
                        # `tl` is 1:k by construction (the position in the
                        # ordered list), not a rank -- the same fake field
                        # the marker node carried until 2026-10-20. The
                        # fraction is real; the index is not recorded.
                        aggregate_values$n_top_labels <- as.integer(sum(!is.na(tl)))
                        # (top_label_ranks removed: see n_top_labels above)
                        aggregate_values$top_label_fractions <- as.numeric(tf)
                        if (!is.null(idents) && length(idents) == n_cells) {
                            tb2 <- table(idents_vec, lab_vec_ch, useNA = "no")
                            cluster_totals <- rowSums(tb2); cluster_totals[cluster_totals == 0L] <- NA_integer_
                            prop <- as.data.frame(tb2 / cluster_totals[row(tb2)])
                            colnames(prop) <- c("cluster", "label", "fraction")
                            prop <- prop[prop$fraction > 0.05 & !is.na(prop$fraction), , drop = FALSE]
                            prop <- prop[order(-prop$fraction), , drop = FALSE]
                            keep_n <- 30L; prop <- utils::head(prop, keep_n)
                            n_rows <- as.integer(nrow(prop))
                            aggregate_values$cluster_summary_rows <- n_rows
                            groups_idx <- integer(keep_n)
                            fractions_vec <- numeric(keep_n)
                            group_levels <- unique(idents_vec)
                            idx_map <- seq_along(group_levels)
                            names(idx_map) <- group_levels
                            if (n_rows > 0L) {
                              ok_clusters <- as.character(prop$cluster[seq_len(n_rows)]) %in% names(idx_map)
                              groups_idx[seq_len(n_rows)] <- ifelse(ok_clusters, idx_map[as.character(prop$cluster[seq_len(n_rows)])], NA_integer_)
                              fractions_vec[seq_len(n_rows)] <- round(as.numeric(prop$fraction[seq_len(n_rows)]), 4L)
                            }
                            if (keep_n > n_rows) {
                              groups_idx[seq(n_rows + 1L, keep_n)] <- NA_integer_
                              fractions_vec[seq(n_rows + 1L, keep_n)] <- NA_real_
                            }
                            aggregate_values$cluster_summary_group_index <- as.integer(groups_idx)
                            aggregate_values$cluster_summary_fractions <- as.numeric(fractions_vec)
                        }
                    }
                    if (!is.null(score_vec)) {
                        score_num <- as.numeric(score_vec)
                        aggregate_values$mean_score <- round(mean(score_num, na.rm = TRUE), 4L)
                        qs <- round(stats::quantile(score_num, c(0.25, 0.5, 0.75), na.rm = TRUE), 4L)
                        aggregate_values$score_quantile_25 <- as.numeric(qs[[1L]])
                        aggregate_values$score_quantile_50 <- as.numeric(qs[[2L]])
                        aggregate_values$score_quantile_75 <- as.numeric(qs[[3L]])
                    }
                    # Method/scope integer codes removed 2026-10-20: opaque
                    # integers with no consumer, the same class of dead
                    # field the marker node carried.
                    aggregate_values$raw_values_included <- FALSE
                    evidence <- list(
                        id = sclet_ai_annotation_evidence_id(name),
                        kind = "deterministic_summary",
                        values = aggregate_values,
                        claim_level = "consistent_with"
                    )
                    # The parent is the annotation record id, not an evidence
                    # id, so it never resolves as an ancestor; it is kept
                    # because it also enters the node's dependency_group,
                    # which keeps two different annotation runs from
                    # collapsing into one line of support.
                    result <- RecordAIEvidence(
                        result,
                        evidence,
                        source = name,
                        parents = if (!is.null(ann_rec$id)) as.character(ann_rec$id) else character(),
                        scope = NULL
                    )
                }
                # The output schema only requires that *some* annotation
                # state exists, which an earlier step could have supplied.
                # Verify this action's own record and its bounded evidence
                # node, so a run that annotated nothing -- or a ledger that
                # refused the node -- cannot be reported as completed.
                if (is.null(ann_rec)) {
                    stop(
                        "annotation_no_result: RunSingleR returned without a completed annotation record.",
                        call. = FALSE
                    )
                }
                recorded <- tryCatch(
                    sclet_ai_evidence_get_all(result)[[sclet_ai_annotation_evidence_id(name)]],
                    error = function(e) NULL
                )
                if (is.null(recorded)) {
                    stop(
                        "annotation_evidence_missing: SingleR completed but its bounded evidence node was not recorded.",
                        call. = FALSE
                    )
                }
                result
            },
            input_schema = list(
                ref = list(type = "object", required = TRUE),
                labels = list(type = "object", required = TRUE),
                name = "character",
                layer = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                ref_arg <- params$ref
                labels_arg <- params$labels
                if (missing(ref_arg) || is.null(ref_arg)) {
                    return("reference_missing: annotation requires an explicit 'ref' parameter: either a celldex dataset name (e.g. HumanPrimaryCellAtlasData) or a SummarizedExperiment reference. Do not rely on the package default of HumanPrimaryCellAtlasData.")
                }
                ref_is_name <- is.character(ref_arg) && length(ref_arg) == 1L && nzchar(ref_arg)
                ref_is_object <- inherits(ref_arg, "SummarizedExperiment") ||
                    (is.matrix(ref_arg) && !is.null(rownames(ref_arg)) && !is.null(colnames(ref_arg)))
                if (!ref_is_name && !ref_is_object) {
                    return("reference_invalid: 'ref' must be either a single non-empty string naming a celldex dataset or a SummarizedExperiment/matrix with dimnames.")
                }
                ref_name <- if (ref_is_name) sclet_celldex_reference_name(ref_arg) else NULL
                if (ref_is_name && is.null(ref_name)) {
                    # A placeholder is not a reference; refuse it at
                    # validation instead of after six upstream steps have
                    # run. It is reported as MISSING, not invalid, because
                    # the choice belongs to the user: the clarification
                    # catalog maps `reference_missing` to the structured
                    # question, while `reference_invalid` has no question
                    # and used to discard the whole plan (2026-10-18).
                    if (sclet_ai_value_is_placeholder(ref_arg)) {
                        return(paste0(
                            "reference_missing: 'ref' is a placeholder ('",
                            ref_arg, "'), not a dataset. Ask the user which ",
                            "reference to use (a celldex dataset such as ",
                            "celldex::MonacoImmuneData, or a ",
                            "SummarizedExperiment/matrix with dimnames)."
                        ))
                    }
                    # A genuine near-miss or typo is still invalid: the
                    # accepted spellings are named in the message.
                    return(paste0(
                        "reference_invalid: 'ref' is a string but not a celldex dataset name: '",
                        ref_arg, "'. Accepted: HumanPrimaryCellAtlasData, ",
                        "celldex::MonacoImmuneData, <celldex::MonacoImmuneData>, ",
                        "or a SummarizedExperiment/matrix with dimnames."
                    ))
                }
                if (missing(labels_arg) || is.null(labels_arg)) {
                    return("labels_missing: annotation requires an explicit 'labels' parameter selecting which labels column of the reference to use (e.g. label.main, label.fine) or a vector of per-reference-sample labels.")
                }
                labels_is_name <- is.character(labels_arg) && length(labels_arg) == 1L && nzchar(labels_arg)
                labels_is_vec <- is.atomic(labels_arg) && length(labels_arg) > 0L &&
                    (is.character(labels_arg) || is.factor(labels_arg))
                if (!labels_is_name && !labels_is_vec) {
                    return("labels_invalid: 'labels' must be either a single non-empty string identifying a reference colData column or one label per reference sample (character or factor vector).")
                }
                if (labels_is_name && sclet_ai_value_is_placeholder(labels_arg)) {
                    return(paste0(
                        "labels_missing: 'labels' is a placeholder ('",
                        labels_arg, "'), not a reference colData column. Ask ",
                        "the user which labels column to use (e.g. ",
                        "label.main, label.fine) or a per-sample label vector."
                    ))
                }
                if (labels_is_vec && ref_is_object) {
                    # One resolver shared with RunSingleR(): a labels string
                    # naming a colData column of the reference is the
                    # documented spelling, and counting it as a one-element
                    # label vector rejected it at validation (measured on
                    # 2026-10-17 with a SummarizedExperiment reference).
                    resolved <- sclet_resolve_singleR_labels(ref_arg, labels_arg)
                    if (!is.null(resolved$problem)) {
                        return(paste0("labels_invalid: ", resolved$problem))
                    }
                    ref_n <- ncol(ref_arg)
                    if (length(resolved$labels) != ref_n) {
                        return(paste0("labels_mismatch: length(labels)=", length(resolved$labels), " does not match ncol(ref)=", ref_n))
                    }
                }
                real_req_ns <- base::requireNamespace
                if (!real_req_ns("SingleR", quietly = TRUE)) {
                    return("optional_package_missing: SingleR is required for run_annotation; install via BiocManager::install('SingleR')")
                }
                if (ref_is_name && !is.null(ref_name)) {
                    if (!real_req_ns("celldex", quietly = TRUE)) {
                        return(paste0("optional_package_missing: celldex is required to resolve ref='", ref_arg, "'; install via BiocManager::install('celldex')"))
                    }
                    # `sclet_celldex_reference_name()` accepts the SHAPE of a
                    # dataset name; only celldex knows which names exist. A
                    # near-miss that ends in "Data" (e.g. BlueprintENCODEData
                    # for BlueprintEncodeData) otherwise passes validation and
                    # fails in the handler -- after the upstream steps have
                    # run, which is the class of failure this benchmark paid
                    # for twice.
                    if (is.null(sclet_celldex_dataset(ref_name))) {
                        return(paste0(
                            "reference_invalid: '", ref_arg, "' is not a dataset exported by ",
                            "celldex. Use a real celldex dataset (e.g. HumanPrimaryCellAtlasData, ",
                            "MonacoImmuneData, BlueprintEncodeData) or pass a ",
                            "SummarizedExperiment/matrix reference directly."
                        ))
                    }
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(annotation = "any"),
                note = "Predicted labels are written into colData columns <name>_labels / <name>_pruned.labels and are NEVER used to overwrite Idents(object) or any existing cluster identity."
            ),
            mutates_object = TRUE,
            allowed_state_types = c("annotation", "mapping", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE
        )
    ))
    }

    # ---- reference mapping, non-SingleR backends --------------------------
    #
    # `run_annotation` wraps RunSingleR only. KNN and Symphony are the other two
    # branches of `RunReferenceMapping`, and both need a SECOND SingleCellExperiment
    # as the reference -- a planner cannot supply an object, so both are
    # `planner_visible = FALSE` and used from the driver path.
    if ("annotation" %in% include) {
    actions <- c(actions, list(
        map_reference_knn = AIAction(
            name = "map_reference_knn",
            description = paste(
                "Transfer cell-type labels to this object from a reference",
                "SingleCellExperiment by k-nearest-neighbour projection, registering one",
                "bounded evidence node (label counts). `ref` is a second object supplied",
                "by the caller; this action never selects a reference."
            ),
            handler = function(object, params) {
                name <- params$name %||% "knn_map"
                args <- list(object = object, ref = params$ref, name = name)
                if (!is.null(params$labels)) args$labels <- as.character(params$labels)
                if (!is.null(params$features)) args$features <- as.character(params$features)
                if (!is.null(params$k)) args$k <- as.integer(params$k)
                result <- do.call(RunKNNPredict, args)
                record <- sclet_ai_mainline_record(result, "mapping", name, "map_reference_knn")
                summary <- record$summary %||% list()
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:mapping_knn_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "measured",
                        values = list(
                            n_labels = sclet_ai_mainline_number(summary$n_labels),
                            n_reference_labels = sclet_ai_mainline_number(
                                summary$n_reference_labels
                            ),
                            n_cells = ncol(result)
                        )
                    ),
                    source = name
                )
            },
            input_schema = list(
                ref = list(type = "object", required = TRUE),
                labels = "character_vector",
                features = "character_vector",
                k = "integer",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!inherits(params$ref, "SingleCellExperiment")) {
                    return("reference_missing: map_reference_knn needs `ref`, a SingleCellExperiment to transfer labels from")
                }
                if (!"counts" %in% SummarizedExperiment::assayNames(params$ref) &&
                    is.null(params$features)) {
                    return("reference_counts_missing: the reference has no 'counts' assay; supply `features` (genes shared with it) instead")
                }
                if (!sclet_ai_mapping_has_labels(params$ref, params$labels)) {
                    return("annotation_levels_missing: the reference carries no 'label'/'label.main' column and no `labels` were supplied")
                }
                if (!base::requireNamespace("BiocNeighbors", quietly = TRUE)) {
                    return("optional_package_missing: BiocNeighbors is required for map_reference_knn")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(mapping = "any", ai_evidence = "any"),
                note = "The transferred labels land in a colData column; the evidence node carries counts only."
            ),
            mutates_object = TRUE,
            allowed_state_types = c("mapping", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE,
            planner_visible = FALSE
        ),
        map_reference_symphony = AIAction(
            name = "map_reference_symphony",
            description = paste(
                "Transfer cell-type labels with Symphony's reference-based embedding,",
                "registering one bounded evidence node (label counts). `ref` is a second",
                "SingleCellExperiment; Symphony repackages it internally. Driver-only: a",
                "personal reference is not something a model can supply."
            ),
            handler = function(object, params) {
                name <- params$name %||% "symphony_map"
                args <- list(object = object, ref = params$ref, name = name)
                if (!is.null(params$labels)) args$labels <- as.character(params$labels)
                if (!is.null(params$vars)) args$vars <- as.character(params$vars)
                if (!is.null(params$k)) args$k <- as.integer(params$k)
                result <- do.call(RunSymphonyMapping, args)
                record <- sclet_ai_mainline_record(result, "mapping", name,
                    "map_reference_symphony")
                summary <- record$summary %||% list()
                # RunSymphonyMapping records n_query_cells / n_ref_cells /
                # n_common_genes (R/annotation.R:660-663); there is no n_labels
                # on the record, so it is read from the predicted column.
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:mapping_symphony_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "measured",
                        values = list(
                            n_query_cells = sclet_ai_mainline_number(summary$n_query_cells),
                            n_ref_cells = sclet_ai_mainline_number(summary$n_ref_cells),
                            n_cells = ncol(result)
                        )
                    ),
                    source = name
                )
            },
            input_schema = list(
                ref = list(type = "object", required = TRUE),
                labels = "character_vector",
                vars = "character_vector",
                k = "integer",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!inherits(params$ref, "SingleCellExperiment")) {
                    return("reference_missing: map_reference_symphony needs `ref`, a SingleCellExperiment to build the reference from")
                }
                if (!sclet_ai_mapping_has_labels(params$ref, params$labels)) {
                    return("annotation_levels_missing: the reference carries no 'label'/'label.main' column and no `labels` were supplied")
                }
                if (!base::requireNamespace("symphony", quietly = TRUE)) {
                    return("optional_package_missing: symphony is required for map_reference_symphony; install via remotes::install_github('immunogenomics/symphony')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(mapping = "any", ai_evidence = "any"),
                note = "Transferred labels land in a colData column; the evidence node carries counts only."
            ),
            mutates_object = TRUE,
            allowed_state_types = c("mapping", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE,
            planner_visible = FALSE
        )
    ))
    }

    actions
}
