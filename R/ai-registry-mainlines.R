# sclet_ai_registry_mainlines(): AI actions for the six mainlines that had none.
#
# A registered action is not the same as a plannable one. Five of these six need
# an input a language model cannot supply -- a gene set, a database path, an
# external phenotype matrix, a base GRN -- so they are marked
# `planner_visible = FALSE` and stay out of the planning catalog. They are still
# validated and executed for the driver path, which is a supported route
# (inst/skills/sclet-driver/SKILL.md). `run_perturbation_priority` is the one
# exception: its external input is a colData column name the planner can read
# from `capabilities`.
#
# Every handler follows the same shape: run the exported function, read the
# record it just wrote, register one evidence node from that record with numeric
# or logical values only (the evidence value whitelist rejects free text), and
# return the object.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R.

# The record one of these functions just wrote, or a typed error saying the
# function returned without one. Silence here would let a run look successful
# while the ledger says nothing happened.
sclet_ai_mainline_record <- function(object, type, id, action) {
    records <- tryCatch(sclet_get_state_records(object, type), error = function(e) list())
    record <- records[[as.character(id)]]
    if (is.null(record)) {
        sclet_ai_error(
            paste0(type, "_no_result"),
            paste0(
                action, " returned without registering a '", type, "' record '",
                id, "'. The evidence node is not registered, and this step ",
                "supports no claim. Check the function's own output for why it ",
                "recorded nothing."
            )
        )
    }
    record
}

# The identity a plan declares it will have created by an earlier step.
#
# Validation runs before execution, so `find_clusters -> run_perturbation_priority`
# has no `ActiveIdent(object)` yet. The planned projection carries the promise
# (`output_schema$active_ident`), and every other gate-consuming action already
# reads it; these two must as well or the only planner-visible new action is
# unusable in the round that creates its prerequisite (measured 2026-10-23).
sclet_ai_mainline_planned_ident <- function(planned) {
    ident <- (planned %||% list())$active_ident
    if (is.null(ident) || !length(ident)) {
        return(NULL)
    }
    as.character(ident)[[1L]]
}

# Does a reference object carry cell-type labels, or did the caller supply them?
sclet_ai_mapping_has_labels <- function(reference, labels = NULL) {
    if (length(labels)) {
        return(TRUE)
    }
    columns <- colnames(SummarizedExperiment::colData(reference))
    any(c("label", "label.main") %in% columns)
}

# The colData column that actually holds the cell identity.
#
# `ActiveIdent()` returns the SENTINEL "colLabels" when the identity is stored
# as colLabels, while `colData(sce)[["colLabels"]]` is empty -- SingleCellExperiment
# keeps them in the "label" column. `RunInSilicoPerturbation()`'s default
# `cluster_key = ActiveIdent(sce)` therefore resolves to a column that does not
# exist on the object this package produces (measured 2026-10-23), and the
# failure surfaces inside Python as a zero-length label vector. The action
# resolves the column itself and passes it explicitly.
sclet_ai_mainline_identity_column <- function(object, key = NULL) {
    columns <- colnames(SummarizedExperiment::colData(object))
    if (!is.null(key)) {
        if (is.character(key) && length(key) == 1L && !is.na(key) && key %in% columns) {
            return(key)
        }
        return(NULL)
    }
    ident <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (is.null(ident) || !length(ident)) {
        return(NULL)
    }
    ident <- as.character(ident)[[1L]]
    if (identical(ident, "colLabels")) {
        labels <- tryCatch(SingleCellExperiment::colLabels(object), error = function(e) NULL)
        if (!length(labels)) {
            return(NULL)
        }
        return(if ("label" %in% columns) "label" else "colLabels")
    }
    if (ident %in% columns) ident else NULL
}

sclet_ai_record_mainline_evidence <- function(object, evidence, source,
                                              replace = FALSE) {
    sclet_ai_record_evidence_for_run(object, evidence, source = source,
        replace = replace)
}

sclet_ai_mainline_number <- function(value, default = NA_real_) {
    if (is.null(value) || !length(value)) {
        return(default)
    }
    out <- suppressWarnings(as.numeric(value)[[1L]])
    if (is.na(out)) default else out
}

sclet_ai_registry_mainlines <- function(object, include) {
    actions <- list()

    # ---- gene-set scoring -------------------------------------------------
    if ("geneset" %in% include) {
    actions <- c(actions, list(
        run_geneset_scoring = AIAction(
            name = "run_geneset_scoring",
            description = paste(
                "Score the cells for user-supplied gene sets with UCell, AUCell or GSVA",
                "and register one bounded evidence node (number of sets, cells and",
                "genes used). The gene sets are the user's modelling choice and are a",
                "required parameter: this action never picks or recommends a set, and",
                "no set, gene or score vector enters the evidence node."
            ),
            handler = function(object, params) {
                gene_sets <- params$gene_sets
                if (!is.list(gene_sets) || !length(gene_sets)) {
                    sclet_ai_error(
                        "gene_sets_empty",
                        "`gene_sets` must be a non-empty named list of character vectors."
                    )
                }
                name <- params$name %||% "Score"
                args <- list(
                    sce = object,
                    gene_sets = gene_sets,
                    method = params$method %||% "UCell",
                    name = name
                )
                if (!is.null(params$assay_use)) args$assay_use <- params$assay_use
                if (!is.null(params$ncores)) args$ncores <- as.integer(params$ncores)
                if (!is.null(params$as_altExp)) args$as_altExp <- isTRUE(params$as_altExp)
                result <- do.call(RunGeneSetScoring, args)
                sclet_ai_mainline_record(result, "geneset_scoring", name, "run_geneset_scoring")
                genes <- unique(unlist(gene_sets, use.names = FALSE))
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:geneset_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "measured",
                        values = list(
                            n_sets = length(gene_sets),
                            n_cells = ncol(result),
                            n_genes_used = sum(genes %in% rownames(object))
                        )
                    ),
                    source = name
                )
            },
            input_schema = list(
                gene_sets = list(type = "list", required = TRUE),
                method = c("UCell", "AUCell", "GSVA"),
                assay_use = "character",
                ncores = "integer",
                name = "character",
                as_altExp = "logical"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!is.list(params$gene_sets) || !length(params$gene_sets)) {
                    return("gene_sets_empty: run_geneset_scoring needs `gene_sets` as a named list of character vectors (the user's gene sets, not a guess)")
                }
                genes <- unique(unlist(params$gene_sets, use.names = FALSE))
                if (!any(genes %in% rownames(object))) {
                    return("gene_sets_not_found: none of the supplied genes is in rownames(object); check the identifier type (symbol vs Ensembl)")
                }
                method <- params$method %||% "UCell"
                package <- switch(method, UCell = "UCell", AUCell = "AUCell", GSVA = "GSVA")
                if (!base::requireNamespace(package, quietly = TRUE)) {
                    return(paste0("geneset_package_missing: ", package, " is required for method '", method, "'; install it before running this action"))
                }
                assay_use <- params$assay_use
                if (!is.null(assay_use) &&
                    !as.character(assay_use)[[1L]] %in% SummarizedExperiment::assayNames(object)) {
                    return(paste0(
                        "assay_missing: `assay_use` names '", assay_use,
                        "', which is not an assay on the object"
                    ))
                }
                if (is.null(assay_use) &&
                    !"counts" %in% SummarizedExperiment::assayNames(object)) {
                    return("counts_assay_missing: run_geneset_scoring scores a counts assay by default and none is present; name an existing assay with assay_use")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(geneset_scoring = "any", ai_evidence = "any"),
                note = paste(
                    "Per-cell scores are added to colData or an altExp. The evidence node",
                    "carries counts only: no gene set name, gene symbol or score vector."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("geneset_scoring", "ai_evidence"),
            estimated_cost = "low",
            idempotent = FALSE,
            planner_visible = FALSE
        )
    ))
    }

    # ---- metacells --------------------------------------------------------
    if ("metacell" %in% include) {
    actions <- c(actions, list(
        run_supercell = AIAction(
            name = "run_supercell",
            description = paste(
                "Aggregate the cells into metacells with SuperCell and register one",
                "bounded evidence node (metacell count, cell count, gamma, k.knn).",
                "IMPORTANT: this REPLACES the object -- the returned SingleCellExperiment",
                "is the metacell object, and the previous analysis records and ledger are",
                "not carried over. Use it at the end of a workflow or on a copy, never as",
                "a mid-plan step."
            ),
            handler = function(object, params) {
                name <- params$name %||% "supercell"
                args <- list(object = object, name = name)
                if (!is.null(params$assay)) args$assay <- params$assay
                if (!is.null(params$layer)) args$layer <- params$layer
                if (!is.null(params$nHVG)) args$nHVG <- as.integer(params$nHVG)
                if (!is.null(params$hvg_method)) args$hvg_method <- params$hvg_method
                if (!is.null(params$gamma)) args$gamma <- as.numeric(params$gamma)
                if (!is.null(params$k.knn)) args$k.knn <- as.integer(params$k.knn)
                result <- do.call(RunSuperCell, args)
                record <- sclet_ai_mainline_record(result, "aggregation", name, "run_supercell")
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:supercell_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "measured",
                        values = list(
                            n_metacells = sclet_ai_mainline_number(record$summary$n_metacells),
                            n_cells = sclet_ai_mainline_number(record$inputs$parent_n_cells),
                            gamma = sclet_ai_mainline_number(record$params$gamma),
                            k_knn = sclet_ai_mainline_number(record$params$k.knn)
                        )
                    ),
                    source = name
                )
            },
            input_schema = list(
                assay = "character",
                layer = "character",
                nHVG = "integer",
                hvg_method = c("scran", "seurat", "scrapper"),
                gamma = "number",
                k.knn = "integer",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                assay <- params$assay %||% "logcounts"
                if (!assay %in% SummarizedExperiment::assayNames(object)) {
                    return(paste0(
                        "logcounts_assay_missing: RunSuperCell expects a normalized '",
                        assay, "' assay and none is present; run NormalizeData() first, ",
                        "or name an existing assay with `assay`"
                    ))
                }
                if (!base::requireNamespace("SuperCell", quietly = TRUE)) {
                    return("optional_package_missing: SuperCell is required for run_supercell; install it from GitHub (supercell::install_supercell())")
                }
                gamma <- sclet_ai_mainline_number(params$gamma, 20)
                if (ncol(object) < 2 * gamma) {
                    return(paste0(
                        "too_few_metacells: ", ncol(object), " cells with gamma = ", gamma,
                        " gives fewer than two metacells; lower gamma or use a larger object"
                    ))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(aggregation = "any", ai_evidence = "any"),
                note = paste(
                    "The returned object is the metacell object. The previous state",
                    "records are not preserved, so run this last or on a copy."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("aggregation", "ai_evidence"),
            estimated_cost = "low",
            idempotent = FALSE,
            planner_visible = FALSE,
            replaces_object = TRUE
        )
    ))
    }

    # ---- perturbation priority (Augur) -----------------------------------
    if ("priority" %in% include) {
    actions <- c(actions, list(
        run_perturbation_priority = AIAction(
            name = "run_perturbation_priority",
            description = paste(
                "Rank cell types by how strongly their transcriptome separates the",
                "conditions (Augur), and register one bounded evidence node (number of",
                "cell types plus AUC quantiles, claim level 'associated'). The condition",
                "column is the user's design and must already exist in colData; the",
                "action does not infer which column it is."
            ),
            handler = function(object, params) {
                name <- params$name %||% "augur"
                args <- list(sce = object, condition_col = params$condition_col, name = name)
                if (!is.null(params$group_col)) args$group_col <- params$group_col
                if (!is.null(params$assay)) args$assay <- params$assay
                if (!is.null(params$layer)) args$layer <- params$layer
                result <- do.call(RunPerturbationPriority, args)
                record <- sclet_ai_mainline_record(result, "priority", name, "run_perturbation_priority")
                auc <- record$artifacts$AUC
                if (!is.data.frame(auc) || !"auc" %in% names(auc) || !nrow(auc)) {
                    sclet_ai_error(
                        "priority_no_result",
                        paste(
                            "run_perturbation_priority wrote a priority record with no",
                            "usable AUC table, so no evidence node is registered and this",
                            "step supports no claim. Check the Augur output in the record."
                        )
                    )
                }
                values <- if (is.data.frame(auc) && "auc" %in% names(auc) && nrow(auc)) {
                    auc_values <- suppressWarnings(as.numeric(auc$auc))
                    auc_values <- auc_values[is.finite(auc_values)]
                    list(
                        n_cell_types = nrow(auc),
                        auc_min = if (length(auc_values)) min(auc_values) else NA_real_,
                        auc_median = if (length(auc_values)) stats::median(auc_values) else NA_real_,
                        auc_max = if (length(auc_values)) max(auc_values) else NA_real_
                    )
                } else {
                    list(n_cell_types = NA_real_, auc_min = NA_real_,
                         auc_median = NA_real_, auc_max = NA_real_)
                }
                # `n_cell_types` is real even when no AUC is finite; an evidence
                # node is not improved by dropping it. No filter here: a node
                # with zero values is a node that says nothing.
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:augur_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "associated",
                        values = values
                    ),
                    source = name
                )
            },
            input_schema = list(
                condition_col = list(type = "character", required = TRUE),
                group_col = "character",
                assay = "character",
                layer = "character",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                condition_col <- params$condition_col
                columns <- colnames(SummarizedExperiment::colData(object))
                if (is.null(condition_col) || !is.character(condition_col) ||
                    length(condition_col) != 1L || !nzchar(condition_col) ||
                    !condition_col %in% columns) {
                    return(paste0(
                        "condition_col_missing: run_perturbation_priority needs `condition_col` naming a colData column with the two conditions; available: ",
                        if (length(columns)) paste(columns, collapse = ", ") else "(none)"
                    ))
                }
                levels <- unique(as.character(SummarizedExperiment::colData(object)[[condition_col]]))
                if (length(levels) < 2L) {
                    return(paste0(
                        "condition_levels_insufficient: '", condition_col, "' has ",
                        length(levels), " distinct value(s); Augur needs at least two conditions"
                    ))
                }
                group_col <- params$group_col
                if (!is.null(group_col)) {
                    if (!is.character(group_col) || length(group_col) != 1L ||
                        !group_col %in% columns) {
                        return(paste0(
                            "group_col_unknown: `group_col` must name a colData column; available: ",
                            paste(columns, collapse = ", ")
                        ))
                    }
                } else if (is.null(sclet_ai_mainline_planned_ident(planned)) &&
                    is.null(tryCatch(ActiveIdent(object), error = function(e) NULL))) {
                    # The callee falls back to ActiveIdent(), so the gate must too.
                    # `planned$active_ident` covers the round that CREATES the
                    # identity (find_clusters -> run_perturbation_priority): the
                    # clustering has not executed at validation time.
                    return("group_column_missing: run_perturbation_priority needs a cell-type grouping; supply `group_col` or run clustering so an identity exists")
                }
                if (!base::requireNamespace("Augur", quietly = TRUE)) {
                    return("optional_package_missing: Augur is required for run_perturbation_priority; install via remotes::install_github('neurorestore/Augur')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(priority = "any", ai_evidence = "any"),
                note = paste(
                    "AUC per cell type is stored in the state record. The evidence node",
                    "carries counts and AUC quantiles at claim level 'associated': a",
                    "separability score is not a measure of perturbation response."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("priority", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE,
            planner_visible = TRUE
        )
    ))
    }

    # ---- phenotype association (scPAS) -----------------------------------
    if ("phenotype" %in% include) {
    actions <- c(actions, list(
        run_phenotype_association = AIAction(
            name = "run_phenotype_association",
            description = paste(
                "Associate cell-type abundance with an external phenotype table through",
                "scPAS, and register one bounded evidence node (number of added columns,",
                "cell count, external_inputs = TRUE at claim level 'observed'). Both",
                "inputs are external data supplied by the user as parameters; this",
                "action cannot produce them from the object."
            ),
            handler = function(object, params) {
                args <- list(
                    sce = object,
                    bulk_matrix = params$bulk_matrix,
                    phenotype_data = params$phenotype_data
                )
                if (!is.null(params$family)) args$family <- params$family
                result <- do.call(RunPhenotypeAssociation, args)
                record <- sclet_ai_mainline_record(result, "phenotype", "scPAS", "run_phenotype_association")
                added <- record$summary$added_columns
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = "ev:phenotype_association",
                        kind = "deterministic_summary",
                        claim_level = "observed",
                        values = list(
                            n_columns_added = length(added %||% character()),
                            n_cells = ncol(result),
                            external_inputs = TRUE
                        )
                    ),
                    source = "scPAS"
                )
            },
            input_schema = list(
                bulk_matrix = list(type = "object", required = TRUE),
                phenotype_data = list(type = "object", required = TRUE),
                family = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!base::requireNamespace("scPAS", quietly = TRUE)) {
                    return("optional_package_missing: scPAS is required for run_phenotype_association; install via remotes::install_github('aiminXie/scPAS')")
                }
                if (!base::requireNamespace("Seurat", quietly = TRUE)) {
                    return("optional_package_missing: Seurat is required because scPAS works on a Seurat object")
                }
                if (is.null(params$bulk_matrix)) {
                    return("bulk_matrix_missing: run_phenotype_association needs `bulk_matrix` (the bulk expression matrix) as a parameter")
                }
                if (is.null(params$phenotype_data)) {
                    return("phenotype_data_missing: run_phenotype_association needs `phenotype_data` (the phenotype table) as a parameter")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(phenotype = "any", ai_evidence = "any"),
                note = paste(
                    "Association results are added to colData. The evidence node records",
                    "only how many columns were added and the cell count, at claim level",
                    "'observed': an association is not a causal effect."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("phenotype", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE,
            planner_visible = FALSE
        )
    ))
    }

    # ---- regulon inference (pySCENIC) ------------------------------------
    if ("regulon" %in% include) {
    actions <- c(actions, list(
        run_scenic = AIAction(
            name = "run_scenic",
            description = paste(
                "Run pySCENIC with user-supplied TF, motif-annotation and database files,",
                "store per-regulon AUC in an altExp, and register one bounded evidence",
                "node (regulon count, cell count, python_backend = TRUE). The three file",
                "inputs are required parameters and are expensive to obtain; this action",
                "runs for hours, so it is a deliberate, user-initiated step."
            ),
            handler = function(object, params) {
                name <- params$name %||% "scenic"
                args <- list(
                    sce = object,
                    tfs_path = params$tfs_path,
                    motif_annotations_path = params$motif_annotations_path,
                    database_paths = params$database_paths,
                    name = name
                )
                if (!is.null(params$assay_use)) args$assay_use <- params$assay_use
                if (!is.null(params$num_workers)) args$num_workers <- as.integer(params$num_workers)
                if (!is.null(params$seed)) args$seed <- as.integer(params$seed)
                if (!is.null(params$output_dir)) args$output_dir <- params$output_dir
                result <- do.call(RunSCENIC, args)
                record <- sclet_ai_mainline_record(result, "scenic", name, "run_scenic")
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:scenic_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "measured",
                        values = list(
                            n_regulons = sclet_ai_mainline_number(record$summary$n_regulons),
                            n_cells = sclet_ai_mainline_number(record$summary$n_cells),
                            python_backend = TRUE
                        )
                    ),
                    source = name
                )
            },
            input_schema = list(
                tfs_path = list(type = "character", required = TRUE),
                motif_annotations_path = list(type = "character", required = TRUE),
                database_paths = list(type = "character_vector", required = TRUE),
                assay_use = "character",
                num_workers = "integer",
                seed = "integer",
                name = "character",
                output_dir = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!base::requireNamespace("basilisk", quietly = TRUE)) {
                    return("optional_package_missing: basilisk is required for run_scenic (the Python environment)")
                }
                assay_use <- params$assay_use %||% "counts"
                if (!as.character(assay_use)[[1L]] %in% SummarizedExperiment::assayNames(object)) {
                    return(paste0(
                        "assay_missing: run_scenic needs assay '", assay_use,
                        "', which is not present on the object"
                    ))
                }
                for (field in c("tfs_path", "motif_annotations_path")) {
                    value <- params[[field]]
                    if (is.null(value) || !is.character(value) || length(value) != 1L ||
                        is.na(value) || !nzchar(value)) {
                        return(paste0("scenic_files_missing: `", field, "` is required and must name an existing file"))
                    }
                    if (!file.exists(value)) {
                        return(paste0("scenic_files_missing: `", field, "` does not exist: ", value))
                    }
                }
                databases <- params$database_paths
                if (!length(databases)) {
                    return("scenic_database_missing: `database_paths` must name at least one motif database file")
                }
                missing <- databases[!file.exists(as.character(databases))]
                if (length(missing)) {
                    return(paste0("scenic_files_missing: database file(s) not found: ", paste(missing, collapse = ", ")))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(scenic = "any", ai_evidence = "any"),
                note = paste(
                    "Per-regulon AUC is stored in an altExp. The evidence node carries",
                    "counts only: no regulon name, TF name or database path."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("scenic", "ai_evidence"),
            estimated_cost = "high",
            idempotent = FALSE,
            planner_visible = FALSE
        )
    ))
    }

    # ---- in-silico perturbation (CellOracle) ------------------------------
    if ("perturbation" %in% include) {
    actions <- c(actions, list(
        run_in_silico_perturbation = AIAction(
            name = "run_in_silico_perturbation",
            description = paste(
                "Simulate knocking down one user-named gene with CellOracle and register",
                "one bounded evidence node (cell count, perturbation value,",
                "simulated = TRUE at claim level 'associated'). The target gene and the",
                "base GRN file are required parameters; this is a simulation, not a",
                "measured perturbation, and the package never proposes a target."
            ),
            handler = function(object, params) {
                args <- list(
                    sce = object,
                    target_gene = params$target_gene,
                    base_grn_path = params$base_grn_path
                )
                if (!is.null(params$perturbation_value)) {
                    args$perturbation_value <- as.numeric(params$perturbation_value)
                }
                # Never rely on the callee default: `ActiveIdent()` may be the
                # "colLabels" sentinel, which is not a colData column here.
                identity <- sclet_ai_mainline_identity_column(object, params$cluster_key)
                if (!is.null(identity)) args$cluster_key <- identity
                if (!is.null(params$reduction)) args$reduction <- params$reduction
                result <- do.call(RunInSilicoPerturbation, args)
                id <- paste0("celloracle_", params$target_gene)
                sclet_ai_mainline_record(result, "perturbation", id, "run_in_silico_perturbation")
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:perturbation_", sclet_ai_evidence_id_token(params$target_gene)),
                        kind = "deterministic_summary",
                        claim_level = "associated",
                        values = list(
                            n_cells = ncol(result),
                            perturbation_value = sclet_ai_mainline_number(params$perturbation_value, 0),
                            simulated = TRUE
                        )
                    ),
                    source = id
                )
            },
            input_schema = list(
                target_gene = list(type = "character", required = TRUE),
                base_grn_path = list(type = "character", required = TRUE),
                perturbation_value = "number",
                cluster_key = "character",
                reduction = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!base::requireNamespace("basilisk", quietly = TRUE)) {
                    return("optional_package_missing: basilisk is required for run_in_silico_perturbation (the Python environment)")
                }
                if (!"counts" %in% SummarizedExperiment::assayNames(object)) {
                    return("counts_assay_missing: run_in_silico_perturbation requires a counts assay")
                }
                target <- params$target_gene
                if (is.null(target) || !is.character(target) || length(target) != 1L ||
                    !nzchar(target) || !target %in% rownames(object)) {
                    return(paste0(
                        "target_gene_not_found: `target_gene` must name one gene in rownames(object)",
                        if (is.character(target) && length(target) == 1L && nzchar(target)) paste0("; '", target, "' is not present") else ""
                    ))
                }
                if (is.null(params$base_grn_path) || !is.character(params$base_grn_path) ||
                    length(params$base_grn_path) != 1L || !file.exists(params$base_grn_path)) {
                    return("base_grn_missing: `base_grn_path` is required and must name an existing base GRN file")
                }
                cluster_key <- params$cluster_key
                if (!is.null(cluster_key)) {
                    if (!is.character(cluster_key) || length(cluster_key) != 1L ||
                        !cluster_key %in% colnames(SummarizedExperiment::colData(object))) {
                        return(paste0(
                            "cluster_key_unknown: `cluster_key` must name a colData column; available: ",
                            paste(colnames(SummarizedExperiment::colData(object)), collapse = ", ")
                        ))
                    }
                } else if (is.null(sclet_ai_mainline_identity_column(object)) &&
                    is.null(sclet_ai_mainline_planned_ident(planned))) {
                    return("cluster_identity_missing: run_in_silico_perturbation needs a cell grouping; run clustering or name `cluster_key`")
                }
                # Validation runs before execution, so a reduction an earlier
                # step of the SAME plan will create is in `planned`, not yet on
                # the object (`run_pca` declares active_reduction "PCA").
                planned_reductions <- as.character(
                    (planned %||% list())$reductions %||% character()
                )
                planned_active <- as.character(
                    (planned %||% list())$active_reduction %||% character()
                )
                reduction <- params$reduction
                if (is.null(reduction) || !length(reduction)) {
                    reduction <- planned_active
                }
                if (is.null(reduction) || !length(reduction)) {
                    reduction <- "PCA"
                }
                available <- unique(c(
                    SingleCellExperiment::reducedDimNames(object),
                    planned_reductions, planned_active
                ))
                if (!reduction %in% available) {
                    return(paste0(
                        "reduction_missing: '", reduction,
                        "' is not a reduction on the object; run the dimensionality reduction or name an existing one"
                    ))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(perturbation = "any", ai_evidence = "any"),
                note = paste(
                    "The simulated shift is stored in the perturbation record. The",
                    "evidence node is claim level 'associated' with simulated = TRUE: a",
                    "model prediction, not an observed perturbation response."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("perturbation", "ai_evidence"),
            estimated_cost = "high",
            idempotent = FALSE,
            planner_visible = FALSE
        )
    ))
    }

    actions
}
