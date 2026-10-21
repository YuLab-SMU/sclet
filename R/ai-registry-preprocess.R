# sclet_ai_registry_preprocess(): built-in AI actions for normalisation, variable-feature and scaling.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# These wrap the basic workflow functions and return the modified object.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_preprocess <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("preprocess" %in% include) {
    actions <- c(actions, list(
        normalize_data = AIAction(
            name = "normalize_data",
            description = "Normalize an assay into the sclet logcounts layer.",
            handler = function(object, params) NormalizeData(
                object,
                scale.factor = params$scale.factor %||% 10000,
                assay = params$assay %||% "counts"
            ),
            input_schema = list(
                scale.factor = "number",
                assay = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                assay <- params$assay %||% "counts"
                available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                if (!assay %in% available) {
                    return(paste0("assay is not available: ", assay))
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(
                assay = "logcounts",
                required_assays = "logcounts",
                required_states = list(preprocess = "normalize_logcounts"),
                active_assay = "logcounts"
            ),
            allowed_state_types = "preprocess",
            estimated_cost = "medium",
            idempotent = TRUE
        ),
        find_variable_features = AIAction(
            name = "find_variable_features",
            description = "Identify highly variable genes using sclet's variance model.",
            handler = function(object, params) FindVariableFeatures(
                object,
                nfeatures = params$nfeatures %||% 2000,
                method = params$method %||% "scran"
            ),
            input_schema = list(
                nfeatures = "integer",
                method = c("scran", "scrapper", "seurat")
            ),
            prerequisites = function(object, params, planned = NULL) {
                available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                if (!"logcounts" %in% available) {
                    return("logcounts assay is required; run normalize_data first")
                }
                TRUE
            },
            returns = "sce",
            output_schema = list(
                required_commands = "FindVariableFeatures",
                required_hvg = TRUE,
                # A referenceable name for the set itself. The static value
                # is a marker, not the genes: `output_schema` is exposed
                # statically through `available_actions` before any step has
                # run, so promising real gene names here would push 2000 of
                # them into every planning context. The real set is written
                # into the execution summary by
                # `sclet_ai_execution_summary()` once it exists, and the
                # projection resolves ${step.output.variable_features} from
                # there. Downstream steps (scale_data$features,
                # run_pca$subset_row) are `character_vector`, so a prose
                # description used to pass the schema check and die in the
                # gene selector -- measured as 4 rejections in 4 spellings
                # on pbmc3k (sclet-benchmark-ai, execute_plan).
                variable_features = "hvg"
            ),
            allowed_state_types = character(),
            estimated_cost = "medium",
            idempotent = TRUE
        ),
        scale_data = AIAction(
            name = "scale_data",
            description = "Scale an assay and register the scaled layer.",
            handler = function(object, params) ScaleData(
                object,
                features = sclet_ai_gene_selector(
                    object,
                    sclet_ai_resolve_deferred_selector(object, params$features, "features"),
                    arg = "features"
                ),
                assay = params$assay %||% "logcounts"
            ),
            input_schema = list(
                features = "character_vector",
                assay = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                assay <- params$assay %||% "logcounts"
                available <- planned$assays %||% SummarizedExperiment::assayNames(object)
                if (!assay %in% available) {
                    return(paste0("assay is not available: ", assay))
                }
                # features is a character_vector, so a prose description
                # passes the schema check. Resolving it here -- the same way
                # the handler will -- means the model is told at validation
                # time and the repair loop can fix it, instead of the plan
                # validating and then dying at execution.
                if (sclet_ai_is_deferred_gene_selector(params$features)) {
                    return(TRUE)
                }
                sclet_ai_selector_problem(object, params$features, "features")
            },
            returns = "sce",
            output_schema = list(
                required_assays = "scaled",
                active_assay = "scaled"
            ),
            allowed_state_types = character(),
            estimated_cost = "medium",
            idempotent = TRUE
        ),
        run_decontx = AIAction(
            name = "run_decontx",
            description = paste(
                "Estimate and subtract ambient RNA contamination with celda::decontX,",
                "writing a cleaned assay (decontXcounts) and per-cell contamination",
                "scores. Labels and cleans; removes no cell."
            ),
            handler = function(object, params) {
                name <- params$name %||% "decontX"
                decontaminated <- RunDecontX(
                    object,
                    assay = params$assay %||% "counts"
                )
                # `RunDecontX()` writes its preprocess record under the fixed id
                # "decontX", and the recorder requires `source` to name a record
                # that actually exists. Passing the evidence-node name here made
                # every run fail with "source must identify a completed
                # analysis" -- the node id and the record id are two different
                # things. Node ids can differ (name =), the source cannot.
                sclet_ai_record_decontx_evidence(decontaminated, name = name,
                    replace = isTRUE(params$replace_evidence))
            },
            input_schema = list(
                assay = "character",
                name = "character",
                # Re-running decontX on the same object reproduces the node id;
                # refresh it instead of refusing, off by default so two runs that
                # disagree stay visible.
                replace_evidence = "logical"
            ),
            prerequisites = function(object, params, planned = NULL) {
                assay <- params$assay %||% "counts"
                available <- planned$assays %||%
                    tryCatch(SummarizedExperiment::assayNames(object),
                        error = function(e) character())
                if (!assay %in% available) {
                    return(paste0(
                        "assay is not available: ", assay,
                        ". decontX estimates the ambient profile from raw counts,",
                        " so it needs a counts assay; available: ",
                        paste(available, collapse = ", ")
                    ))
                }
                # A log-normalised assay stored under `counts` would produce a
                # contamination estimate about a matrix decontX was not designed
                # for, and the number reads the same either way.
                values <- tryCatch(SummarizedExperiment::assay(object, assay),
                    error = function(e) NULL)
                if (!is.null(values)) {
                    sample_values <- utils::head(as.numeric(
                        if (inherits(values, "dgCMatrix") || inherits(values, "dgTMatrix")) {
                            values@x
                        } else {
                            as.vector(values)
                        }
                    ), 200L)
                    sample_values <- sample_values[!is.na(sample_values)]
                    if (length(sample_values) &&
                        !all(sample_values >= 0) ) {
                        return(paste0(
                            "counts_assay_not_count: assay '", assay,
                            "' holds negative values, so it is not a raw count",
                            " matrix; decontX needs non-negative integers"
                        ))
                    }
                    if (length(sample_values) &&
                        any(abs(sample_values - round(sample_values)) > 1e-8)) {
                        return(paste0(
                            "counts_assay_not_count: assay '", assay,
                            "' holds non-integer values; decontX estimates the",
                            " ambient RNA profile from raw integer counts"
                        ))
                    }
                }
                if (!base::requireNamespace("celda", quietly = TRUE)) {
                    return("optional_package_missing: celda is required for decontX; install via BiocManager::install('celda')")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_assays = "decontXcounts",
                required_columns = "decontX_contamination",
                note = paste(
                    "Adds a cleaned assay (decontXcounts) and per-cell contamination",
                    "scores (colData decontX_contamination). The original counts assay",
                    "is preserved and NO cell is removed, filtered or relabeled.",
                    "Contamination is an EM estimate of where counts came from, not a",
                    "measured property of the cell, so evidence is 'associated'.",
                    "Downstream steps that should use the cleaned matrix must name",
                    "decontXcounts explicitly: the active assay is NOT switched."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("preprocess", "ai_evidence"),
            estimated_cost = "medium",
            idempotent = FALSE
        )
    ))
    }

    actions
}
