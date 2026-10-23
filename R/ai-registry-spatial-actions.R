# sclet_ai_registry_spatial_actions(): EXECUTING spatial actions.
#
# Separate from `ai-registry-spatial-multimodal.R` on purpose: that file is the
# read-only catalog (readiness / evidence summaries / comparisons, none of them
# mutating). Mixing a mutating action in would destroy the property that a test
# asserts about that group. `include = "spatial"` now yields both.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R.
sclet_ai_registry_spatial_actions <- function(object, include) {
    actions <- list()

    # ---- cell2location deconvolution --------------------------------------
    #
    # Needs a SECOND object (the scRNA-seq reference), which a planner cannot
    # supply: `planner_visible = FALSE`, driver path only. The spatial object keeps
    # its identity; this action does not replace it.
    if ("spatial" %in% include) {
    actions <- c(actions, list(
        run_spatial_deconvolution = AIAction(
            name = "run_spatial_deconvolution",
            description = paste(
                "Estimate per-spot cell-type abundance with cell2location against a",
                "scRNA-seq reference SingleCellExperiment, registering one bounded",
                "evidence node (common genes, abundance columns, spot count). The",
                "reference and the reference's group column are supplied by the caller;",
                "this action never chooses a reference. It runs a Python model through",
                "basilisk and it is computationally intensive."
            ),
            handler = function(object, params) {
                name <- params$name %||% "spatial_deconv"
                args <- list(
                    sce_spatial = object,
                    sce_ref = params$sce_ref,
                    ref_group_key = params$ref_group_key,
                    name = name
                )
                if (!is.null(params$ref_batch_key)) args$ref_batch_key <- params$ref_batch_key
                if (!is.null(params$max_epochs)) args$max_epochs <- as.integer(params$max_epochs)
                result <- do.call(RunSpatialDeconvolution, args)
                record <- sclet_ai_mainline_record(result, "spatial", name,
                    "run_spatial_deconvolution")
                columns <- record$artifacts$columns %||% character()
                sclet_ai_record_mainline_evidence(
                    result,
                    list(
                        id = paste0("ev:spatial_deconv_", sclet_ai_evidence_id_token(name)),
                        kind = "deterministic_summary",
                        claim_level = "measured",
                        values = list(
                            n_common_genes = sclet_ai_mainline_number(
                                record$inputs$n_common_genes
                            ),
                            n_columns = length(columns),
                            n_cells = ncol(result)
                        )
                    ),
                    source = name
                )
            },
            input_schema = list(
                sce_ref = list(type = "object", required = TRUE),
                ref_group_key = list(type = "character", required = TRUE),
                ref_batch_key = "character",
                max_epochs = "integer",
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                if (!base::requireNamespace("basilisk", quietly = TRUE)) {
                    return("optional_package_missing: basilisk is required for run_spatial_deconvolution (the Python environment)")
                }
                if (!inherits(params$sce_ref, "SingleCellExperiment")) {
                    return("reference_missing: run_spatial_deconvolution needs `sce_ref`, a scRNA-seq reference SingleCellExperiment")
                }
                for (side in c("spatial", "reference")) {
                    target <- if (identical(side, "spatial")) object else params$sce_ref
                    if (!"counts" %in% SummarizedExperiment::assayNames(target)) {
                        return(paste0(
                            "counts_assay_missing: run_spatial_deconvolution needs a 'counts' assay in the ",
                            side, " object"
                        ))
                    }
                }
                key <- params$ref_group_key
                if (is.null(key) || !is.character(key) || length(key) != 1L ||
                    !key %in% colnames(SummarizedExperiment::colData(params$sce_ref))) {
                    return(paste0(
                        "reference_group_missing: `ref_group_key` must name a colData column of the reference carrying cell types; available: ",
                        paste(colnames(SummarizedExperiment::colData(params$sce_ref)),
                            collapse = ", ")
                    ))
                }
                if (!length(intersect(rownames(object), rownames(params$sce_ref)))) {
                    return("no_common_genes: the spatial and reference objects share no gene identifiers")
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(spatial = "any", ai_evidence = "any"),
                note = paste(
                    "Per-spot abundances land in colData. The evidence node carries",
                    "counts of common genes and columns, never an abundance value."
                )
            ),
            mutates_object = TRUE,
            allowed_state_types = c("spatial", "ai_evidence"),
            estimated_cost = "high",
            idempotent = FALSE,
            planner_visible = FALSE
        )
    ))
    }

    actions
}
