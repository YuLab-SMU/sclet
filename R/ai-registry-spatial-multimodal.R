# sclet_ai_registry_spatial_multimodal(): built-in AI actions for spatial and multimodal analysis.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# Spatial precedes multimodal, as before.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_spatial_multimodal <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("spatial" %in% include) {
    actions <- c(actions, list(
        check_spatial_readiness = AIAction(
            name = "check_spatial_readiness",
            description = paste(
                "Read-only readiness report for the spatial analyses recorded on this object.",
                "It reports which spatial state records exist and what their stored model estimates are,",
                "and it never runs a spatial algorithm, never ranks or recommends a spatial method, and",
                "never writes an evidence node or mutates the object."
            ),
            handler = function(object, params) {
                check_spatial_readiness(object, analysis = params$analysis)
            },
            input_schema = list(
                analysis = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                records <- sclet_get_state_records(object, "spatial")
                if (!length(records)) {
                    return(paste0(
                        "no_spatial_record: no spatial state record is available on this object; ",
                        "no spatial analysis has been recorded yet."
                    ))
                }
                analysis <- params$analysis
                if (!is.null(analysis) && !analysis %in% names(records)) {
                    return(paste0(
                        "analysis_not_available: spatial analysis '", analysis,
                        "' is not a recorded spatial analysis. Available spatial records: ",
                        paste(names(records), collapse = ", ")
                    ))
                }
                TRUE
            },
            returns = "value",
            requires_confirmation = FALSE,
            output_schema = list(
                type = "object",
                required_states = list(spatial = "any"),
                note = paste(
                    "Read-only: reports spatial model estimates and record counts only. It does not rank",
                    "or recommend a spatial method, and no per-cell or per-spot value is included. No",
                    "evidence node is written and the object is not mutated."
                )
            ),
            mutates_object = FALSE,
            allowed_state_types = c("spatial"),
            estimated_cost = "low",
            idempotent = TRUE
        ),
        summarize_spatial_evidence = AIAction(
            name = "summarize_spatial_evidence",
            description = paste(
                "Read-only bounded summary of the stored evidence for one recorded spatial analysis.",
                "It reports aggregate model estimates and record counts only; it never ranks or recommends",
                "a spatial method, and it never writes an evidence node or mutates the object."
            ),
            handler = function(object, params) {
                summarize_spatial_evidence(object, analysis = params$analysis)
            },
            input_schema = list(
                analysis = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                records <- sclet_get_state_records(object, "spatial")
                if (!length(records)) {
                    return(paste0(
                        "no_spatial_record: no spatial state record is available on this object; ",
                        "no spatial analysis has been recorded yet."
                    ))
                }
                analysis <- params$analysis
                if (!is.null(analysis) && !analysis %in% names(records)) {
                    return(paste0(
                        "analysis_not_available: spatial analysis '", analysis,
                        "' is not a recorded spatial analysis. Available spatial records: ",
                        paste(names(records), collapse = ", ")
                    ))
                }
                TRUE
            },
            returns = "value",
            requires_confirmation = FALSE,
            output_schema = list(
                type = "object",
                required_states = list(spatial = "any"),
                note = paste(
                    "Read-only: reports spatial model estimates and record counts only. It does not rank",
                    "or recommend a spatial method, and no per-cell or per-spot value is included. No",
                    "evidence node is written and the object is not mutated."
                )
            ),
            mutates_object = FALSE,
            allowed_state_types = c("spatial"),
            estimated_cost = "low",
            idempotent = TRUE
        ),
        compare_spatial_analyses = AIAction(
            name = "compare_spatial_analyses",
            description = paste(
                "Read-only comparison of the spatial analyses recorded on this object, side by side.",
                "It reports the stored estimates of each record so a difference can be seen; it does not",
                "rank or recommend a spatial method, never nominates a \"best\" analysis, and never writes",
                "an evidence node or mutates the object."
            ),
            handler = function(object, params) {
                compare_spatial_analyses(object)
            },
            input_schema = list(),
            prerequisites = function(object, params, planned = NULL) {
                records <- sclet_get_state_records(object, "spatial")
                if (!length(records)) {
                    return(paste0(
                        "no_spatial_record: no spatial state record is available on this object; ",
                        "no spatial analysis has been recorded yet."
                    ))
                }
                TRUE
            },
            returns = "value",
            requires_confirmation = FALSE,
            output_schema = list(
                type = "object",
                required_states = list(spatial = "any"),
                note = paste(
                    "Read-only: reports spatial model estimates and record counts only. It does not rank",
                    "or recommend a spatial method, and no per-cell or per-spot value is included. No",
                    "evidence node is written and the object is not mutated."
                )
            ),
            mutates_object = FALSE,
            allowed_state_types = c("spatial"),
            estimated_cost = "low",
            idempotent = TRUE
        )
    ))
    }

    if ("multimodal" %in% include) {
    actions <- c(actions, list(
        check_multimodal_readiness = AIAction(
            name = "check_multimodal_readiness",
            description = paste(
                "Read-only readiness report for the multimodal (assay-platform) records on this object.",
                "It reports which multimodal state records exist and what their stored estimates are, and",
                "it never registers a modality, never runs a multimodal workflow, never ranks or recommends",
                "a modality or method, and never writes an evidence node or mutates the object."
            ),
            handler = function(object, params) {
                check_multimodal_readiness(object, modalities = params$modalities)
            },
            input_schema = list(
                modalities = "character_vector"
            ),
            prerequisites = function(object, params, planned = NULL) {
                records <- sclet_get_state_records(object, "multimodal")
                if (!length(records)) {
                    return(paste0(
                        "no_multimodal_record: no multimodal state record is available on this object; ",
                        "no multimodal analysis has been recorded yet."
                    ))
                }
                TRUE
            },
            returns = "value",
            requires_confirmation = FALSE,
            output_schema = list(
                type = "object",
                required_states = list(multimodal = "any"),
                note = paste(
                    "Read-only: reports modality-platform labels, multimodal model estimates and record",
                    "counts only. It does not rank or recommend a modality or method, and no per-cell value,",
                    "barcode or feature name is included. No evidence node is written and the object is not",
                    "mutated."
                )
            ),
            mutates_object = FALSE,
            allowed_state_types = c("multimodal"),
            estimated_cost = "low",
            idempotent = TRUE
        ),
        summarize_multimodal_evidence = AIAction(
            name = "summarize_multimodal_evidence",
            description = paste(
                "Read-only bounded summary of the stored evidence for one recorded multimodal analysis.",
                "It reports aggregate model estimates and record counts only; it never ranks or recommends",
                "a modality or method, and it never writes an evidence node or mutates the object."
            ),
            handler = function(object, params) {
                summarize_multimodal_evidence(object, name = params$name)
            },
            input_schema = list(
                name = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                records <- sclet_get_state_records(object, "multimodal")
                if (!length(records)) {
                    return(paste0(
                        "no_multimodal_record: no multimodal state record is available on this object; ",
                        "no multimodal analysis has been recorded yet."
                    ))
                }
                name <- params$name
                if (!is.null(name) && !name %in% names(records)) {
                    return(paste0(
                        "analysis_not_available: multimodal analysis '", name,
                        "' is not a recorded multimodal analysis. Available multimodal records: ",
                        paste(names(records), collapse = ", ")
                    ))
                }
                TRUE
            },
            returns = "value",
            requires_confirmation = FALSE,
            output_schema = list(
                type = "object",
                required_states = list(multimodal = "any"),
                note = paste(
                    "Read-only: reports modality-platform labels, multimodal model estimates and record",
                    "counts only. It does not rank or recommend a modality or method, and no per-cell value,",
                    "barcode or feature name is included. No evidence node is written and the object is not",
                    "mutated."
                )
            ),
            mutates_object = FALSE,
            allowed_state_types = c("multimodal"),
            estimated_cost = "low",
            idempotent = TRUE
        ),
        compare_multimodal_evidence = AIAction(
            name = "compare_multimodal_evidence",
            description = paste(
                "Read-only comparison of the multimodal analyses recorded on this object, side by side.",
                "It reports the stored estimates of each record so a difference can be seen; it does not",
                "rank or recommend a modality or method, never nominates a \"best\" analysis, and never",
                "writes an evidence node or mutates the object."
            ),
            handler = function(object, params) {
                compare_multimodal_evidence(object)
            },
            input_schema = list(),
            prerequisites = function(object, params, planned = NULL) {
                records <- sclet_get_state_records(object, "multimodal")
                if (!length(records)) {
                    return(paste0(
                        "no_multimodal_record: no multimodal state record is available on this object; ",
                        "no multimodal analysis has been recorded yet."
                    ))
                }
                TRUE
            },
            returns = "value",
            requires_confirmation = FALSE,
            output_schema = list(
                type = "object",
                required_states = list(multimodal = "any"),
                note = paste(
                    "Read-only: reports modality-platform labels, multimodal model estimates and record",
                    "counts only. It does not rank or recommend a modality or method, and no per-cell value,",
                    "barcode or feature name is included. No evidence node is written and the object is not",
                    "mutated."
                )
            ),
            mutates_object = FALSE,
            allowed_state_types = c("multimodal"),
            estimated_cost = "low",
            idempotent = TRUE
        )
    ))
    }

    actions
}
