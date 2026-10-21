# sclet_ai_registry_integration(): built-in AI actions for batch integration and integration-route comparison.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# Batch keys are design semantics, so their levels stay available to planning.
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_integration <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("integration" %in% include) {
    actions <- c(actions, list(
        run_integration = AIAction(
            name = "run_integration",
            description = "Run batch-correction integration via fastMNN, Harmony, or scVI and register the result.",
            handler = function(object, params) {
                method <- match.arg(params$method %||% "fastMNN", c("fastMNN", "Harmony", "scVI"))
                batch <- params$batch
                features <- sclet_ai_gene_selector(object, params$features, arg = "features")
                # RunIntegration() chooses the dimensions internally and
                # has no `dims` argument. Reject a supplied value during
                # plan validation rather than silently accepting an inert
                # parameter that gives the user/model a false impression.
                if (!is.null(params$dims)) {
                    stop("`run_integration` does not support `dims`; omit it so RunIntegration() can choose internally.", call. = FALSE)
                }
                route_name <- params$name %||% tolower(method)
                started <- Sys.time()
                result <- RunIntegration(
                    object,
                    method = method,
                    batch = batch,
                    features = features,
                    layer = params$layer %||% NULL,
                    reduction = params$reduction %||% NULL,
                    name = route_name
                )
                # Record the route's OWN metrics, on the route's own
                # embedding. Without this the plan-driven path produced an
                # integration record with no numbers at all, so a route
                # comparison could only be made through the separate
                # RunIntegrationRoutes() entry point -- and that one
                # measured the raw PCA for every route (see
                # sclet_ai_route_reduction()).
                elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
                result <- sclet_ai_record_route_metrics(
                    result, route_name = route_name, batch = batch,
                    runtime_sec = elapsed)
                result
            },
            input_schema = list(
                batch = list(type = "character", required = TRUE),
                method = c("fastMNN", "Harmony", "scVI"),
                name = "character",
                # Intentionally no `dims`: RunIntegration chooses this
                # internally, so the plan must not advertise a no-op input.
                features = "character_vector",
                layer = "character",
                reduction = "character"
            ),
            prerequisites = function(object, params, planned = NULL) {
                method <- match.arg(params$method %||% "fastMNN", c("fastMNN", "Harmony", "scVI"))
                batch <- params$batch
                if (is.null(batch) || !is.character(batch) || length(batch) != 1L || !nzchar(batch)) {
                    return("missing required parameter: batch (a single colData column name)")
                }
                cd_names <- planned$colData_columns %||% colnames(SummarizedExperiment::colData(object))
                if (!batch %in% cd_names) {
                    return(paste0("batch column is not available in colData: ", batch))
                }
                confirmation <- sclet_ai_find_design_confirmation(object, batch = batch)
                # A `confirm_design` step earlier in the SAME plan will
                # write the confirmation at execution time; validation must
                # accept that promise (deferred projection), otherwise an
                # integration plan can never validate on an unconfirmed
                # object. Execution re-checks against the real state.
                planned_design <- isTRUE(planned$design_confirmed)
                if (is.null(confirmation) && !planned_design) {
                    return(paste0(
                        "design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, design = list(batch = '",
                        batch, "', ...)) before running integration"
                    ))
                }
                if (identical(method, "Harmony")) {
                    if (!requireNamespace("harmony", quietly = TRUE)) {
                        return("optional_package_missing: harmony is required for method Harmony; install.packages('harmony')")
                    }
                }
                if (identical(method, "scVI")) {
                    if (!requireNamespace("basilisk", quietly = TRUE)) {
                        return("optional_package_missing: basilisk is required for method scVI")
                    }
                    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
                        return("optional_package_missing: zellkonverter is required for method scVI")
                    }
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(integration = "any_qualified"),
                # Measured per method on this machine (2026-10-14, offline):
                # fastMNN registers "corrected", Harmony registers
                # "HARMONY". scVI (unused here) registers "scVI". fastMNN
                # ignores the `name` parameter and keeps "corrected".
                allowed_reductions = c("corrected", "HARMONY", "scVI")
            ),
            mutates_object = TRUE,
            allowed_state_types = c("integration", "reduction"),
            estimated_cost = "medium",
            idempotent = FALSE
        ),
        compare_integration_routes = AIAction(
            name = "compare_integration_routes",
            description = paste(
                "Compare the integration routes recorded on this object and",
                "register ONE bounded evidence node for the comparison.",
                "Reads the metrics each route recorded on its OWN embedding",
                "(batch mixing, biological preservation, runtime) and reports",
                "the trade-off: a route that removes batch structure at the cost",
                "of biological signal is a trade-off, not a win. Route names",
                "stay in the state records; the evidence carries anonymized",
                "route codes and numbers only."
            ),
            handler = function(object, params) {
                sclet_ai_record_integration_comparison(object,
                    name = params$name %||% "integration_comparison")
            },
            input_schema = list(name = "character"),
            prerequisites = function(object, params, planned = NULL) {
                records <- tryCatch(
                    sclet_get_state_records(object, "integration"),
                    error = function(e) list())
                measured <- Filter(function(rec) {
                    length(rec$summary$metrics %||% list()) > 0L
                }, records)
                if (!length(measured)) {
                    return(paste0(
                        "no_integration_route_recorded: run_integration (or register the raw ",
                        "baseline) first -- there is nothing to compare without a recorded route."
                    ))
                }
                TRUE
            },
            returns = "sce",
            requires_confirmation = TRUE,
            output_schema = list(
                required_states = list(ai_evidence = "ev:integration_comparison"),
                active_evidence = "ev:integration_comparison"
            ),
            mutates_object = TRUE,
            allowed_state_types = c("integration", "ai_evidence"),
            estimated_cost = "low",
            idempotent = TRUE
        )
    ))
    }

    actions
}
