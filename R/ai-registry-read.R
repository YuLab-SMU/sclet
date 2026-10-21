# sclet_ai_registry_read(): built-in AI actions for the read-only status, ledger and QC views.
#
# Split out of AIDefaultExecutionRegistry() as a pure move: the action
# definitions below are unchanged from the original registry body. Nothing
# was reworded, reordered or re-parameterised.
# The read group never mutates the object; the registry default is "read".
#
# Internal: called by AIDefaultExecutionRegistry() in ai-execution.R, never
# exported. `include` is already validated and normalised by the caller, and
# each block keeps its original `%in% include` guard, so group selection
# (including the composite confirm_design predicate) behaves as before.
sclet_ai_registry_read <- function(object, include) {
    actions <- list()
    # Blocks are mutually exclusive by construction, so appending to an
    # initially empty list reproduces the old `actions <- c(actions, ...)`
    # and keeps registry order identical.
    if ("read" %in% include) {
    actions <- c(actions, list(
        inspect_status = AIAction(
            name = "inspect_status",
            description = "Inspect deterministic sclet status without modifying the object.",
            handler = function(object, params) Status(object),
            returns = "value",
            requires_confirmation = FALSE,
            mutates_object = FALSE,
            output_schema = list(type = "object"),
            estimated_cost = "low",
            idempotent = TRUE
        ),
        inspect_ledger = AIAction(
            name = "inspect_ledger",
            description = "Inspect the bounded AI-facing analysis ledger.",
            handler = function(object, params) {
                detail <- params$detail %||% "summary"
                target <- params$target %||% NULL
                GetAnalysisLedger(object, detail = detail, target = target)
            },
            input_schema = list(
                detail = c("summary", "full"),
                target = "character"
            ),
            returns = "value",
            requires_confirmation = FALSE,
            mutates_object = FALSE,
            output_schema = list(type = "object"),
            estimated_cost = "low",
            idempotent = TRUE
        ),
        check_qc = AIAction(
            name = "check_qc",
            description = "Inspect deterministic quality checks without modifying the object.",
            handler = function(object, params) GetAnalysisLedger(object)$quality_checks,
            returns = "value",
            requires_confirmation = FALSE,
            mutates_object = FALSE,
            output_schema = list(type = "object"),
            estimated_cost = "low",
            idempotent = TRUE
        )
    ))
    }

    actions
}
