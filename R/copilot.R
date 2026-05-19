#' Summarize SingleCellExperiment Context for LLM
#'
#' This function extracts the full analysis-state provenance (bloodline) from a 
#' SingleCellExperiment object and formats it into a highly structured JSON-like 
#' string. This output is specifically designed to be consumed as the `system_prompt`
#' or `context` for an LLM via the `aisdk` package.
#'
#' @param sce A SingleCellExperiment object.
#' @return A character string containing the structured context.
#' @export
SummarizeContextForLLM <- function(sce) {
    
    # 1. Get basic metadata
    n_cells <- ncol(sce)
    n_genes <- nrow(sce)
    assays_avail <- SummarizedExperiment::assayNames(sce)
    reductions_avail <- SingleCellExperiment::reducedDimNames(sce)
    
    # 2. Get full analysis state via get_analysis_context
    ctx <- get_analysis_context(sce)
    
    # 3. Format active states
    active_view <- paste0(
        "Active Assay: ", ctx$active$assay, "\n",
        "Active Layer: ", ctx$active$layer, "\n",
        "Active Reduction: ", ctx$active$reduction, "\n",
        "Active Ident (Clustering): ", ctx$active$ident
    )
    
    # 4. Format provenance records (the DAG)
    provenance <- ""
    if (length(ctx$records) > 0) {
        prov_list <- lapply(names(ctx$records), function(rec_name) {
            rec <- ctx$records[[rec_name]]
            paste0(
                "  - ", toupper(rec_name), ": ", 
                "Method=[", rec$method, "], ",
                "Inputs=[", paste(names(rec$inputs), unlist(rec$inputs), sep="=", collapse="; "), "], ",
                "Summary=[", paste(names(rec$summary), unlist(rec$summary), sep="=", collapse="; "), "]"
            )
        })
        provenance <- paste(prov_list, collapse = "\n")
    } else {
        provenance <- "  - No advanced analysis records found."
    }
    
    # 5. Build the final prompt context
    llm_context <- sprintf(
        "You are 'sclet_copilot', an expert single-cell RNA-seq computational biologist and code reviewer.\n\n=== CURRENT DATASET STATUS ===\nCells: %d\nGenes: %d\nAvailable Assays: %s\nAvailable Reductions: %s\n\n=== ACTIVE VIEW ===\n%s\n\n=== ANALYSIS PROVENANCE (DAG) ===\n%s\n\nBased on the above strict computational provenance, answer the user's question, diagnose potential errors, or evaluate causality and artifacts in the analysis chain.",
        n_cells,
        n_genes,
        paste(assays_avail, collapse=", "),
        paste(reductions_avail, collapse=", "),
        active_view,
        provenance
    )

    return(llm_context)
}

#' sclet AI Copilot
#'
#' An intelligent agent built on `aisdk` that uses the `analysis-state` provenance 
#' of the SingleCellExperiment object to diagnose analysis chains, interpret biological 
#' results, and detect potential computational artifacts.
#'
#' @param sce A SingleCellExperiment object.
#' @param question Character string. The question or task for the AI Copilot.
#' @param model Optional. An `aisdk` model object. If NULL, the function will attempt 
#'   to initialize a default model using the environment variable OPENAI_MODEL.
#' @return A character string containing the AI's response.
#' @export
sclet_copilot <- function(sce, question, model = NULL) {
    if (!requireNamespace("aisdk", quietly = TRUE)) {
        stop("Please install 'aisdk' to use sclet_copilot: remotes::install_github('YuLab-SMU/aisdk')")
    }
    
    # 1. Extract context
    context <- SummarizeContextForLLM(sce)
    
    # 2. Setup AI Model
    if (is.null(model)) {
        # Prioritize getting the default model from aisdk
        try({
            if (exists("get_model", envir = asNamespace("aisdk"))) {
                model <- aisdk::get_model()
            }
        }, silent = TRUE)
        
        # Fallback if no default model is found
        if (is.null(model)) {
            model_name <- Sys.getenv("OPENAI_MODEL", unset = NA)
            if (is.na(model_name)) {
                stop("No 'model' provided, no default model set in aisdk, and 'OPENAI_MODEL' is not set.")
            }
            provider <- aisdk::create_openai()
            model <- provider$language_model(model_name)
        }
    }
    
    # 3. Create Copilot Agent
    # We use create_agent from aisdk to set up the persona and inject the provenance context
    copilot_agent <- aisdk::create_agent(
        name = "sclet_copilot",
        description = "A Single-Cell Analysis Diagnostician and Biological Interpreter",
        system_prompt = context
    )
    
    # 4. Run Agent
    message("Consulting sclet_copilot via aisdk...")
    response <- copilot_agent$run(task = question, model = model)
    
    return(response$text)
}

#' Audit Analysis Chain
#'
#' Evaluates whether a set of candidate features (e.g., DE genes) is robust across 
#' the analysis chain (e.g., before and after batch correction/integration) and queries 
#' `sclet_copilot` to assign a State-Dependency Confidence Score.
#'
#' @param sce A SingleCellExperiment object.
#' @param features A character vector of features to audit.
#' @param raw_layer String. The name of the uncorrected layer. Defaults to "counts".
#' @param integrated_layer String. The name of the integrated/corrected layer. If NULL, automatically inferred from integration state.
#' @param ... Additional arguments passed to `sclet_copilot`.
#' @return A character string containing the audit report from the AI.
#' @export
AuditAnalysisChain <- function(sce, features, raw_layer = "counts", integrated_layer = NULL, ...) {
    
    # Check if features exist
    features <- intersect(features, rownames(sce))
    if (length(features) == 0) stop("None of the specified features are found in the SCE object.")
    
    # Identify integrated layer if not provided
    if (is.null(integrated_layer)) {
        integ_state <- get_integration(sce)
        if (is.null(integ_state)) {
            stop("No integration state found and 'integrated_layer' not specified. Cannot audit batch effects.")
        }
        # Look for the corrected layer from integration state
        if (!is.null(integ_state$artifacts$layer)) {
            integrated_layer <- integ_state$artifacts$layer
        } else {
            stop("The active integration method did not output a corrected layer (e.g., it output a reduction like Harmony). Expression-level audit requires a corrected layer like fastMNN's output.")
        }
    }
    
    if (!raw_layer %in% SummarizedExperiment::assayNames(sce)) stop(sprintf("Raw layer '%s' not found.", raw_layer))
    if (!integrated_layer %in% SummarizedExperiment::assayNames(sce)) stop(sprintf("Integrated layer '%s' not found.", integrated_layer))
    
    # Simple statistical heuristic: mean fold-change variance across layers (Optional computation to feed AI)
    # Here we just feed the fact that we want the AI to reason about these features based on the provenance.
    
    prompt <- sprintf(
        "I have a set of candidate features (e.g., DE genes): %s. \n\nI want to audit if these features are true biological signals or artifacts generated by the integration step (Raw layer: '%s', Integrated layer: '%s'). \n\nPlease trace the provenance DAG provided in your system context. Evaluate the assumptions of the integration method used, and calculate a conceptual 'State-Dependency Confidence Score' (Low/Medium/High) for these markers. Explain your reasoning based on Fisher's design foresight and Box's model humility.",
        paste(features, collapse = ", "),
        raw_layer,
        integrated_layer
    )
    
    message("Running AuditAnalysisChain via sclet_copilot...")
    report <- sclet_copilot(sce, question = prompt, ...)
    return(report)
}
