#' Print a textual summary of the analysis pipeline
#' 
#' @title PipelineSummary
#' @param object a SingleCellExperiment object
#' @return a textual summary of the pipeline
#' @export
PipelineSummary <- function(object) {
    commands <- CommandLog(object, details = TRUE)
    if (nrow(commands) == 0) {
        cat("No sclet analysis pipeline recorded.\n")
        return(invisible(NULL))
    }
    
    cat("========================================\n")
    cat("sclet Analysis Pipeline Summary\n")
    cat("========================================\n\n")
    
    cat("1. Preprocessing & Feature Selection:\n")
    prep_cmds <- c("NormalizeData", "ScaleData", "FindVariableFeatures", "BatchRemover")
    prep <- commands[commands$command %in% prep_cmds, ]
    if (nrow(prep) == 0) {
        cat("  - Not recorded.\n")
    } else {
        for (i in seq_len(nrow(prep))) {
            cat(sprintf("  - [%s] %s (%s)\n", prep$timestamp[i], prep$command[i], prep$params_summary[i]))
        }
    }
    cat("\n")
    
    cat("2. Dimensional Reduction & Clustering:\n")
    dim_cmds <- c("RunPCA", "RunUMAP", "FindNeighbors", "FindClusters")
    dim_steps <- commands[commands$command %in% dim_cmds, ]
    if (nrow(dim_steps) == 0) {
        cat("  - Not recorded.\n")
    } else {
        for (i in seq_len(nrow(dim_steps))) {
            cat(sprintf("  - [%s] %s (%s)\n", dim_steps$timestamp[i], dim_steps$command[i], dim_steps$params_summary[i]))
        }
    }
    cat("\n")
    
    cat("3. Downstream Analysis:\n")
    downstream <- commands[!commands$command %in% c(prep_cmds, dim_cmds, "sce_merge"), ]
    if (nrow(downstream) == 0) {
        cat("  - Not recorded.\n")
    } else {
        for (i in seq_len(nrow(downstream))) {
            cat(sprintf("  - [%s] %s (%s)\n", downstream$timestamp[i], downstream$command[i], downstream$params_summary[i]))
        }
    }
    
    cat("\n========================================\n")
    invisible(commands)
}
