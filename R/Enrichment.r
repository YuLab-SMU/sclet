#' Run Functional Enrichment Analysis
#' 
#' Wrapper for clusterProfiler enrichment analysis on marker genes.
#' 
#' @title RunEnrichment
#' @param object A SingleCellExperiment object
#' @param id ID of the detest state to use for enrichment. If NULL, uses the active detest state.
#' @param p_val_cutoff p-value threshold for marker selection (default 0.05)
#' @param logfc_cutoff logFC threshold for marker selection (default 0.25)
#' @param db Database to use ("GO" or "KEGG")
#' @param orgDb OrgDb object (e.g., org.Hs.eg.db)
#' @param keyType Key type of gene IDs (e.g., "SYMBOL", "ENTREZID")
#' @param ont Ontology for GO ("BP", "MF", "CC", "ALL")
#' @param name Analysis record id. Defaults to "enrichment".
#' @param ... Additional arguments passed to enrichGO or enrichKEGG
#' @return A SingleCellExperiment object with enrichment state added
#' @export
RunEnrichment <- function(object, id = NULL, p_val_cutoff = 0.05, logfc_cutoff = 0.25, 
                          db = c("GO", "KEGG"), orgDb = "org.Hs.eg.db", keyType = "SYMBOL", 
                          ont = "BP", name = "enrichment", ...) {
    db <- match.arg(db)
    
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is needed for this function to work.")
    }
    
    detest_state <- get_detest(object, id)
    if (is.null(detest_state)) {
        stop("No detest state found. Please run RunDEtest() first.")
    }
    
    markers <- detest_state$artifacts$result
    
    if (!"pval" %in% colnames(markers)) {
        if ("p_val" %in% colnames(markers)) colnames(markers)[colnames(markers) == "p_val"] <- "pval"
    }
    
    fc_col <- grep("log.*FC", colnames(markers), value = TRUE, ignore.case = TRUE)
    if (length(fc_col) == 0) {
        stop("Could not find logFC column in detest results.")
    }
    fc_col <- fc_col[1]
    
    # Filter for significant up-regulated markers for enrichment
    sig_markers <- markers[markers$pval < p_val_cutoff & markers[[fc_col]] > logfc_cutoff, ]
    
    if (nrow(sig_markers) == 0) {
        warning("No significant up-regulated markers found with given cutoffs. Returning original object.")
        return(object)
    }
    
    if (!"gene" %in% colnames(sig_markers)) {
        sig_markers$gene <- rownames(sig_markers)
    }
    if (!"cluster" %in% colnames(sig_markers)) {
        sig_markers$cluster <- "Cluster1"
    }
    
    if (db == "GO") {
        if (is.character(orgDb)) {
            if (!requireNamespace(orgDb, quietly = TRUE)) {
                stop(paste0("Package '", orgDb, "' is needed. Please install it."))
            }
        }
        
        res <- clusterProfiler::compareCluster(
            gene ~ cluster, 
            data = sig_markers, 
            fun = "enrichGO", 
            OrgDb = orgDb,
            keyType = keyType,
            ont = ont,
            ...
        )
    } else if (db == "KEGG") {
        if (keyType != "ENTREZID") {
            message("Note: KEGG analysis usually requires ENTREZID. If your genes are SYMBOLS, this might fail or return empty.")
        }
        
        res <- clusterProfiler::compareCluster(
            gene ~ cluster,
            data = sig_markers,
            fun = "enrichKEGG",
            ...
        )
    }
    
    state_inputs <- list(
        detest_id = detest_state$id,
        p_val_cutoff = p_val_cutoff,
        logfc_cutoff = logfc_cutoff,
        db = db,
        orgDb = orgDb
    )
    
    object <- sclet_set_analysis_state(
        object = object,
        type = "enrichment",
        id = name,
        method = "clusterProfiler::compareCluster",
        inputs = state_inputs,
        artifacts = list(result = res),
        summary = list(n_enriched = nrow(as.data.frame(res))),
        active = TRUE
    )
    
    object <- sclet_log_command(
        object,
        "RunEnrichment",
        params = list(db = db, orgDb = orgDb, name = name),
        outputs = list(enrichment = name)
    )
    
    return(object)
}

#' Get Enrichment Analysis Results
#'
#' @title get_enrichment
#' @param object A SingleCellExperiment object
#' @param id Analysis record id. If NULL, uses the active enrichment state.
#' @return The enrichment state record
#' @export
get_enrichment <- function(object, id = NULL) {
    if (is.null(id)) id <- sclet_get_active_state(object, "enrichment")
    if (is.null(id)) return(NULL)
    sclet_get_state_record(object, "enrichment", id)
}

#' Check if Enrichment Analysis Results exist
#'
#' @title has_enrichment
#' @param object A SingleCellExperiment object
#' @param id Analysis record id. If NULL, uses the active enrichment state.
#' @return Logical indicating if the enrichment state exists
#' @export
has_enrichment <- function(object, id = NULL) {
    !is.null(get_enrichment(object, id))
}
