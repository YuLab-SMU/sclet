#' Run Functional Enrichment Analysis
#' 
#' Wrapper for clusterProfiler enrichment analysis on marker genes.
#' 
#' @title RunEnrichment
#' @param markers Data frame returned by FindAllMarkers or FindMarkers
#' @param db Database to use ("GO" or "KEGG")
#' @param orgDb OrgDb object (e.g., org.Hs.eg.db)
#' @param keyType Key type of gene IDs (e.g., "SYMBOL", "ENTREZID")
#' @param ont Ontology for GO ("BP", "MF", "CC", "ALL")
#' @param ... Additional arguments passed to enrichGO or enrichKEGG
#' @return A clusterProfiler result object (compareClusterResult)
#' @export
RunEnrichment <- function(markers, db = c("GO", "KEGG"), orgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", ...) {
    db <- match.arg(db)
    
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is needed for this function to work.")
    }
    
    # Check if markers is a data frame with 'gene' and 'cluster' columns
    if (!all(c("gene", "cluster") %in% colnames(markers))) {
        # Try to infer if it's from FindMarkers (single cluster)
        if ("gene" %in% colnames(markers) && !"cluster" %in% colnames(markers)) {
             markers$cluster <- "Cluster1"
        } else if (!"gene" %in% colnames(markers)) {
            # Maybe rownames are genes
            markers$gene <- rownames(markers)
            # Check if we still don't have gene column (e.g. empty dataframe or rownames numeric)
            if (is.null(markers$gene) || length(markers$gene) == 0) {
                 stop("Input markers must have 'gene' and 'cluster' columns.")
            }
            if (!"cluster" %in% colnames(markers)) markers$cluster <- "Cluster1"
        } else {
             stop("Input markers must have 'gene' and 'cluster' columns.")
        }
    }
    
    # Filter for significant markers? Usually user does this before.
    
    if (db == "GO") {
        # Check OrgDb
        if (is.character(orgDb)) {
            if (!requireNamespace(orgDb, quietly = TRUE)) {
                stop(paste0("Package '", orgDb, "' is needed. Please install it."))
            }
            # Note: clusterProfiler::enrichGO accepts character string for OrgDb.
            # We do NOT need to load it with get() or library().
        }
        
        # Use compareCluster for multiple clusters
        res <- clusterProfiler::compareCluster(
            gene ~ cluster, 
            data = markers, 
            fun = "enrichGO", 
            OrgDb = orgDb,
            keyType = keyType,
            ont = ont,
            ...
        )
    } else if (db == "KEGG") {
        # For KEGG, we usually need Entrez IDs.
        # If keyType is SYMBOL, we might need to convert.
        # But clusterProfiler::enrichKEGG takes gene input.
        
        # If input is SYMBOL, warn or convert?
        # clusterProfiler::bitr can be used.
        
        if (keyType != "ENTREZID") {
            # Attempt conversion
            if (is.character(orgDb)) {
                 if (!requireNamespace(orgDb, quietly = TRUE)) stop(paste0("Package '", orgDb, "' needed for ID conversion."))
                 odb <- get(orgDb) # Load package? 
                 # Safer to use clusterProfiler::bitr
            }
            
            # This is complex to automate robustly without loading the OrgDb.
            # Let's assume user provides ENTREZID for KEGG or we rely on them to do it.
            # Or we try to convert if orgDb is standard.
            
            message("Note: KEGG analysis usually requires ENTREZID. If your genes are SYMBOLS, this might fail or return empty.")
        }
        
        res <- clusterProfiler::compareCluster(
            gene ~ cluster,
            data = markers,
            fun = "enrichKEGG",
            ...
        )
    }
    
    return(res)
}
