#' Find Marker Pathways
#'
#' This function identifies cell-type specific pathways by performing differential testing
#' on pathway scores stored as an alternative experiment (altExp).
#'
#' @param sce A SingleCellExperiment object.
#' @param altExp_name Character. The name of the altExp containing pathway scores. Defaults to "PathwayScores".
#' @param group.by Character. The column in colData to use for grouping. If NULL, uses ActiveIdent.
#' @param assay Character. The assay name inside the altExp to use. Defaults to "scores".
#' @param ... Additional arguments passed to `FindAllMarkers()`.
#'
#' @return A data frame of significant marker pathways.
#' @export
FindMarkerPathways <- function(sce, altExp_name = "PathwayScores", group.by = NULL, assay = "scores", ...) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("SingleCellExperiment is required.")
    }
    
    if (!altExp_name %in% SingleCellExperiment::altExpNames(sce)) {
        stop(sprintf("altExp '%s' not found. Please run RunGeneSetScoring with as_altExp=TRUE first.", altExp_name))
    }
    
    # Extract altExp
    alt_sce <- SingleCellExperiment::altExp(sce, altExp_name)
    
    # Transfer groups
    if (is.null(group.by)) {
        idents <- Idents(sce)
        if (is.null(idents)) stop("No active identity found. Please specify group.by.")
        SingleCellExperiment::colLabels(alt_sce) <- idents
    } else {
        if (!group.by %in% colnames(SummarizedExperiment::colData(sce))) {
            stop(sprintf("Column '%s' not found in colData.", group.by))
        }
        SingleCellExperiment::colLabels(alt_sce) <- factor(SummarizedExperiment::colData(sce)[[group.by]])
    }
    
    message(sprintf("Running differential testing on pathway scores in altExp '%s'...", altExp_name))
    
    # Run FindAllMarkers on altExp
    res <- FindAllMarkers(alt_sce, assay = assay, ...)
    
    # Clean up column names for pathway context
    if ("gene" %in% colnames(res)) {
        colnames(res)[colnames(res) == "gene"] <- "pathway"
    }
    
    return(res)
}

#' Run KEGG Pathway Scoring and Marker Detection
#'
#' A convenient wrapper to download KEGG pathways for a specific species, score cells,
#' and optionally identify cell-type specific marker pathways.
#'
#' @param sce A SingleCellExperiment object.
#' @param species Character. Species code for KEGG (e.g., "hsa" for human, "mmu" for mouse).
#' @param method Character. Scoring method (e.g., "UCell", "AUCell", "GSVA"). Defaults to "UCell".
#' @param find_markers Logical. Whether to run differential pathway analysis after scoring. Defaults to TRUE.
#' @param group.by Character. Grouping column for marker detection. If NULL, uses ActiveIdent.
#' @param name Character. Name of the altExp to store scores. Defaults to "KEGG".
#' @param ... Additional arguments passed to `RunGeneSetScoring()`.
#'
#' @return A list containing the updated `SingleCellExperiment` and the marker `data.frame` (if find_markers=TRUE).
#' @export
RunKEGG <- function(sce, species = "hsa", method = "UCell", find_markers = TRUE, group.by = NULL, name = "KEGG", ...) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Please install 'clusterProfiler' to download KEGG pathways.")
    }
    
    message(sprintf("Downloading KEGG pathways for species '%s'...", species))
    # We use internal function from clusterProfiler to get KEGG data, same as sc2p
    kk <- clusterProfiler::download_KEGG(species)
    gs <- kk$KEGGPATHID2EXTID
    
    # Map ENTREZ to SYMBOL since most scRNA-seq uses SYMBOL
    message("Mapping ENTREZID to SYMBOL...")
    
    # Determine the OrgDb based on species
    orgdb_name <- switch(species,
        "hsa" = "org.Hs.eg.db",
        "mmu" = "org.Mm.eg.db",
        "rno" = "org.Rn.eg.db",
        NULL
    )
    
    if (is.null(orgdb_name)) {
        stop("Currently only 'hsa', 'mmu', and 'rno' are automatically mapped. Please map genes manually.")
    }
    
    if (!requireNamespace(orgdb_name, quietly = TRUE)) {
        stop(sprintf("Please install '%s' via BiocManager.", orgdb_name))
    }
    
    gene <- clusterProfiler::bitr(gs$to, 'ENTREZID', 'SYMBOL', orgdb_name)
    gs2 <- merge(gs, gene, by.x = 'to', by.y = 'ENTREZID')
    gset <- split(gs2$SYMBOL, gs2$from)
    
    message(sprintf("Found %d KEGG pathways.", length(gset)))
    
    # Score cells
    sce <- RunGeneSetScoring(sce, gene_sets = gset, method = method, name = name, as_altExp = TRUE, ...)
    
    res <- list(sce = sce)
    
    # Differential analysis
    if (find_markers) {
        markers <- FindMarkerPathways(sce, altExp_name = name, group.by = group.by)
        
        # Merge with pathway names
        term2name <- kk$KEGGPATHID2NAME
        colnames(term2name) <- c("pathway", "Description")
        
        markers <- merge(markers, term2name, by = "pathway", all.x = TRUE)
        # Reorder columns
        if ("cluster" %in% colnames(markers)) {
            markers <- markers[order(markers$cluster, markers$pval), ]
        }
        res$markers <- markers
    }
    
    return(res)
}
