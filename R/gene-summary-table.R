
#' @title gene_summary_table
#' This function append gene summary information to the input gene data frame and
#' open the result in various formats.
#' @param gene a data frame of genes (e.g., marker genes, DE genes)
#' @param gene_col column of the 'gene' that contains gene IDs
#' @param keyType ID type of the gene IDs
#' @param cluster_col column to cluster genes (optional)
#' @param output one of 'html', 'pdf', 'png', 'md', 'csv' or 'data.frame'
#' @param browse whether open the output file
#' @return output file path
#' @importFrom fanyi gene_summary
#' @importFrom tinytable tt
#' @importFrom tinytable group_tt
#' @importFrom tinytable style_tt
#' @importFrom tinytable save_tt
#' @importFrom yulab.utils show_in_excel
#' @importFrom yulab.utils o
#' @examples 
#' \dontrun{
#' n <- 3
#' topn <- pbmc.markers |>
#'     group_by(cluster) |>
#'     dplyr::filter(avg_log2FC > 1) |>
#'     slice_head(n = n) |>
#'     ungroup()
#' 
#' gene_summary_table(topn, gene_col = 'gene', keyType='SYMBOL')
#' 
#' gene_summary_table(topn, gene_col = 'gene', cluster_col = 'cluster', keyType='SYMBOL')
#' }
#' @export
gene_summary_table <- function(gene, gene_col = "gene", keyType = "SYMBOL", cluster_col, output = "html", browse = TRUE) {
    
    message(yulab.utils:::pkg_ref("fanyi"))

    output <- match.arg(output, c("html", "csv", "pdf", "png", "md", "data.frame"))

    gene[[gene_col]] <- sub("\\.\\d", "", gene[[gene_col]])

    if (keyType != "ENTREZID") {
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
             stop("Package 'org.Hs.eg.db' is needed for this function to work. Please install it.")
        }
        g <- clusterProfiler::bitr(
            gene[[gene_col]], 
            fromType = keyType, 
            toType = 'ENTREZID', 
            OrgDb = 'org.Hs.eg.db'
        )
        gid <- g$ENTREZID
        gene <- merge(gene, g, by.x = gene_col, by.y=keyType, all.x=TRUE)
        eid <- 'ENTREZID'
    } else {
        gid <- gene[[gene_col]]
        eid <- gene_col
    }
    
    gs <- fanyi::gene_summary(gid)  

    g2 <- merge(gene, gs, by.x=eid, by.y='uid', all.x=TRUE)

    if (output == "data.frame") {
        return(g2)
    }

    if (output == "csv") {
        return(show_in_excel(g2))
    }

    outfile <- tempfile(fileext = sprintf(".%s", output))

    width <- rep(1, ncol(g2))
    width[length(width)] <- 4

    res <- tt(g2, width=width) 
    
    if (!missing(cluster_col)) {
        cls <- g2[[cluster_col]]
        g2 <- g2[order(cls),]
        rownames(g2) <- NULL
    
        cls <- g2[[cluster_col]]
        idx <- which(!duplicated(cls)) |>
            setNames(cls[!duplicated(cls)])

        group_idx <- idx + cumsum(rep(1, length(idx))) -1 

        j <- which(names(g2) == cluster_col)
        res <- tt(g2[, -j], width=width[-1]) |>
            group_tt(i = as.list(idx)) |> 
            style_tt(i = group_idx, background='grey', color='white', fontsize=1.2)
    }

    save_tt(res, outfile, overwrite=TRUE)

    o(outfile)
    invisible(outfile)
}

