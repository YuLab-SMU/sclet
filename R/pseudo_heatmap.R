#' Plot gene expression heatmap along pseudotime trajectory
#'
#' This function generates a heatmap of gene expression along a given pseudotime trajectory, 
#' with options for scaling, sorting, clustering decision, and lineage selection. It also handles genes with all 
#' `NaN` expression values by issuing a warning and removing them from the analysis.
#'
#' @param sce A `SingleCellExperiment` object containing gene expression data and pseudotime information.
#' @param features A character vector of gene names to be visualized in the heatmap.
#' @param assay_name A string specifying the assay to use for gene expression data. Default is `"logcounts"`.
#' @param pseudotime.data A string specifying the pseudotime data column in `colData(sce)`. Default is `"slingPseudotime"`.
#' @param lineage A string specifying the lineage (cell fate) for the pseudotime trajectory. Default is `NULL`, in which case the first column of the pseudotime data is used.
#' @param scale A logical indicating whether to scale the gene expression values. Default is `TRUE`.
#' @param sort A logical indicating whether to sort genes by gradient of expression. Default is `TRUE`.
#' @param cluster_columns A logical indicating whether to perform hierarchical clustering on cells' pseudotime.  Default is `FALSE`.
#' @param cluster_rows A logical indicating whether to perform hierarchical clustering on gene expressions.  Default is `FALSE`.
#'
#'
#' @return A `ComplexHeatmap` object representing the gene expression heatmap.
#'
#' @import ggplot2
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom stats na.omit
#' @importFrom stats ave
#' @importFrom methods is
#' @export
pseudo_heatmap <- function(
    sce, 
    features,
    assay_name = "logcounts",
    pseudotime.data = "slingPseudotime",
    lineage = NULL,
    scale = TRUE,
    sort = TRUE,
    cluster_columns = FALSE, 
    cluster_rows = FALSE) {  
  
 
  if (!is(sce, "SingleCellExperiment")){
    stop("The input must be a SingleCellExperiment object.")
  }
  
  if (!pseudotime.data %in% colnames(colData(sce))){
    stop("SCE lacks pseudotime data. Please run 'runSlingshot' before this function.")
  }
  
  if(! "logcounts" %in% assayNames(sce)){
    sce <- NormalizeData(sce)
  }
  
  if (cluster_rows && sort) {
    # stop("Both `cluster_rows = TRUE` and `sort = TRUE` are set. These options modify row (gene) order differently. Please choose one.")
    message("As 'cluster_rows = TRUE', 'sort' will set to 'FALSE' automatically.")
  }
  
  
  ## extract pseudotime data
  if (is.null(features) || length(features) == 0) {
    stop("Please provide at least one gene in the 'features' argument.")
  }
  
  if (length(lineage) > 1){
    stop("Please provide only one lineage at a time")
  }
  
  pseudo <- colData(sce)[[pseudotime.data]]
  pseudo <- as.data.frame(pseudo)
  fate_names <- colnames(pseudo)
  
  ## extract gene expression data
  gene_expression <- t(assay(sce, assay_name))
  gene_expression <- as.data.frame(gene_expression[, features])
  colnames(gene_expression) <- features
  
  ## prepare lineage input
  if (is.null(lineage)) {
    lineage <- fate_names[1]
  }
  
  if (!lineage %in% fate_names) {
    stop(paste("The provided lineage '", lineage, "' is not valid. Please choose one from the following cell fates:", paste(fate_names, collapse = ", "), "."))
  }
  
  pseudo_new <- pseudo[, lineage, drop = FALSE]
  df <- cbind(pseudo_new, gene_expression)
  
  ## melt df to the format which is prepared for plotting heatmap
  df_melt_gene <- reshape2::melt(
    df, 
    id.vars = colnames(pseudo_new),
    measure.vars = features,
    variable.name = "Gene",
    value.name = "Expression"
  )
  
  df_melt <- reshape2::melt(
    df_melt_gene, 
    id.vars = c("Gene", "Expression"),
    measure.vars = lineage,
    variable.name = "cell_fate",
    value.name = "Pseudotime"
  )
  
  df_melt$Gene <- factor(df_melt$Gene, levels = unique(df_melt$Gene))
  df_melt <- na.omit(df_melt)
  
  ## predict gene expression
  predictions <- gam_analysis_SCE(df_melt = df_melt, features = features, fate_names = lineage)
  
  if (scale){
    predictions$Predicted <- ave(predictions$Predicted, predictions$Gene, FUN = rescale)
  }
  
  ## prepare heatmap matrix
  heatmap_matrix <- reshape2::dcast(predictions, Gene ~ Pseudotime, value.var = "Predicted")
  
  heatmap_matrix <- as.matrix(heatmap_matrix[, -1])
  rownames(heatmap_matrix) <- levels(predictions$Gene)
  
  ## check if some all-NaN genes exist
  nan_genes <- apply(heatmap_matrix, 1, function(row) all(is.na(row)))
  if (any(nan_genes)) {
    warning("The following genes have all NaN expression values and will be removed: ",
            paste(rownames(heatmap_matrix)[nan_genes], collapse = ", "))
    heatmap_matrix <- heatmap_matrix[!nan_genes, , drop = FALSE]
  }
  
  ## choose whether to sort the matrix
  if (sort && nrow(heatmap_matrix) > 1){
    gradient_scores <- apply(heatmap_matrix, 1, compute_gradient)
    
    sort_df <- data.frame(
      row = rownames(heatmap_matrix),
      max_position = apply(heatmap_matrix, 1, which.max),
      gradient = gradient_scores
    ) 
    
    # Arrange by gradient first and then by max position second
    sort_df <- sort_df %>%
      arrange(.data$gradient) %>%
      arrange(.data$max_position)
    
    if (nrow(sort_df) > 0) {
      heatmap_matrix <- heatmap_matrix[sort_df$row, , drop = FALSE]
    }
  }
  
  
  ## gene expression heatmap
  title_fill <- ifelse(scale, "Relative\nExpression", "Expression")
  
  p <- ComplexHeatmap::Heatmap(heatmap_matrix, 
                               name = "Expression", 
                               col = viridis::viridis(100),  
                               show_column_names = FALSE, 
                               show_row_names = TRUE,
                               row_dend_width = unit(2, "mm"),
                               height = unit(6, "cm"),  
                               column_title = lineage,  
                               column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
                               row_title = "Gene",
                               row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                               cluster_columns = cluster_columns,
                               cluster_rows = cluster_rows,
                               heatmap_legend_param = list(title = title_fill))
  
  return(p)
}
