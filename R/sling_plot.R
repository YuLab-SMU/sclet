#' Visualize Single-Cell Trajectories on the basis of ggsc package
#' 
#' This function visualizes single-cell trajectories using UMAP (or other specified reduction method) and overlays the trajectory curves. 
#' The cells are colored by their respective lineage trajectories and the lineages are drawn as curves based on Slingshot's inferred trajectories.
#' 
#' @param sce A SingleCellExperiment object containing the results from Slingshot.
#' @param dims A vector specifying the dimensions of the reduction result to plot, default the top 2 dimensions are used.
#' @param group The name of the column in colData(sce) containing cluster labels.
#' @param reduction The name of the dimensionality reduction to plot. Default is `"UMAP"`.
#' @param lineages The lineages plotted in this function. Default all lineages are used.
#' @param cells A vector of cell names (IDs) to be included in the plot. If `NULL`, all cells are plotted.
#' @param slot The assay slot to use for dimensionality reduction, default is `"data"`. Alternatively, can be set to `"logcounts"`, `"counts"`, etc.
#' @param mapping A mapping for aesthetics such as `color`, `size`, etc., passed to `ggplot2`. By default, `NULL` and `color` will be mapped automatically.
#' @param geom A ggplot2 geometry to use for plotting. Default is `sc_geom_point`, which generates scatter points.
#' @param linewidth A numeric value indicating the line width for the trajectory paths. Default is 0.9.
#' @param ... Additional parameters passed to `ggplot2` plotting functions or the geometry.
#' 
#' @return A ggplot object displaying the UMAP plot with overlaid trajectories.
#' 
#' @examples
#' \dontrun{
#‘ pancreas <- readRDS("./pancreas_sub_sce.rds")
#’ pancreas <- RunSlingshot(sce = pancreas, group = "celltype", reduction = "UMAP")
#' p <- lineage_plot(sce = pancreas, group = "celltype", reduction = "UMAP")
#' print(p)
#' }
#' 
#' @import SingleCellExperiment
#' @import ggplot2
#' @import ggsc
#' 
#' @export
lineage_plot <- function(
    sce,
    dims = c(1, 2),
    group,
    reduction = "UMAP",
    lineages = NULL,
    cells = NULL,
    slot = "data",
    mapping = NULL,
    geom = sc_geom_point,
    linewidth = 0.9, ...){
  
  if (!has_trajectory(sce)) {
      stop("No slingshot results found. Please run RunSlingshot first.")
  }
  sds <- get_trajectory(sce, "dataset")
  
  # Use accessor to get curves safely
  lineage_curves <- slingshot::slingCurves(sds)
  
  if(is.null(lineages)){
    lineages <- names(lineage_curves)
  }
  
  ## set the lineage colors (adjust to ensure enough colors for the number of lineages)
  ## set the lineage data for the input of ggplot2
  lineage_data <- lapply(seq_len(length(lineage_curves)), function(i) {
    temp <- as.data.frame(lineage_curves[[i]]$s) |> setNames(c("x", "y"))
    temp$lineage <- sprintf("Lineage%d", i)
    return(temp)
  }) |> yulab.utils::rbindlist()
  
  # Subset lineage_data based on the input lineage
  lineage_data <- subset(lineage_data, lineage_data$lineage %in% lineages)

  # Existing UMAP plot with labels
  p1 <- ggsc::sc_dim(object = sce, reduction = reduction, cells = cells, slot = slot,
                     mapping = mapping, geom = geom)  +
    geom_path(data = lineage_data, 
            aes(x = .data$x, y = .data$y, 
                color = .data$lineage, group = .data$lineage),
            arrow = arrow(type = "open", length = unit(0.1, "inches")),
            linewidth = linewidth) 
  
  return(p1)
}


#' Visualize pseudotime data on reduction plots
#'
#' @param sce A SingleCellExperiment object after `RunSlingshot()` which contains the pseudotime and reduction data.
#' @param reduction A string specifying the reduction method to be used for plotting. Default is `NULL`.
#' @param pseudotime.data The name of the column in the colData of `sce` that contains the pseudotime data. Default is `"slingPseudotime"`.
#' @param lineage selected lineage to be displayed, all lineages will be used by default.
#' @param ... Additional parameters to be passed to `ggplot2` and other plot settings.
#'
#' @return A `ggplot` object representing the pseudotime visualization.
#'
#'
#' @import ggplot2
#' @import SingleCellExperiment
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 scale_color_viridis_c
#'
#' @export
pseudo_plot <- function(
    sce,
    reduction = NULL,
    pseudotime.data = "slingPseudotime",
    lineage = NULL,
    ...){
  
  reduction <- match.arg(reduction, c("UMAP", "PCA", "tSNE"))
  
  if (!is(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  if (!pseudotime.data %in% colnames(colData(sce))) {
    stop("SCE lacks pseudotime data. Please run 'RunSlingshot' before this function.")
  }
  
  ## Extract pseudotime data
  sling <- sce[[pseudotime.data]]
  # Ensure it is a data frame for manipulation
  if (is.matrix(sling)) {
      sling <- as.data.frame(sling)
  }
  sling$cell <- rownames(sling)
  if (is.null(sling$cell)) {
      sling$cell <- colnames(sce)
  }
  d2 <- tidyr::pivot_longer(sling, 
                            cols = starts_with("Lineage"),
                            names_to = "Lineage", 
                            values_to = "Pseudotime")
  
  ## Extract reduction data
  reduction_data <- as.data.frame(reducedDims(sce)[[reduction]]) |>
    setNames(c("x", "y"))
  reduction_data$cell <- rownames(reduction_data)
  
  ## Merge pseudotime data and reduction data
  merged <- merge(d2, reduction_data, by = "cell")
  
  ## If a specific lineage is provided, filter the data and set title
  if (!is.null(lineage)) {
    if (!lineage %in% unique(d2$Lineage)) {
      stop(paste("The provided lineage '", lineage, "' is not valid. Please choose one from:", 
                 paste(unique(d2$Lineage), collapse = ", "), "."))
    }
    merged <- merged[merged$Lineage == lineage, ]
    plot_title <- paste(lineage)
  } else {
    plot_title <- NULL
  }
  
  ## Use ggplot2 to plot pseudotime data of cells
  p <- ggplot(merged, aes(x = .data$x, y = .data$y, color = .data$Pseudotime)) +
    geom_point() +
    theme_bw() +
    scale_color_viridis_c(option = "C") +
    theme(
      axis.text = element_blank(),    
      axis.title = element_blank(),   
      axis.ticks = element_blank(),   
      panel.grid = element_blank(),   
      strip.text = element_text(size = 12, face = "bold"),  
      legend.position = "right",      
      legend.title = element_text(size = 12),  
      legend.text = element_text(size = 10),   
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14) 
    ) +
    labs(
      title = plot_title,  
      x = NULL,      
      y = NULL       
    )
  
  ## If a specific lineage is selected, remove facet_wrap
  if (is.null(lineage)) {
    p <- p + facet_wrap(~Lineage, ncol = 3)
  }
  
  return(p)
}


#' Generate gene expression curves across pseudotime
#'
#' This function visualizes the expression of specified genes across pseudotime using either 
#' smoothed line plots or GAM (Generalized Additive Models) predictions. It generates a 
#' ggplot showing the expression of selected genes along pseudotime for different cell fates.
#'
#' @param sce A SingleCellExperiment object containing the expression data and pseudotime.
#' @param assay_name The assay to use for extracting gene expression data (default: "logcounts").
#' @param features A vector of gene names for which to plot expression curves.
#' @param lineaegs A vector of cell fates to visualize (optional). If NULL, all available fates will be used.
#' @param pseudotime.data The column name in the SingleCellExperiment object containing pseudotime data (default: "slingPseudotime").
#' @param method A character vector specifying the method for curve fitting. Options are "smooth" (default) or "gam" (Generalized Additive Model).
#' @param point A logical value indicating whether to include points for individual cells (default: FALSE).
#' @param line.size Line thickness for the plot (default: 0.8).
#' @param alpha Transparency for the points in the plot (default: 0.5).
#' @param size Size of the points in the plot (default: 1).
#' @param ncol Number of columns to use in facet_wrap (default: NULL).
#' @param se A logical value indicating whether to include the standard error for smoothing (default: TRUE).
#'
#' @return A ggplot object displaying gene expression curves across pseudotime.
#' @export
genecurve_plot <- function(
    sce,
    assay_name = "logcounts",
    features,
    lineaegs = NULL,
    pseudotime.data = "slingPseudotime",
    method = c("smooth","gam"),
    point = FALSE,
    line.size = 0.8,
    alpha = 0.5,
    size = 1,
    ncol = NULL,
    se = TRUE){
  
  
  if( !"logcounts" %in% assayNames(sce)){
    sce <- NormalizeData(sce)
  }
  
  # Ensure that pseudotime data is available
  if (!pseudotime.data %in% colnames(colData(sce))) {
    stop("The specified pseudotime data is not available in the SCE object.")
  }
  
  ## extract pseudotime data
  pseudo <- colData(sce)[[pseudotime.data]]
  pseudo <- as.data.frame(pseudo)
  fate_names <- colnames(pseudo) 
  
  ## extract gene expression data
  gdf <- t(assay(sce, assay_name))
  gdf <- as.data.frame(gdf[, features])
  colnames(gdf) <- features
  df <- cbind(pseudo, gdf)
  
  ## reshape data
  d2_gene <- reshape2::melt(df, id.vars = fate_names, 
                            variable.name = "Gene", 
                            value.name = "Expression")
  
  d2 <- reshape2::melt(d2_gene, 
                       id.vars = c("Gene", "Expression"), 
                       measure.vars = fate_names, 
                       variable.name = "cell_fate", 
                       value.name = "Pseudotime")
  
  d2$Gene <- factor(d2$Gene, levels = unique(d2$Gene))
  d2 <- na.omit(d2)
  
  ## baseplot
  p <- base_plot(point, df_melt = d2, ncol)
  
  ## predict gene expression among lineages
  predictions <- data.frame()
  
  if( !is.null(lineaegs)){
    fate_names <- lineaegs
  }
  else{
    fate_names <- colnames(pseudo)
  }
  
  ## perform GAM analysis
  if (method[1] == "gam") {
    predictions <- gam_analysis_SCE(d2, features, fate_names)
  } 
  else {
    predictions <- d2
  }
  
  p <- add_plot_layer(p, method, line.size, se = se, predictions)
  
  return(p)
}



#' Add plot layers based on method choice
#'
#' This function adds the appropriate plot layer to the base plot depending on the selected method. 
#' It either adds a GAM-based line plot or a smoothed curve for the gene expression across pseudotime.
#'
#' @param p The base ggplot object.
#' @param method A character vector specifying the plot method, either "smooth" or "gam".
#' @param line.size Size of the plot line (default: NULL).
#' @param se A logical value indicating whether to display the standard error (default: TRUE).
#' @param plot_data The data used for plotting, either predictions or raw data.
#'
#' @return A ggplot object with the added plot layer.
# @export
#' @noRd
add_plot_layer <- function(p, method, line.size, se, plot_data) {
  
  if (method[1] == "gam") {
    ## using prediction result to plot
    p <- p + geom_line(
      data = plot_data,
      aes(x = .data$Pseudotime, y = .data$Predicted, color = .data$cell_fate),
      size = line.size)  
  } 
  else if (method[1] == "smooth") {
    # using actual gene expression data to plot
    p <- p + geom_smooth(
      data = plot_data,
      aes(x = .data$Pseudotime, y = .data$Expression, color = .data$cell_fate),
      size = line.size, se = se)  
  }
  
  return(p)
}


#' Create the base plot for gene expression over pseudotime
#'
#' This function generates the base plot framework, which includes labels, faceting by gene, and 
#' basic styling. It can optionally add points to represent individual cell data.
#'
#' @param point A logical value indicating whether to add points for individual cells (default: FALSE).
#' @param df_melt A data frame containing the melted version of the gene expression data.
#' @param alpha The transparency level for points (default: 0.5).
#' @param size The size of the points (default: 1).
#' @param ncol The number of columns for faceting (default: NULL).
#'
#' @return A ggplot object representing the base plot.
# @export
#' @noRd
base_plot <- function(
    point = FALSE, 
    df_melt, 
    alpha, 
    size, 
    ncol = NULL) {
  p <- ggplot() +
    labs(x = "Pseudotime", y = "Gene Expression", color = "Cell Fate") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(face = "bold", size = 10)) +
    facet_wrap(~ Gene, scales = "free", ncol = ncol) +
    scale_x_continuous(expand = c(0, 0))
  
  if (point) {
    p <- p + geom_point(
      data = df_melt,
      aes(.data$Pseudotime, .data$Expression, color = .data$cell_fate),
      alpha = alpha, size = size,
      shape = 15)
  }
  
  return(p)
}


#' Perform GAM analysis on gene expression data
#'
#' This function fits a Generalized Additive Model (GAM) to the gene expression data for each
#' gene and lineage, predicting gene expression across pseudotime.
#'
#' @param df_melt A data frame containing melted gene expression data, including pseudotime and cell fate.
#' @param features A vector of gene names to analyze.
#' @param fate_names A vector of cell fate names to consider for each gene.
#'
#' @return A data frame with predicted gene expression values based on the fitted GAM models.
#' @importFrom stats predict
# @export
#' @noRd
gam_analysis_SCE <- function(df_melt, features, fate_names) {
  

  predictions <- data.frame()
  
  for (gene in features){
    for (fate in fate_names){
      
      # subset specific gene and lineages data from d2
      subset_df <- df_melt[df_melt$Gene == gene & df_melt$cell_fate == fate, ]
  
      # using GAM model to predict pseudotime expression
      model <- mgcv::gam(Expression ~ s(Pseudotime), data = subset_df)
  
      pseudo_seq <- seq(min(subset_df$Pseudotime), max(subset_df$Pseudotime), length.out = 200)
      grid_df <- data.frame(Pseudotime = pseudo_seq, Gene = gene, cell_fate = fate, Expression = NA)
  
      # add prediction result to final dataframe
      grid_df$Predicted <- predict(model, newdata = grid_df)
      predictions <- rbind(predictions, grid_df)
    }
  }
  predictions$Gene <- factor(predictions$Gene, levels = unique(predictions$Gene))
  
  return(predictions)
}



# Rescale values to the range [0, 1]
rescale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}


#' Compute gradient of a vector around the peak position
#'
#' This function calculates the gradient score of a numeric vector `row` by first identifying the position
#' of the peak (maximum value) in the vector. It then computes the average difference between subsequent values
#' in a small window around the peak. If all values in the vector are `NaN` or all zeros, the function returns `NA`.
#'
#' @param row A numeric vector representing the expression values of a gene across pseudotime.
#' @param n The window size (in number of elements) around the peak to calculate the gradient. Default is 10.
#'
#' @return A numeric value representing the average gradient score around the peak, or `NA` if the row is 
#'         all `NaN` or all zeros.
#' 
# @export
#' @noRd
compute_gradient <- function(row, n=10) {
  if (all(is.nan(row)) || all(row == 0)) {
    return(NA)  
  }
  
  peak_position <- which.max(row)
  left_bound <- max(1, peak_position - n)
  right_bound <- min(length(row), peak_position + n)
  
  # Extract the sub-section of the row around the peak
  sub_section <- row[left_bound:right_bound]
  
  # Calculate the gradient score as the average difference between subsequent positions
  gradient_score <- mean(diff(sub_section))
  
  return(gradient_score)
}

#' @name lineage_plot
#' @aliases plot_lineage
#' @export
plot_lineage <- lineage_plot

#' @name pseudo_plot
#' @aliases plot_pseudotime
#' @export
plot_pseudotime <- pseudo_plot

#' @name genecurve_plot
#' @aliases plot_genecurve
#' @export
plot_genecurve <- genecurve_plot
