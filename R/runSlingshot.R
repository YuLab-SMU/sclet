#' Run Slingshot Trajectory Analysis
#' 
#' Runs Slingshot on a SingleCellExperiment object with specified clustering and dimensionality reduction.
#' 
#' @param sce A SingleCellExperiment object.
#' @param group The name of the column in colData(sce) containing cluster labels.
#' @param reduction The name of the reduced dimension method in reducedDims(sce) (e.g., "UMAP" or "PCA").
#' @param start_cluster Optional. The starting cluster for trajectory inference. Default is NULL.
#' @param end_cluster Optional. The ending cluster for trajectory inference. Default is NULL.
#' @param reverse A logical value indicating whether to reverse the pseudotime. Default is FALSE.
#' @param align_start A logical value indicating whether to align the starting pseudotime values at the maximum pseudotime. Default is FALSE.
#' @param seed Random seed for reproducibility.
#' 
#' @importFrom SingleCellExperiment colData reducedDims
#' @importFrom graphics plot lines
#' @return A SingleCellExperiment object with added pseudotime information in colData.
#' @export

RunSlingshot <- function(
    sce,
    group,
    reduction = "UMAP",
    start_cluster = NULL,
    end_cluster = NULL,
    reverse = FALSE,
    align_start = FALSE,
    seed = 2025){
  
  if (!requireNamespace("slingshot", quietly = TRUE)) {
      stop("Package 'slingshot' is needed for this function to work. Please install it.")
  }

  reduction <- match.arg(reduction, c("UMAP", "PCA", "tSNE"))
  
  ## check the type of data
  if (!is(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  ## check the cluster result in sce obj
  if (!group %in% colnames(colData(sce))) {
    stop("The specified cluster column does not exist in colData(sce).")
  }
  
  ## check the reduction result in sce obj
  if (!reduction %in% names(reducedDims(sce))) {
    stop("The specified reduction method does not exist in reducedDims(sce).")
  }

  ## run slingshot
  set.seed(seed)
  
  sce <- slingshot::slingshot(data = sce, clusterLabels = group, 
                              reducedDim = reduction, start.clus = start_cluster, 
                              end.clus = end_cluster)
  
  ## reverse the pseudotime and align the start_cell to the maximum if necessary
  slingpseudotime_cols <- grep("^slingPseudotime", colnames(colData(sce)), value = TRUE)
  
  df <- as.data.frame(slingshot::slingPseudotime(sce))
  colnames(df) <- slingpseudotime_cols
  
  if (isTRUE(reverse)){
    if (isTRUE(align_start)){
      df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    }
    else{
      df <- max(df, na.rm = TRUE) - df
    }
  }
  colData(sce)[, slingpseudotime_cols] <- df
  
  ssds <- slingshot::SlingshotDataSet(sce) 
  
  # Store Slingshot results
  colData(sce)$slingPseudotime <- as.data.frame(slingshot::slingPseudotime(ssds))
  colData(sce)$slingBranchID <- slingshot::slingBranchID(ssds)
  # SummarizedExperiment::metadata(sce)$slingshot_info <- ssds
  sce@metadata$slingshot_info <- ssds

  # Store PseudotimeOrdering S4 object in metadata
  # SummarizedExperiment::metadata(sce)$slingshot <- sce$slingshot
  sce@metadata$slingshot <- sce$slingshot
  sce$slingshot <- NULL
  
  return(sce)
}

