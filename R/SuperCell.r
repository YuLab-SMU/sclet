#' wraper function to run SuperCell
#' 
#' @title RunSuperCell
#' @param object SingleCellExperiment object
#' @param assay selected assay used to build KNN graph
#' @param nHVG number of HVGs to use
#' @param hvg_method one of 'seurat' or 'scran', see also `FindVariableFeatures()`
#' @param cellname IDs to label the cells, use 'Barcode' by default
#' @param gamma graining level, that is number_of_cells / number_of_metacells
#' @param k.knn number of nearest neighbors to build KNN graph
#' @return SingleCellExperiment object
#' @export
RunSuperCell <- function(object, assay = "logcounts", nHVG = 2000, hvg_method = "scran",
                        cellname = "Barcode", gamma = 20, k.knn = 5) {

  if (!requireNamespace("SuperCell", quietly = TRUE)) {
      stop("Package 'SuperCell' is needed for this function to work. Please install it.")
  }
  
  # Ensure VariableFeatures are computed
  if (is.null(object@metadata$nVariableFeatures)) {
      object <- FindVariableFeatures(object, nfeatures = nHVG, method = hvg_method)
  }

  # gene expression matrix
  GE <- SummarizedExperiment::assay(object, assay)
  if (is.null(colnames(GE))) {
      colnames(GE) <- paste0("Cell", seq_len(ncol(GE)))
  }
  
  if (cellname == "Barcode") {
      # If Barcode is requested but not in colData, check if we can use colnames or add it
      if (!"Barcode" %in% names(colData(object))) {
          # Use colnames as Barcode
          colData(object)$Barcode <- colnames(object)
      }
  }

  colnames(GE) <- colData(object)[, cellname]
  
  hvgs <- VariableFeatures(object, method = hvg_method)

  SC <- SuperCell::SCimplify(
    X = GE, 
    genes.use = hvgs, 
    gamma = gamma,
    k.knn = k.knn
  )

  mdata <- S4Vectors::metadata(object)
  mdata$nVariableFeatures <- NULL
  mdata$SuperCell <- SC

  assay_names <- names(SummarizedExperiment::assays(object))
  asy <- lapply(assay_names, \(nn)  SuperCell::supercell_GE(assay(object, nn), SC$membership))
  names(asy) <- assay_names

  cdata <- data.frame(size = SC$supercell_size, row.names = colnames(asy[[1]]), 
                     stringsAsFactors = FALSE)

  sce <- SingleCellExperiment(
    assays = asy,
    rowData = rowData(object),
    colData = cdata,
    metadata = mdata
  )
  return(sce)
}
