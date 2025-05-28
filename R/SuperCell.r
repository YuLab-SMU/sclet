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
RunSuperCell <- function(object, assay = "logcounts", nHVG = 2000, hvg_method = "seurat",
                        cellname = "Barcode", gamma = 20, k.knn = 5) {

  # gene expression matrix
  GE <- SummarizedExperiment::assay(object, assay)
  colnames(GE) <- colData(object)[, cellname]

  hvgs <- VariableFeatures(object, nfeatures = nHVG, method = hvg_method)

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
