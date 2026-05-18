#' wraper function to run SuperCell
#' 
#' @title RunSuperCell
#' @param object SingleCellExperiment object
#' @param assay selected assay used to build KNN graph. If `layer` is provided,
#' this is treated as a compatibility alias.
#' @param layer layer used to build the metacell representation. If `NULL`,
#' sclet resolves it from `DefaultLayer(object)`.
#' @param nHVG number of HVGs to use
#' @param hvg_method one of 'seurat' or 'scran', see also `FindVariableFeatures()`
#' @param cellname IDs to label the cells, use 'Barcode' by default
#' @param gamma graining level, that is number_of_cells / number_of_metacells
#' @param k.knn number of nearest neighbors to build KNN graph
#' @param name aggregation record id. Defaults to `"supercell"`.
#' @return SingleCellExperiment object
#' @export
RunSuperCell <- function(object, assay = "logcounts", layer = NULL, nHVG = 2000, hvg_method = "scran",
                        cellname = "Barcode", gamma = 20, k.knn = 5, name = "supercell") {

  if (!requireNamespace("SuperCell", quietly = TRUE)) {
      stop("Package 'SuperCell' is needed for this function to work. Please install it.")
  }
  
  # Ensure VariableFeatures are computed
  if (is.null(sclet_get_hvg_nfeatures(object))) {
      object <- FindVariableFeatures(object, nfeatures = nHVG, method = hvg_method)
  }

  # gene expression matrix
  source <- sclet_resolve_expression_source(
      object = object,
      layer = layer,
      assay = assay,
      prefer_nonscaled = TRUE,
      context = "SuperCell aggregation"
  )
  assay <- source$assay
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
  state <- sclet_get_state(object)
  if (!is.null(state$features$hvg)) {
      state$features$hvg$n <- NULL
  }
  mdata$sclet <- state

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
  sce <- sclet_set_active_assay(sce, assay)
  if (!is.null(source$layer) && source$layer %in% Layers(sce)) {
      sce <- sclet_set_active_layer(sce, source$layer)
  }
  sce <- sclet_set_analysis(
      sce,
      "supercell",
      list(
          method = "SuperCell::SCimplify",
          id = name,
          gamma = gamma,
          k.knn = k.knn,
          assay = assay,
          layer = source$layer,
          cellname = cellname,
          parent = list(
              assay = assay,
              layer = source$layer,
              n_cells = ncol(object),
              n_genes = nrow(object)
          ),
          child = list(
              n_cells = length(SC$supercell_size)
          ),
          object = SC
      )
  )
  sce <- sclet_set_analysis_state(
      object = sce,
      type = "aggregation",
      id = name,
      method = "SuperCell::SCimplify",
      inputs = list(
          assay = assay,
          layer = source$layer,
          cellname = cellname,
          n_hvg = length(hvgs),
          parent_n_cells = ncol(object),
          parent_n_genes = nrow(object)
      ),
      artifacts = list(
          analysis_key = "supercell",
          size_col = "size"
      ),
      params = list(
          gamma = gamma,
          k.knn = k.knn,
          hvg_method = hvg_method
      ),
      summary = list(
          n_metacells = length(SC$supercell_size),
          mean_metacell_size = mean(SC$supercell_size)
      )
  )
  sce <- sclet_log_command(
      sce,
      "RunSuperCell",
      params = list(
          assay = assay,
          layer = source$layer,
          nHVG = nHVG,
          hvg_method = hvg_method,
          cellname = cellname,
          gamma = gamma,
          k.knn = k.knn,
          name = name
      ),
      outputs = list(
          analysis = "supercell",
          aggregation = name
      )
  )
  return(sce)
}
