#' wrapper function for running SCENIC on default parameter. 
#' 
#' https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
#' 
#' @param sce A singlecellExperiment object.
#' @param species The species of the data, one of human, mouse or drosophila_melanogaster.
#' @param nCores The number of cores to use for parallel processing.
#' @param dbDir The directory path to the cisTarget databases files.
#' @param assay_name the name of assay.
#' 
#' @return The class contains the options/settings for a run of SCENIC. 
#' 
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom utils data
#' @export 
runSCENIC <- function(sce,species,nCores,dbDir,
                      assay_name = "counts"){
  # Create directories to store process data and output results
  if (!dir.exists("int")) {
    dir.create("int")
  }
  if (!dir.exists("output")) {
    dir.create("output")
  }

  # Initialize settings
  org <- switch(species, 
                "human" = "hgnc", 
                "mouse" = "mgi", 
                "drosophila_melanogaster" = "dmel", 
                stop("Unsupported species"))
  data(defaultDbNames,envir = environment())
  dbs <- SCENIC::defaultDbNames[[org]]
  
  if (org == "mgi") {
    data(list = "motifAnnotations_mgi",envir = environment())
  } else if (org == "hgnc") {
    data(list = "motifAnnotations_hgnc",envir = environment())
  } else if (org == "dmel") {
    data(list = "motifAnnotations_dmel",envir = environment())
  } else {
    stop("Unsupported species. Please choose 'mgi', 'hgnc', or 'dmel'.")
  }

  scenicOptions <- SCENIC::initializeScenic(org = org, dbDir = dbDir, dbs = dbs, nCores = nCores)

  # Process the cell annotation information of the sce object for plot later
  if ("CellType" %in% colnames(SingleCellExperiment::colData(sce))) {
    cellInfo <- as.data.frame(SingleCellExperiment::colData(sce))
    saveRDS(cellInfo, file="int/cellInfo.Rds") 
    scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
} else {
  warning("Warning: 'cellInfo' not found, proceeding without cell annotations.")
}

  if ("CellTypeColor" %in% colnames(SingleCellExperiment::colData(sce))) {
    colVars <- as.list(sce$CellTypeColor)
    saveRDS(colVars, file="int/colVars.Rds")
    scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
} else {
  warning("Warning: 'colVars' not found, proceeding without color annotations.")
}
  
  # Co-expression network 
  exprMat <- SummarizedExperiment::assay(sce, assay_name)
  genesKept <- SCENIC::geneFiltering(exprMat, scenicOptions = scenicOptions,
                                     minCountsPerGene = 3 * 0.01 * ncol(exprMat),
                                     minSamples = ncol(exprMat) * 0.01)
  exprMat_filtered <- exprMat[genesKept, ]
  SCENIC::runCorrelation(exprMat_filtered, scenicOptions)
  
  set.seed(123)
  exprMat_filtered <- log2(exprMat_filtered+1) 
  SCENIC::runGenie3(exprMat_filtered, scenicOptions)

  # Build and score the GRN
  exprMat_log <- log2(exprMat+1)
  scenicOptions <- SCENIC::runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- SCENIC::runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- SCENIC::runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
  
  saveRDS(scenicOptions,file="int/scenicOptions.Rds")
  return(scenicOptions)
}

