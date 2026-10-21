#' @keywords internal
#' @aliases package-sclet
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- NormalizeData(sce)
#' sce <- RunPCA(sce, ncomponents = 5)
#' Status(sce)$available$reductions
"_PACKAGE"
