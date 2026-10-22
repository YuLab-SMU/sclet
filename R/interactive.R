#' Run Interactive Explorer
#' 
#' Launch an iSEE interactive application for exploring the SingleCellExperiment object.
#' 
#' @title RunExplorer
#' @param object A SingleCellExperiment object
#' @param ... Additional arguments passed to iSEE::iSEE
#' @return A Shiny app object
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' \dontrun{
#' # needs the 'iSEE' Suggests package
#' app <- RunExplorer(sce)
#' }
RunExplorer <- function(object, ...) {
    if (!requireNamespace("iSEE", quietly = TRUE)) {
        stop("Package 'iSEE' is needed for this function to work. Please install it.")
    }
    
    # iSEE works directly on SCE objects
    iSEE::iSEE(object, ...)
}
