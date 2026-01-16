#' Read H5AD file
#' 
#' Wrapper for zellkonverter::readH5AD to import AnnData objects as SingleCellExperiment.
#' 
#' @title ReadH5AD
#' @param file Path to the .h5ad file
#' @param ... Additional arguments passed to zellkonverter::readH5AD
#' @return A SingleCellExperiment object
#' @export
ReadH5AD <- function(file, ...) {
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
        stop("Package 'zellkonverter' is needed for this function to work. Please install it.")
    }
    
    zellkonverter::readH5AD(file, ...)
}

#' Write H5AD file
#' 
#' Wrapper for zellkonverter::writeH5AD to export SingleCellExperiment objects as .h5ad.
#' 
#' @title WriteH5AD
#' @param object A SingleCellExperiment object
#' @param file Path to the output .h5ad file
#' @param ... Additional arguments passed to zellkonverter::writeH5AD
#' @return Invisible NULL
#' @export
WriteH5AD <- function(object, file, ...) {
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
        stop("Package 'zellkonverter' is needed for this function to work. Please install it.")
    }
    
    zellkonverter::writeH5AD(object, file, ...)
}

#' Convert to Seurat object
#' 
#' Convert a SingleCellExperiment object to a Seurat object.
#' 
#' @title as.Seurat
#' @param object A SingleCellExperiment object
#' @param ... Additional arguments passed to Seurat::as.Seurat
#' @return A Seurat object
#' @export
as.Seurat <- function(object, ...) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is needed for this function to work. Please install it.")
    }
    
    # Ensure logcounts/counts are present and named correctly for Seurat conversion
    # Seurat typically expects 'counts' and 'logcounts'
    
    Seurat::as.Seurat(object, ...)
}

#' Convert to SingleCellExperiment object
#' 
#' Convert a Seurat object to a SingleCellExperiment object.
#' 
#' @title as.SCE
#' @param object A Seurat object
#' @param ... Additional arguments passed to Seurat::as.SingleCellExperiment
#' @return A SingleCellExperiment object
#' @export
as.SCE <- function(object, ...) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is needed for this function to work. Please install it.")
    }
    
    Seurat::as.SingleCellExperiment(object, ...)
}
