#' read 10x single cell expression matrix folder
#' 
#' @title Read10X
#' @param data.dir directory of the data
#' @param backend backend of the data, default is "memory", support "HDF5"
#' @param ... additional parameters passed to 'DropletUtils::read10xCounts()'
#' @return a 'SingleCellExperiment' object
#' @importFrom DropletUtils read10xCounts
#' @importFrom scuttle uniquifyFeatureNames
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @export
Read10X <- function(data.dir, backend = c("memory", "HDF5"), ...) {
    backend <- match.arg(backend)
    
    if (backend == "HDF5") {
        # Currently DropletUtils::read10xCounts returns a DelayedArray if type="HDF5" is specified,
        # but that usually requires an HDF5 file, not a directory of MTX.
        # However, we can read into memory and then realize to HDF5, 
        # OR if the input is already an HDF5 file (10x .h5), we can read it directly.
        
        # Check if data.dir is a directory or a file
        if (dir.exists(data.dir)) {
            # It's a directory (MTX), read as sparse then write to HDF5
            # Since user wants HDF5 for big data, we should try to avoid full memory load if possible.
            # But DropletUtils::read10xCounts for directory loads into memory (sparse).
            # We will convert it immediately to HDF5-backed.
            sce <- DropletUtils::read10xCounts(samples = data.dir, ...)
            
            # Convert to HDF5Backed
            sce <- HDF5Array::saveHDF5SummarizedExperiment(sce, dir = tempfile())
            
        } else if (file.exists(data.dir) && grepl("\\.h5$", data.dir)) {
            # It's an HDF5 file
            # DropletUtils can read 10x h5 directly, but returns sparse matrix by default unless specified?
            # Actually DropletUtils::read10xCounts supports .h5
            # To get on-disk, we might need to use HDF5Array::HDF5Array or similar,
            # but standard workflow is often read -> saveHDF5.
            
            sce <- DropletUtils::read10xCounts(samples = data.dir, ...)
            sce <- HDF5Array::saveHDF5SummarizedExperiment(sce, dir = tempfile())
        } else {
             sce <- DropletUtils::read10xCounts(samples = data.dir, ...)
             sce <- HDF5Array::saveHDF5SummarizedExperiment(sce, dir = tempfile())
        }
        
    } else {
        sce <- DropletUtils::read10xCounts(samples = data.dir, ...)
    }

    rownames(sce) <- scuttle::uniquifyFeatureNames(
        rowData(sce)$ID, rowData(sce)$Symbol)

    sce <- QCMetrics(sce)

    return(sce)
}
