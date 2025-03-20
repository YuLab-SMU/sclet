
IntersectGenes <- function(..., verbose = TRUE) {
    # only common genes were kept
    all.batches <- scuttle::.unpackLists(...)  
    all.genes <- lapply(all.batches, rownames)
    common_genes <- Reduce(intersect, all.genes)
    
    if (length(common_genes) == 0) stop("there are no common genes")
    
    if (verbose) {
      r <- round(length(common_genes) / sapply(all.genes, length) * 100)
      message("Common gene ratio:")
      for (i in seq_along(r)) {
        message(sprintf("\t%s\t%g%s", names(r)[i], r[i], '%'))
      }
    }

    all.batches <- lapply(all.batches, function(batch) {
                              batch[common_genes, , drop = FALSE]
                            })
    return(all.batches)
}

get_hvg_method <- function(object) {
    return(object@metadata$hvgmethod)
}

get_hvginfo <- function(object) {
    hvgcols <- object@metadata$hvgcols
    return(rowData(object)[, hvgcols])    
}

CombineVariableFeatures <- function(all.batches, ...) {
    method <- sapply(all.batches, get_hvg_method)
    no_hvg <- sapply(method, is.null)
    
    if (any(no_hvg)) {
        msg <- sprintf("%s%s%s.",
            'Forget to run `FindVariableFeatures()`?\n', 
            'Please check the object with the index of ', 
            paste(which(no_hvg), collapse=",")
        )
        stop(msg)
    }
    method <- unique(method)
    if (length(method) > 1) stop("Please make sure HVGs were identified with the same method.")

    # all.batches <- IntersectGenes(all.batches)

    hvg.info <- lapply(all.batches, get_hvginfo)

    if (method == "seurat") {
        hvg.info <- lapply(hvg.info, function(x) {
            x$p.value = 1
            return(x)
        })
    }
    res <- scran::combineVar(hvg.info, ...)
    if (method == "seurat") {
        res$p.value <- NULL
        res$FDR <- NULL
    }

    return(res)
}

#' Batch correction
#' 
#' This is a wrapper function of `batchelor::batchCorrect()`.
#' 
#' @title BatchRemover
#' @param ... One or more or a list of SingleCellExperiment objects 
#' @param batch A factor specifying the batch of origin for each cell if only one batch is supplied in `...`. This will be ignored if two or more batches are supplied.
#' @param HVG user specified high variable genes. If NULL, it will be extracted automatically from the HVG supplied in `...`.
#' @param nHVG number of HVG. Only works for `HVG = NULL` and will be extracted from the info supplied in `...` by default. 
#' @param assay.type specify the assay used in batch correction. If `assay.type` is 'counts', it will try to use 'logcounts' or perform log transformation.
#' @param PARAM batch correction method, see also `batchelor::batchCorrect()`.
#' @param restrict A list of length equal to the number of objects in `...`. Each entry of the list corresponds to one batch and specifies the cells to use when computing the correction.
#' @param correct.all A logical scalar indicating whether to return corrected expression values for all genes, even if `subset.row` is set. Used to ensure that the output is of the same dimensionality as the input.
#' @param combineVarParams parameters passed to `scran::combineVar()` for combining HVG information when extracting HVG supplied in `...`.
#' @return A SingleCellExperiment object
#' @export
BatchRemover <- function (..., batch = NULL, HVG = NULL, nHVG = NULL, 
          assay.type = "logcounts", PARAM = batchelor::FastMnnParam(), 
          restrict = NULL, correct.all = FALSE, 
          combineVarParams = list(equiweight = TRUE, ncells = NULL)
        ) { 
    all.batches <- IntersectGenes(...) 
    if (assay.type == "counts") {
        assay_nm <- lapply(all.batches, \(x) names(SummarizedExperiment::assays(x)))
        has_logcounts <- sapply(assay_nm, \(x) "logcounts" %in% x)
        if (all(has_logcounts)) {
            message("'assay.type' will be set to 'logcounts' automatically.")
        } else {
            message("'assay.type' will be log transformed using 'batchelor::multiBatchNor()")
            all.batches <- do.call(batchelor::multiBatchNorm, 
                            c(all.batches, list(batch = batch, assay.type = assay.type, preserve.single = TRUE), 
                            multi.norm.args = list())
                        )
        }
    } 

    if (is.null(nHVG)) {
        nHVG <- sapply(all.batches, \(x) x@metadata$nVariableFeatures) |> max()
        message("nHVG is automatically set to ", nHVG, '.')
    }
    combineVarParams$all.batches = all.batches
    hvginfo <- do.call(CombineVariableFeatures, combineVarParams)
    method <- get_hvg_method(all.batches[[1]])

    if (is.null(HVG)) {    
        if (method == "scran") {
            HVG <- scran::getTopHVGs(hvginfo, nHVG)
        } else {
            HVG <- getTopHVGs_seurat(hvginfo, nHVG)
        }
    }

    # this is actually equivalent to 
    # merge the all.batches to a single SingleCellExperiemnt object
    # and then apply fastMNN to peform the batch correction

    corrected <- batchelor::batchCorrect(all.batches, batch = batch, 
                            subset.row = HVG, PARAM = PARAM, 
                            restrict = restrict, correct.all = correct.all) 


    
    corrected@metadata$nVariableFeatures <- nHVG
    corrected@metadata$hvgmethod <- method

    i <- setdiff(names(rowData(all.batches[[1]])), names(hvginfo))
    xx <- rowData(all.batches[[1]])[,i]
    yy <- cbind(xx, hvginfo[rownames(xx), ])
    rd <- rowData(corrected)
    yy <- cbind(yy[rownames(rd), ], rd)
    rowData(corrected) <- yy

    return(corrected)
} 
