
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
#' @param sce A SingleCellExperiment object
#' @param batch A factor specifying the batch of origin for each cell if only one batch is supplied in `...`. This will be ignored if two or more batches are supplied.
#' @param HVG user specified high variable genes. If NULL, it will be extracted automatically from the HVG supplied in `...`.
#' @param nHVG number of HVG. Only works for `HVG = NULL` and will be extracted from the info supplied in `...` by default. 
#' @param assay.type specify the assay used in batch correction. If `assay.type` is 'counts', it will try to use 'logcounts' or perform log transformation.
#' @param PARAM batch correction method, see also `batchelor::batchCorrect()`.
#' @param restrict A list of length equal to the number of objects in `...`. Each entry of the list corresponds to one batch and specifies the cells to use when computing the correction.
#' @param correct.all A logical scalar indicating whether to return corrected expression values for all genes, even if `subset.row` is set. Used to ensure that the output is of the same dimensionality as the input.
#' @return A SingleCellExperiment object
#' @export
BatchRemover <- function (sce, batch = NULL, HVG = NULL, nHVG = 5000, 
                          assay.type = "logcounts", PARAM = batchelor::FastMnnParam(), 
                          restrict = NULL, correct.all = FALSE) { 
  
  if(!is(sce,"SingleCellExperiment")){
    stop("The input object needs to be a SingleCellExperiment object.")
  }
  
  if(is.null(batch)) {  
    batch <- sce$batch
  }
  
  if (assay.type == "counts") {
    assay_all <- names(SummarizedExperiment::assays(sce))
    has_logcounts <-  "logcounts" %in% assay_all
    if (all(has_logcounts)) {
      message("'assay.type' will be set to 'logcounts' automatically.")
    } else {
      message("'assay.type' will be log transformed using 'batchelor::multiBatchNor()")
      sce <- batchelor::multiBatchNorm(sce, batch = batch, assay.type = "counts", preserve.single = TRUE)
    }
  } 
  
  #默认设置5000，结果图会好看一点不会是密密的团
  if (is.null(nHVG)) {
    nHVG <- sce@metadata$nVariableFeatures
    message("nHVG is automatically set to ", nHVG, '.')
  }
  
  hvginfo <- rowData(sce)
  
  #提取高度可变基因
  method <- get_hvg_method(sce)
  
  if (is.null(HVG)) {    
    if (method == "scran") {
      HVG <- scran::getTopHVGs(hvginfo, nHVG)
    } else {
      HVG <- getTopHVGs_seurat(hvginfo, nHVG)
    }
  }
  
  #校正
  corrected <- batchelor::batchCorrect(sce, batch = batch, 
                                       subset.row = HVG, PARAM = PARAM, 
                                       restrict = restrict, correct.all = correct.all) 
  
  
  corrected@metadata <- sce@metadata
  corrected@metadata$nVariableFeatures <- nHVG
  subset.rowdata <- rowData(sce)[HVG,]
  new.rowdata <- rowData(corrected)
  
  #为防行名的顺序出错多做一步（行数和行肯定是对应的）
  subset.rowdata <- subset.rowdata[rownames(new.rowdata),]
  rowData(corrected) <- cbind(subset.rowdata, new.rowdata)
  return(corrected)
} 



overlap <- function(x) {
  if (!is.list(x)) return(x)
  if (length(x) <= 1) return(x)

  res <- x[[1]]
  for (i in 2:length(x)) {
    res <- intersect(res, x[[i]])
  }
  return(res)
}

#' merge a list of SingleCellExperiment objects to a single SingleCellExperiment object
#' 
#' 
#' 
#' @title sce_merge 
#' @param sce_list a named list of SingleCellExperiment objects 
#' @param combineVarParams parameters that passed to `CombineVariableFeatures()`
#' @return A SingleCellExperiment object
#' @importFrom utils modifyList
#' @export
sce_merge <- function(sce_list, combineVarParams = list(equiweight = TRUE, ncells = NULL)){

  if (length(sce_list) < 2) {
    stop("At least two objects are required for merging.")
  }
  
  is_sce <- sapply(sce_list, function(x) is(x, "SingleCellExperiment"))
  
  if (!all(is_sce)) {
    bad_indices <- which(!is_sce)
    stop_message <- paste(
      "Object",paste(bad_indices, collapse = ","),"is/are not SingleCellExperiment objects."
    )
    stop(stop_message)
  }
  
  sce_list <- IntersectGenes(sce_list)
  
  asys <- lapply(sce_list, \(x) names(SummarizedExperiment::assays(x))) 
  asys <- overlap(asys)

  cbn_assays <- list()
  for (i in seq_along(asys)) {
    lasy <- lapply(sce_list, \(x) assay(x, asys[i]))
    cbn_assays[[i]] <- do.call(cbind, lasy)
  }
  names(cbn_assays) <- asys

  combineVarParams$all.batches <- sce_list
  combined_hvginfo <- do.call(CombineVariableFeatures, combineVarParams)
  
  hvg_temp <- setdiff(colnames(rowData(sce_list[[1]])),colnames(combined_hvginfo))
  rowdata_combined <- cbind(rowData(sce_list[[1]])[,hvg_temp],combined_hvginfo)

  for (i in 2:length(sce_list)) {
    idx <- setdiff(colnames(rowData(sce_list[[i]])),colnames(rowdata_combined))
    rowdata_combined <- cbind(rowData(sce_list[[i]])[, idx], rowdata_combined)
  }
  
  coldata_combined <- do.call(rbind, lapply(sce_list, colData))
  
  metadata_combined <- S4Vectors::metadata(sce_list[[1]])
  
  for (i in 2:length(sce_list)) {
    metadata_combined <- modifyList(metadata_combined, S4Vectors::metadata(sce_list[[i]]))
  }
  
  
  combined_sce <- SingleCellExperiment(
    assays = cbn_assays,
    rowData = rowdata_combined,
    colData = coldata_combined,
    metadata = metadata_combined
  )
  
  combined_sce$batch <- factor(rep(seq_along(sce_list), sapply(sce_list, ncol)))
  
  return(combined_sce)
}
