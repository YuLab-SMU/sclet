
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

CombineVariableFeatures <- function(all.batches, equiweight = TRUE, ncells = NULL) {
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

    all.batches <- IntersectGenes(all.batches)

    hvg.info <- lapply(all.batches, function(sce) {
        hvgcols <- sce@metadata$hvgcols
        return(rowData(sce)[, hvgcols])
    })

    if (method == "scran") {
        res <- scran::combineVar(hvg.info, equiweight = equiweight, ncells = ncells)
    } else {
        res <- combine_blocked(hvg.info, equiweight = equiweight, ncells = ncells)
    }

    return(res)
}


extract_hvgs <- function(all.batches, equiweight = TRUE, ncells = NULL, nHVG = NULL){
    hvginfo <- CombineVariableFeatures(all.batches = all.batches, equiweight = equiweight, ncells = ncells)
    
    if (is.null(nHVG)) {
        nHVG <- sapply(all.batches, \(x) x@metadata$nVariableFeatures) |> max()
        message("nHVG is automatically set to ", nHVG, '.')
    }

    method <- get_hvg_method(all.batches[[1]])
    if (method == "scran") {
        return(scran::getTopHVGs(hvginfo, nHVG))
    }

    hvg <- rownames(hvginfo)[order(hvginfo$variance.standardized, decreasing = TRUE)]
    return(hvg[1:nHVG])
}


BatchRemover <- function (..., batch = NULL, HVG = NULL, nHVG = NULL, restrict = NULL, correct.all = TRUE, 
          assay.type = "counts", PARAM = batchelor::FastMnnParam(), multi.norm.args = list()) 
{ 
    all.batches <- scuttle::.unpackLists(...)  
    #重新计算了scale.factor,并进行LogNorm
    all.batches <- do.call(batchelor::multiBatchNorm, 
                        c(all.batches, list(batch = batch, assay.type = assay.type, preserve.single = TRUE
                            ), multi.norm.args)
                    ) #进行多批次的修正

    if (is.null(HVG)) {
       HVG <- extract_hvgs(all.batches, equiweight = TRUE, ncells = NULL, nHVG = nHVG)
    }

    # batch correction
    corrected <- batchelor:::batchCorrect(all.batches, batch = batch, restrict = restrict, 
                              correct.all = correct.all, subset.row = HVG, PARAM = PARAM) 

    return(corrected)
} 


combine_blocked <- function (blocks, equiweight=TRUE, ncells=ncells) { # <-- 这个ncells是没有赋值的
  valid = ncells >= 2L # <-- 为什么是2
  weights = ncells #比重
  #length(blocks) == 1L时，不需要合并
  if (length(blocks) == 1L) {
    return(blocks[[1]])
  }
  rn <- unique(lapply(blocks, rownames))
  if (length(rn) != 1L) {
    stop("gene should be the same")
  } #检查各批次的基因是否一致

  if (equiweight) {
    weights <- rep(1, length(blocks))
    #若设置等比重，则比重为length(blocks)
  } else if (is.null(weights)) {
    stop("'weights must be specified")
    #若没有默认等比重，且没有设置比重，则报错
  }
  
  original <- blocks
  
  if (!any(valid)) {
    stop("no entry of 'blocks' has positive weights")
  } #若有valid为FALSE，则存在部分数据并没有以一定的比重进行合并
  
  blocks <- blocks[valid]
  weights <- weights[valid]
  combined <- list()
  for (i in ave.fields) {  
    extracted <- lapply(blocks, "[[", i = i)
    extracted <- mapply("*", extracted, weights, SIMPLIFY = FALSE, 
                        USE.NAMES = FALSE)
    averaged <- Reduce("+", extracted)/sum(weights)
    combined[[i]] <- averaged
  }#加权平均
  
  #输出结果
  output <- S4Vectors::DataFrame(combined, row.names = rn[[1]])
  # <-- 会报错
  # output$per.block <- do.call(S4Vectors::DataFrame, c(lapply(original, I),
  #                                         list(check.names = FALSE))) 
  output
} 

