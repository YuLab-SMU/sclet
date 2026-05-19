
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
    return(get_hvg(object, "method"))
}

get_hvginfo <- function(object) {
    hvginfo <- get_hvg(object, "rowData")
    if (is.null(hvginfo)) {
        stop("HVG statistics not found. Please run `FindVariableFeatures()` first.")
    }
    return(hvginfo)    
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
    if (method == "scrapper") {
        res <- combine_variable_features_scrapper(hvg.info)
    } else {
        res <- sclet_combine_var(hvg.info, ...)
    }
    if (method == "seurat") {
        res$p.value <- NULL
        res$FDR <- NULL
    }

    return(res)
}

combine_variable_features_scrapper <- function(hvg.info) {
    suffixes <- c("means", "variances", "fitted", "residuals")
    output <- S4Vectors::DataFrame(row.names = rownames(hvg.info[[1]]))
    for (suffix in suffixes) {
        cur_cols <- vapply(hvg.info, sclet_find_hvg_column, suffix = suffix, hvg_cols = NULL, FUN.VALUE = character(1), USE.NAMES = FALSE)
        if (!length(cur_cols) || any(cur_cols == "")) {
            next
        }
        cur_mat <- do.call(cbind, Map(function(info, col) info[[col]], hvg.info, cur_cols))
        output[[suffix]] <- rowMeans(cur_mat, na.rm = TRUE)
    }
    hvg_cols <- vapply(hvg.info, sclet_find_hvg_column, suffix = "hvg", hvg_cols = NULL, FUN.VALUE = character(1), USE.NAMES = FALSE)
    if (!length(hvg_cols) || any(hvg_cols == "")) {
        stop("scrapper HVG indicator column not found when combining batches.")
    }
    hvg_mat <- do.call(cbind, Map(function(info, col) as.numeric(info[[col]]), hvg.info, hvg_cols))
    output$hvg_frequency <- rowMeans(hvg_mat, na.rm = TRUE)
    output$hvg <- output$hvg_frequency > 0
    output
}

#' Batch correction
#' 
#' State-aware batch correction built on `batchelor::batchCorrect()`.
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
#' @param name Integration record id. Defaults to "batchcorrect".
#' @return A SingleCellExperiment object
#' @export
BatchRemover <- function (sce, batch = NULL, HVG = NULL, nHVG = 5000, 
                          assay.type = "logcounts", PARAM = NULL, 
                          restrict = NULL, correct.all = FALSE, name = "batchcorrect") { 
  
  if (!requireNamespace("batchelor", quietly = TRUE)) {
      stop("Package 'batchelor' is needed for this function to work. Please install it.")
  }

  if (is.null(PARAM)) {
      PARAM <- batchelor::FastMnnParam()
  }
  
  if(!is(sce,"SingleCellExperiment")){
    stop("The input object needs to be a SingleCellExperiment object.")
  }
  
  batch_supplied <- !is.null(batch)

  if(is.null(batch)) {  
    batch <- sce$batch
  }
  
  if (assay.type == "counts") {
    assay_all <- names(SummarizedExperiment::assays(sce))
    has_logcounts <-  "logcounts" %in% assay_all
    if (all(has_logcounts)) {
      message("'assay.type' will be set to 'logcounts' automatically.")
      assay.type <- "logcounts"
    } else {
      message("'assay.type' will be log transformed using 'batchelor::multiBatchNor()")
      sce <- batchelor::multiBatchNorm(sce, batch = batch, assay.type = "counts", preserve.single = TRUE)
      assay.type <- "logcounts"
    }
  } 
  
  #默认设置5000，结果图会好看一点不会是密密的团
  if (is.null(nHVG)) {
    nHVG <- sclet_get_hvg_nfeatures(sce)
  }
  
  hvginfo <- rowData(sce)
  
  # Extract HVGs
  method <- get_hvg_method(sce)
  
  if (is.null(HVG)) {    
    # Fix: Check if method is NULL first to avoid comparison error
    if (is.null(method)) {
      HVG <- getTopHVGs_seurat(hvginfo, nHVG)
    } else if (method == "scran") {
      var_field <- "bio"
      if (!"bio" %in% colnames(hvginfo)) {
          if ("variance.standardized" %in% colnames(hvginfo)) {
              HVG <- getTopHVGs_seurat(hvginfo, nHVG)
              # avoid calling scran::getTopHVGs below
              var_field <- NULL 
          } else if ("total" %in% colnames(hvginfo)) {
              var_field <- "total"
          }
      }
      
      if (!is.null(var_field)) {
          HVG <- choose_hvgs_by_variance(hvginfo, n = nHVG, var.field = var_field)
      }
    } else if (method == "scrapper") {
      HVG <- getTopHVGs_scrapper(hvginfo, nHVG)
    } else {
      # seurat
      HVG <- getTopHVGs_seurat(hvginfo, nHVG)
    }
  }
  
  # Batch correction
  source_state <- sclet_get_state(sce)
  source_commands <- sclet_get_commands(sce)

  corrected <- sclet_muffle_known_warnings(
      batchelor::batchCorrect(
          sce,
          batch = batch,
          subset.row = HVG,
          PARAM = PARAM,
          restrict = restrict,
          correct.all = correct.all
      ),
      patterns = c(
          "'normalizeCounts' is deprecated",
          "You're computing too large a percentage of total singular values, use a standard svd instead.",
          "more singular values/vectors requested than available"
      )
  )
  SingleCellExperiment::reducedDims(corrected) <- S4Vectors::SimpleList()
  SingleCellExperiment::colLabels(corrected) <- NULL
  S4Vectors::metadata(corrected) <- sclet_merge_external_metadata(sce, corrected)
  
  corrected <- sclet_set_hvg_state(
    corrected,
    nfeatures = nHVG,
    method = method,
    hvgcols = sclet_get_hvg_cols(sce),
    selected = HVG
  )
  
  batch_record <- list(
      method = "batchelor::batchCorrect",
      param = PARAM,
      hvg_n = nHVG,
      assay.type = assay.type,
      batch_var = if(is.null(batch)) "internal_batch" else "user_provided",
      timestamp = Sys.time()
  )
  corrected <- sclet_set_analysis(corrected, "batch", batch_record)
  corrected_assay <- if ("corrected" %in% SummarizedExperiment::assayNames(corrected)) {
    "corrected"
  } else if ("reconstructed" %in% SummarizedExperiment::assayNames(corrected)) {
    "reconstructed"
  } else if (assay.type %in% SummarizedExperiment::assayNames(corrected)) {
    assay.type
  } else {
    SummarizedExperiment::assayNames(corrected)[[1]]
  }
  corrected <- sclet_rebuild_internal_state(
    corrected,
    hvg = source_state$features$hvg,
    commands = source_commands,
    active_assay = corrected_assay
  )
  if (assay.type %in% SummarizedExperiment::assayNames(corrected)) {
    corrected <- sclet_set_layer(
      corrected,
      name = assay.type,
      assay = assay.type,
      role = assay.type,
      active = FALSE
    )
  }
  corrected <- sclet_set_layer(
    corrected,
    name = "corrected",
    assay = corrected_assay,
    role = "corrected",
    params = list(
      method = "batchelor::batchCorrect",
      source = assay.type
    )
  )
  corrected <- sclet_set_analysis_state(
    object = corrected,
    type = "integration",
    id = name,
    method = "batchelor::batchCorrect",
    inputs = list(
      assay = assay.type,
      layer = if (assay.type %in% Layers(corrected)) assay.type else NULL,
      batch = batch,
      hvg = HVG,
      merge = {
        merge_state <- sclet_get_state_record(sce, "integration", "merged_inputs")
        if (is.null(merge_state)) {
          NULL
        } else {
          list(
            type = "integration",
            id = "merged_inputs",
            method = merge_state$method
          )
        }
      }
    ),
    artifacts = list(
      layer = "corrected",
      assay = corrected_assay,
      analysis_key = "batch"
    ),
    params = list(
      nHVG = nHVG,
      correct.all = correct.all,
      restrict = restrict
    ),
    summary = list(
      n_hvg = length(HVG),
      batch_source = if (batch_supplied) "supplied" else "colData$batch"
    )
  )
  corrected <- sclet_set_active_assay(corrected, corrected_assay)
  corrected <- sclet_log_command(
    corrected,
    "BatchRemover",
    params = list(
      nHVG = nHVG,
      assay.type = assay.type,
      correct.all = correct.all,
      name = name
    ),
    outputs = list(analysis = "batch", hvg = "selected", integration = name)
  )

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
  
  metadata_combined <- do.call(sclet_merge_external_metadata, sce_list)
  active_assay <- sclet_get_active_assay(sce_list[[1]])
  if (!active_assay %in% asys) {
    active_assay <- NULL
  }
  active_layer <- sclet_get_active_layer(sce_list[[1]])

  
  combined_sce <- SingleCellExperiment(
    assays = cbn_assays,
    rowData = rowdata_combined,
    colData = coldata_combined,
    metadata = metadata_combined
  )
  SingleCellExperiment::colLabels(combined_sce) <- NULL
  
  combined_sce <- sclet_rebuild_internal_state(
    combined_sce,
    commands = list(),
    active_assay = active_assay
  )
  combined_sce <- sclet_set_hvg_state(
    combined_sce,
    nfeatures = sclet_get_hvg_nfeatures(sce_list[[1]]),
    method = get_hvg_method(sce_list[[1]]),
    hvgcols = colnames(combined_hvginfo)
  )
  if (!is.null(active_layer) && active_layer %in% Layers(combined_sce)) {
    combined_sce <- sclet_set_active_layer(combined_sce, active_layer)
  }
  
  combined_sce$batch <- factor(rep(seq_along(sce_list), sapply(sce_list, ncol)))
  combined_sce <- sclet_set_analysis_state(
    object = combined_sce,
    type = "integration",
    id = "merged_inputs",
    method = "sce_merge",
    inputs = list(
      objects = names(sce_list),
      assays = asys
    ),
    artifacts = list(
      batch_col = "batch",
      assays = asys
    ),
    summary = list(
      n_inputs = length(sce_list),
      n_cells = ncol(combined_sce)
    ),
    active = FALSE
  )
  combined_sce <- sclet_log_command(
    combined_sce,
    "sce_merge",
    params = list(
      n_inputs = length(sce_list),
      assays = asys
    ),
    outputs = list(
      state = "integration",
      integration = "merged_inputs"
    )
  )
  
  return(combined_sce)
}
