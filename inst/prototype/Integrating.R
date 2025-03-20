Integrating <- function(sce_list, 
                        n_pcs = 30, 
                        check_preprocessing = TRUE, 
                        chosen_method = NULL) {
  # 1. 输入数据检查
  if (!all(sapply(sce_list, function(x) is(x, "SingleCellExperiment")))) {
    stop("所有输入对象必须为 SingleCellExperiment 类型。")
  }
  
  # 检查数据是否经过必要的预处理步骤
  if (check_preprocessing) {
    for (i in seq_along(sce_list)) {
      sce <- sce_list[[i]]
      if (!("counts" %in% assayNames(sce) && "logcounts" %in% assayNames(sce))) {
        stop(paste("输入数据", i, "未经过标准化，请确保数据已完成质控和标准化步骤。"))
      }
      if (is.null(rowData(sce)$highly_variable)) {
        stop(paste("输入数据", i, "未进行特征选择，请确保数据包含高度可变基因的信息。"))
      }
    }
  }
  
  # 2. 合并数据集
  # 找到基因的交集
  universe <- Reduce(intersect, lapply(sce_list, rownames))
  
  # 仅保留交集中的基因
  sce_list <- lapply(sce_list, function(sce) sce[universe, ])
  
  # 合并 counts 和 logcounts
  counts_combined <- do.call(cbind, lapply(sce_list, counts))
  logcounts_combined <- do.call(cbind, lapply(sce_list, logcounts))
  row_data_combined <- rowData(sce_list[[1]]) 
  col_data_combined <- do.call(rbind, lapply(sce_list, colData))
  
  # 创建合并后的 SingleCellExperiment 对象
  combined_sce <- SingleCellExperiment(
    assays = list(counts = counts_combined, logcounts = logcounts_combined),
    rowData = row_data_combined,
    colData = col_data_combined
  )
  # 添加批次信息
  combined_sce$batch <- factor(rep(seq_along(sce_list), sapply(sce_list, ncol)))
  # UMAP 可视化检查批次效应
  combined_sce <- scater::runPCA(combined_sce, exprs_values = "logcounts", ncomponents = n_pcs)
  combined_sce <- scater::runUMAP(combined_sce, dimred = "PCA")
  p1 <- plotUMAP(combined_sce, colour_by = "batch") + ggtitle("Before Batch Correction")
  
  print(p1)  # 确保绘图被显示
  
  if (is.null(chosen_method)) {
    stop("未指定去批次方法，请在调用时通过 chosen_method 参数指定 ('harmony', 'Combat', 'limma')。")
  }
  
  if (chosen_method == "fastMNN") {
    set.seed(1000101001)
    highly_variable <- rownames(sce_list)[rowData(combined_sce)$highly_variable]
    
    combined_corrected <- fastMNN(sce_list, d=50, k=20, subset.row=highly_variable,
                                  BSPARAM=BiocSingular::RandomParam(deferred=TRUE))}
  
  combined_corrected$batch <- factor(rep(seq_along(sce_list), sapply(sce_list, ncol)))
  combined_corrected <- scater::runPCA(combined_corrected, exprs_values = "reconstructed", ncomponents = n_pcs)
  combined_corrected <- scater::runUMAP(combined_corrected, dimred = "PCA")
  p2 <-plotUMAP(combined_corrected, colour_by = "batch") + ggtitle("After Batch Correction")
  print(p1+p2) 
  return(combined_corrected)
}




