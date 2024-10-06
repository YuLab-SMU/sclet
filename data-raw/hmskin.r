# make human skin in the form of singlecellExperiment object.
# raw data from https://figshare.com/articles/dataset/Example_data_for_running_CellChat_s_tutorial/13520015?file=25950872
# script https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html

library(SingleCellExperiment)

# load("./data_humanSkin_CellChat.rda")
data.input = data_humanSkin$data # normalized data matrix
meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
# Subset the input data for CelChat analysis
cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
data.input = data.input[, cell.use]
meta = meta[cell.use, ]

hmskin_sce <- SingleCellExperiment(assay = list(counts = data.input))
identical(colnames(hmskin_sce),rownames(meta))
colData(hmskin_sce) <- cbind(colData(hmskin_sce), meta)

usethis::use_data(hmskin_sce, overwrite = TRUE)
