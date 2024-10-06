#' wrapper function for running cellchat on default parameter. 
#' 
#' https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#c-starting-from-a-singlecellexperiment-object
#' 
#' @param sce singlecellExperiment object.
#' @param group parameter to group data, must be one column in colData(obj).
#' @param assay_name the name of assay
#' @param species one of human, mouse
#' @param db_item one of  "all", "except Non-protein","Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact" and "Non-protein Signaling"
#' @param type one of  "triMean", "truncatedMean", "thresholdedMean", "median"
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed
#' @param min.cells min cells to filter
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @export 
runCellChat <- function(sce, group = NULL, 
                        assay_name = "logcounts",
                        species = "human",
                        db_item = c("Secreted Signaling"),
                        type = "triMean",
                        trim = 0.1,
                        min.cells = 10){
    
    species <- match.arg(species, c("human", "mouse", "zebrafish"))
    db_item <- match.arg(db_item, c("all", "except Non-protein",
                                    "Secreted Signaling", "ECM-Receptor", 
                                    "Cell-Cell Contact" , "Non-protein Signaling"))
    type <- match.arg(type, c("triMean", "truncatedMean", "thresholdedMean", "median"))

    if(! "logcounts" %in% assayNames(sce)){
        sce <- NormalizeData(sce)
    }

    if(is.null(group)){
        cat(">> Grouping data by default\t", format(Sys.time(), "%Y-%m-%d %X"), "\n")
        
        sce <- FindVariableFeatures(sce)
        sce <- ScaleData(sce)
        sce <- runPCA(sce, subset_row = VariableFeatures(sce), 
                      exprs_values = "scaled")
        sce <- FindNeighbors(sce, dims = 1:10)
        sce <- FindClusters(sce)
        group <- "label"
    }

    data <- assay(sce, assay_name)
    meta <- as.data.frame(colData(sce)) 

    cellchat_obj <- CellChat::createCellChat(object = data, 
                                             meta = meta, 
                                             group.by = group)
    
    

    if(species == "human"){
        db <- CellChat::CellChatDB.human
    }

    if(species == "mouse"){
        db <- CellChat::CellChatDB.mouse
    }

    if(species == "zebrafish"){
        db <- CellChat::CellChatDB.zebrafish
    }

    

    if(db_item %in% c("Secreted Signaling","ECM-Receptor", 
                      "Cell-Cell Contact")){
        db.use <- CellChat::subsetDB(db, search = db_item, key = "annotation")
    }

    if(db_item == "except Non-protein"){
        db.use <- CellChat::subsetDB(db)
    }

    if(db_item == "all"){
        db.use <- db
    }

    cellchat_obj@DB <- db.use
    cellchat_obj <- CellChat::subsetData(cellchat_obj)

    cellchat_obj <- CellChat::identifyOverExpressedGenes(cellchat_obj)
    cellchat_obj <- CellChat::identifyOverExpressedInteractions(cellchat_obj)

    
    cellchat_obj <- CellChat::computeCommunProb(cellchat_obj, type = type)
    cellchat_obj <- CellChat::filterCommunication(cellchat_obj, min.cells = min.cells)

    cellchat_obj <- CellChat::computeCommunProbPathway(cellchat_obj)
    cellchat_obj <- CellChat::aggregateNet(cellchat_obj)

    return(cellchat_obj)
}