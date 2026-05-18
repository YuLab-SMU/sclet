#' wrapper function for running cellchat on default parameter.
#' 
#' https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#c-starting-from-a-singlecellexperiment-object
#' 
#' @title RunCellChat
#' @param sce singlecellExperiment object.
#' @param group parameter to group data, must be one column in colData(obj).
#' @param assay_name the name of assay. If NULL, sclet resolves it from the
#' selected `layer`.
#' @param layer layer used for communication analysis. If NULL, use
#' `DefaultLayer(sce)`.
#' @param species one of human, mouse
#' @param db_item one of  "all", "except Non-protein","Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact" and "Non-protein Signaling"
#' @param type one of  "triMean", "truncatedMean", "thresholdedMean", "median"
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed
#' @param min.cells min cells to filter
#' @param name communication record id. Defaults to `"cellchat"`.
#' @param return one of `"cellchat"`, `"sce"` or `"both"`
#' @return by default a CellChat object; if `return = "sce"`, returns the
#' updated SingleCellExperiment; if `return = "both"`, returns a list with both
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @export
RunCellChat <- function(sce, group = "label",
                        assay_name = NULL,
                        layer = NULL,
                        species = "human",
                        db_item = c("Secreted Signaling"),
                        type = "triMean",
                        trim = 0.1,
                        min.cells = 10,
                        name = "cellchat",
                        return = c("cellchat", "sce", "both")) {

    if (!requireNamespace("CellChat", quietly = TRUE)) {
        stop("Package 'CellChat' is needed for this function to work. Please install it.")
    }

    species <- match.arg(species, c("human", "mouse", "zebrafish"))
    db_item <- match.arg(db_item, c("all", "except Non-protein",
                                    "Secreted Signaling", "ECM-Receptor", 
                                    "Cell-Cell Contact" , "Non-protein Signaling"))
    return <- match.arg(return)
    type <- match.arg(type, c("triMean", "truncatedMean", "thresholdedMean", "median"))

    if(! "logcounts" %in% assayNames(sce)){
        sce <- NormalizeData(sce)
    }

    if( !group %in% colnames(colData(sce))){
        stop("group must be in the colData of sce")
    }

    source <- sclet_resolve_expression_source(
        object = sce,
        layer = layer,
        assay = assay_name,
        prefer_nonscaled = TRUE,
        context = "CellChat analysis"
    )
    assay_name <- source$assay
    data <- assay(sce, assay_name)
    meta <- as.data.frame(colData(sce)) 

    cellchat_obj <- CellChat::createCellChat(object = data, 
                                             meta = meta, 
                                             group.by = group)


    db <- switch(species,
        human = CellChat::CellChatDB.human,
        mouse = CellChat::CellChatDB.mouse,
        zebrafish = CellChat::CellChatDB.zebrafish
    )

    if (db_item == "all"){
        db.use <- db
    } else if(db_item == "except Non-protein"){
        db.use <- CellChat::subsetDB(db)
    } else {
        db.use <- CellChat::subsetDB(db, search = db_item, key = "annotation")
    }


    cellchat_obj@DB <- db.use
    cellchat_obj <- CellChat::subsetData(cellchat_obj)

    cellchat_obj <- CellChat::identifyOverExpressedGenes(cellchat_obj)
    cellchat_obj <- CellChat::identifyOverExpressedInteractions(cellchat_obj)


    cellchat_obj <- CellChat::computeCommunProb(cellchat_obj, type = type)
    cellchat_obj <- CellChat::filterCommunication(cellchat_obj, min.cells = min.cells)

    cellchat_obj <- CellChat::computeCommunProbPathway(cellchat_obj)
    cellchat_obj <- CellChat::aggregateNet(cellchat_obj)

    sce <- sclet_set_analysis(
        sce,
        "cellchat",
        list(
            id = name,
            method = "CellChat",
            object = cellchat_obj,
            group = group,
            assay = assay_name,
            layer = source$layer,
            species = species,
            db_item = db_item,
            params = list(
                type = type,
                trim = trim,
                min.cells = min.cells
            )
        )
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "communication",
        id = name,
        method = "CellChat",
        inputs = list(
            assay = assay_name,
            layer = source$layer,
            group = group
        ),
        artifacts = list(
            analysis_key = "cellchat"
        ),
        params = list(
            species = species,
            db_item = db_item,
            type = type,
            trim = trim,
            min.cells = min.cells
        ),
        summary = list(
            n_groups = length(unique(meta[[group]]))
        ),
        active = TRUE
    )
    sce <- sclet_log_command(
        sce,
        "RunCellChat",
        params = list(
            group = group,
            assay_name = assay_name,
            layer = source$layer,
            species = species,
            db_item = db_item,
            type = type,
            trim = trim,
            min.cells = min.cells,
            name = name
        ),
        outputs = list(
            analysis = "communication",
            communication = name
        )
    )

    if (identical(return, "sce")) {
        return(sce)
    }
    if (identical(return, "both")) {
        return(list(
            sce = sce,
            cellchat = cellchat_obj
        ))
    }

    attr(cellchat_obj, "sce") <- sce
    return(cellchat_obj)
}

#' @rdname RunCellChat
#' @export
runCellChat <- function(sce, group = "label",
                        assay_name = NULL,
                        layer = NULL,
                        species = "human",
                        db_item = c("Secreted Signaling"),
                        type = "triMean",
                        trim = 0.1,
                        min.cells = 10,
                        name = "cellchat",
                        return = c("cellchat", "sce", "both")) {
    RunCellChat(
        sce = sce,
        group = group,
        assay_name = assay_name,
        layer = layer,
        species = species,
        db_item = db_item,
        type = type,
        trim = trim,
        min.cells = min.cells,
        name = name,
        return = return
    )
}
