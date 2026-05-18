#' list available layers
#'
#' @title Layers
#' @param object a SingleCellExperiment object
#' @return character vector of available layer names
#' @export
Layers <- function(object) {
    names(sclet_get_layer_registry(object))
}

#' access layer data
#'
#' @title LayerData
#' @param object a SingleCellExperiment object
#' @param layer layer name. If NULL, use `DefaultLayer(object)`.
#' @return assay matrix associated with the selected layer
#' @importFrom SummarizedExperiment assay
#' @export
LayerData <- function(object, layer = NULL) {
    resolved <- sclet_resolve_layer(object, layer = layer)
    if (is.null(resolved)) {
        stop("No valid layer found in object.")
    }
    SummarizedExperiment::assay(object, resolved$assay)
}

#' get or set the default layer
#'
#' @title DefaultLayer
#' @param object a SingleCellExperiment object
#' @param value layer name
#' @return current default layer or updated object
#' @export
DefaultLayer <- function(object) {
    sclet_get_active_layer(object)
}

#' @rdname DefaultLayer
#' @export
"DefaultLayer<-" <- function(object, value) {
    if (!value %in% Layers(object)) {
        stop("Layer '", value, "' not found in sclet layer registry.")
    }
    sclet_set_active_layer(object, value)
}
