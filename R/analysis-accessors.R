#' Access trajectory analysis results
#'
#' @title get_trajectory
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the stored trajectory record
#' @return the full trajectory analysis record, a selected element, or `NULL`
#' @export
get_trajectory <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "trajectory")
    if (is.null(result)) {
        legacy <- S4Vectors::metadata(object)$slingshot_info
        if (!is.null(legacy)) {
            result <- list(
                method = "slingshot",
                dataset = legacy
            )
        }
    }
    extract_analysis_element(result, element)
}

#' Check whether trajectory results are available
#'
#' @title has_trajectory
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_trajectory <- function(object) {
    !is.null(get_trajectory(object))
}

#' Access milo analysis results
#'
#' @title get_milo
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the stored milo record
#' @return the full milo analysis record, a selected element, or `NULL`
#' @export
get_milo <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "milo")
    if (is.null(result)) {
        result <- S4Vectors::metadata(object)$milo_results
    }
    extract_analysis_element(result, element)
}

#' Check whether milo results are available
#'
#' @title has_milo
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_milo <- function(object) {
    !is.null(get_milo(object))
}

#' Access SuperCell analysis results
#'
#' @title get_supercell
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the stored SuperCell record
#' @return the full SuperCell analysis record, a selected element, or `NULL`
#' @export
get_supercell <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "supercell")
    if (is.null(result)) {
        legacy <- S4Vectors::metadata(object)$SuperCell
        if (!is.null(legacy)) {
            result <- list(
                object = legacy
            )
        }
    }
    extract_analysis_element(result, element)
}

#' Check whether SuperCell results are available
#'
#' @title has_supercell
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_supercell <- function(object) {
    !is.null(get_supercell(object))
}

#' Access CellChat analysis results
#'
#' @title get_cellchat
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the stored CellChat record
#' @return the full CellChat record, a selected element, or `NULL`
#' @export
get_cellchat <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "cellchat")
    extract_analysis_element(result, element)
}

#' Check whether CellChat results are available
#'
#' @title has_cellchat
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_cellchat <- function(object) {
    !is.null(get_cellchat(object))
}

extract_analysis_element <- function(result, element = NULL) {
    if (is.null(result)) {
        return(NULL)
    }
    if (is.null(element)) {
        return(result)
    }
    result[[element]]
}

#' Access batch correction results
#'
#' @title get_batch
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the stored batch record
#' @return the full batch record, a selected element, or `NULL`
#' @export
get_batch <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "batch")
    if (is.null(result)) {
        result <- S4Vectors::metadata(object)$batch_correction
    }
    extract_analysis_element(result, element)
}

#' Check whether batch correction results are available
#'
#' @title has_batch
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_batch <- function(object) {
    !is.null(get_batch(object))
}

#' Access highly variable feature metadata
#'
#' @title get_hvg
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the stored HVG record
#' @return the full HVG record, a selected element, or `NULL`
#' @export
get_hvg <- function(object, element = NULL) {
    rd <- SummarizedExperiment::rowData(object)
    result <- list(
        n = sclet_get_hvg_nfeatures(object),
        method = sclet_get_hvg_method(object),
        rowData_cols = sclet_get_hvg_cols(object)
    )
    if (!is.null(result$rowData_cols)) {
        result$rowData <- rd[, result$rowData_cols, drop = FALSE]
    }
    state <- sclet_get_state(object)
    if (!is.null(state$features$hvg$selected)) {
        result$selected <- state$features$hvg$selected
    }
    if (all(vapply(result, is.null, logical(1)))) {
        return(NULL)
    }
    extract_analysis_element(result, element)
}

#' Check whether HVG metadata are available
#'
#' @title has_hvg
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_hvg <- function(object) {
    !is.null(get_hvg(object))
}

#' Access graph registry entries
#'
#' @title get_graph
#' @param object a SingleCellExperiment object
#' @param name graph name; if NULL, use the active graph
#' @param element optional element name to extract from the stored graph record
#' @return the full graph record, a selected element, or `NULL`
#' @export
get_graph <- function(object, name = NULL, element = NULL) {
    state <- sclet_get_state(object)
    if (is.null(name)) {
        name <- sclet_get_active_graph(object)
    }
    result <- NULL
    if (!is.null(name)) {
        result <- state$graphs[[name]]
    }
    if (is.null(result) && identical(name, "knn_graph")) {
        legacy <- S4Vectors::metadata(object)$knn_graph
        if (!is.null(legacy)) {
            result <- list(object = legacy)
        }
    }
    if (!is.null(result) && is.null(result$name)) {
        result$name <- name
    }
    extract_analysis_element(result, element)
}

#' Check whether graph results are available
#'
#' @title has_graph
#' @param object a SingleCellExperiment object
#' @param name graph name; if NULL, use the active graph
#' @return logical
#' @export
has_graph <- function(object, name = NULL) {
    !is.null(get_graph(object, name = name))
}
