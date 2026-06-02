#' Access trajectory analysis results
#'
#' @title get_trajectory
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected trajectory record
#' @param id optional trajectory record id. If NULL, use the active trajectory record.
#' @return the full trajectory analysis record, a selected element, or `NULL`
#' @export
get_trajectory <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "trajectory")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "trajectory", id = id)
    }
    result <- sclet_normalize_analysis_record(
        type = "trajectory",
        id = if (is.null(id)) "trajectory" else id,
        record = result,
        analysis = sclet_get_analysis(object, "trajectory"),
        default_method = "slingshot"
    )
    extract_analysis_element(result, element)
}

#' Check whether trajectory results are available
#'
#' @title has_trajectory
#' @param object a SingleCellExperiment object
#' @param id optional trajectory record id
#' @return logical
#' @export
has_trajectory <- function(object, id = NULL) {
    !is.null(get_trajectory(object, id = id))
}

#' Access milo analysis results
#'
#' @title get_milo
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected milo record
#' @param id optional milo record id. If NULL, use the active milo record.
#' @return the full milo analysis record, a selected element, or `NULL`
#' @export
get_milo <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "milo")
    }
    result <- NULL
    if (!is.null(id)) {
        candidate <- sclet_get_state_record(object, "milo", id = id)
        if (!is.null(candidate) && (identical(candidate$method, "miloR") ||
            identical(candidate$artifacts$analysis_key, "milo"))) {
            result <- candidate
        }
    }
    result <- sclet_normalize_analysis_record(
        type = "milo",
        id = if (is.null(id)) "milo" else id,
        record = result,
        analysis = sclet_get_analysis(object, "milo"),
        default_method = "miloR"
    )
    extract_analysis_element(result, element)
}

#' Check whether milo results are available
#'
#' @title has_milo
#' @param object a SingleCellExperiment object
#' @param id optional milo record id
#' @return logical
#' @export
has_milo <- function(object, id = NULL) {
    !is.null(get_milo(object, id = id))
}

#' Access SuperCell analysis results
#'
#' @title get_supercell
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected SuperCell record
#' @param id optional aggregation record id. If NULL, use the active aggregation record.
#' @return the full SuperCell analysis record, a selected element, or `NULL`
#' @export
get_supercell <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "aggregation")
    }
    result <- NULL
    if (!is.null(id)) {
        candidate <- sclet_get_state_record(object, "aggregation", id = id)
        if (!is.null(candidate) && (identical(candidate$method, "SuperCell::SCimplify") ||
            identical(candidate$artifacts$analysis_key, "supercell"))) {
            result <- candidate
        }
    }
    result <- sclet_normalize_analysis_record(
        type = "aggregation",
        id = if (is.null(id)) "supercell" else id,
        record = result,
        analysis = sclet_get_analysis(object, "supercell"),
        default_method = "SuperCell::SCimplify"
    )
    extract_analysis_element(result, element)
}

#' Check whether SuperCell results are available
#'
#' @title has_supercell
#' @param object a SingleCellExperiment object
#' @param id optional aggregation record id
#' @return logical
#' @export
has_supercell <- function(object, id = NULL) {
    !is.null(get_supercell(object, id = id))
}

#' Access CellChat analysis results
#'
#' @title get_cellchat
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected CellChat record
#' @param id optional communication record id. If NULL, use the active communication record.
#' @return the full CellChat record, a selected element, or `NULL`
#' @export
get_cellchat <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "communication")
    }
    result <- NULL
    if (!is.null(id)) {
        candidate <- sclet_get_state_record(object, "communication", id = id)
        if (!is.null(candidate) && (identical(candidate$method, "CellChat") ||
            identical(candidate$artifacts$analysis_key, "cellchat"))) {
            result <- candidate
        }
    }
    result <- sclet_normalize_analysis_record(
        type = "communication",
        id = if (is.null(id)) "cellchat" else id,
        record = result,
        analysis = sclet_get_analysis(object, "cellchat"),
        default_method = "CellChat"
    )
    extract_analysis_element(result, element)
}

#' Check whether CellChat results are available
#'
#' @title has_cellchat
#' @param object a SingleCellExperiment object
#' @param id optional communication record id
#' @return logical
#' @export
has_cellchat <- function(object, id = NULL) {
    !is.null(get_cellchat(object, id = id))
}

#' Access CellPhoneDB analysis results
#'
#' @title get_cellphonedb
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected CellPhoneDB record
#' @param id optional communication record id. If NULL, use the active communication record.
#' @return the full CellPhoneDB record, a selected element, or `NULL`
#' @export
get_cellphonedb <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "communication")
    }
    result <- NULL
    if (!is.null(id)) {
        candidate <- sclet_get_state_record(object, "communication", id = id)
        if (!is.null(candidate) && (identical(candidate$method, "CellPhoneDB") ||
            identical(candidate$artifacts$analysis_key, "cellphonedb"))) {
            result <- candidate
        }
    }
    result <- sclet_normalize_analysis_record(
        type = "communication",
        id = if (is.null(id)) "cellphonedb" else id,
        record = result,
        analysis = sclet_get_analysis(object, "cellphonedb"),
        default_method = "CellPhoneDB"
    )
    extract_analysis_element(result, element)
}

#' Check whether CellPhoneDB results are available
#'
#' @title has_cellphonedb
#' @param object a SingleCellExperiment object
#' @param id optional communication record id
#' @return logical
#' @export
has_cellphonedb <- function(object, id = NULL) {
    !is.null(get_cellphonedb(object, id = id))
}

#' Access NicheNet analysis results
#'
#' @title get_nichenet
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected NicheNet record
#' @param id optional communication record id. If NULL, use the active communication record.
#' @return the full NicheNet record, a selected element, or `NULL`
#' @export
get_nichenet <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "communication")
    }
    result <- NULL
    if (!is.null(id)) {
        candidate <- sclet_get_state_record(object, "communication", id = id)
        if (!is.null(candidate) && (identical(candidate$method, "NicheNet") ||
            identical(candidate$artifacts$analysis_key, "nichenet"))) {
            result <- candidate
        }
    }
    result <- sclet_normalize_analysis_record(
        type = "communication",
        id = if (is.null(id)) "nichenet" else id,
        record = result,
        analysis = sclet_get_analysis(object, "nichenet"),
        default_method = "NicheNet"
    )
    extract_analysis_element(result, element)
}

#' Check whether NicheNet results are available
#'
#' @title has_nichenet
#' @param object a SingleCellExperiment object
#' @param id optional communication record id
#' @return logical
#' @export
has_nichenet <- function(object, id = NULL) {
    !is.null(get_nichenet(object, id = id))
}

#' Access standard Cell-Cell Communication (CCI) interactions
#'
#' @title get_cci
#' @param object a SingleCellExperiment object
#' @param id optional communication record id. If NULL, use the active communication record.
#' @return a data.frame containing standardized CCI interactions, or `NULL`
#' @export
get_cci <- function(object, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "communication")
    }
    if (is.null(id)) return(NULL)
    
    record <- sclet_get_state_record(object, "communication", id = id)
    if (is.null(record)) return(NULL)
    
    if (identical(record$method, "CellChat") || identical(record$artifacts$analysis_key, "cellchat")) {
        return(get_cellchat(object, element = "interactions", id = id))
    }
    
    if (identical(record$method, "CellPhoneDB") || identical(record$artifacts$analysis_key, "cellphonedb")) {
        return(get_cellphonedb(object, element = "interactions", id = id))
    }

    if (identical(record$method, "NicheNet") || identical(record$artifacts$analysis_key, "nichenet")) {
        return(get_nichenet(object, element = "interactions", id = id))
    }
    
    # Placeholder for other methods
    return(NULL)
}

#' Check whether standardized CCI results are available
#'
#' @title has_cci
#' @param object a SingleCellExperiment object
#' @param id optional communication record id
#' @return logical
#' @export
has_cci <- function(object, id = NULL) {
    !is.null(get_cci(object, id = id))
}

#' Access annotation state records
#'
#' @title get_annotation
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected annotation record
#' @param id optional annotation record id. If NULL, use the active annotation record.
#' @return the full annotation state record, a selected element, or `NULL`
#' @export
get_annotation <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "annotation")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "annotation", id = id)
    }
    if (is.null(result)) {
        cd <- SummarizedExperiment::colData(object)
        if ("SingleR_labels" %in% colnames(cd)) {
            result <- list(
                id = "singler",
                type = "annotation",
                method = "SingleR",
                artifacts = list(
                    labels_col = "SingleR_labels",
                    pruned_labels_col = if ("SingleR_pruned.labels" %in% colnames(cd)) "SingleR_pruned.labels" else NULL
                ),
                summary = list(
                    n_labels = length(unique(cd$SingleR_labels))
                )
            )
        }
    }
    extract_analysis_element(result, element)
}

#' Check whether annotation results are available
#'
#' @title has_annotation
#' @param object a SingleCellExperiment object
#' @param id optional annotation record id
#' @return logical
#' @export
has_annotation <- function(object, id = NULL) {
    !is.null(get_annotation(object, id = id))
}

#' Access mapping state records
#'
#' @title get_mapping
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected mapping record
#' @param id optional mapping record id. If NULL, use the active mapping record.
#' @return the full mapping state record, a selected element, or `NULL`
#' @export
get_mapping <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "mapping")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "mapping", id = id)
    }
    if (is.null(result)) {
        cd <- SummarizedExperiment::colData(object)
        if ("SingleR_labels" %in% colnames(cd)) {
            result <- list(
                id = "singler",
                type = "mapping",
                method = "SingleR",
                inputs = list(
                    mode = "label_transfer"
                ),
                artifacts = list(
                    labels_col = "SingleR_labels",
                    pruned_labels_col = if ("SingleR_pruned.labels" %in% colnames(cd)) "SingleR_pruned.labels" else NULL,
                    mapping_type = "reference_mapping"
                ),
                summary = list(
                    n_labels = length(unique(cd$SingleR_labels))
                )
            )
        }
    }
    extract_analysis_element(result, element)
}

#' Check whether mapping results are available
#'
#' @title has_mapping
#' @param object a SingleCellExperiment object
#' @param id optional mapping record id
#' @return logical
#' @export
has_mapping <- function(object, id = NULL) {
    !is.null(get_mapping(object, id = id))
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

sclet_normalize_analysis_record <- function(type, id, record = NULL, analysis = NULL, default_method = NULL) {
    if (is.null(record) && is.null(analysis)) {
        return(NULL)
    }

    result <- record
    if (is.null(result)) {
        result <- list()
    }

    if (!is.null(analysis)) {
        for (nm in names(analysis)) {
            if (is.null(result[[nm]])) {
                result[[nm]] <- analysis[[nm]]
            }
        }
    }

    if (is.null(result$id)) {
        result$id <- id
    }
    if (is.null(result$type)) {
        result$type <- type
    }
    if (is.null(result$method) && !is.null(default_method)) {
        result$method <- default_method
    }
    if (is.null(result$inputs)) {
        result$inputs <- list()
    }
    if (is.null(result$artifacts)) {
        result$artifacts <- list()
    }
    if (is.null(result$summary)) {
        result$summary <- list()
    }

    result
}

#' Access batch correction results
#'
#' @title get_batch
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the stored batch record
#' @return the full batch record, a selected element, or `NULL`
#' @export
get_batch <- function(object, element = NULL) {
    result <- sclet_normalize_analysis_record(
        type = "batch",
        id = "batchcorrect",
        analysis = sclet_get_analysis(object, "batch")
    )
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

#' Access integration state records
#'
#' @title get_integration
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected integration record
#' @param id optional integration record id. If NULL, use the active integration record.
#' @return the full integration state record, a selected element, or `NULL`
#' @export
#' @examples
#' # get_integration(sce)
#' # get_integration(sce, id = "merged_inputs")
get_integration <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "integration")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "integration", id = id)
    }
    if (is.null(result)) {
        batch <- get_batch(object)
        if (is.null(id) && !is.null(batch)) {
            result <- sclet_normalize_analysis_record(
                type = "integration",
                id = "batchcorrect",
                analysis = batch,
                default_method = batch$method
            )
            result$type <- "integration"
            result$artifacts <- utils::modifyList(
                list(layer = "corrected", analysis_key = "batch"),
                result$artifacts
            )
            if (is.null(result$summary$n_hvg) && !is.null(batch$hvg_n)) {
                result$summary$n_hvg <- batch$hvg_n
            }
        }
    }
    extract_analysis_element(result, element)
}

#' Check whether integration results are available
#'
#' @title has_integration
#' @param object a SingleCellExperiment object
#' @param id optional integration record id
#' @return logical
#' @export
has_integration <- function(object, id = NULL) {
    !is.null(get_integration(object, id = id))
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
        result <- sclet_get_legacy_graph_entry(object, name)
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


summarize_analysis_context_record <- function(record) {
    if (is.null(record)) {
        return(NULL)
    }
    keep <- intersect(c("id", "type", "method", "inputs", "artifacts", "summary"), names(record))
    record[keep]
}

#' Access the current analysis context
#'
#' @title get_analysis_context
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the context summary
#' @return a lightweight summary of the current active view and active analysis records
#' @export
get_analysis_context <- function(object, element = NULL) {
    records <- list(
        integration = summarize_analysis_context_record(get_integration(object)),
        annotation = summarize_analysis_context_record(get_annotation(object)),
        mapping = summarize_analysis_context_record(get_mapping(object)),
        trajectory = summarize_analysis_context_record(get_trajectory(object)),
        communication = summarize_analysis_context_record(get_cellchat(object)),
        milo = summarize_analysis_context_record(get_milo(object)),
        aggregation = summarize_analysis_context_record(get_supercell(object)),
        velocity = summarize_analysis_context_record(get_velocity(object)),
        scenic = summarize_analysis_context_record(get_scenic(object)),
        geneset_scoring = summarize_analysis_context_record(get_geneset_scoring(object)),
        cellrank = summarize_analysis_context_record(get_cellrank(object)),
        spatial = summarize_analysis_context_record(get_spatial(object)),
        perturbation = summarize_analysis_context_record(get_perturbation(object))
    )
    records <- records[!vapply(records, is.null, logical(1))]

    result <- list(
        active = list(
            assay = DefaultAssay(object),
            layer = DefaultLayer(object),
            reduction = DefaultReduction(object),
            graph = DefaultGraph(object),
            ident = ActiveIdent(object)
        ),
        records = records
    )

    extract_analysis_element(result, element)
}

#' Access RNA velocity analysis results
#'
#' @title get_velocity
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected velocity record
#' @param id optional velocity record id. If NULL, use the active velocity record.
#' @return the full velocity analysis record, a selected element, or `NULL`
#' @export
get_velocity <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "velocity")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "velocity", id = id)
    }
    analysis <- NULL
    if (is.null(result) || isTRUE(is.null(requested_id) || identical(id, sclet_get_active_state(object, "velocity")))) {
        analysis <- sclet_get_analysis(object, "velocity")
    }
    result <- sclet_normalize_analysis_record(
        type = "velocity",
        id = if (is.null(id)) "velocity" else id,
        record = result,
        analysis = analysis,
        default_method = "velociraptor::scvelo"
    )
    extract_analysis_element(result, element)
}

#' Check whether velocity results are available
#'
#' @title has_velocity
#' @param object a SingleCellExperiment object
#' @param id optional velocity record id
#' @return logical
#' @export
has_velocity <- function(object, id = NULL) {
    !is.null(get_velocity(object, id = id))
}

#' Access pySCENIC analysis results
#'
#' @title get_scenic
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected scenic record
#' @param id optional scenic record id. If NULL, use the active scenic record.
#' @return the full scenic analysis record, a selected element, or `NULL`
#' @export
get_scenic <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "scenic")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "scenic", id = id)
    }
    analysis <- NULL
    if (is.null(result) || isTRUE(is.null(requested_id) || identical(id, sclet_get_active_state(object, "scenic")))) {
        analysis <- sclet_get_analysis(object, "scenic")
    }
    result <- sclet_normalize_analysis_record(
        type = "scenic",
        id = if (is.null(id)) "scenic" else id,
        record = result,
        analysis = analysis,
        default_method = "pySCENIC"
    )
    extract_analysis_element(result, element)
}

#' Check whether SCENIC results are available
#'
#' @title has_scenic
#' @param object a SingleCellExperiment object
#' @param id optional scenic record id
#' @return logical
#' @export
has_scenic <- function(object, id = NULL) {
    !is.null(get_scenic(object, id = id))
}

#' Access gene set scoring analysis results
#'
#' @title get_geneset_scoring
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected geneset_scoring record
#' @param id optional geneset_scoring record id. If NULL, use the active geneset_scoring record.
#' @return the full geneset_scoring analysis record, a selected element, or `NULL`
#' @export
get_geneset_scoring <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "geneset_scoring")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "geneset_scoring", id = id)
    }
    analysis <- NULL
    if (is.null(result) || isTRUE(is.null(requested_id) || identical(id, sclet_get_active_state(object, "geneset_scoring")))) {
        analysis <- sclet_get_analysis(object, "geneset_scoring")
    }
    result <- sclet_normalize_analysis_record(
        type = "geneset_scoring",
        id = if (is.null(id)) "geneset_scoring" else id,
        record = result,
        analysis = analysis,
        default_method = if (!is.null(result$method)) result$method else NULL
    )
    extract_analysis_element(result, element)
}

#' Check whether gene set scoring results are available
#'
#' @title has_geneset_scoring
#' @param object a SingleCellExperiment object
#' @param id optional geneset_scoring record id
#' @return logical
#' @export
has_geneset_scoring <- function(object, id = NULL) {
    !is.null(get_geneset_scoring(object, id = id))
}

#' Access CellRank analysis results
#'
#' @title get_cellrank
#' @param object a SingleCellExperiment object
#' @param element optional element name
#' @param id optional cellrank record id
#' @return the full cellrank analysis record, a selected element, or `NULL`
#' @export
get_cellrank <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "trajectory")
    }
    result <- NULL
    if (!is.null(id)) {
        candidate <- sclet_get_state_record(object, "trajectory", id = id)
        if (!is.null(candidate) && (
            identical(candidate$method, "CellRank") ||
            identical(candidate$artifacts$analysis_key, "cellrank")
        )) {
            result <- candidate
        }
    }
    analysis <- NULL
    if (is.null(result) || isTRUE(is.null(requested_id) || identical(id, sclet_get_active_state(object, "trajectory")))) {
        analysis <- sclet_get_analysis(object, "cellrank")
    }
    result <- sclet_normalize_analysis_record(
        type = "trajectory",
        id = if (is.null(id)) "cellrank" else id,
        record = result,
        analysis = analysis,
        default_method = "CellRank"
    )
    extract_analysis_element(result, element)
}

#' Check whether CellRank results are available
#'
#' @title has_cellrank
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_cellrank <- function(object) {
    !is.null(get_cellrank(object))
}

#' Access Spatial Deconvolution analysis results
#'
#' @title get_spatial
#' @param object a SingleCellExperiment object
#' @param element optional element name
#' @param id optional spatial record id
#' @return the full spatial analysis record, a selected element, or `NULL`
#' @export
get_spatial <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "spatial")
    }
    result <- NULL
    if (!is.null(id)) {
        candidate <- sclet_get_state_record(object, "spatial", id = id)
        if (!is.null(candidate)) {
            result <- candidate
        } else if (id == "spatial_deconv") {
            result <- sclet_get_analysis(object, "spatial_deconv")
        } else if (id == "colocalization") {
            result <- sclet_get_analysis(object, "colocalization")
        } else if (id == "spatial_niche") {
            result <- sclet_get_analysis(object, "spatial_niche")
        }
    } else {
        result <- sclet_get_analysis(object, "spatial_deconv")
    }
    extract_analysis_element(result, element)
}

#' Check whether Spatial Deconvolution results are available
#'
#' @title has_spatial
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_spatial <- function(object) {
    !is.null(get_spatial(object))
}

#' Access spatial colocalization results
#'
#' @title get_colocalization
#' @param object a SingleCellExperiment object
#' @param element optional element name
#' @return the colocalization record, a selected element, or `NULL`
#' @export
get_colocalization <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "colocalization")
    extract_analysis_element(result, element)
}

#' Check whether colocalization results are available
#'
#' @title has_colocalization
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_colocalization <- function(object) {
    !is.null(get_colocalization(object))
}

#' Access spatial niche results
#'
#' @title get_spatial_niche
#' @param object a SingleCellExperiment object
#' @param element optional element name
#' @return the niche record, a selected element, or `NULL`
#' @export
get_spatial_niche <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "spatial_niche")
    extract_analysis_element(result, element)
}

#' Check whether spatial niche results are available
#'
#' @title has_spatial_niche
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_spatial_niche <- function(object) {
    !is.null(get_spatial_niche(object))
}

#' Access CellOracle Perturbation analysis results
#'
#' @title get_perturbation
#' @param object a SingleCellExperiment object
#' @param element optional element name
#' @param id optional perturbation record id
#' @return the full perturbation analysis record, a selected element, or `NULL`
#' @export
get_perturbation <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "perturbation")
    }
    result <- NULL
    if (!is.null(id) && id == "celloracle") {
        result <- sclet_get_analysis(object, "celloracle")
    } else if (is.null(id)) {
        result <- sclet_get_analysis(object, "celloracle")
    }
    extract_analysis_element(result, element)
}

#' Check whether CellOracle Perturbation results are available
#'
#' @title has_perturbation
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_perturbation <- function(object) {
    !is.null(get_perturbation(object))
}

#' Access perturbation priority analysis results (Augur)
#'
#' @title get_perturbation_priority
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract (e.g. `"AUC"`, `"result"`)
#' @return the full priority record, a selected element, or `NULL`
#' @export
get_perturbation_priority <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "augur")
    extract_analysis_element(result, element)
}

#' Check whether perturbation priority results are available
#'
#' @title has_perturbation_priority
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_perturbation_priority <- function(object) {
    !is.null(get_perturbation_priority(object))
}

#' Access rare cell detection results
#'
#' @title get_rare_cells
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract
#' @return the full rare cell record, a selected element, or `NULL`
#' @export
get_rare_cells <- function(object, element = NULL) {
    result <- sclet_get_analysis(object, "rare_cells")
    extract_analysis_element(result, element)
}

#' Check whether rare cell detection results are available
#'
#' @title has_rare_cells
#' @param object a SingleCellExperiment object
#' @return logical
#' @export
has_rare_cells <- function(object) {
    !is.null(get_rare_cells(object))
}

#' Access unified program activity data
#'
#' Resolves a program (signature, pathway, or regulon) from the gene set scoring
#' or SCENIC state records and returns the per-cell activity vector.
#'
#' @title get_program
#' @param object A SingleCellExperiment object.
#' @param program Character. Program name to resolve.
#' @param source One of `"auto"`, `"geneset_scoring"`, or `"scenic"`.
#' @param id Optional analysis record id.
#' @param assay Optional assay name when reading from an altExp.
#' @return A named numeric vector of program activity values.
#' @export
get_program <- function(object, program, source = c("auto", "geneset_scoring", "scenic"),
                        id = NULL, assay = NULL) {
    source <- match.arg(source)
    resolved <- sclet_resolve_program_activity(
        object, program = program, source = source, id = id, assay = assay
    )
    resolved$activity
}

#' Check whether a program result is available
#'
#' @title has_program
#' @param object A SingleCellExperiment object.
#' @param program Character. Program name to check.
#' @param source One of `"auto"`, `"geneset_scoring"`, or `"scenic"`.
#' @param id Optional analysis record id.
#' @return logical.
#' @export
has_program <- function(object, program, source = c("auto", "geneset_scoring", "scenic"),
                        id = NULL) {
    activity <- tryCatch(
        get_program(object, program = program, source = source, id = id),
        error = function(e) NULL
    )
    !is.null(activity) && length(activity) > 0
}
