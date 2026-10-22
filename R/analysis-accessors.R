#' Access trajectory analysis results
#'
#' @title get_trajectory
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected trajectory record
#' @param id optional trajectory record id. If NULL, use the active trajectory record.
#' @return the full trajectory analysis record, a selected element, or `NULL`
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "trajectory", method = "slingshot")
#' state <- list(states = list(active = list(trajectory = "run1"),
#'     records = list(trajectory = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_trajectory(sce)$method
get_trajectory <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "trajectory")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "trajectory", id = id)
    }
    analysis <- sclet_legacy_analysis(object, "trajectory",
        "trajectory", requested_id, id)
    result <- sclet_normalize_analysis_record(
        type = "trajectory",
        id = if (is.null(id)) "trajectory" else id,
        record = result,
        analysis = analysis,
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "trajectory", method = "slingshot")
#' state <- list(states = list(active = list(trajectory = "run1"),
#'     records = list(trajectory = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_trajectory(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "milo", method = "miloR")
#' state <- list(states = list(active = list(milo = "run1"),
#'     records = list(milo = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_milo(sce)$method
get_milo <- function(object, element = NULL, id = NULL) {
    requested_id <- id
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
    analysis <- sclet_legacy_analysis(object, "milo",
        "milo", requested_id, id)
    result <- sclet_normalize_analysis_record(
        type = "milo",
        id = if (is.null(id)) "milo" else id,
        record = result,
        analysis = analysis,
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "milo", method = "miloR")
#' state <- list(states = list(active = list(milo = "run1"),
#'     records = list(milo = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_milo(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "aggregation", method = "SuperCell::SCimplify")
#' state <- list(states = list(active = list(aggregation = "run1"),
#'     records = list(aggregation = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_supercell(sce)$method
get_supercell <- function(object, element = NULL, id = NULL) {
    requested_id <- id
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
    analysis <- sclet_legacy_analysis(object, "supercell",
        "supercell", requested_id, id)
    result <- sclet_normalize_analysis_record(
        type = "aggregation",
        id = if (is.null(id)) "supercell" else id,
        record = result,
        analysis = analysis,
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "aggregation", method = "SuperCell::SCimplify")
#' state <- list(states = list(active = list(aggregation = "run1"),
#'     records = list(aggregation = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_supercell(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "communication", method = "CellChat")
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_cellchat(sce)$method
get_cellchat <- function(object, element = NULL, id = NULL) {
    requested_id <- id
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
    analysis <- sclet_legacy_analysis(object, "cellchat",
        "cellchat", requested_id, id)
    result <- sclet_normalize_analysis_record(
        type = "communication",
        id = if (is.null(id)) "cellchat" else id,
        record = result,
        analysis = analysis,
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "communication", method = "CellChat")
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_cellchat(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "communication", method = "CellPhoneDB")
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_cellphonedb(sce)$method
get_cellphonedb <- function(object, element = NULL, id = NULL) {
    requested_id <- id
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
    analysis <- sclet_legacy_analysis(object, "cellphonedb",
        "cellphonedb", requested_id, id)
    result <- sclet_normalize_analysis_record(
        type = "communication",
        id = if (is.null(id)) "cellphonedb" else id,
        record = result,
        analysis = analysis,
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "communication", method = "CellPhoneDB")
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_cellphonedb(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "communication", method = "NicheNet")
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_nichenet(sce)$method
get_nichenet <- function(object, element = NULL, id = NULL) {
    requested_id <- id
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
    analysis <- sclet_legacy_analysis(object, "nichenet",
        "nichenet", requested_id, id)
    result <- sclet_normalize_analysis_record(
        type = "communication",
        id = if (is.null(id)) "nichenet" else id,
        record = result,
        analysis = analysis,
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "communication", method = "NicheNet")
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_nichenet(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' interactions <- data.frame(source = "alpha", target = "beta",
#'     ligand = "L1", receptor = "R1")
#' record <- list(id = "run1", type = "communication", method = "CellChat",
#'     interactions = interactions)
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' nrow(get_cci(sce))
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' interactions <- data.frame(source = "alpha", target = "beta",
#'     ligand = "L1", receptor = "R1")
#' record <- list(id = "run1", type = "communication", method = "CellChat",
#'     interactions = interactions)
#' state <- list(states = list(active = list(communication = "run1"),
#'     records = list(communication = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_cci(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$SingleR_labels <- rep(c("alpha", "beta"), length.out = 20)
#' get_annotation(sce)$summary$n_labels
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

#' Access malignant-cell call records
#'
#' @title get_malignancy
#' @param object a SingleCellExperiment object
#' @param element optional element name to extract from the selected record
#' @param id optional record id. If NULL, the most recent one.
#' @return the full record, a selected element, or `NULL`
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "malignant", method = "four-state malignant-cell call")
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(malignancy = list(malignant = record)))
#' get_malignancy(sce)$method
get_malignancy <- function(object, element = NULL, id = NULL) {
    records <- sclet_get_analysis(object, "malignancy")
    if (is.null(records)) {
        return(NULL)
    }
    result <- if (is.null(id)) {
        records[[length(records)]]
    } else {
        records[[id]]
    }
    extract_analysis_element(result, element)
}

#' Check whether a malignant-cell call is available
#'
#' @title has_malignancy
#' @param object a SingleCellExperiment object
#' @param id optional record id
#' @return logical
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' has_malignancy(sce)
has_malignancy <- function(object, id = NULL) {
    !is.null(get_malignancy(object, id = id))
}

#' Access annotation state records
#'
#' @title has_annotation
#' @param object a SingleCellExperiment object
#' @param id optional annotation record id
#' @return logical
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$SingleR_labels <- rep(c("alpha", "beta"), length.out = 20)
#' has_annotation(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$SingleR_labels <- rep(c("alpha", "beta"), length.out = 20)
#' get_mapping(sce)$artifacts$mapping_type
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$SingleR_labels <- rep(c("alpha", "beta"), length.out = 20)
#' has_mapping(sce)
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

# The legacy `analyses` record to stand in for the one the caller asked for.
#
# Returns the legacy record when standing in is honest, and NULL when it is not.
# Standing in is honest when:
#   * the caller named no id (they want "the" record for this analysis), or
#   * the id they named is the active one, or
#   * the id they named IS the legacy record's own id -- the legacy slot holds a
#     single record and does not index by id, so a caller who read the id off
#     that record and passed it back must still get it.
#
# It is NOT honest when they named some other id and no record by that id
# exists: the answer is "there is no such record", and returning a different one
# is worse than returning NULL, because `has_<type>(sce, id = "typo")` then
# reports TRUE and the caller reads a record it never asked for. Measured
# 2026-10-22: `get_velocity(sce, id = "typo")` returned the record whose id is
# "velocity", and `has_velocity(sce, id = "typo")` returned TRUE.
sclet_legacy_analysis <- function(object, state_type, legacy_key, requested_id, resolved_id) {
    analysis <- tryCatch(sclet_get_analysis(object, legacy_key), error = function(e) NULL)
    if (is.null(analysis)) {
        return(NULL)
    }
    if (is.null(requested_id)) {
        return(analysis)
    }
    active <- tryCatch(sclet_get_active_state(object, state_type), error = function(e) NULL)
    if (!is.null(active) && identical(as.character(resolved_id), as.character(active))) {
        return(analysis)
    }
    legacy_id <- tryCatch(analysis$id, error = function(e) NULL)
    if (!is.null(legacy_id) && length(legacy_id) == 1L &&
        identical(as.character(legacy_id), as.character(requested_id))) {
        return(analysis)
    }
    NULL
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "batchcorrect", type = "batch", method = "BatchRemover")
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(batch = record))
#' get_batch(sce)$method
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "batchcorrect", type = "batch", method = "BatchRemover")
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(batch = record))
#' has_batch(sce)
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
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "batchcorrect", type = "batch", method = "BatchRemover",
#'     artifacts = list(layer = "corrected"))
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(batch = record))
#' get_integration(sce)$artifacts$layer
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "batchcorrect", type = "batch", method = "BatchRemover",
#'     artifacts = list(layer = "corrected"))
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(batch = record))
#' has_integration(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- NormalizeData(sce)
#' sce <- FindVariableFeatures(sce, nfeatures = 20)
#' get_hvg(sce)$n
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- NormalizeData(sce)
#' sce <- FindVariableFeatures(sce, nfeatures = 20)
#' has_hvg(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- RunPCA(sce, ncomponents = 5)
#' sce <- FindNeighbors(sce, dims = 1:5)
#' get_graph(sce)$name
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- RunPCA(sce, ncomponents = 5)
#' sce <- FindNeighbors(sce, dims = 1:5)
#' has_graph(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$SingleR_labels <- rep(c("alpha", "beta"), length.out = 20)
#' names(get_analysis_context(sce)$records)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "velocity", method = "velociraptor::scvelo")
#' state <- list(states = list(active = list(velocity = "run1"),
#'     records = list(velocity = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_velocity(sce)$method
get_velocity <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "velocity")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "velocity", id = id)
    }
    analysis <- sclet_legacy_analysis(object, "velocity",
        "velocity", requested_id, id)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "velocity", method = "velociraptor::scvelo")
#' state <- list(states = list(active = list(velocity = "run1"),
#'     records = list(velocity = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_velocity(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "scenic", method = "pySCENIC")
#' state <- list(states = list(active = list(scenic = "run1"),
#'     records = list(scenic = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_scenic(sce)$method
get_scenic <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "scenic")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "scenic", id = id)
    }
    analysis <- sclet_legacy_analysis(object, "scenic",
        "scenic", requested_id, id)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "scenic", method = "pySCENIC")
#' state <- list(states = list(active = list(scenic = "run1"),
#'     records = list(scenic = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_scenic(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$Score_sigA <- seq_len(20)
#' record <- list(id = "Score", type = "geneset_scoring", method = "UCell",
#'     artifacts = list(storage = "colData", score_columns = "Score_sigA"))
#' state <- list(states = list(active = list(geneset_scoring = "Score"),
#'     records = list(geneset_scoring = list(Score = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_geneset_scoring(sce)$method
get_geneset_scoring <- function(object, element = NULL, id = NULL) {
    requested_id <- id
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "geneset_scoring")
    }
    result <- NULL
    if (!is.null(id)) {
        result <- sclet_get_state_record(object, "geneset_scoring", id = id)
    }
    analysis <- sclet_legacy_analysis(object, "geneset_scoring",
        "geneset_scoring", requested_id, id)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$Score_sigA <- seq_len(20)
#' record <- list(id = "Score", type = "geneset_scoring", method = "UCell",
#'     artifacts = list(storage = "colData", score_columns = "Score_sigA"))
#' state <- list(states = list(active = list(geneset_scoring = "Score"),
#'     records = list(geneset_scoring = list(Score = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_geneset_scoring(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "trajectory", method = "CellRank")
#' state <- list(states = list(active = list(trajectory = "run1"),
#'     records = list(trajectory = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_cellrank(sce)$method
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
    analysis <- sclet_legacy_analysis(object, "trajectory",
        "cellrank", requested_id, id)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "trajectory", method = "CellRank")
#' state <- list(states = list(active = list(trajectory = "run1"),
#'     records = list(trajectory = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_cellrank(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "spatial", method = "cell2location")
#' state <- list(states = list(active = list(spatial = "run1"),
#'     records = list(spatial = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_spatial(sce)$method
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "run1", type = "spatial", method = "cell2location")
#' state <- list(states = list(active = list(spatial = "run1"),
#'     records = list(spatial = list(run1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_spatial(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "colocalization", type = "spatial", method = "SVP")
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(colocalization = record))
#' get_colocalization(sce)$method
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "colocalization", type = "spatial", method = "SVP")
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(colocalization = record))
#' has_colocalization(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "spatial_niche", type = "spatial", method = "SVP")
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(spatial_niche = record))
#' get_spatial_niche(sce)$method
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "spatial_niche", type = "spatial", method = "SVP")
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(spatial_niche = record))
#' has_spatial_niche(sce)
has_spatial_niche <- function(object) {
    !is.null(get_spatial_niche(object))
}

#' Access CellOracle Perturbation analysis results
#'
#' @title get_perturbation
#' @param object a SingleCellExperiment object
#' @param element optional element name
#' @param id optional perturbation record id (e.g. `"celloracle_<gene>"`)
#' @return the full perturbation analysis record, a selected element, or `NULL`
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "celloracle_GENE1", type = "perturbation", method = "CellOracle")
#' state <- list(states = list(active = list(perturbation = "celloracle_GENE1"),
#'     records = list(perturbation = list(celloracle_GENE1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' get_perturbation(sce)$method
get_perturbation <- function(object, element = NULL, id = NULL) {
    if (is.null(id)) {
        id <- sclet_get_active_state(object, "perturbation")
    }
    result <- sclet_get_state_record(object, "perturbation", id)
    if (is.null(result)) {
        # Fallback for objects created before the state contract
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "celloracle_GENE1", type = "perturbation", method = "CellOracle")
#' state <- list(states = list(active = list(perturbation = "celloracle_GENE1"),
#'     records = list(perturbation = list(celloracle_GENE1 = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_perturbation(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "augur", type = "priority", method = "Augur",
#'     AUC = data.frame(cell_type = c("alpha", "beta"), auc = c(0.7, 0.4)))
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(augur = record))
#' get_perturbation_priority(sce)$AUC$auc
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "augur", type = "priority", method = "Augur",
#'     AUC = data.frame(cell_type = c("alpha", "beta"), auc = c(0.7, 0.4)))
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(augur = record))
#' has_perturbation_priority(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "rareq", type = "rare_cells", method = "density",
#'     artifacts = list(label_col = "rare_cluster"),
#'     summary = list(n_clusters = 3, n_rare_clusters = 1))
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(rare_cells = record))
#' get_rare_cells(sce)$summary$n_clusters
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' record <- list(id = "rareq", type = "rare_cells", method = "density",
#'     artifacts = list(label_col = "rare_cluster"),
#'     summary = list(n_clusters = 3, n_rare_clusters = 1))
#' S4Vectors::metadata(sce)$sclet <- list(analyses = list(rare_cells = record))
#' has_rare_cells(sce)
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$Score_sigA <- seq_len(20)
#' record <- list(id = "Score", type = "geneset_scoring", method = "UCell",
#'     artifacts = list(storage = "colData", score_columns = "Score_sigA"))
#' state <- list(states = list(active = list(geneset_scoring = "Score"),
#'     records = list(geneset_scoring = list(Score = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' head(get_program(sce, "sigA"))
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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' SummarizedExperiment::colData(sce)$Score_sigA <- seq_len(20)
#' record <- list(id = "Score", type = "geneset_scoring", method = "UCell",
#'     artifacts = list(storage = "colData", score_columns = "Score_sigA"))
#' state <- list(states = list(active = list(geneset_scoring = "Score"),
#'     records = list(geneset_scoring = list(Score = record))))
#' S4Vectors::metadata(sce)$sclet <- state
#' has_program(sce, "sigA")
has_program <- function(object, program, source = c("auto", "geneset_scoring", "scenic"),
                        id = NULL) {
    activity <- tryCatch(
        get_program(object, program = program, source = source, id = id),
        error = function(e) NULL
    )
    !is.null(activity) && length(activity) > 0
}
