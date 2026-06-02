#' Check whether an ADT modality is present
#'
#' @param sce A SingleCellExperiment object.
#' @return logical. TRUE if `altExp(sce, "ADT")` exists.
#' @export
has_multimodal_adt <- function(sce) {
    "ADT" %in% SingleCellExperiment::altExpNames(sce)
}

#' Check whether an ATAC modality is present
#'
#' @param sce A SingleCellExperiment object.
#' @return logical. TRUE if `altExp(sce, "ATAC")` exists.
#' @export
has_multimodal_atac <- function(sce) {
    "ATAC" %in% SingleCellExperiment::altExpNames(sce)
}

#' Register ADT modality as a multimodal state record
#'
#' Records the presence and basic metadata of the ADT altExp in the
#' analysis-state layer so that downstream code can reason about
#' multimodal data without inspecting object internals directly.
#'
#' @param sce A SingleCellExperiment object with `altExp(sce, "ADT")`.
#' @param name Character. Multimodal record id. Defaults to `"ADT"`.
#'
#' @return Updated SingleCellExperiment with multimodal state record.
#' @export
register_multimodal_adt <- function(sce, name = "ADT") {
    stopifnot(is(sce, "SingleCellExperiment"))
    if (!has_multimodal_adt(sce)) {
        stop('altExp(sce, "ADT") not found. Load ADT data into the object first.')
    }
    adt <- SingleCellExperiment::altExp(sce, "ADT")
    adt_dims <- dim(adt)

    sce <- sclet_set_analysis(sce, "multimodal_adt", list(
        id = name,
        modality = "ADT",
        n_features = adt_dims[1],
        n_cells = adt_dims[2],
        timestamp = Sys.time()
    ))
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "multimodal",
        id = name,
        method = "altExp_ADT",
        inputs = list(
            altExp_name = "ADT"
        ),
        artifacts = list(
            analysis_key = "multimodal_adt",
            n_features = adt_dims[1],
            n_cells = adt_dims[2]
        ),
        summary = list(
            modality = "ADT"
        )
    )
    sce <- sclet_log_command(
        sce,
        "register_multimodal_adt",
        params = list(name = name),
        outputs = list(
            analysis = "multimodal_adt",
            state = "multimodal"
        )
    )

    sce
}

#' Register ATAC modality as a multimodal state record
#'
#' Records the presence and basic metadata of the ATAC altExp in the
#' analysis-state layer.
#'
#' @param sce A SingleCellExperiment object with `altExp(sce, "ATAC")`.
#' @param name Character. Multimodal record id. Defaults to `"ATAC"`.
#'
#' @return Updated SingleCellExperiment with multimodal state record.
#' @export
register_multimodal_atac <- function(sce, name = "ATAC") {
    stopifnot(is(sce, "SingleCellExperiment"))
    if (!has_multimodal_atac(sce)) {
        stop('altExp(sce, "ATAC") not found. Load ATAC data into the object first.')
    }
    atac <- SingleCellExperiment::altExp(sce, "ATAC")
    atac_dims <- dim(atac)

    sce <- sclet_set_analysis(sce, "multimodal_atac", list(
        id = name,
        modality = "ATAC",
        n_features = atac_dims[1],
        n_cells = atac_dims[2],
        timestamp = Sys.time()
    ))
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "multimodal",
        id = name,
        method = "altExp_ATAC",
        inputs = list(
            altExp_name = "ATAC"
        ),
        artifacts = list(
            analysis_key = "multimodal_atac",
            n_features = atac_dims[1],
            n_cells = atac_dims[2]
        ),
        summary = list(
            modality = "ATAC"
        )
    )
    sce <- sclet_log_command(
        sce,
        "register_multimodal_atac",
        params = list(name = name),
        outputs = list(
            analysis = "multimodal_atac",
            state = "multimodal"
        )
    )

    sce
}

#' Run a multimodal expansion workflow
#'
#' This is the reserved semantic shell over the multimodal mainline. It
#' currently detects and registers ADT and ATAC altExps but does not
#' perform modality-specific analysis. The point is to establish the
#' object contract now so that analysis backends can be plugged in later
#' without breaking the navigation model.
#'
#' @param sce A SingleCellExperiment object.
#' @param register_adt Logical. Whether to register ADT altExp if present.
#' @param register_atac Logical. Whether to register ATAC altExp if present.
#' @param name Workflow analysis id. Defaults to `"multimodal"`.
#'
#' @return Updated SingleCellExperiment with multimodal workflow record.
#' @export
RunMultimodalWorkflow <- function(sce, register_adt = TRUE,
    register_atac = TRUE, name = "multimodal") {

    component_ids <- list()

    if (isTRUE(register_adt) && has_multimodal_adt(sce)) {
        adt_id <- paste0(name, "_ADT")
        sce <- register_multimodal_adt(sce, name = adt_id)
        component_ids$ADT <- adt_id
    }

    if (isTRUE(register_atac) && has_multimodal_atac(sce)) {
        atac_id <- paste0(name, "_ATAC")
        sce <- register_multimodal_atac(sce, name = atac_id)
        component_ids$ATAC <- atac_id
    }

    if (!length(component_ids)) {
        message(
            "RunMultimodalWorkflow(): no altExp ADT or ATAC found. ",
            "Nothing was registered. Use altExp(sce, \"ADT\") <- adt_sce ",
            "to attach a modality first."
        )
        return(sce)
    }

    workflow_analysis <- list(
        id = name,
        method = "multimodal_mainline",
        inputs = list(),
        artifacts = component_ids,
        summary = list(
            has_ADT = !is.null(component_ids$ADT),
            has_ATAC = !is.null(component_ids$ATAC)
        ),
        created_at = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "multimodal_workflow", workflow_analysis)
    sce <- sclet_log_command(
        sce,
        "RunMultimodalWorkflow",
        params = list(
            register_adt = register_adt,
            register_atac = register_atac,
            name = name
        ),
        outputs = list(
            workflow = name,
            ADT = component_ids$ADT,
            ATAC = component_ids$ATAC
        )
    )

    sce
}
