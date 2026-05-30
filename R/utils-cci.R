#' Resolve Group Column
#' 
#' Internal helper to resolve the grouping column for analyses.
#' If `group` is NULL, it falls back to `ActiveIdent(object)`.
#' If the active identity is `colLabels`, it creates a temporary column in the metadata.
#' 
#' @param object A SingleCellExperiment object
#' @param meta A data.frame of metadata (usually `colData(object)`)
#' @param group The group column name or NULL
#' @return A list with `meta` (updated if necessary) and `group` (the resolved column name).
#' @noRd
sclet_resolve_group_column <- function(object, meta, group = NULL) {
    if (is.null(group)) {
        group <- ActiveIdent(object)
    }
    
    if (identical(group, "colLabels")) {
        if (!is.null(SingleCellExperiment::colLabels(object))) {
            meta$sclet_active_ident <- SingleCellExperiment::colLabels(object)
            group <- "sclet_active_ident"
        } else {
            stop("Active identity is 'colLabels', but colLabels(object) is NULL.")
        }
    } else if (!group %in% colnames(meta)) {
        stop(sprintf("Group column '%s' not found in colData.", group))
    }
    
    list(meta = meta, group = group)
}