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
    # An object with no active identity used to fall through to
    # `!group %in% colnames(meta)` with `group = NULL`, which is
    # `logical(0)`, and `if (logical(0))` reports base R's "argument is of
    # length zero" -- an error that names neither the object nor the missing
    # value. Hit by RunCellChat(), RunCellPhoneDB() and RunNicheNet(). Measured
    # 2026-10-22.
    if (is.null(group) || !nzchar(group)) {
        stop(
            "No cell grouping is available: the object has no active identity ",
            "and no `group` / `group_col` was supplied. Cluster the object ",
            "(FindClusters()) or name a colData column explicitly.",
            call. = FALSE
        )
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