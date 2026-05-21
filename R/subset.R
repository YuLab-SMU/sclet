#' @importFrom dplyr filter
#' @importFrom rlang enquo
#' @importFrom methods setMethod
setMethod("subset", "SingleCellExperiment",
    function(x, subset, select, ...) {
        subset <- rlang::enquo(arg = subset)
        cd <- as.data.frame(colData(x), row.names = colnames(x), optional = TRUE)
        cd$.sclet_cell <- colnames(x)
        y <- dplyr::filter(cd, !!subset)
        x[, y$.sclet_cell, drop = FALSE]
    }
)


#' @title subset_feature
#' @rdname subset-feature
#' @param x a SingleCellExperiment object
#' @param mincell feature should be detected with selected minimal cells, otherwise will be filter out
#' @param peek logical. If TRUE, print information of the feature and not performing subsetting
#' @return updated object if peek is FALSE
#' @importFrom SingleCellExperiment counts
#' @export
subset_feature <- function(x, mincell, peek = TRUE) {
    if (!inherits(x, "SingleCellExperiment")) {
        stop("Only SingleCellExperiment objects are supported.")
    }
    y <- counts(x)
    
    if (!peek) {
        return(x[Matrix::rowSums(y > 0) >= mincell, ])
    }

    n <- 0:mincell
    yy <- sapply(n, \(i) sum(Matrix::rowSums(y > 0) >= i))
    r <- yy/yy[1] * 100
    j <- round(r/2)
    for (i in seq_along(j)) {
        cat(rep("*", j[i]), sep="")
        lab <- sprintf(" >=%s | %s (%s%s)\n", 
                    n[i], yy[i], round(r[i]), '%')
        cat(lab)
    }
}

##' @title subset_cell
##' @rdname subset-cell
##' @param x a SingleCellExperiment object
##' @param feature selected features to filter out cells
##' @param method selected method to detect outliers to be filtered. One of 'mad', 'sd', and 'quantile'
##' @param n if method='quantile' and e.g. n=.98, then the lower and upper side of 1% data will be filter out; 
##'     if method is 'mad' or 'sd' (data will be transform by `log1p``), only data within the interval `median (mean) +/- n * mad (sd)` will be retained.
##' @return updated object
##' @export
subset_cell <- function(x, feature = "nFeature_RNA", method = 'mad', n = 3) {
    if (!inherits(x, "SingleCellExperiment")) {
        stop("Only SingleCellExperiment objects are supported.")
    }

    if (length(feature) == 1) {
        res <- subset_cell_internal(x, feature, method = method, n)
        return(x[, res$idx])
    }
    
    if (length(n) == 1) n <- rep(n, length(feature))
    res <- lapply(seq_along(feature), function(i) {
        subset_cell_internal(x, feature = feature[i], method = method, n = n[i])
    })

    criteria <- paste(sapply(res, function(y) sprintf("(%s)", y$criteria)), collapse = ' & ')

    for (i in seq_along(res)) {
        y <- res[[i]]
        if (i == 1) {
            idx <- y$idx
        } else {
            idx <- idx & y$idx
        }
    }
 
    msg <- sprintf("\nFilter out %s (%s/%s) percent of cells with:\n\t%s", 
            round(1-mean(idx), 3)*100,
            length(idx) - sum(idx),
            length(idx),
            criteria)
   
    cat(msg, "\n\n")
    x[, idx]
}

#' @importFrom stats mad median quantile sd
subset_cell_internal <- function(x, feature = "nFeature_RNA", method = 'mad', n = 3) {
    md <- colData(x)
    y <- md[[feature]]
    method <- tolower(method)
    if (method == 'quantile') {
        cutoff_l <- quantile(y, (1-n)/2)
        cutoff_r <- quantile(y, 1-(1-n)/2)
    } else {
        if (method == 'mad') {
            m <- mad(y)
        } else if (method == 'sd') {
            y <- log1p(y)
            m <- sd(y)
        }
        cutoff_l <- median(y) - n * m
        cutoff_r <- median(y) + n * m
    } 

    if (cutoff_l < 0) cutoff_l <- 0
    cutoff_l <- round(cutoff_l)
    cutoff_r <- round(cutoff_r)

    criteria <- sprintf("%s > %s & %s < %s",
                    feature,
                    cutoff_l,
                    feature,
                    cutoff_r
                )

    idx2 <- y > cutoff_r | y < cutoff_l
    msg <- sprintf("Filter out %s (%s/%s) percent of cells with:\n\t%s", 
            round(mean(idx2), 3)*100,
            sum(idx2),
            length(idx2),
            criteria)

    cat(msg, "\n")

    invisible(list(criteria = criteria, 
                   idx = y >= cutoff_l & y <= cutoff_r
            ))
}

