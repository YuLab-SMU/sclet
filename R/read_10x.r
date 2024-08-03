library(R6)


SCE <- R6Class("SCE", list(
    Count = NA,
    initialize = function(Count) {
        self$Count <- Count
    },
    get_Count = function() {
        self$Count
    }
))

read_10x <- function(dir, gene_column = 2, cell_column = 1) {
    ff <- list.files(dir)
    rf <- c("matrix.mtx", "barcodes.tsv", "genes.tsv")
    idx <- match(rf, ff)
    if (any(is.na(idx))) {
        mf <- paste(rf[which(is.na(idx))], collapse=", ")
        stop("required file(s) not found: ", mf)
    }

    rf <- file.path(dir, rf)    
    mat <- Matrix::readMM(rf[1])
    cell <- read.delim(rf[2], header=FALSE)
    feature <- read.delim(rf[3], header=FALSE)
    rownames(mat) <- feature[, gene_column]
    colnames(mat) <- cell[, cell_column]
    requireNamespace("HDF5Array", quietly = TRUE)
    Count <- as(mat, "HDF5Matrix")
    SCE$new(Count)
}




x <- read_10x('hg19')
obj_size(xx)

# seurat
x2 <- Read10X('hg19')
obj_size(yy)

y = xx$get_Count()

y[1:3, 1:3]

ncount <- colSums(y)
ncount2 <- colSums(x2)

identical(ncount, ncount2)
