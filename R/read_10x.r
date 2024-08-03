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

    # 这里还是全读进内存，然后转对象。
    # https://arrow.apache.org/cookbook/r/index.html
    # https://www.r-bloggers.com/2023/02/large-matrix-multiplication-duckdb-vs-sqlite/
    # 后续改进方案，应该是用arrow读，直接就是on disk的,
    # 然后对接duckdb,
    # 然后对接DelayedArray
    
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

class(xx)

subset.SCE <- function(x) {
    x$get_Count()[1:3, 1:3]
}

subset(xx)


