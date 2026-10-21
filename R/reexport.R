#' @importFrom aplot plot_list
#' @export
aplot::plot_list

#' @importFrom SummarizedExperiment colData
#' @export
SummarizedExperiment::colData

#' @importFrom SummarizedExperiment rowData
#' @export
SummarizedExperiment::rowData

#' @importFrom SummarizedExperiment 'colData<-'
#' @export
SummarizedExperiment::'colData<-'

#' @importFrom SummarizedExperiment 'rowData<-'
#' @export
SummarizedExperiment::'rowData<-'

#' @importFrom SummarizedExperiment assay
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' dim(assay(sce))
SummarizedExperiment::assay

#' @importFrom scater plotColData
#' @export
scater::plotColData

#' @importFrom scuttle logNormCounts
#' @export
scuttle::logNormCounts

