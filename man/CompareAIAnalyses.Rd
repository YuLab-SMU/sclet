\name{CompareAIAnalyses}
\alias{CompareAIAnalyses}
\title{Compare recorded analysis routes without executing a route}
\usage{
CompareAIAnalyses(object, ids = NULL,
    criterion = c("evidence", "stability", "biological_preservation", "cost"))
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
\item{ids}{Optional analysis ids to compare.}
\item{criterion}{Comparison criterion.}
}
\value{
A read-only structured comparison. If no comparable integration records are
present, the result has \code{status = "not_available"}. No route is selected
automatically and no integration action is executed.
}
\description{
\code{CompareAIAnalyses()} compares already recorded integration or batch
analysis summaries, including metrics, trade-offs, fingerprints, and warnings.
It does not run Harmony, fastMNN, scVI, or any other modifying analysis.
}
