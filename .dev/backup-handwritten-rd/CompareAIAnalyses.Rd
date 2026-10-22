\name{CompareAIAnalyses}
\alias{CompareAIAnalyses}
\title{Compare recorded analysis routes without executing a route}
\usage{
CompareAIAnalyses(
  object,
  ids = NULL,
  criterion = c("evidence", "stability", "biological_preservation", "cost")
)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{ids}{Optional analysis ids to compare.}

\item{criterion}{Comparison criterion used to order routes. Supported values
are \code{"evidence"}, \code{"stability"},
\code{"biological_preservation"}, and \code{"cost"}.}
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
\details{
Routes are ordered according to \code{criterion} when the required metric is
available. For example, \code{"biological_preservation"} orders routes by the
recorded preservation metric, while \code{"cost"} orders by runtime. If a
metric is missing or explicitly unavailable, the comparison normalizes it to
\code{value = NA_real_} with \code{uncertainty$status = "not_available"}
instead of inventing a score.

Trade-offs are reported explicitly. For example, if a route improves batch
mixing relative to baseline while degrading biological preservation beyond the
configured threshold, the result records a \code{batch_vs_bio} trade-off and
still leaves \code{recommendation = NULL}. The function is observational only:
it compares recorded evidence and never auto-selects a winner.
}
