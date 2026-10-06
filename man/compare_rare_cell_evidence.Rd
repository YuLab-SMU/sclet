\name{compare_rare_cell_evidence}
\alias{compare_rare_cell_evidence}
\title{Compare rare-cell evidence summaries}
\usage{
compare_rare_cell_evidence(object, ids = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
\item{ids}{Optional character vector of completed rare-cell analysis ids to compare. If \code{NULL}, all completed rare-cell runs with evidence are used.}
}
\value{
A bounded list with \code{status}, run summaries, recurring population summaries,
and an explicit matching limitation. Fewer than two usable runs return
\code{status = "not_available"}. The result contains no per-cell values and
always sets \code{raw_values_included = FALSE}.
}
\description{
Compares aggregate evidence nodes recorded by multiple completed rare-cell
runs. Populations are matched only by their recorded label within the same
object. The function reports recurrence, population-size summaries, signal
counts, and claim-level summaries, but it does not infer cell-level overlap,
biological identity, or a preferred rare-cell run. Populations emitted by one
run are not treated as independent evidence lines.
}
