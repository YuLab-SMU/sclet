\name{check_de_readiness}
\alias{check_de_readiness}
\title{Check whether an object is ready for differential-expression testing}
\usage{
check_de_readiness(object, group = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{group}{Optional group column name. When \code{NULL} it resolves from the
active identity.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks} (counts
assay availability and the available assays, group column resolvability and the
resolved column, \code{n_groups}, \code{smallest_group_n}, existing detest
records and their ids, availability of the differential-expression backend, and
\code{raw_values_included = FALSE}), \code{blocked_actions} (\code{"de"} when
not ready), \code{questions} and \code{notes}.
}
\description{
Reports whether an active group identity resolves to at least two populations
with enough cells to test, and whether a detest record already exists.
}
\details{
The diagnostic guards the silent hazard of this stage: \code{RunDEtest()} with
one population, or with a group too small to have within-group variance, returns
an empty marker table rather than failing, and that emptiness then surfaces far
away as an enrichment run that produced nothing. A single-cell count below two is
reported for the same reason.

No group column, clustering resolution or significance threshold is chosen or
recommended. Marker results are statistical associations under the selected test
and thresholds, not measured quantities, and they are the input to enrichment
rather than a second independent line of evidence for it. No cell is removed or
filtered.
}
\seealso{\code{\link{RunDEtest}}, \code{\link{check_enrichment_readiness}}}
