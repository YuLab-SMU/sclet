\name{compare_spatial_analyses}
\alias{compare_spatial_analyses}
\title{Compare the available spatial analyses}
\usage{
compare_spatial_analyses(object)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}
}
\value{
A read-only comparison list with elements \code{status}
(\code{"available"} or \code{"not_available"}), \code{reason}
(\code{"at_least_two_analyses_required"} when fewer than two spatial
analyses are stored), \code{n_analyses}, \code{analysis_codes} (anonymized
codes matching \code{^(group|cluster|condition|sample|route)_[0-9]+$},
derived from a stable ordering of the analysis keys), \code{analyses}
(per-analysis bounded counts keyed by those codes), \code{claim_level},
\code{caveats}, \code{recommendation = NULL} and
\code{raw_values_included = FALSE}.
}
\description{
Observational comparison across the spatial analyses stored on the object.
The function reports how many analyses are present and how many records and
feature columns each carries, and it never selects a best analysis, never
ranks one above another and never reports which analysis is right:
deconvolution, colocalization and niche detection are different algorithms
answering different questions, so their outputs are not competing answers
to one question. No analysis is preferred, and \code{n_analyses} alone says
nothing about the quality of any of them. The function is read-only, records
no AI evidence, and emits no per-cell or per-spot value.
}
