\name{check_spatial_readiness}
\alias{check_spatial_readiness}
\title{Check whether an object is ready for spatial diagnostics}
\usage{
check_spatial_readiness(object, analysis = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{analysis}{Optional analysis key to report on. One of
\code{"spatial_deconv"}, \code{"colocalization"} or
\code{"spatial_niche"}. When \code{NULL}, every analysis present is
reported.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks}
(spatial record availability, the available analysis keys, the number of
spatial state records, whether spatial coordinates were found in the
modality, whether a \code{counts} assay is present, the requested analysis
and whether it is available, and \code{raw_values_included = FALSE}),
\code{blocked_actions} (\code{"spatial"} when not ready), \code{questions}
(human-readable descriptions of missing prerequisites, including the
\code{RunSpatialDeconvolution()}, \code{RunSpatialColocalization()} and
\code{RunSpatialNiche()} entry points when nothing is stored), and
\code{notes} (an explicit statement that no spatial method is recommended,
ranked or preferred, and that deconvolution proportions and niche
assignments are model estimates rather than measured cell-type composition
or ground-truth tissue architecture).
}
\description{
Reports what spatial analyses are stored on the object and whether the
modality exposes the prerequisites a spatial record needs. The result
contains counts, availability flags and analysis keys only: no spatial
method is recommended, ranked or preferred, because deconvolution,
colocalization and niche detection answer different questions and the
choice belongs to the user. The diagnostic is read-only, records no AI
evidence, and removes or filters no cell.
}
