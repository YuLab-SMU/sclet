\name{check_velocity_readiness}
\alias{check_velocity_readiness}
\title{Check whether an object is ready for RNA velocity analysis}
\usage{
check_velocity_readiness(object, reduction = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{reduction}{Optional preferred reduction name. Defaults to UMAP when
available, then any usable embedding.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks}
(spliced/unspliced assay availability, resolved reduction, available
reductions), \code{blocked_actions} (\code{"velocity"} when not ready),
\code{questions} (human-readable descriptions of missing prerequisites),
and \code{notes} (an explicit statement that no velocity mode or design
confirmation is recommended at the readiness stage).
}
\description{
Reports whether the \code{"spliced"} and \code{"unspliced"} assays and a
usable embedding are present for \code{RunVelocity()}. No velocity mode
(deterministic/stochastic/dynamical) is recommended here: that is a
numerical-method choice for \code{run_velocity}'s own params, not a
readiness gate. No design confirmation is required for this domain,
because velocity has no sample/batch/root semantic parameter the way
trajectory or integration do.
}
