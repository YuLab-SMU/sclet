\name{check_fate_readiness}
\alias{check_fate_readiness}
\title{Check whether an object is ready for cell fate analysis}
\usage{
check_fate_readiness(object, reduction = NULL, velocity_id = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{reduction}{Optional preferred reduction name. Defaults to UMAP when
available, then any usable embedding.}

\item{velocity_id}{Optional velocity record id. When supplied, that exact
completed velocity record must be present.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks}
(velocity record availability, resolved velocity id, spliced/unspliced or
velocity assay availability, resolved reduction, cluster key resolvability,
available reductions, available velocity ids, and
\code{raw_values_included = FALSE}), \code{blocked_actions} (\code{"fate"}
when not ready), \code{questions} (human-readable descriptions of missing
prerequisites), and \code{notes} (an explicit statement that no terminal
state, fate or lineage is suggested, ranked or scored here, and that fate
probabilities are CellRank model estimates rather than measured cell fate).
}
\description{
Reports whether a completed velocity record, a usable embedding and a cluster
assignment are present for \code{RunCellFate()}. No terminal state, fate or
lineage is proposed, ranked or scored: which states count as terminal, and what
a fate probability means biologically, is a modeling assumption only the user
can make. Fate (absorption) probabilities are CellRank model estimates derived
from the velocity field, not measured cell fate. The diagnostic is read-only
and removes or filters no cell.
}
