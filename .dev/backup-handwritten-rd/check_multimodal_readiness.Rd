\name{check_multimodal_readiness}
\alias{check_multimodal_readiness}
\title{Check whether an object is ready for multimodal diagnostics}
\usage{
check_multimodal_readiness(object, modalities = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{modalities}{Optional character vector of altExp names the caller is
asking about. Used only to report which of the requested names are present;
no requested name is ever reported as missing-in-error, and no name outside
the caller's own request enters the result.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks},
\code{blocked_actions} (\code{"multimodal"} when not ready), \code{questions}
(human-readable descriptions that distinguish an unattached modality from an
attached-but-unregistered one) and \code{notes}.

\code{checks} contains \code{alt_exp_names} (the altExps that actually exist),
\code{requested_modalities} and \code{requested_modalities_present},
\code{adt_available} and \code{atac_available} (from
\code{has_multimodal_adt()} and \code{has_multimodal_atac()}),
\code{n_modalities_available}, \code{adt_registered}, \code{atac_registered},
\code{workflow_record_available} and \code{raw_values_included = FALSE}.

\code{status} is \code{"not_ready"} only when no modality is available
\strong{and} no multimodal record exists.
}
\description{
Reports which multimodal modalities are attached to the object as altExps,
which of them carry a multimodal analysis-state record, and whether a
multimodal workflow record is present.

Attached and registered are different states. A modality is \emph{attached}
when an altExp such as \code{"ADT"} or \code{"ATAC"} exists on the object; it is
\emph{registered} when a \code{"multimodal"} analysis-state record exists for
it. Registration makes the modality trackable by the mainline and computes
nothing, so the two branches of \code{questions} read differently: an
unattached modality is reported as needing
\code{altExp(sce, "ADT") <- adt_sce}, while an attached but unregistered
modality is reported as needing \code{register_multimodal_adt(sce)} and
explicitly not as needing data that is already present to be re-attached.

The diagnostic never says which modalities an object \emph{should} have, never
recommends an integration method, never ranks or scores a modality, and removes
or filters no cell. Which modalities are expected, and how they should be
combined, is a modeling assumption only the user can make.
}
