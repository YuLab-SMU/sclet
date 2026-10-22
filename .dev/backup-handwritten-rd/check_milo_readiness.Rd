\name{check_milo_readiness}
\alias{check_milo_readiness}
\title{Check whether an object is ready for differential-abundance testing}
\usage{
check_milo_readiness(object, sample_col = NULL, condition_col = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{sample_col}{Optional column name the caller intends as the sample.}

\item{condition_col}{Optional column name the caller intends as the condition.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks} (miloR
availability, PCA availability and the available reductions, the resolved
sample and condition columns, \code{design_roles_confirmed},
\code{design_roles_supplied_by_caller}, \code{n_samples},
\code{n_conditions}, existing Milo records and their ids, and
\code{raw_values_included = FALSE}), \code{blocked_actions} (\code{"milo"} when
not ready), \code{questions} and \code{notes}.
}
\description{
Reports whether miloR is installed, whether a PCA embedding exists, and whether
the sample and condition design roles are both known and confirmed.
}
\details{
Neither design role is ever guessed. A role counts as known only when the
caller supplies it or when an earlier
\code{\link{ConfirmAIDesignSemantics}} call declared it, because which colData
column "is the sample" is the user's own statement about the experiment.
Readiness deliberately agrees with the \code{run_milo} gate: columns that
resolve to real, testable data but were never confirmed produce
\code{"not_ready"}, since a diagnostic green-lighting an action the gate would
block is worse than one that blocks -- a planning agent acts on the green
light.

Differential abundance results are associations between neighbourhoods and
conditions under a chosen graph and test, not causal or measured effects. Only
counts are reported here; no per-neighbourhood, per-cell or sample-label value
is returned.
}
\seealso{\code{\link{ConfirmAIDesignSemantics}}, \code{\link{RunMilo}}}
