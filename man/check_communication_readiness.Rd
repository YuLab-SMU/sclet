\name{check_communication_readiness}
\alias{check_communication_readiness}
\title{Check whether an object is ready for cell-cell interaction analysis}
\usage{
check_communication_readiness(object, method = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{method}{Optional engine name, one of \code{"CellChat"},
\code{"CellPhoneDB"}, \code{"NicheNet"}. When supplied, that engine's package
must be available.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks} (group
column resolvability and the resolved column, expression source availability,
available assays, the engines whose packages are installed, existing
communication records and their ids, and \code{raw_values_included = FALSE}),
\code{blocked_actions} (\code{"communication"} when not ready),
\code{questions}, and \code{notes}.
}
\description{
Reports whether a resolvable group column, an expression source and at least one
installed CCI engine are present for \code{RunCCI()}.
}
\details{
No engine is recommended or ranked. CellChat, CellPhoneDB and NicheNet differ in
their ligand-receptor database, their statistical test and their sender/receiver
framing, so the choice is a modelling decision reserved for the user; when
several engines are possible and none was named, the questions report only the
count of available engines. Inferred interactions are model estimates derived
from expression plus a curated database, not observed signalling events. Two
engines run over the same expression matrix and the same group labels are not
independent measurements, and the evidence layer enforces that by sharing a
dependency group across communication runs.
}
\seealso{\code{\link{summarize_communication_evidence}}, \code{\link{RunCCI}}}
