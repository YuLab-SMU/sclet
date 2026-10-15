\name{check_fate_drivers_readiness}
\alias{check_fate_drivers_readiness}
\title{Check whether an object is ready for lineage driver summary}
\usage{
check_fate_drivers_readiness(object, id = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{id}{Optional CellRank trajectory record id. When supplied, that exact
completed record must be present.}
}
\value{
A read-only readiness summary list with elements \code{status}
(\code{"ready_for_diagnostic"} or \code{"not_ready"}), \code{checks}
(CellRank record availability, the resolved record id, lineage driver table
availability, whether the table carries a \code{lineage} column, the resolved
driver metric column name or \code{NULL}, \code{n_lineages}, \code{n_drivers},
and \code{raw_values_included = FALSE}), \code{blocked_actions}
(\code{"fate_drivers"} when not ready), \code{questions} (human-readable
descriptions of missing prerequisites), and \code{notes} (an explicit statement
that no gene name, driver ranking or nominated fate is reported here). No gene
symbol appears anywhere in the result: only counts and column names are
returned.
}
\description{
Reports whether a completed CellRank trajectory record with a lineage driver
table is present. No gene name, driver ranking or nominated fate is reported:
driver identity and the choice of lineage are modeling questions reserved for
the user. Lineage drivers are CellRank model estimates derived from the velocity
kernel, not measured biology. The diagnostic is read-only and removes or
filters no cell.
}
