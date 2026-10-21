\name{summarize_fate_lineage_drivers}
\alias{summarize_fate_lineage_drivers}
\title{Summarize the lineage drivers of one user-specified lineage}
\usage{
summarize_fate_lineage_drivers(object, id = NULL, lineage = NULL, top_n = 10L)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{id}{Optional CellRank trajectory record id.}

\item{lineage}{The lineage to summarize. When \code{NULL} and several lineages
are available, no lineage is chosen for the user and the result asks for an
explicit choice.}

\item{top_n}{Maximum number of drivers to summarize. Defaults to 10.}
}
\value{
A read-only driver summary list. When drivers are available it contains
\code{status = "available"}, \code{lineage_code} (a single anonymized
\code{route_N} code), \code{n_drivers_total}, \code{n_drivers_summarized},
\code{driver_metric} (the metric column name as a plain identifier such as
\code{"corr"}), \code{driver_rank_codes} (anonymized \code{route_N} codes for
the summarized drivers, ordered by absolute metric descending),
\code{driver_abs_metric_summary} (bounded \code{n_observed}, \code{mean},
\code{quantile_25}, \code{median}, \code{quantile_75} and \code{max} of the
absolute metric over the summarized drivers), \code{claim_level = "associated"},
\code{caveats}, and \code{raw_values_included = FALSE}.

When the lineage is ambiguous the result has
\code{status = "clarification_required"} and a \code{question} naming only the
count of available lineages. When no record or no driver table is available the
result has \code{status = "not_available"} with \code{reason} one of
\code{"no_cellrank_record"}, \code{"no_lineage_driver_table"} or
\code{"lineage_unknown"}.

The function is read-only: it returns a summary list and does not modify the
input object. No gene symbol and no per-cell value is ever returned.
}
\description{
Reads the stored CellRank lineage driver table for a single lineage named by the
user, ranks its drivers by the absolute value of the resolved driver metric, and
returns a bounded summary. Gene symbols exist in the underlying table but never
leave this function: only counts, the metric column name, anonymized
\code{route_N} codes and aggregate statistics of the absolute metric are
returned. The function never selects a lineage on the user's behalf, and records
no evidence node.
}
