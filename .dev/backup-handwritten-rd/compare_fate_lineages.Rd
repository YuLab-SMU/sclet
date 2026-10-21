\name{compare_fate_lineages}
\alias{compare_fate_lineages}
\title{Compare lineage drivers across lineages}
\usage{
compare_fate_lineages(object, id = NULL)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{id}{Optional CellRank trajectory record id.}
}
\value{
A read-only comparison list with elements \code{status}, \code{n_lineages},
\code{lineage_codes} (anonymized \code{route_N} codes only), \code{lineages}
(per-lineage summaries keyed by \code{route_N} code, each with \code{n_drivers},
\code{n_drivers_summarized} and the \code{metric_mean}, \code{metric_median} and
\code{metric_max} of the absolute metric), \code{driver_metric},
\code{distinctness} (how much the lineages' summarized driver sets differ,
expressed only as numeric Jaccard overlap fractions), \code{claim_level},
\code{caveats}, \code{recommendation = NULL} and
\code{raw_values_included = FALSE}.

Fewer than two lineages with a driver table return
\code{status = "not_available"} with \code{reason =
"at_least_two_lineages_required"}. No gene symbol, no per-cell value and no
\code{"best lineage"} field is ever emitted.
}
\description{
Observationally compares the stored CellRank driver table across at least two
lineages: how many drivers each anonymized lineage carries, how large their
absolute metric values are, and how much the lineages' driver sets overlap.
The comparison never nominates, ranks or describes any lineage as \dQuote{the}
fate of the dataset, never reads per-cell values, and writes no evidence node.
Overlapping driver sets are expected between related lineages, so a high overlap
fraction is not evidence that the lineages are biologically identical.
}
